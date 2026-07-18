"""Experimental aether grid built on the public ``fdtd`` package.

``AetherGrid`` is an opt-in companion to :class:`fdtd.grid.Grid`. It keeps the
same educational package style:

- the grid owns the timestep loop;
- sources, boundaries, detectors, and objects register through ``grid[...] =``;
- and all numerics go through the shared backend abstraction.

The theory path is different from Maxwell, however. Instead of updating only
``E`` and ``H`` with first-order curl equations, the aether path evolves a
velocity-like linear state together with a separate angular sector and derived
scalar, electromagnetic, acceleration, and jerk fields.

This module is still experimental. It currently uses collocated vector fields
for clarity and package compatibility even though the longer-term STPT program
likely requires a more explicit split between linear and angular sectors. The
current refactor makes that split visible in the code structure while still
using the legacy bridge where angular quantities are derived from the linear
sector.
"""

## Imports

# standard library
import os
from os import path, makedirs, chdir, remove
from subprocess import check_call, CalledProcessError
from glob import glob
from datetime import datetime

# 3rd party
from tqdm import tqdm
from numpy import savez

# typing
from .typing_ import Tuple, Number, Tensorlike

# relative
from .backend import backend as bd

from .operators import (
    angular_clock_eigenvalue,
    angular_to_linear_bridge,
    div,
    grad,
    linear_to_angular_bridge,
)


from math import pi

## Constants

# base constants

c       = 299792458.0       # speed of light                [m/s]

eta     = 1/(4*pi*1e-7)     # viscosity (1/mu_0)            [kg/m-s],   [Pa-s]

h       = 6.62607015e-34    # Planck's constant             [kg-m^2/s], [J-s]

e       = 1.602176634e-19   # elementary charge             [kg/s]

# derived constants

k       = c**2              # quantum circulation constant
                            # 8.987551787368176e+16         [m^2/s]

rho     = eta/k             # mass density (eps_0) 
                            # 8.85418781762039e-12          [kg/m^3]

m       = h/k               # elementary mass 
                            # 7.372497323812708e-51         [kg]
                      
rho_q0  = e/m * rho         # vacuum charge density 
                            # 1.9241747011042014e+20        [kg/m^3-s]

eta_e   = eta/e
e_eta   = e/eta
inv_rho = 1/rho

inv_rho_q0 = 1/rho_q0      


## FDTD Grid Class
class AetherGrid:
    """Experimental aether grid that mirrors the public ``fdtd.Grid`` API."""

    from .visualization import visualize

    def __init__(
        self,
        shape: Tuple[Number, Number, Number],
        grid_spacing: float = 155e-9,
        permittivity: float = 1.0,
        permeability: float = 1.0,
        courant_number: float = None,
        native_angular_transport: bool = False,
        native_angular_collect_sources: bool = True,
        native_angular_update_clocks: bool = False,
        native_angular_feedback_mode=None,
        native_angular_feedback_scale: float = 1.0,
    ):
        """
        Args:
            shape: shape of the FDTD grid.
            grid_spacing: distance between the grid cells.
            permittivity: the relative permittivity of the background.
            permeability: the relative permeability of the background.
            courant_number: the courant number of the FDTD simulation.
                Defaults to the inverse of the square root of the number of
                dimensions > 1 (optimal value). The timestep of the simulation
                will be derived from this number using the CFL-condition.
            native_angular_transport: opt-in native angular momentum transport
                inside ``update()``. Defaults to false.
            native_angular_collect_sources: collect opt-in native angular
                source and boundary hooks when native transport runs.
            native_angular_update_clocks: update ``omega_t`` and ``omega_p``
                from positive inertia after native transport runs.
            native_angular_feedback_mode: optional opt-in feedback mode for
                applying ``native_angular_linear_response`` to ``linear_a``.
                Use ``"add"``, ``"replace"``, or ``None``.
            native_angular_feedback_scale: scale used when projecting native
                angular stress for opt-in feedback.
        """
        if native_angular_feedback_mode not in (None, "add", "replace"):
            raise ValueError(
                "native_angular_feedback_mode must be None, 'add', or 'replace'"
            )

        # save the grid spacing
        self.grid_spacing = float(grid_spacing)
        self.native_angular_transport_enabled = bool(native_angular_transport)
        self.native_angular_collect_sources = bool(native_angular_collect_sources)
        self.native_angular_update_clocks = bool(native_angular_update_clocks)
        self.native_angular_feedback_mode = native_angular_feedback_mode
        self.native_angular_feedback_scale = float(native_angular_feedback_scale)

        # save grid shape as integers
        self.Nx, self.Ny, self.Nz = self._handle_tuple(shape)

        # dimension of the simulation:
        self.D = int(self.Nx > 1) + int(self.Ny > 1) + int(self.Nz > 1)

        # courant number of the simulation (optimal value)
        max_courant_number = float(self.D) ** (-0.5)
        if courant_number is None:
            # slight stability factor added
            self.courant_number = 0.99 * max_courant_number
        elif courant_number > max_courant_number:
            raise ValueError(
                f"courant_number {courant_number} too high for "
                f"a {self.D}D simulation"
            )
        else:
            self.courant_number = float(courant_number)
            
        """
        For now, we assume our fields to propagate with a maximum speed of pi/2
        * c, since we assume Tesla's longitudinal sound-like waves to also exist
        and to propagate at that speed, even though Dr. Steffen Kuehn has
        demonstrated the transmission of information in electrically short
        coaxial cables at speeds of up to 3c, which is a strong argument to be
        made that the speed of light is not the maximum speed of propagation of
        dielectric wave phenomena nor information. 
        
        It is these dielectric wave phenomena that are distincly different from
        the electromagnetic waves we are familiar with that appear to have been
        overlooked by science, even though Tesla demonstrated the transmission
        and reception of telluric currents with a speed of pi/2*c, even
        wirelessly powering light bulbs at a distance of over half a mile from
        his laboratory in Colorado Springs, according to Hugo Gernsback.
        """
        # timestep of the simulation
        self.time_step = self.courant_number * self.grid_spacing / ((pi/2) * c)
        
        # The linear and angular sectors currently share the same sampled grid
        # spacing, but we keep both handles explicit because the angular branch
        # carries its own lever-arm length scale in the intended geometry.
        self.linear_grid_spacing = self.grid_spacing
        self.angular_grid_spacing = self.grid_spacing
        self.angular_metric_length = (
            bd.ones((self.Nx, self.Ny, self.Nz, 1), dtype=bd.float)
            * self.angular_grid_spacing
        )

        # define linear-sector fields
        self.linear_v = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.linear_f = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_E = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_a = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_dpdt = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.linear_yank = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_dEdt = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_j = bd.zeros((self.Nx, self.Ny, self.Nz, 3))

        # define angular-sector fields
        self.angular_omega = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_tau = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_H = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_A = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_alpha = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_dtau_dt = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_dHdt = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_dAdt = bd.zeros((self.Nx, self.Ny, self.Nz, 3))

        # Native angular-clock placeholders. These are not yet coupled into
        # the bridge dynamics; they expose the intended two-clock state for
        # analytic benchmarks and future opt-in angular-sector work.
        self.omega_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.omega_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_gamma = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.theta_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.theta_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_chi = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_clock_lambda = bd.zeros((self.Nx, self.Ny, self.Nz, 1))

        # Native angular geometry placeholders. The default local angular
        # frame uses x/y host-grid axes as a neutral basis until a true
        # cell-local orientation update is derived.
        self.angular_e_t = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_e_p = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.angular_e_t[:, :, :, 0] = 1.0
        self.angular_e_p[:, :, :, 1] = 1.0
        self.ell_t = (
            bd.ones((self.Nx, self.Ny, self.Nz, 1), dtype=bd.float)
            * self.angular_grid_spacing
        )
        self.ell_p = (
            bd.ones((self.Nx, self.Ny, self.Nz, 1), dtype=bd.float)
            * self.angular_grid_spacing
        )
        self.angular_inertia_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_inertia_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_momentum_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_momentum_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_torque_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.angular_torque_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_tau = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.native_angular_tau_divergence = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_tau_metric_divergence_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_tau_metric_divergence_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_transport_residual_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_transport_residual_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_source_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_source_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_exchange_power_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_exchange_power_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_momentum_rhs_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_momentum_rhs_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_transport_power_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_transport_power_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_momentum_candidate_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_momentum_candidate_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_linear_response = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.linear_charge_flux_candidate = bd.zeros((self.Nx, self.Ny, self.Nz, 3))
        self.native_angular_charge_flux_t = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_charge_flux_p = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_charge_reduction = bd.zeros((self.Nx, self.Ny, self.Nz, 1))
        self.native_angular_boundary_normal_flux = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_incident_flux = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_outgoing_flux = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_incident_flux_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_incident_flux_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_outgoing_flux_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_outgoing_flux_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_direct_normal_flux_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_direct_normal_flux_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_direct_incident_flux_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_direct_incident_flux_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_direct_outgoing_flux_t = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )
        self.native_angular_boundary_direct_outgoing_flux_p = bd.zeros(
            (self.Nx, self.Ny, self.Nz, 1)
        )

        self._sync_public_aliases()
        
        # save the inverse of the relative permittiviy and the relative permeability
        # these tensors can be anisotropic!

        if bd.is_array(permittivity) and len(permittivity.shape) == 3:
            permittivity = permittivity[:, :, :, None]
        self.inverse_permittivity = bd.ones((self.Nx, self.Ny, self.Nz, 3)) / bd.array(
            permittivity, dtype=bd.float
        )

        if bd.is_array(permeability) and len(permeability.shape) == 3:
            permeability = permeability[:, :, :, None]
        self.inverse_permeability = bd.ones((self.Nx, self.Ny, self.Nz, 3)) / bd.array(
            permeability, dtype=bd.float
        )

        # save current time index
        self.time_steps_passed = 0

        # dictionary containing the sources:
        self.sources = []

        # dictionary containing the boundaries
        self.boundaries = []

        # dictionary containing the detectors
        self.detectors = []

        # dictionary containing the objects in the grid
        self.objects = []

        # folder path to store the simulation
        self.folder = None

    def _handle_distance(self, distance: Number) -> int:
        """transform a distance to an integer number of gridpoints"""
        if not isinstance(distance, int):
            return int(float(distance) / self.grid_spacing + 0.5)
        return distance

    def _handle_time(self, time: Number) -> int:
        """transform a time value to an integer number of timesteps"""
        if not isinstance(time, int):
            return int(float(time) / self.time_step + 0.5)
        return time

    def _handle_tuple(
        self, shape: Tuple[Number, Number, Number]
    ) -> Tuple[int, int, int]:
        """validate the grid shape and transform to a length-3 tuple of ints"""
        if len(shape) != 3:
            raise ValueError(
                f"invalid grid shape {shape}\n"
                f"grid shape should be a 3D tuple containing floats or ints"
            )
        x, y, z = shape
        x = self._handle_distance(x)
        y = self._handle_distance(y)
        z = self._handle_distance(z)
        return x, y, z

    def _handle_slice(self, s: slice) -> slice:
        """validate the slice and transform possibly float values to ints"""
        start = (
            s.start
            if not isinstance(s.start, float)
            else self._handle_distance(s.start)
        )
        stop = (
            s.stop if not isinstance(s.stop, float) else self._handle_distance(s.stop)
        )
        step = (
            s.step if not isinstance(s.step, float) else self._handle_distance(s.step)
        )
        return slice(start, stop, step)

    def _handle_single_key(self, key):
        """transform a single index key to a slice or list"""
        try:
            len(key)
            return [self._handle_distance(k) for k in key]
        except TypeError:
            if isinstance(key, slice):
                return self._handle_slice(key)
            else:
                return [self._handle_distance(key)]
        return key

    @property
    def x(self) -> int:
        """get the number of grid cells in the x-direction"""
        return self.Nx * self.grid_spacing

    @property
    def y(self) -> int:
        """get the number of grid cells in the y-direction"""
        return self.Ny * self.grid_spacing

    @property
    def z(self) -> int:
        """get the number of grid cells in the y-direction"""
        return self.Nz * self.grid_spacing

    @property
    def shape(self) -> Tuple[int, int, int]:
        """get the shape of the FDTD grid"""
        return (self.Nx, self.Ny, self.Nz)

    @property
    def time_passed(self) -> float:
        """get the total time passed"""
        return self.time_steps_passed * self.time_step

    def _sync_public_aliases(self):
        """Expose legacy field names as views onto explicit sector state."""

        self.v = self.linear_v
        self.p = self.linear_p
        self.f = self.linear_f
        self.E = self.linear_E
        self.a = self.linear_a
        self.dpdt = self.linear_dpdt
        self.yank = self.linear_yank
        self.dEdt = self.linear_dEdt
        self.j = self.linear_j

        self.omega = self.angular_omega
        self.tau = self.angular_tau
        self.H = self.angular_H
        self.A = self.angular_A
        self.alpha = self.angular_alpha
        self.dtau_dt = self.angular_dtau_dt
        self.dHdt = self.angular_dHdt
        self.dAdt = self.angular_dAdt

    def run(self, total_time: Number, progress_bar: bool = True):
        """run an FDTD simulation.

        Args:
            total_time: the total time for the simulation to run.
            progress_bar: choose to show a progress bar during
                simulation

        """
        if isinstance(total_time, float):
            total_time /= self.time_step
        time = range(0, int(total_time), 1)
        if progress_bar:
            time = tqdm(time)
        for _ in time:
            self.step()
    
    def step(self):
        """Advance the experimental aether system by one timestep."""

        self.update()
        
        self.time_steps_passed += 1
        
    def update(self):
        """Update all experimental aether fields once.

        The current chain is intentionally explicit:

        1. apply any native aether sources to the linear sector;
        2. derive first-order linear fields from the translational state;
        3. derive the current temporary angular bridge from the linear sector;
        4. expose ``E`` and ``H``-like observables and run package-style scene
           hooks;
        5. couple the angular sector back into linear acceleration;
        6. optionally run explicit native-angular transport / feedback;
        7. repeat the same pattern one derivative level higher to obtain jerk;
        8. advance the linear state with the resulting Taylor step.
        """

        self.updateBoundaries()

        self.apply_native_sources()
        self.update_linear_sector()
        self.update_angular_sector()
        self._sync_public_aliases()

        # Existing package hooks still operate on E/H-like observables.
        self.updateEH()
        self.update_linear_angular_coupling()
        self.update_native_angular_opt_in_step()
        self.detect_native_angular()
        self.update_second_order_linear_sector()
        self.update_second_order_angular_sector()
        self.update_second_order_coupling()
        self._sync_public_aliases()
        self.advance_linear_sector()

    def apply_native_sources(self):
        """Inject any aether-native sources into the linear translational sector."""

        for src in self.sources:
            update_v = getattr(src, "update_v", None)
            if update_v is not None:
                update_v()

    def update_linear_sector(self):
        """Update first-order linear quantities from the linear state."""

        self.linear_p = eta * div(self.linear_v)
        self.linear_f = -grad(self.linear_p)
        self.linear_E = inv_rho_q0 * self.linear_f

    def update_angular_sector(self):
        """Update the angular sector from the current temporary bridge.

        The long-term design goal is a native angular state with its own
        constitutive update. For now this sector is still derived from the
        linear branch so we can refactor the grid structure without discarding
        the current runnable baseline.
        """

        self.angular_omega = linear_to_angular_bridge(self.linear_v)
        self.angular_tau = eta * self.angular_omega
        self.angular_H = eta_e * self.angular_tau

    def update_linear_angular_coupling(self):
        """Couple first-order angular content back into linear acceleration."""

        self.angular_A = angular_to_linear_bridge(e_eta * self.angular_H)
        self.linear_a = rho_q0 * self.linear_E + inv_rho * self.angular_A

    def update_native_angular_opt_in_step(self):
        """Run optional native-angular transport and feedback inside update().

        This hook is inert unless ``AetherGrid`` was constructed with
        ``native_angular_transport=True`` or a non-``None``
        ``native_angular_feedback_mode``. It is the first production-lifecycle
        integration point for the native-angular staging helpers, while still
        keeping the default aether bridge unchanged.
        """

        if self.native_angular_transport_enabled:
            self.advance_native_angular_momentum_transport(
                delta=self.time_step,
                collect_sources=self.native_angular_collect_sources,
                update_clocks=self.native_angular_update_clocks,
            )

        if self.native_angular_feedback_mode is not None:
            self.apply_native_angular_linear_response(
                metric_length=self.angular_metric_length,
                scale=self.native_angular_feedback_scale,
                mode=self.native_angular_feedback_mode,
            )

    def update_second_order_linear_sector(self):
        """Update second-order linear quantities from the linear acceleration."""

        self.linear_dpdt = eta * div(self.linear_a)
        self.linear_yank = -grad(self.linear_dpdt)
        self.linear_dEdt = inv_rho_q0 * self.linear_yank

    def update_second_order_angular_sector(self):
        """Update second-order angular quantities from the current bridge."""

        self.angular_alpha = linear_to_angular_bridge(self.linear_a)
        self.angular_dtau_dt = eta * self.angular_alpha
        self.angular_dHdt = eta_e * self.angular_dtau_dt

    def update_second_order_coupling(self):
        """Couple second-order angular content back into linear jerk."""

        self.angular_dAdt = angular_to_linear_bridge(e_eta * self.angular_dHdt)
        self.linear_j = rho_q0 * self.linear_dEdt + inv_rho * self.angular_dAdt

    def evaluate_angular_clock_benchmark(self, delta=None):
        """Evaluate the native angular-clock eigenvalue benchmark.

        The current aether bridge still derives ``angular_omega`` from the
        linear state. This method is deliberately separate: it lets tests and
        future native angular-sector experiments populate ``omega_t``,
        ``omega_p``, and an independent ``angular_gamma`` benchmark channel
        and evaluate the benchmark without altering the timestep update or the
        classic Maxwell path.
        """

        if delta is None:
            delta = self.time_step

        self.theta_t = self.omega_t * delta
        self.theta_p = self.omega_p * delta
        self.angular_chi = self.angular_gamma * delta
        self.angular_clock_lambda = angular_clock_eigenvalue(
            self.omega_t,
            self.omega_p,
            self.angular_gamma,
            delta,
        )
        return self.angular_clock_lambda

    def update_native_angular_inertia(self, density=None):
        """Populate the passive native angular inertia placeholders.

        Args:
            density: Optional scalar or array-like mass density ``[kg/m^3]``.
                Defaults to the package-level reference density used by the
                current aether bridge.

        Returns:
            ``(angular_inertia_t, angular_inertia_p)``.

        Notes:
            The provisional closure is ``I = rho * ell**2``. Its units are
            ``[kg/m]`` for a per-volume moment-density style quantity. This
            gives the angular metric lengths a testable constitutive meaning
            without coupling them into the timestep update yet.
        """

        density = rho if density is None else bd.asarray(density)
        self.angular_inertia_t = density * self.ell_t * self.ell_t
        self.angular_inertia_p = density * self.ell_p * self.ell_p
        return self.angular_inertia_t, self.angular_inertia_p

    def update_native_angular_momentum(self):
        """Populate passive native angular momentum placeholders.

        Returns:
            ``(angular_momentum_t, angular_momentum_p)``.

        Notes:
            The provisional closure is ``L = I * omega`` on the two native
            angular clocks. It is intentionally not a torque law; torque would
            require a time update or flux law for this momentum state.
        """

        self.angular_momentum_t = self.angular_inertia_t * self.omega_t
        self.angular_momentum_p = self.angular_inertia_p * self.omega_p
        return self.angular_momentum_t, self.angular_momentum_p

    def evaluate_native_angular_kinetic_energy(
        self,
        momentum_t=None,
        momentum_p=None,
        inertia_t=None,
        inertia_p=None,
    ):
        """Evaluate passive native angular kinetic-energy diagnostics.

        Args:
            momentum_t: Optional toroidal momentum field. Defaults to
                ``angular_momentum_t``.
            momentum_p: Optional poloidal momentum field. Defaults to
                ``angular_momentum_p``.
            inertia_t: Optional toroidal inertia field. Defaults to
                ``angular_inertia_t``.
            inertia_p: Optional poloidal inertia field. Defaults to
                ``angular_inertia_p``.

        Returns:
            ``(energy_t, energy_p, total_energy)`` using
            ``E_L = 0.5*L**2/I`` per native angular channel.

        Notes:
            This is diagnostic bookkeeping for candidate tests. It is not a
            Hamiltonian, not a detector-facing observable, and not coupled into
            ``step()``.
        """

        momentum_t = (
            self.angular_momentum_t if momentum_t is None else bd.asarray(momentum_t)
        )
        momentum_p = (
            self.angular_momentum_p if momentum_p is None else bd.asarray(momentum_p)
        )
        inertia_t = self.angular_inertia_t if inertia_t is None else bd.asarray(inertia_t)
        inertia_p = self.angular_inertia_p if inertia_p is None else bd.asarray(inertia_p)

        if bd.max(inertia_t <= 0) or bd.max(inertia_p <= 0):
            raise ValueError("native angular inertia must be positive")

        energy_t = 0.5 * momentum_t * momentum_t / inertia_t
        energy_p = 0.5 * momentum_p * momentum_p / inertia_p
        total_energy = bd.sum(energy_t + energy_p)
        return energy_t, energy_p, total_energy

    def evaluate_native_angular_exchange_power(
        self,
        source_t=None,
        source_p=None,
        momentum_t=None,
        momentum_p=None,
        inertia_t=None,
        inertia_p=None,
    ):
        """Evaluate native angular source power diagnostics.

        Args:
            source_t: Optional toroidal source/exchange term. Defaults to
                ``native_angular_source_t``.
            source_p: Optional poloidal source/exchange term. Defaults to
                ``native_angular_source_p``.
            momentum_t: Optional toroidal momentum field. Defaults to
                ``angular_momentum_t``.
            momentum_p: Optional poloidal momentum field. Defaults to
                ``angular_momentum_p``.
            inertia_t: Optional toroidal inertia field. Defaults to
                ``angular_inertia_t``.
            inertia_p: Optional poloidal inertia field. Defaults to
                ``angular_inertia_p``.

        Returns:
            ``(power_t, power_p, total_power)`` using
            ``P_L = S_L * L/I`` in each native angular channel.

        Notes:
            This diagnostic classifies explicit source or boundary exchange as
            energy injecting, neutral, or dissipative. It does not advance
            momentum, define a Hamiltonian, or make boundary hooks physical.
        """

        if source_t is None:
            source_t = self.native_angular_source_t
        if source_p is None:
            source_p = self.native_angular_source_p
        source_t, source_p = self._validate_native_angular_source_terms(
            source_t,
            source_p,
        )
        momentum_t = (
            self.angular_momentum_t if momentum_t is None else bd.asarray(momentum_t)
        )
        momentum_p = (
            self.angular_momentum_p if momentum_p is None else bd.asarray(momentum_p)
        )
        inertia_t = self.angular_inertia_t if inertia_t is None else bd.asarray(inertia_t)
        inertia_p = self.angular_inertia_p if inertia_p is None else bd.asarray(inertia_p)

        if bd.max(inertia_t <= 0) or bd.max(inertia_p <= 0):
            raise ValueError("native angular inertia must be positive")

        self.native_angular_exchange_power_t = source_t * momentum_t / inertia_t
        self.native_angular_exchange_power_p = source_p * momentum_p / inertia_p
        total_power = bd.sum(
            self.native_angular_exchange_power_t
            + self.native_angular_exchange_power_p
        )
        return (
            self.native_angular_exchange_power_t,
            self.native_angular_exchange_power_p,
            total_power,
        )

    def update_native_angular_torque(
        self,
        previous_momentum_t=None,
        previous_momentum_p=None,
        delta=None,
    ):
        """Populate passive native torque placeholders from momentum changes.

        Args:
            previous_momentum_t: Previous toroidal angular momentum state.
                Defaults to zero.
            previous_momentum_p: Previous poloidal angular momentum state.
                Defaults to zero.
            delta: Time interval for the finite difference. Defaults to the
                grid timestep.

        Returns:
            ``(angular_torque_t, angular_torque_p)``.

        Notes:
            This is the time-update interpretation of native torque,
            ``tau = dL/dt``. It is not yet a spatial flux law or a coupling law
            back into the linear branch.
        """

        if delta is None:
            delta = self.time_step
        if delta == 0:
            raise ValueError("delta must be non-zero")

        if previous_momentum_t is None:
            previous_momentum_t = bd.zeros_like(self.angular_momentum_t)
        else:
            previous_momentum_t = bd.asarray(previous_momentum_t)
        if previous_momentum_p is None:
            previous_momentum_p = bd.zeros_like(self.angular_momentum_p)
        else:
            previous_momentum_p = bd.asarray(previous_momentum_p)

        self.angular_torque_t = (
            self.angular_momentum_t - previous_momentum_t
        ) / delta
        self.angular_torque_p = (
            self.angular_momentum_p - previous_momentum_p
        ) / delta
        return self.angular_torque_t, self.angular_torque_p

    def configure_native_angular_frame(
        self,
        e_t,
        e_p,
        require_orthogonal=True,
    ):
        """Set the native angular frame explicitly.

        Args:
            e_t: Toroidal frame vector, either shape ``(3,)`` or
                ``(Nx, Ny, Nz, 3)``.
            e_p: Poloidal frame vector, either shape ``(3,)`` or
                ``(Nx, Ny, Nz, 3)``.
            require_orthogonal: When true, reject frames where ``e_t`` and
                ``e_p`` have a nonzero local dot product.

        Returns:
            ``(angular_e_t, angular_e_p)``.

        Notes:
            This only configures the bridge projection frame. It does not
            derive a physical local angular orientation law.
        """

        expected_shape = self.angular_e_t.shape

        def _expand_frame(name, value):
            value = bd.asarray(value)
            if value.shape == (3,):
                value = bd.ones(expected_shape) * value
            if value.shape != expected_shape:
                raise ValueError(
                    f"{name} must have shape (3,) or {expected_shape}, "
                    f"got {value.shape}"
                )
            return value

        e_t = _expand_frame("e_t", e_t)
        e_p = _expand_frame("e_p", e_p)
        e_t_np = bd.numpy(e_t)
        e_p_np = bd.numpy(e_p)
        norm_t = (e_t_np * e_t_np).sum(axis=-1)
        norm_p = (e_p_np * e_p_np).sum(axis=-1)
        if (norm_t <= 0.0).any() or (norm_p <= 0.0).any():
            raise ValueError("native angular frame vectors must be nonzero")
        if require_orthogonal:
            dot = (e_t_np * e_p_np).sum(axis=-1)
            if (abs(dot) > 1.0e-12).any():
                raise ValueError("native angular frame vectors must be orthogonal")

        self.angular_e_t = e_t
        self.angular_e_p = e_p
        return self.angular_e_t, self.angular_e_p

    def project_native_angular_torque(self):
        """Project the two native torque channels onto the local angular frame.

        Returns:
            A vector field ``native_angular_tau`` with shape ``(Nx, Ny, Nz, 3)``.

        Notes:
            This is only a projection from native scalar channels to the host
            vector layout:

            ``tau_native = tau_t * e_t + tau_p * e_p``.

            It does not replace the existing bridge ``angular_tau`` and does
            not couple torque back into the linear branch.
        """

        self.native_angular_tau = (
            self.angular_torque_t * self.angular_e_t
            + self.angular_torque_p * self.angular_e_p
        )
        return self.native_angular_tau

    def evaluate_native_angular_torque_divergence(self):
        """Evaluate a passive spatial diagnostic for native angular torque.

        Returns:
            Scalar field ``native_angular_tau_divergence`` with shape
            ``(Nx, Ny, Nz, 1)``.

        Notes:
            This contracts the already projected ``native_angular_tau`` field
            with the current bridge ``div`` operator. It is a diagnostic for
            spatial imbalance in the staged native torque field, not a torque
            flux law, angular update, or feedback term into ``linear_a``.
        """

        self.native_angular_tau_divergence = div(self.native_angular_tau)
        return self.native_angular_tau_divergence

    def evaluate_native_angular_metric_divergence(self):
        """Evaluate metric-weighted torque-divergence candidates.

        Returns:
            ``(native_angular_tau_metric_divergence_t,
            native_angular_tau_metric_divergence_p)``.

        Notes:
            ``div(native_angular_tau)`` has one extra inverse-length factor
            compared with the staged torque channels. Multiplying by the local
            lever arms ``ell_t`` and ``ell_p`` produces passive candidates that
            can be compared with ``angular_torque_t`` and
            ``angular_torque_p`` in a later residual test. This still does not
            define a transport law, source term, angular update, or feedback
            into the linear branch.
        """

        self.evaluate_native_angular_torque_divergence()
        self.native_angular_tau_metric_divergence_t = (
            self.ell_t * self.native_angular_tau_divergence
        )
        self.native_angular_tau_metric_divergence_p = (
            self.ell_p * self.native_angular_tau_divergence
        )
        return (
            self.native_angular_tau_metric_divergence_t,
            self.native_angular_tau_metric_divergence_p,
        )

    def evaluate_native_angular_transport_residual(
        self,
        source_t=None,
        source_p=None,
    ):
        """Evaluate a passive native angular transport residual candidate.

        Args:
            source_t: Optional explicit toroidal source term with the same
                units and shape as ``angular_torque_t``. Defaults to zero.
            source_p: Optional explicit poloidal source term with the same
                units and shape as ``angular_torque_p``. Defaults to zero.

        Returns:
            ``(native_angular_transport_residual_t,
            native_angular_transport_residual_p)``.

        Notes:
            The residual candidate is
            ``R_L = dL/dt + ell*div(tau_native) - S_L``. The ``dL/dt`` part is
            represented by the staged ``angular_torque_t`` and
            ``angular_torque_p`` fields. This is only a diagnostic residual; it
            does not advance angular momentum, call boundary hooks, define
            boundary exchange, or feed back into the linear branch. Boundary or
            external exchange must be supplied explicitly through this method's
            source arrays until native angular boundary semantics are derived;
            previously collected source buffers are not read implicitly.
        """

        rhs_t, rhs_p = self.evaluate_native_angular_momentum_rhs(
            source_t=source_t,
            source_p=source_p,
        )

        self.native_angular_transport_residual_t = self.angular_torque_t - rhs_t
        self.native_angular_transport_residual_p = self.angular_torque_p - rhs_p
        return (
            self.native_angular_transport_residual_t,
            self.native_angular_transport_residual_p,
        )

    def _validate_native_angular_source_terms(self, source_t, source_p):
        """Validate explicit native angular source/exchange terms."""

        source_t = bd.asarray(source_t)
        source_p = bd.asarray(source_p)
        expected_shape = self.angular_torque_t.shape
        if source_t.shape != expected_shape:
            raise ValueError(
                "native angular toroidal source must have shape "
                f"{expected_shape}, got {source_t.shape}"
            )
        if source_p.shape != expected_shape:
            raise ValueError(
                "native angular poloidal source must have shape "
                f"{expected_shape}, got {source_p.shape}"
            )
        return source_t, source_p

    def collect_native_angular_source_terms(
        self,
        include_sources=True,
        include_boundaries=True,
    ):
        """Collect explicit native angular source/exchange terms.

        Scene elements may opt in by exposing
        ``native_angular_source_terms()``, returning ``(source_t, source_p)``
        arrays with the same shape and units as ``angular_torque_t`` and
        ``angular_torque_p``. This method only sums those explicit terms into
        ``native_angular_source_t`` and ``native_angular_source_p``. It is
        called by ``step()`` only when opt-in native angular transport and
        source collection are enabled, does not itself update angular momentum,
        and does not make boundary exchange implicit in the residual helper.
        The source buffers are reset on every call, so include flags select the
        current accounting view rather than accumulating previous collections.
        Current boundary hooks are source-buffer contracts only: they may read
        native angular momentum, but they must not mutate clocks, moments,
        torque fields, residuals, candidates, or linear-sector state.
        """

        self.native_angular_source_t *= 0.0
        self.native_angular_source_p *= 0.0

        actors = []
        if include_sources:
            actors.extend(self.sources)
        if include_boundaries:
            actors.extend(self.boundaries)

        for actor in actors:
            hook = getattr(actor, "native_angular_source_terms", None)
            if hook is None:
                continue
            terms = hook()
            if terms is None:
                continue
            try:
                source_t, source_p = terms
            except (TypeError, ValueError):
                raise ValueError(
                    "native_angular_source_terms() must return "
                    "(source_t, source_p)"
                )
            source_t, source_p = self._validate_native_angular_source_terms(
                source_t,
                source_p,
            )
            self.native_angular_source_t += source_t
            self.native_angular_source_p += source_p

        return self.native_angular_source_t, self.native_angular_source_p

    def collect_native_angular_boundary_contracts(self):
        """Return metadata contracts for native-angular boundary hooks.

        Boundaries may opt in by exposing
        ``native_angular_boundary_contract()``. The returned dictionaries are
        copied and annotated with ``boundary_type`` and ``boundary_name`` so
        examples and tests can report the current boundary semantics without
        inferring physical absorber or reflector behavior from source hooks.
        """

        contracts = []
        for boundary in self.boundaries:
            hook = getattr(boundary, "native_angular_boundary_contract", None)
            if hook is None:
                continue
            contract = dict(hook())
            contract["boundary_type"] = boundary.__class__.__name__
            contract["boundary_name"] = getattr(boundary, "name", None)
            contracts.append(contract)
        return contracts

    def validate_native_angular_boundary_contracts(
        self,
        require_physical_laws=False,
    ):
        """Validate advertised native-angular boundary semantics.

        The current accepted contract is intentionally narrow:
        source-accounting hooks may return explicit source buffers and may not
        claim to be physical boundary laws. Passing
        ``require_physical_laws=True`` turns this into the acceptance gate for
        future absorbers or reflectors: every boundary contract must explicitly
        claim a physical law and must name the extra geometry, flux, and
        balance fields needed to make that claim reviewable.
        """

        contracts = self.collect_native_angular_boundary_contracts()
        required_keys = {
            "scope",
            "returns",
            "reads",
            "mutates",
            "physical_boundary_law",
        }
        physical_law_keys = {
            "boundary_slots",
            "flux_split",
            "metric_frame",
            "energy_balance",
            "acceptance_test",
        }

        if require_physical_laws and not contracts:
            raise ValueError("native angular physical boundary law required")

        for contract in contracts:
            label = contract.get("boundary_name") or contract.get("boundary_type")
            missing = sorted(required_keys - set(contract))
            if missing:
                raise ValueError(
                    f"native angular boundary contract {label} is missing "
                    f"required keys: {missing}"
                )

            if contract["scope"] == "source_accounting":
                if contract["physical_boundary_law"]:
                    raise ValueError(
                        "source-accounting native angular boundary contracts "
                        "must not claim physical_boundary_law=True"
                    )
                if tuple(contract["returns"]) != (
                    "native_angular_source_t",
                    "native_angular_source_p",
                ):
                    raise ValueError(
                        "source-accounting native angular boundary contracts "
                        "must return native_angular_source_t/p"
                    )
                if tuple(contract["mutates"]) != ():
                    raise ValueError(
                        "source-accounting native angular boundary contracts "
                        "must declare no direct mutation targets"
                    )
            elif contract["physical_boundary_law"]:
                missing = sorted(physical_law_keys - set(contract))
                if missing:
                    raise ValueError(
                        f"physical native angular boundary contract {label} "
                        f"is missing required keys: {missing}"
                    )
            else:
                raise ValueError(
                    "native angular boundary contract must either be "
                    "source_accounting or declare physical_boundary_law=True"
                )

            if require_physical_laws and not contract["physical_boundary_law"]:
                raise ValueError("native angular physical boundary law required")

        return contracts

    def _native_angular_boundary_face(self, axis, side):
        """Return axis metadata and a single outer-face slice."""

        axis_map = {"x": 0, "y": 1, "z": 2, 0: 0, 1: 1, 2: 2}
        if axis not in axis_map:
            raise ValueError("axis must be 'x', 'y', 'z', 0, 1, or 2")
        axis_index = axis_map[axis]
        axis_label = ("x", "y", "z")[axis_index]
        if side not in ("low", "high"):
            raise ValueError("side must be 'low' or 'high'")

        face_index = 0 if side == "low" else -1
        normal_sign = -1.0 if side == "low" else 1.0
        face = [slice(None), slice(None), slice(None), slice(None)]
        face[axis_index] = face_index
        return axis_index, axis_label, normal_sign, tuple(face)

    def evaluate_native_angular_boundary_flux(self, axis, side, field=None):
        """Split projected native-angular boundary flux into incident/outgoing.

        Args:
            axis: Boundary-normal axis, either ``"x"``, ``"y"``, ``"z"``, or
                integer ``0`` / ``1`` / ``2``.
            side: ``"low"`` or ``"high"`` outer face on that axis.
            field: Optional projected angular torque/flux vector. Defaults to
                ``native_angular_tau``.

        Returns:
            A dictionary containing the boundary face, outward-normal flux,
            incident flux, outgoing flux, and scalar totals.

        Notes:
            This is an observability diagnostic for the current collocated
            bridge. It separates signs of the projected ``native_angular_tau``
            normal component, but it is not yet a native t/p flux law and does
            not define absorption or reflection.
        """

        if field is None:
            field = self.native_angular_tau
        field = bd.asarray(field)
        if field.shape != self.native_angular_tau.shape:
            raise ValueError(
                "native angular boundary flux field must have shape "
                f"{self.native_angular_tau.shape}, got {field.shape}"
            )

        axis_index, axis_label, normal_sign, face = self._native_angular_boundary_face(
            axis,
            side,
        )
        component_slice = list(face)
        component_slice[-1] = slice(axis_index, axis_index + 1)
        component_slice = tuple(component_slice)
        normal_flux = normal_sign * field[component_slice]
        magnitude = abs(normal_flux)
        outgoing_flux = 0.5 * (normal_flux + magnitude)
        incident_flux = 0.5 * (magnitude - normal_flux)

        self.native_angular_boundary_normal_flux *= 0.0
        self.native_angular_boundary_incident_flux *= 0.0
        self.native_angular_boundary_outgoing_flux *= 0.0
        self.native_angular_boundary_normal_flux[face] = normal_flux
        self.native_angular_boundary_incident_flux[face] = incident_flux
        self.native_angular_boundary_outgoing_flux[face] = outgoing_flux

        return {
            "axis": axis_label,
            "side": side,
            "normal_sign": normal_sign,
            "face": face,
            "normal_flux": normal_flux,
            "incident_flux": incident_flux,
            "outgoing_flux": outgoing_flux,
            "net_flux": bd.sum(normal_flux),
            "incident_total": bd.sum(incident_flux),
            "outgoing_total": bd.sum(outgoing_flux),
        }

    def collect_native_angular_boundary_fluxes(self, contracts=None, field=None):
        """Evaluate boundary-flux diagnostics for registered boundary faces.

        Args:
            contracts: Optional boundary-contract dictionaries. Defaults to
                ``validate_native_angular_boundary_contracts()``.
            field: Optional projected angular torque/flux vector passed to
                ``evaluate_native_angular_boundary_flux()``.

        Returns:
            A list of per-boundary summary dictionaries containing contract
            identity metadata and scalar incident/outgoing/net flux totals.

        Notes:
            This connects boundary-contract metadata to the projected
            incident/outgoing diagnostic. It does not collect source terms,
            advance momentum, or make the boundary contract physical.
        """

        if contracts is None:
            contracts = self.validate_native_angular_boundary_contracts()

        fluxes = []
        for contract in contracts:
            if "boundary_axis" not in contract or "boundary_side" not in contract:
                continue
            flux = self.evaluate_native_angular_boundary_flux(
                axis=contract["boundary_axis"],
                side=contract["boundary_side"],
                field=field,
            )
            fluxes.append(
                {
                    "boundary_type": contract.get("boundary_type"),
                    "boundary_name": contract.get("boundary_name"),
                    "exchange": contract.get("exchange"),
                    "physical_boundary_law": contract.get(
                        "physical_boundary_law",
                    ),
                    "axis": flux["axis"],
                    "side": flux["side"],
                    "normal_sign": flux["normal_sign"],
                    "incident_total": flux["incident_total"],
                    "outgoing_total": flux["outgoing_total"],
                    "net_flux": flux["net_flux"],
                }
            )
        return fluxes

    def evaluate_native_angular_boundary_channel_flux(self, axis, side, field=None):
        """Split projected boundary flux into native t/p channel diagnostics.

        This first channel diagnostic weights the incident/outgoing projected
        normal flux by the absolute projection of the local native angular
        frame vectors ``angular_e_t`` and ``angular_e_p`` on the boundary
        normal. It makes the matched-flux candidate measurable per channel, but
        it is still a bridge diagnostic rather than a derived native flux law.
        """

        axis_index, _, normal_sign, face = self._native_angular_boundary_face(
            axis,
            side,
        )
        flux = self.evaluate_native_angular_boundary_flux(
            axis=axis,
            side=side,
            field=field,
        )
        component = (slice(None), slice(None), slice(None), axis_index)
        weight_t = abs((normal_sign * self.angular_e_t[component])[:, :, :, None])
        weight_p = abs((normal_sign * self.angular_e_p[component])[:, :, :, None])
        incident_t = flux["incident_flux"] * weight_t[face]
        incident_p = flux["incident_flux"] * weight_p[face]
        outgoing_t = flux["outgoing_flux"] * weight_t[face]
        outgoing_p = flux["outgoing_flux"] * weight_p[face]

        self.native_angular_boundary_incident_flux_t *= 0.0
        self.native_angular_boundary_incident_flux_p *= 0.0
        self.native_angular_boundary_outgoing_flux_t *= 0.0
        self.native_angular_boundary_outgoing_flux_p *= 0.0
        self.native_angular_boundary_incident_flux_t[face] = incident_t
        self.native_angular_boundary_incident_flux_p[face] = incident_p
        self.native_angular_boundary_outgoing_flux_t[face] = outgoing_t
        self.native_angular_boundary_outgoing_flux_p[face] = outgoing_p

        channel_flux = dict(flux)
        channel_flux.update(
            {
                "incident_flux_t": incident_t,
                "incident_flux_p": incident_p,
                "outgoing_flux_t": outgoing_t,
                "outgoing_flux_p": outgoing_p,
                "incident_total_t": bd.sum(incident_t),
                "incident_total_p": bd.sum(incident_p),
                "outgoing_total_t": bd.sum(outgoing_t),
                "outgoing_total_p": bd.sum(outgoing_p),
            }
        )
        return channel_flux

    def collect_native_angular_boundary_channel_fluxes(
        self,
        contracts=None,
        field=None,
    ):
        """Evaluate t/p channel boundary fluxes for registered contracts."""

        if contracts is None:
            contracts = self.validate_native_angular_boundary_contracts()

        fluxes = []
        for contract in contracts:
            if "boundary_axis" not in contract or "boundary_side" not in contract:
                continue
            flux = self.evaluate_native_angular_boundary_channel_flux(
                axis=contract["boundary_axis"],
                side=contract["boundary_side"],
                field=field,
            )
            fluxes.append(
                {
                    "boundary_type": contract.get("boundary_type"),
                    "boundary_name": contract.get("boundary_name"),
                    "exchange": contract.get("exchange"),
                    "physical_boundary_law": contract.get(
                        "physical_boundary_law",
                    ),
                    "axis": flux["axis"],
                    "side": flux["side"],
                    "incident_total": flux["incident_total"],
                    "outgoing_total": flux["outgoing_total"],
                    "net_flux": flux["net_flux"],
                    "incident_total_t": flux["incident_total_t"],
                    "incident_total_p": flux["incident_total_p"],
                    "outgoing_total_t": flux["outgoing_total_t"],
                    "outgoing_total_p": flux["outgoing_total_p"],
                }
            )
        return fluxes

    def evaluate_native_angular_boundary_direct_channel_flux(
        self,
        axis,
        side,
        torque_t=None,
        torque_p=None,
    ):
        """Split native t/p boundary flux before summing projected channels.

        This comparator projects ``angular_torque_t`` and ``angular_torque_p``
        separately through ``angular_e_t/p`` and then sign-splits each channel.
        It is useful for falsifying the current matched-flux candidate, whose
        first diagnostic splits the summed projected normal flux and then
        weights that total by channel visibility.
        """

        if torque_t is None:
            torque_t = self.angular_torque_t
        else:
            torque_t = bd.asarray(torque_t)
        if torque_p is None:
            torque_p = self.angular_torque_p
        else:
            torque_p = bd.asarray(torque_p)
        if torque_t.shape != self.angular_torque_t.shape:
            raise ValueError(
                "native angular direct channel flux torque_t must have shape "
                f"{self.angular_torque_t.shape}, got {torque_t.shape}"
            )
        if torque_p.shape != self.angular_torque_p.shape:
            raise ValueError(
                "native angular direct channel flux torque_p must have shape "
                f"{self.angular_torque_p.shape}, got {torque_p.shape}"
            )

        axis_index, axis_label, normal_sign, face = self._native_angular_boundary_face(
            axis,
            side,
        )
        component = (slice(None), slice(None), slice(None), axis_index)
        frame_t = (normal_sign * self.angular_e_t[component])[:, :, :, None]
        frame_p = (normal_sign * self.angular_e_p[component])[:, :, :, None]
        normal_t = torque_t[face] * frame_t[face]
        normal_p = torque_p[face] * frame_p[face]
        magnitude_t = abs(normal_t)
        magnitude_p = abs(normal_p)
        outgoing_t = 0.5 * (normal_t + magnitude_t)
        outgoing_p = 0.5 * (normal_p + magnitude_p)
        incident_t = 0.5 * (magnitude_t - normal_t)
        incident_p = 0.5 * (magnitude_p - normal_p)

        self.native_angular_boundary_direct_normal_flux_t *= 0.0
        self.native_angular_boundary_direct_normal_flux_p *= 0.0
        self.native_angular_boundary_direct_incident_flux_t *= 0.0
        self.native_angular_boundary_direct_incident_flux_p *= 0.0
        self.native_angular_boundary_direct_outgoing_flux_t *= 0.0
        self.native_angular_boundary_direct_outgoing_flux_p *= 0.0
        self.native_angular_boundary_direct_normal_flux_t[face] = normal_t
        self.native_angular_boundary_direct_normal_flux_p[face] = normal_p
        self.native_angular_boundary_direct_incident_flux_t[face] = incident_t
        self.native_angular_boundary_direct_incident_flux_p[face] = incident_p
        self.native_angular_boundary_direct_outgoing_flux_t[face] = outgoing_t
        self.native_angular_boundary_direct_outgoing_flux_p[face] = outgoing_p

        return {
            "axis": axis_label,
            "side": side,
            "normal_sign": normal_sign,
            "face": face,
            "normal_flux_t": normal_t,
            "normal_flux_p": normal_p,
            "incident_flux_t": incident_t,
            "incident_flux_p": incident_p,
            "outgoing_flux_t": outgoing_t,
            "outgoing_flux_p": outgoing_p,
            "net_flux_t": bd.sum(normal_t),
            "net_flux_p": bd.sum(normal_p),
            "net_flux": bd.sum(normal_t + normal_p),
            "incident_total_t": bd.sum(incident_t),
            "incident_total_p": bd.sum(incident_p),
            "incident_total": bd.sum(incident_t + incident_p),
            "outgoing_total_t": bd.sum(outgoing_t),
            "outgoing_total_p": bd.sum(outgoing_p),
            "outgoing_total": bd.sum(outgoing_t + outgoing_p),
        }

    def collect_native_angular_boundary_direct_channel_fluxes(
        self,
        contracts=None,
        torque_t=None,
        torque_p=None,
    ):
        """Evaluate direct native t/p channel fluxes for registered contracts."""

        if contracts is None:
            contracts = self.validate_native_angular_boundary_contracts()

        fluxes = []
        for contract in contracts:
            if "boundary_axis" not in contract or "boundary_side" not in contract:
                continue
            flux = self.evaluate_native_angular_boundary_direct_channel_flux(
                axis=contract["boundary_axis"],
                side=contract["boundary_side"],
                torque_t=torque_t,
                torque_p=torque_p,
            )
            fluxes.append(
                {
                    "boundary_type": contract.get("boundary_type"),
                    "boundary_name": contract.get("boundary_name"),
                    "exchange": contract.get("exchange"),
                    "physical_boundary_law": contract.get(
                        "physical_boundary_law",
                    ),
                    "axis": flux["axis"],
                    "side": flux["side"],
                    "net_flux_t": flux["net_flux_t"],
                    "net_flux_p": flux["net_flux_p"],
                    "net_flux": flux["net_flux"],
                    "incident_total_t": flux["incident_total_t"],
                    "incident_total_p": flux["incident_total_p"],
                    "incident_total": flux["incident_total"],
                    "outgoing_total_t": flux["outgoing_total_t"],
                    "outgoing_total_p": flux["outgoing_total_p"],
                    "outgoing_total": flux["outgoing_total"],
                }
            )
        return fluxes

    def evaluate_native_angular_momentum_rhs(
        self,
        source_t=None,
        source_p=None,
    ):
        """Evaluate the passive native angular momentum transport RHS.

        The staged transport equation is
        ``dL/dt = S_L - ell*div(tau_native)``. This method stores that right
        hand side in ``native_angular_momentum_rhs_t`` and
        ``native_angular_momentum_rhs_p``. It does not advance
        ``angular_momentum_t`` or ``angular_momentum_p``. Previously collected
        source buffers are not read unless passed through ``source_t`` and
        ``source_p``.
        """

        metric_t, metric_p = self.evaluate_native_angular_metric_divergence()
        if source_t is None:
            source_t = bd.zeros_like(self.angular_torque_t)
        else:
            source_t = bd.asarray(source_t)
        if source_p is None:
            source_p = bd.zeros_like(self.angular_torque_p)
        else:
            source_p = bd.asarray(source_p)
        source_t, source_p = self._validate_native_angular_source_terms(
            source_t,
            source_p,
        )

        self.native_angular_momentum_rhs_t = source_t - metric_t
        self.native_angular_momentum_rhs_p = source_p - metric_p
        return (
            self.native_angular_momentum_rhs_t,
            self.native_angular_momentum_rhs_p,
        )

    def evaluate_native_angular_transport_power(
        self,
        source_t=None,
        source_p=None,
        momentum_t=None,
        momentum_p=None,
        inertia_t=None,
        inertia_p=None,
    ):
        """Evaluate native angular transport-power diagnostics.

        Args:
            source_t: Optional toroidal source/exchange term for the transport
                RHS. Defaults to zero.
            source_p: Optional poloidal source/exchange term for the transport
                RHS. Defaults to zero.
            momentum_t: Optional toroidal momentum field. Defaults to
                ``angular_momentum_t``.
            momentum_p: Optional poloidal momentum field. Defaults to
                ``angular_momentum_p``.
            inertia_t: Optional toroidal inertia field. Defaults to
                ``angular_inertia_t``.
            inertia_p: Optional poloidal inertia field. Defaults to
                ``angular_inertia_p``.

        Returns:
            ``(power_t, power_p, total_power)`` using
            ``P_rhs = (S_L - ell*div(tau_native))*L/I``.

        Notes:
            This is the full-RHS companion to
            ``evaluate_native_angular_exchange_power()``. It is useful for
            finite-step energy accounting of the staged transport law, but it
            does not define a Hamiltonian or make the native angular transport
            law physical.
        """

        rhs_t, rhs_p = self.evaluate_native_angular_momentum_rhs(
            source_t=source_t,
            source_p=source_p,
        )
        momentum_t = (
            self.angular_momentum_t if momentum_t is None else bd.asarray(momentum_t)
        )
        momentum_p = (
            self.angular_momentum_p if momentum_p is None else bd.asarray(momentum_p)
        )
        inertia_t = self.angular_inertia_t if inertia_t is None else bd.asarray(inertia_t)
        inertia_p = self.angular_inertia_p if inertia_p is None else bd.asarray(inertia_p)

        if bd.max(inertia_t <= 0) or bd.max(inertia_p <= 0):
            raise ValueError("native angular inertia must be positive")

        self.native_angular_transport_power_t = rhs_t * momentum_t / inertia_t
        self.native_angular_transport_power_p = rhs_p * momentum_p / inertia_p
        total_power = bd.sum(
            self.native_angular_transport_power_t
            + self.native_angular_transport_power_p
        )
        return (
            self.native_angular_transport_power_t,
            self.native_angular_transport_power_p,
            total_power,
        )

    def predict_native_angular_momentum_step(
        self,
        delta=None,
        source_t=None,
        source_p=None,
    ):
        """Predict one passive native angular momentum transport step.

        Returns candidate next momentum fields computed as
        ``L_next = L + delta*(S_L - ell*div(tau_native))``. The candidate is
        stored separately and is not promoted into ``angular_momentum_t`` or
        ``angular_momentum_p``.
        """

        if delta is None:
            delta = self.time_step

        rhs_t, rhs_p = self.evaluate_native_angular_momentum_rhs(
            source_t=source_t,
            source_p=source_p,
        )
        self.native_angular_momentum_candidate_t = (
            self.angular_momentum_t + delta * rhs_t
        )
        self.native_angular_momentum_candidate_p = (
            self.angular_momentum_p + delta * rhs_p
        )
        return (
            self.native_angular_momentum_candidate_t,
            self.native_angular_momentum_candidate_p,
        )

    def advance_native_angular_momentum_transport(
        self,
        delta=None,
        source_t=None,
        source_p=None,
        collect_sources=False,
        include_sources=True,
        include_boundaries=True,
        update_clocks=False,
    ):
        """Advance native angular momentum through the passive transport law.

        Args:
            delta: Time interval for the update. Defaults to the grid timestep.
            source_t: Optional explicit toroidal source term.
            source_p: Optional explicit poloidal source term.
            collect_sources: If true, add current opt-in scene source and
                boundary hooks from ``collect_native_angular_source_terms()``.
            include_sources: Passed to ``collect_native_angular_source_terms``
                when ``collect_sources`` is true.
            include_boundaries: Passed to ``collect_native_angular_source_terms``
                when ``collect_sources`` is true.
            update_clocks: If true, update ``omega_t`` and ``omega_p`` from the
                promoted momentum and positive native inertia.

        Returns:
            ``(angular_momentum_t, angular_momentum_p)`` after promotion.

        Notes:
            This is the first explicit opt-in transport advance for the native
            angular staging layer:

            ``L_next = L + delta*(S_L - ell*div(tau_native))``.

            It runs from ``step()`` only when ``native_angular_transport`` is
            enabled, does not update the legacy bridge fields, and does not
            define physical boundary semantics.
        """

        if update_clocks:
            if bd.max(self.angular_inertia_t <= 0) or bd.max(
                self.angular_inertia_p <= 0
            ):
                raise ValueError(
                    "native angular inertia must be positive to update clocks"
                )

        if collect_sources:
            collected_t, collected_p = self.collect_native_angular_source_terms(
                include_sources=include_sources,
                include_boundaries=include_boundaries,
            )
            if source_t is None:
                source_t = collected_t
            else:
                source_t = bd.asarray(source_t) + collected_t
            if source_p is None:
                source_p = collected_p
            else:
                source_p = bd.asarray(source_p) + collected_p

        candidate_t, candidate_p = self.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        self.angular_momentum_t = candidate_t
        self.angular_momentum_p = candidate_p

        if update_clocks:
            self.omega_t = self.angular_momentum_t / self.angular_inertia_t
            self.omega_p = self.angular_momentum_p / self.angular_inertia_p

        return self.angular_momentum_t, self.angular_momentum_p

    def evaluate_native_angular_linear_response(
        self,
        metric_length=None,
        scale=1.0,
    ):
        """Project native angular stress into a passive linear response.

        Args:
            metric_length: Optional angular metric length for the bridge.
                Defaults to no additional weighting.
            scale: Optional scalar multiplier for candidate comparisons.

        Returns:
            ``native_angular_linear_response`` with shape ``(Nx, Ny, Nz, 3)``.

        Notes:
            This is a coupling diagnostic only. It projects
            ``native_angular_tau`` through the current ``angular_to_linear``
            bridge but does not write ``angular_A``, ``linear_a``, or any
            production update state.
        """

        self.native_angular_linear_response = scale * angular_to_linear_bridge(
            self.native_angular_tau,
            metric_length=metric_length,
        )
        return self.native_angular_linear_response

    def apply_native_angular_linear_response(
        self,
        response=None,
        metric_length=None,
        scale=1.0,
        mode="add",
    ):
        """Apply a native-angular linear response to ``linear_a`` explicitly.

        Args:
            response: Optional precomputed response field. If omitted, the
                current ``native_angular_tau`` is projected with
                ``evaluate_native_angular_linear_response()``.
            metric_length: Optional metric length used only when ``response``
                is omitted.
            scale: Scalar multiplier used only when ``response`` is omitted.
            mode: ``"add"`` adds the response to ``linear_a``;
                ``"replace"`` replaces ``linear_a`` with the response.

        Returns:
            The updated ``linear_a`` field.

        Notes:
            This is an opt-in coupling experiment. It does not run from
            ``step()``, does not write the legacy ``angular_A`` bridge field,
            and does not define the physical feedback law.
        """

        if mode not in ("add", "replace"):
            raise ValueError("mode must be 'add' or 'replace'")

        if response is None:
            response = self.evaluate_native_angular_linear_response(
                metric_length=metric_length,
                scale=scale,
            )
        else:
            response = bd.asarray(response)
            expected_shape = self.linear_a.shape
            if response.shape != expected_shape:
                raise ValueError(
                    "native angular linear response must have shape "
                    f"{expected_shape}, got {response.shape}"
            )
            self.native_angular_linear_response = response

        if mode == "add":
            self.linear_a = self.linear_a + response
        else:
            self.linear_a = response

        self._sync_public_aliases()
        return self.linear_a

    def evaluate_linear_charge_flux_candidate(
        self,
        delta_t=None,
        area_measure=1.0,
        normalization_area=1.0,
        sign=1.0,
        field=None,
    ):
        """Evaluate the local linear charge-flux candidate.

        Args:
            delta_t: Time scale used in ``Q_q = delta_t*eta*v``. Defaults to
                the grid timestep.
            area_measure: Optional dimensionless ``dA``-like weight or array.
                Use together with ``normalization_area`` to represent
                ``dA/A0``.
            normalization_area: Reference area ``A0``. It must be chosen by
                geometry, not by fitting measured charge.
            sign: Explicit linear-channel polarity convention.
            field: Optional velocity-like vector field. Defaults to
                ``linear_v``.

        Returns:
            ``linear_charge_flux_candidate`` with shape ``(Nx, Ny, Nz, 3)``.

        Notes:
            This is only the simulator-facing candidate for the integrand of
            ``q_l``. It does not integrate over a surface, choose ``A0``, or
            derive a measured scalar charge.
        """

        if delta_t is None:
            delta_t = self.time_step
        if normalization_area == 0:
            raise ValueError("normalization_area must be non-zero")

        if field is None:
            field = self.linear_v
        else:
            field = bd.asarray(field)
        area_measure = bd.asarray(area_measure)

        self.linear_charge_flux_candidate = (
            sign * delta_t * eta * field * area_measure / normalization_area
        )
        return self.linear_charge_flux_candidate

    def evaluate_native_angular_charge_flux_candidates(
        self,
        delta_t=None,
        loop_measure_t=1.0,
        loop_measure_p=1.0,
        angular_normalization=1.0,
        sign_t=1.0,
        sign_p=1.0,
    ):
        """Evaluate native two-clock angular charge-flux candidates.

        Args:
            delta_t: Time scale used in ``Q_q = delta_t*eta*ell*omega``.
                Defaults to the grid timestep.
            loop_measure_t: Dimensionless toroidal loop measure, such as
                ``dtheta_t/(2*pi)``.
            loop_measure_p: Dimensionless poloidal loop measure, such as
                ``dtheta_p/(2*pi)``.
            angular_normalization: Candidate angular normalization ``N_a``.
                It must be chosen independently of measured charge.
            sign_t: Explicit toroidal polarity convention.
            sign_p: Explicit poloidal polarity convention.

        Returns:
            ``(native_angular_charge_flux_t, native_angular_charge_flux_p)``.

        Notes:
            This exposes the two native angular channels only. It does not
            choose the final scalar angular reduction and is not coupled into
            ``step()``.
        """

        if delta_t is None:
            delta_t = self.time_step
        if angular_normalization == 0:
            raise ValueError("angular_normalization must be non-zero")

        loop_measure_t = bd.asarray(loop_measure_t)
        loop_measure_p = bd.asarray(loop_measure_p)
        self.native_angular_charge_flux_t = (
            sign_t
            * delta_t
            * eta
            * self.ell_t
            * self.omega_t
            * loop_measure_t
            / angular_normalization
        )
        self.native_angular_charge_flux_p = (
            sign_p
            * delta_t
            * eta
            * self.ell_p
            * self.omega_p
            * loop_measure_p
            / angular_normalization
        )
        return (
            self.native_angular_charge_flux_t,
            self.native_angular_charge_flux_p,
        )

    def reduce_native_angular_charge_flux(self, mode="additive"):
        """Reduce the staged angular charge candidates to one scalar field.

        Args:
            mode: ``"additive"`` for ``q_t+q_p`` or ``"geometric_mean"`` for
                ``sqrt(q_t*q_p)``.

        Returns:
            ``native_angular_charge_reduction``.

        Notes:
            This is an explicit comparison helper for candidate observables.
            It does not decide which reduction is physical. The geometric-mean
            candidate requires a non-negative channel product after the caller
            has chosen signs.
        """

        if mode == "additive":
            self.native_angular_charge_reduction = (
                self.native_angular_charge_flux_t + self.native_angular_charge_flux_p
            )
        elif mode == "geometric_mean":
            product = (
                self.native_angular_charge_flux_t * self.native_angular_charge_flux_p
            )
            if bd.max(product < 0):
                raise ValueError(
                    "geometric_mean reduction requires non-negative q_t*q_p"
                )
            self.native_angular_charge_reduction = product ** 0.5
        else:
            raise ValueError(
                "mode must be 'additive' or 'geometric_mean'"
            )
        return self.native_angular_charge_reduction

    def advance_linear_sector(self):
        """Advance the primary linear state with a short Taylor step."""

        self.linear_v += self.time_step * self.linear_a
        self.linear_v += 0.5 * (self.time_step ** 2) * self.linear_j


    def updateBoundaries(self):
        """Run any pre-update boundary hooks supported by the scene."""

        for boundary in self.boundaries:
            update_phi_E = getattr(boundary, "update_phi_E", None)
            if update_phi_E is not None:
                update_phi_E()
            update_phi_H = getattr(boundary, "update_phi_H", None)
            if update_phi_H is not None:
                update_phi_H()

        return


    def updateEH(self):    
        """Run package-style hooks around the derived ``E`` / ``H`` fields.

        Existing Maxwell-style sources and detectors can therefore still be
        reused as comparison tools, while aether-native sources may inject into
        ``v`` earlier in the update chain via ``update_v``.
        """
        for boundary in self.boundaries:
            update_E = getattr(boundary, "update_E", None)
            if update_E is not None:
                update_E()
            update_H = getattr(boundary, "update_H", None)
            if update_H is not None:
                update_H()
           
        for src in self.sources:
            update_E = getattr(src, "update_E", None)
            if update_E is not None:
                update_E()
            update_H = getattr(src, "update_H", None)
            if update_H is not None:
                update_H()
           
           
        for det in self.detectors:
            detect_E = getattr(det, "detect_E", None)
            if detect_E is not None:
                detect_E()
            detect_H = getattr(det, "detect_H", None)
            if detect_H is not None:
                detect_H()

        return


    def detect_native_angular(self):
        """Run optional detector hooks for native-angular observables."""

        for det in self.detectors:
            detect = getattr(det, "detect_native_angular", None)
            if detect is not None:
                detect()

        return



    def reset(self):
        """reset the grid by setting all fields to zero"""
        self.linear_v *= 0.0
        self.linear_p *= 0.0
        self.linear_f *= 0.0
        self.linear_E *= 0.0
        self.linear_a *= 0.0
        self.linear_dpdt *= 0.0
        self.linear_yank *= 0.0
        self.linear_dEdt *= 0.0
        self.linear_j *= 0.0

        self.angular_omega *= 0.0
        self.angular_tau *= 0.0
        self.angular_H *= 0.0
        self.angular_A *= 0.0
        self.angular_alpha *= 0.0
        self.angular_dtau_dt *= 0.0
        self.angular_dHdt *= 0.0
        self.angular_dAdt *= 0.0
        self.omega_t *= 0.0
        self.omega_p *= 0.0
        self.angular_gamma *= 0.0
        self.theta_t *= 0.0
        self.theta_p *= 0.0
        self.angular_chi *= 0.0
        self.angular_clock_lambda *= 0.0
        self.angular_momentum_t *= 0.0
        self.angular_momentum_p *= 0.0
        self.angular_torque_t *= 0.0
        self.angular_torque_p *= 0.0
        self.native_angular_tau *= 0.0
        self.native_angular_tau_divergence *= 0.0
        self.native_angular_tau_metric_divergence_t *= 0.0
        self.native_angular_tau_metric_divergence_p *= 0.0
        self.native_angular_transport_residual_t *= 0.0
        self.native_angular_transport_residual_p *= 0.0
        self.native_angular_source_t *= 0.0
        self.native_angular_source_p *= 0.0
        self.native_angular_exchange_power_t *= 0.0
        self.native_angular_exchange_power_p *= 0.0
        self.native_angular_momentum_rhs_t *= 0.0
        self.native_angular_momentum_rhs_p *= 0.0
        self.native_angular_transport_power_t *= 0.0
        self.native_angular_transport_power_p *= 0.0
        self.native_angular_momentum_candidate_t *= 0.0
        self.native_angular_momentum_candidate_p *= 0.0
        self.native_angular_linear_response *= 0.0
        self.linear_charge_flux_candidate *= 0.0
        self.native_angular_charge_flux_t *= 0.0
        self.native_angular_charge_flux_p *= 0.0
        self.native_angular_charge_reduction *= 0.0
        self.native_angular_boundary_normal_flux *= 0.0
        self.native_angular_boundary_incident_flux *= 0.0
        self.native_angular_boundary_outgoing_flux *= 0.0
        self.native_angular_boundary_incident_flux_t *= 0.0
        self.native_angular_boundary_incident_flux_p *= 0.0
        self.native_angular_boundary_outgoing_flux_t *= 0.0
        self.native_angular_boundary_outgoing_flux_p *= 0.0
        self.native_angular_boundary_direct_normal_flux_t *= 0.0
        self.native_angular_boundary_direct_normal_flux_p *= 0.0
        self.native_angular_boundary_direct_incident_flux_t *= 0.0
        self.native_angular_boundary_direct_incident_flux_p *= 0.0
        self.native_angular_boundary_direct_outgoing_flux_t *= 0.0
        self.native_angular_boundary_direct_outgoing_flux_p *= 0.0

        self._sync_public_aliases()
        self.time_steps_passed = 0

    def add_source(self, name, source):
        """add a source to the grid"""
        source._register_grid(self)
        self.sources[name] = source

    def add_boundary(self, name, boundary):
        """add a boundary to the grid"""
        boundary._register_grid(self)
        self.boundaries[name] = boundary

    def add_detector(self, name, detector):
        """add a detector to the grid"""
        detector._register_grid(self)
        self.detectors[name] = detector

    def add_object(self, name, obj):
        """add an object to the grid"""
        obj._register_grid(self)
        self.objects[name] = obj
    
    def promote_dtypes_to_complex(self):
        self.linear_E = self.linear_E.astype(bd.complex)
        self.angular_H = self.angular_H.astype(bd.complex)
        self._sync_public_aliases()
        [boundary.promote_dtypes_to_complex() for boundary in self.boundaries]

    def __setitem__(self, key, attr):
        if not isinstance(key, tuple):
            x, y, z = key, slice(None), slice(None)
        elif len(key) == 1:
            x, y, z = key[0], slice(None), slice(None)
        elif len(key) == 2:
            x, y, z = key[0], key[1], slice(None)
        elif len(key) == 3:
            x, y, z = key
        else:
            raise KeyError("maximum number of indices for the grid is 3")

        attr._register_grid(
            grid=self,
            x=self._handle_single_key(x),
            y=self._handle_single_key(y),
            z=self._handle_single_key(z),
        )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(shape=({self.Nx},{self.Ny},{self.Nz}), "
            f"grid_spacing={self.grid_spacing:.2e}, courant_number={self.courant_number:.2f})"
        )

    def __str__(self):
        """string representation of the grid

        lists all the components and their locations in the grid.
        """
        s = repr(self) + "\n"
        if self.sources:
            s = s + "\nsources:\n"
            for src in self.sources:
                s += str(src)
        if self.detectors:
            s = s + "\ndetectors:\n"
            for det in self.detectors:
                s += str(det)
        if self.boundaries:
            s = s + "\nboundaries:\n"
            for bnd in self.boundaries:
                s += str(bnd)
        if self.objects:
            s = s + "\nobjects:\n"
            for obj in self.objects:
                s += str(obj)
        return s

    def save_simulation(self, sim_name=None):
        """
        Creates a folder and initializes environment to store simulation or related details.
        saveSimulation() needs to be run before running any function that stores data (generate_video(), save_data()).

        Parameters:-
            (optional) sim_name (string): Preferred name for simulation
        """
        makedirs("fdtd_output", exist_ok=True)  # Output master folder declaration
        # making full_sim_name with timestamp
        full_sim_name = (
            str(datetime.now().year)
            + "-"
            + str(datetime.now().month)
            + "-"
            + str(datetime.now().day)
            + "-"
            + str(datetime.now().hour)
            + "-"
            + str(datetime.now().minute)
            + "-"
            + str(datetime.now().second)
        )
        # Simulation name (optional)
        if sim_name is not None:
            full_sim_name = full_sim_name + " (" + sim_name + ")"
        folder = "fdtd_output_" + full_sim_name
        # storing folder path for saving simulation
        self.folder = os.path.abspath(path.join("fdtd_output", folder))
        # storing timestamp title for self.generate_video
        self.full_sim_name = full_sim_name
        makedirs(self.folder, exist_ok=True)
        return self.folder

    def generate_video(self, delete_frames=False):
        """Compiles frames into a video

        These framed should be saved through ``fdtd.Grid.visualize(save=True)`` while having ``fdtd.Grid.save_simulation()`` enabled.

        Args:
            delete_frames (optional, bool): delete stored frames after conversion to video.

        Returns:
            the filename of the generated video.

        Note:
            this function requires ``ffmpeg`` to be available in your path.
        """
        if self.folder is None:
            raise Exception(
                "Save location not initialized. Please read about 'fdtd.Grid.saveSimulation()' or try running 'grid.saveSimulation()'."
            )
        cwd = path.abspath(os.getcwd())
        chdir(self.folder)
        try:
            check_call(
                [
                    "ffmpeg",
                    "-y",
                    "-framerate",
                    "8",
                    "-i",
                    "file%04d.png",
                    "-r",
                    "30",
                    "-pix_fmt",
                    "yuv420p",
                    "fdtd_sim_video_" + self.full_sim_name + ".mp4",
                ]
            )
        except (FileNotFoundError, CalledProcessError):
            raise CalledProcessError(
                "Error when calling ffmpeg. Is ffmpeg installed and available in your path?"
            )
        if delete_frames:  # delete frames
            for file_name in glob("*.png"):
                remove(file_name)
        video_path = path.abspath(
            path.join(self.folder, f"fdtd_sim_video_{self.full_sim_name}.mp4")
        )
        chdir(cwd)
        return video_path

    def save_data(self):
        """
        Saves readings from all detectors in the grid into a numpy zip file. Each detector is stored in separate arrays. Electric and magnetic field field readings of each detector are also stored separately with suffix " (E)" and " (H)" (Example: ['detector0 (E)', 'detector0 (H)']). Therefore, the numpy zip file contains arrays twice the number of detectors.
        REQUIRES 'fdtd.Grid.save_simulation()' to be run before this function.

        Parameters: None
        """
        def _numpyfy(item):
            if isinstance(item, list):
                return [_numpyfy(el) for el in item]
            elif bd.is_array(item):
                return bd.numpy(item)
            else:
                return item
                
        if self.folder is None:
            raise Exception(
                "Save location not initialized. Please read about 'fdtd.Grid.saveSimulation()' or try running 'grid.saveSimulation()'."
            )
        dic = {}
        for detector in self.detectors:
            values = detector.detector_values()
            for field_name, field_values in values.items():
                dic[detector.name + f" ({field_name})"] = _numpyfy(field_values)
        savez(path.join(self.folder, "detector_readings"), **dic)
