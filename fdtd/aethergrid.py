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
        """
        # save the grid spacing
        self.grid_spacing = float(grid_spacing)

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
        6. repeat the same pattern one derivative level higher to obtain jerk;
        7. advance the linear state with the resulting Taylor step.
        """

        self.updateBoundaries()

        self.apply_native_sources()
        self.update_linear_sector()
        self.update_angular_sector()
        self._sync_public_aliases()

        # Existing package hooks still operate on E/H-like observables.
        self.updateEH()
        self.update_linear_angular_coupling()
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
        ``omega_p``, and ``angular_gamma`` and evaluate the benchmark without
        altering the timestep update or the classic Maxwell path.
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
            dic[detector.name + " (E)"] = _numpyfy(values['E'])
            dic[detector.name + " (H)"] = _numpyfy(values['H'])
        savez(path.join(self.folder, "detector_readings"), **dic)
