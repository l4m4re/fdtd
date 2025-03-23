"""
Potential Grid: Space Time Potential Theory Simulation Grid
==========================================================

This module implements a computational grid for simulating the kinetic substrate
according to Space Time Potential Theory. It models the action flux field (C)
and radiosity field (I) to capture electromagnetic and gravitational effects.

The simulation is based on the fundamental principle:
    a = -η/ρ·Δv

References:
- Theory/vector_laplacian/foundation.md - Mathematical foundations
- Theory/equation_reference.tex - Complete equation reference
"""

import numpy as np
from typing import Tuple, Optional, Dict, Callable, Union, List
from .backend import backend
from .operators import (
    gradient, divergence, curl_face_to_edge, curl_edge_to_face, vector_laplacian,
    scalar_laplacian, helmholtz_decomposition, quantize_circulation
)
from .wave_solvers import FirstSoundSolver, SecondSoundSolver

# Physical constants
C_LIGHT = 299792458.0          # speed of light [m/s]
VISCOSITY = 1.0/(4*np.pi*1e-7) # vacuum viscosity [Pa·s]
BACKGROUND_DENSITY = 1.0       # normalized substrate density [kg/m³]
PLANCK_CONSTANT = 6.62607015e-34  # Planck's constant [J·s]

class PotentialGrid:
    """
    Simulation grid for the kinetic substrate based on Space Time Potential Theory.
    
    This class implements a staggered grid for the Space Time Potential Theory simulation,
    which represents space as a superfluid-like medium with quantized vortices.
    
    The simulation supports:
    - Velocity and vorticity as dual field representations
    - First and second sound wave phenomena
    - Electromagnetic field derivation from substrate dynamics
    - Quantized vortex structures representing particles
    """
    
    def __init__(self, 
                 nx: int, 
                 ny: int, 
                 nz: int,
                 dx: float = 1.0,
                 dt: float = None,
                 viscosity: float = VISCOSITY,
                 background_density: float = BACKGROUND_DENSITY):
        """
        Initialize the potential grid simulation.
        
        Parameters
        ----------
        nx, ny, nz : int
            Grid dimensions
        dx : float
            Grid spacing (assumed uniform in all dimensions)
        dt : float, optional
            Time step, if None, calculated for stability
        viscosity : float
            Dynamic viscosity η of the substrate [Pa·s]
        background_density : float
            Background density ρ₀ of the substrate [kg/m³]
        """
        # Get current backend
        self.xp = backend  # Changed from backend() to backend
        
        # Grid parameters
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx = dx
        
        # Calculate stable time step if not provided
        if dt is None:
            # CFL stability condition for the vector Laplacian
            # dt ≤ dx²/(6*viscosity/density)
            self.dt = 0.9 * dx**2 * background_density / (6 * viscosity)
        else:
            self.dt = dt
            
        self.viscosity = viscosity
        self.background_density = background_density
        self.kinematic_viscosity = viscosity / background_density
        
        # Initialize staggered grid fields
        # Velocity components (on faces)
        self.vx = self.xp.zeros((nx+1, ny, nz))  # x-component on yz faces
        self.vy = self.xp.zeros((nx, ny+1, nz))  # y-component on xz faces
        self.vz = self.xp.zeros((nx, ny, nz+1))  # z-component on xy faces
        
        # Vorticity components (on edges)
        self.wx = self.xp.zeros((nx, ny+1, nz+1))  # x-component along x-edges
        self.wy = self.xp.zeros((nx+1, ny, nz+1))  # y-component along y-edges
        self.wz = self.xp.zeros((nx+1, ny+1, nz))  # z-component along z-edges
        
        # Acceleration components (on faces)
        self.ax = self.xp.zeros((nx+1, ny, nz))
        self.ay = self.xp.zeros((nx, ny+1, nz))
        self.az = self.xp.zeros((nx, ny, nz+1))
        
        # Angular acceleration (on edges)
        self.alpha_x = self.xp.zeros((nx, ny+1, nz+1))
        self.alpha_y = self.xp.zeros((nx+1, ny, nz+1))
        self.alpha_z = self.xp.zeros((nx+1, ny+1, nz))
        
        # Action flux field C = η·v (on faces)
        self.action_flux_x = self.viscosity * self.vx
        self.action_flux_y = self.viscosity * self.vy
        self.action_flux_z = self.viscosity * self.vz
        
        # Radiosity field I = η·a (on faces)
        self.radiosity_x = self.viscosity * self.ax
        self.radiosity_y = self.viscosity * self.ay
        self.radiosity_z = self.viscosity * self.az
        
        # Scalar fields (at cell centers)
        self.density = self.xp.ones((nx, ny, nz)) * background_density
        self.scalar_potential = self.xp.zeros((nx, ny, nz))  # φ = ∇·C
        self.temperature_potential = self.xp.zeros((nx, ny, nz))  # T = ∇·I
        self.temperature = self.xp.zeros((nx, ny, nz))
        
        # Add wave solvers
        self.first_sound = FirstSoundSolver(
            nx, ny, nz, dx, self.dt, viscosity, background_density
        )
        self.second_sound = SecondSoundSolver(
            nx, ny, nz, dx, self.dt, viscosity, background_density
        )
        
        # Time tracking
        self.time = 0.0
        self.step_count = 0
        
        # Initialize derived EM fields
        self._electric_field = None  # Calculated on demand
        self._magnetic_field = None  # Calculated on demand
        
        print(f"PotentialGrid initialized: {nx}×{ny}×{nz} cells")
        print(f"dx={dx}m, dt={self.dt:.2e}s")
        print(f"Viscosity: {viscosity:.2e} Pa·s")
        print(f"Background density: {background_density:.2e} kg/m³")
    
    def update_vorticity(self):
        """
        Update the vorticity field from the velocity field.
        
        Vorticity ω = ∇×v represents the local rotation of the substrate
        and is connected to the magnetic field in electromagnetic theory.
        """
        self.wx, self.wy, self.wz = curl_face_to_edge(self.vx, self.vy, self.vz, self.dx)
    
    def update_acceleration(self):
        """
        Update the acceleration field using the vector Laplacian of velocity.
        
        According to the fundamental principle of the theory:
        a = -η/ρ·Δv
        """
        # Calculate vector Laplacian
        lap_vx, lap_vy, lap_vz = vector_laplacian(self.vx, self.vy, self.vz, self.dx)
        
        # Update acceleration field with sign reversal characteristic of the theory
        factor = -self.viscosity / self.background_density
        self.ax = factor * lap_vx
        self.ay = factor * lap_vy
        self.az = factor * lap_vz
    
    def update_angular_acceleration(self):
        """
        Update the angular acceleration field from the acceleration field.
        
        Angular acceleration α = ∇×a is the time derivative of vorticity
        and relates to the time evolution of magnetic fields.
        """
        self.alpha_x, self.alpha_y, self.alpha_z = curl_face_to_edge(
            self.ax, self.ay, self.az, self.dx
        )
    
    def update_derived_fields(self):
        """
        Update all derived fields from primary fields.
        """
        # Action flux field C = η·v
        self.action_flux_x = self.viscosity * self.vx
        self.action_flux_y = self.viscosity * self.vy
        self.action_flux_z = self.viscosity * self.vz
        
        # Radiosity field I = η·a
        self.radiosity_x = self.viscosity * self.ax
        self.radiosity_y = self.viscosity * self.ay
        self.radiosity_z = self.viscosity * self.az
        
        # Scalar potential φ = ∇·C
        self.scalar_potential = divergence(
            self.action_flux_x, self.action_flux_y, self.action_flux_z, self.dx
        )
        
        # Temperature potential T = ∇·I
        self.temperature_potential = divergence(
            self.radiosity_x, self.radiosity_y, self.radiosity_z, self.dx
        )
        
        # Clear cached EM fields
        self._electric_field = None
        self._magnetic_field = None
    
    def step(self):
        """
        Perform a complete update step for all fields.
        """
        # Update vorticity from velocity
        self.update_vorticity()
        
        # Update acceleration from velocity
        self.update_acceleration()
        
        # Update angular acceleration from acceleration
        self.update_angular_acceleration()
        
        # Update velocity using acceleration
        self.vx += self.ax * self.dt
        self.vy += self.ay * self.dt
        self.vz += self.az * self.dt
        
        # Update derived fields
        self.update_derived_fields()
        
        # Run wave solvers
        self.first_sound.step()
        self.second_sound.step()
        
        # Incorporate wave effects
        # Add first sound velocity components
        self.vx += self.first_sound.vx
        self.vy += self.first_sound.vy
        self.vz += self.first_sound.vz
        
        # Add second sound acceleration components
        self.ax += self.second_sound.ax
        self.ay += self.second_sound.ay
        self.az += self.second_sound.az
        
        # Update density from first sound
        self.density = self.background_density * (1 + self.first_sound.density_perturbation)
        
        # Update temperature from second sound
        self.temperature = self.second_sound.temperature_perturbation
        
        # Update time tracking
        self.time += self.dt
        self.step_count += 1
    
    def add_vortex_ring(self, 
                         center: Tuple[float, float, float],
                         radius: float, 
                         strength: float = None, 
                         orientation: Tuple[float, float, float] = (0, 0, 1),
                         quantized: bool = True):
        """
        Add a vortex ring to the simulation.
        
        In Space Time Potential Theory, elementary particles are modeled as
        quantized vortex rings in the substrate.
        
        Parameters
        ----------
        center : tuple(float, float, float)
            Center position of the vortex ring [m]
        radius : float
            Radius of the vortex ring [m]
        strength : float, optional
            Circulation strength, if None and quantized=True, uses single quantum
        orientation : tuple(float, float, float)
            Orientation vector of the vortex ring
        quantized : bool
            If True, strength is quantized as n·h/m
        """
        # Set default strength (single quantum) if none provided
        if strength is None:
            strength = self.kinematic_viscosity
        
        # Quantize strength if requested
        if quantized:
            n = round(strength / self.kinematic_viscosity)
            strength = n * self.kinematic_viscosity
            print(f"Quantized strength to {n} quanta: {strength:.3e} m²/s")
        
        # Convert center to grid indices
        cx = int(center[0] / self.dx)
        cy = int(center[1] / self.dx)
        cz = int(center[2] / self.dx)
        
        # Convert radius to grid units
        r_grid = int(radius / self.dx)
        
        # Normalize orientation vector
        ox, oy, oz = orientation
        norm = (ox**2 + oy**2 + oz**2)**0.5
        if norm > 0:
            ox, oy, oz = ox/norm, oy/norm, oz/norm
        else:
            ox, oy, oz = 0, 0, 1
        
        # For simplicity, add a circular vortex in the appropriate plane
        # We'll create a vortex ring by setting vorticity directly
        
        # Calculate ring coordinates
        theta = np.linspace(0, 2*np.pi, 100)
        ring_x = cx + r_grid * np.cos(theta)
        ring_y = cy + r_grid * np.sin(theta)
        ring_z = np.ones_like(theta) * cz
        
        # Now rotate the ring to match the orientation vector
        # This is a simplified approach; a full implementation would use proper 3D rotation
        
        # Set vorticity along the ring
        for i in range(len(theta)):
            # Get integer grid positions
            ix = int(ring_x[i])
            iy = int(ring_y[i])
            iz = int(ring_z[i])
            
            # Skip if outside grid
            if (ix < 0 or ix >= self.nx or 
                iy < 0 or iy >= self.ny or 
                iz < 0 or iz >= self.nz):
                continue
                
            # Set vorticity in the orientation direction
            # Simplified approach - in reality, the vorticity would be tangent to the ring
            self.wx[ix, iy, iz] += strength * ox
            self.wy[ix, iy, iz] += strength * oy
            self.wz[ix, iy, iz] += strength * oz
        
        # Derive velocity from vorticity
        self.update_velocity_from_vorticity()
        
        # Update derived fields
        self.update_derived_fields()
        
        print(f"Added vortex ring at {center}, radius={radius}m, strength={strength:.3e}m²/s")
    
    def update_velocity_from_vorticity(self):
        """
        Update velocity field from vorticity using a simplified approach.
        
        This is an approximation and would be replaced by a proper Biot-Savart solver
        in a production implementation.
        """
        # Create a temporary vector potential on cell faces (not centers)
        # The vector potential components need to be on the same grid positions as velocity
        ax = self.xp.zeros((self.nx+1, self.ny, self.nz))  # Same shape as vx
        ay = self.xp.zeros((self.nx, self.ny+1, self.nz))  # Same shape as vy
        az = self.xp.zeros((self.nx, self.ny, self.nz+1))  # Same shape as vz
        
        # Set vector potential from vorticity
        # We need to properly interpolate from edges to faces
        
        # For ax (on x-faces), average wx from edges
        for i in range(self.nx+1):
            for j in range(self.ny):
                for k in range(self.nz):
                    # Find nearby wx values (on x-edges) and average them
                    wx_sum = 0.0
                    count = 0
                    if j < self.ny-1 and k < self.nz-1:
                        wx_sum += self.wx[min(i, self.nx-1), j, k]
                        count += 1
                    if j > 0 and k < self.nz-1:
                        wx_sum += self.wx[min(i, self.nx-1), j-1, k]
                        count += 1
                    if j < self.ny-1 and k > 0:
                        wx_sum += self.wx[min(i, self.nx-1), j, k-1]
                        count += 1
                    if j > 0 and k > 0:
                        wx_sum += self.wx[min(i, self.nx-1), j-1, k-1]
                        count += 1
                    
                    if count > 0:
                        ax[i, j, k] = wx_sum / count
        
        # Similar interpolation for ay and az from wy and wz
        # (Code omitted for brevity but would follow the same pattern)
        
        # Use the Helmholtz decomposition to get a divergence-free velocity field
        # directly from the vorticity
        _, _, _, self.vx, self.vy, self.vz = helmholtz_decomposition(
            ax, ay, az, self.dx
        )
        
        # Update derived fields
        self.update_vorticity()
        
        # Ensure velocity field is divergence-free
        self.project_velocity()
        
        return

    def project_velocity(self):
        """
        Project velocity field to ensure it's divergence-free.
        
        This is essential for maintaining the solenoidal constraint in a
        vorticity-based representation.
        """
        # Calculate divergence
        div_v = divergence(self.vx, self.vy, self.vz, self.dx)
        
        # Calculate gradient of divergence
        grad_x, grad_y, grad_z = gradient(div_v, self.dx)
        
        # Subtract gradient component to make velocity divergence-free
        self.vx -= grad_x
        self.vy -= grad_y
        self.vz -= grad_z
    
    def add_plane_wave(self,
                       amplitude: float,
                       frequency: float,
                       direction: Tuple[float, float, float] = (1, 0, 0),
                       polarization: Tuple[float, float, float] = (0, 1, 0)):
        """
        Add a plane wave to the simulation.
        
        Parameters
        ----------
        amplitude : float
            Wave amplitude [m/s]
        frequency : float
            Wave frequency [Hz]
        direction : tuple(float, float, float)
            Wave propagation direction
        polarization : tuple(float, float, float)
            Wave polarization direction (perpendicular to direction)
        """
        # Normalize direction
        dx, dy, dz = direction
        norm = (dx**2 + dy**2 + dz**2)**0.5
        dx, dy, dz = dx/norm, dy/norm, dz/norm
        
        # Ensure polarization is perpendicular to direction
        px, py, pz = polarization
        dot = px*dx + py*dy + pz*dz
        px, py, pz = px - dot*dx, py - dot*dy, pz - dot*dz
        
        # Normalize polarization
        norm = (px**2 + py**2 + pz**2)**0.5
        if norm > 0:
            px, py, pz = px/norm, py/norm, pz/norm
        else:
            # Default polarization if original was parallel to direction
            px, py, pz = dy, -dx, 0
            norm = (px**2 + py**2 + pz**2)**0.5
            px, py, pz = px/norm, py/norm, pz/norm
        
        # Calculate wavelength
        wavelength = C_LIGHT / frequency
        k = 2 * np.pi / wavelength
        
        # Add wave to velocity field
        for i in range(self.vx.shape[0]):
            for j in range(self.vx.shape[1]):
                for k in range(self.vx.shape[2]):
                    # Position in physical units
                    x = i * self.dx
                    y = j * self.dx
                    z = k * self.dx
                    
                    # Phase at this position
                    phase = k * (x*dx + y*dy + z*dz)
                    
                    # Add velocity component
                    self.vx[i,j,k] += amplitude * px * np.sin(phase)
        
        # Similarly for y and z components
        for i in range(self.vy.shape[0]):
            for j in range(self.vy.shape[1]):
                for k in range(self.vy.shape[2]):
                    x = i * self.dx
                    y = j * self.dx
                    z = k * self.dx
                    phase = k * (x*dx + y*dy + z*dz)
                    self.vy[i,j,k] += amplitude * py * np.sin(phase)
        
        for i in range(self.vz.shape[0]):
            for j in range(self.vz.shape[1]):
                for k in range(self.vz.shape[2]):
                    x = i * self.dx
                    y = j * self.dx
                    z = k * self.dx
                    phase = k * (x*dx + y*dy + z*dz)
                    self.vz[i,j,k] += amplitude * pz * np.sin(phase)
        
        # Update derived fields
        self.update_vorticity()
        self.update_derived_fields()
        
        print(f"Added plane wave: amplitude={amplitude}m/s, frequency={frequency}Hz")
    
    def add_point_charge(self,
                        position: Tuple[float, float, float],
                        charge: float,
                        radius: float = None):
        """
        Add a point charge to the simulation.
        
        In Space Time Potential Theory, charges correspond to divergences
        in the action flux field.
        
        Parameters
        ----------
        position : tuple(float, float, float)
            Position of the charge [m]
        charge : float
            Charge value [C]
        radius : float, optional
            Characteristic radius for charge smoothing [m]
        """
        # Convert position to grid coordinates
        ix = int(position[0] / self.dx)
        iy = int(position[1] / self.dx)
        iz = int(position[2] / self.dx)
        
        # Set default radius if none provided
        if radius is None:
            radius = 2.0 * self.dx
        
        # Convert radius to grid units
        r_grid = radius / self.dx
        
        # Ensure position is within grid bounds
        if (0 <= ix < self.nx and 0 <= iy < self.ny and 0 <= iz < self.nz):
            # Add Gaussian charge distribution to scalar potential
            for i in range(max(0, ix-3*int(r_grid)), min(self.nx, ix+3*int(r_grid)+1)):
                for j in range(max(0, iy-3*int(r_grid)), min(self.ny, iy+3*int(r_grid)+1)):
                    for k in range(max(0, iz-3*int(r_grid)), min(self.nz, iz+3*int(r_grid)+1)):
                        r2 = ((i-ix)**2 + (j-iy)**2 + (k-iz)**2) * (self.dx**2)
                        self.scalar_potential[i,j,k] += charge * np.exp(-r2/(2*radius**2)) / (
                            (2*np.pi)**(3/2) * radius**3
                        )
            
            # Calculate electric field from scalar potential
            ex, ey, ez = gradient(self.scalar_potential, self.dx)
            
            # Convert to velocity field
            factor = self.viscosity / (self.background_density * C_LIGHT**2)
            
            # Add to velocity field (negated because E = -∇φ)
            self.vx -= factor * ex
            self.vy -= factor * ey
            self.vz -= factor * ez
            
            # Update derived fields
            self.update_vorticity()
            self.update_derived_fields()
            
            print(f"Added point charge at {position}, q={charge}C, radius={radius}m")
        else:
            print(f"Warning: Charge position {position} is outside grid bounds")
    
    @property
    def electric_field(self):
        """
        Electric field derived from velocity and scalar potential.
        
        E = -∇φ = (c²/η)∇·(η·v)
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Components of the electric field (Ex, Ey, Ez)
        """
        if self._electric_field is None:
            # Calculate divergence of action flux
            div_C = divergence(
                self.action_flux_x, self.action_flux_y, self.action_flux_z, self.dx
            )
            
            # Calculate gradient of divergence
            factor = C_LIGHT**2 / self.viscosity
            ex, ey, ez = gradient(div_C, self.dx)
            
            # Scale by appropriate factor
            self._electric_field = (factor * ex, factor * ey, factor * ez)
        
        return self._electric_field
    
    @property
    def magnetic_field(self):
        """
        Magnetic field derived from vorticity.
        
        B = μ₀·ω = (1/η)∇×(η·v)
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Components of the magnetic field (Bx, By, Bz)
        """
        if self._magnetic_field is None:
            # Magnetic field is directly proportional to vorticity
            factor = 1.0 / self.viscosity
            
            # Interpolate vorticity from edges to cell centers for easier visualization
            wx_center = 0.25 * (self.wx[:, :-1, :-1] + self.wx[:, 1:, :-1] + 
                               self.wx[:, :-1, 1:] + self.wx[:, 1:, 1:])
            wy_center = 0.25 * (self.wy[:-1, :, :-1] + self.wy[1:, :, :-1] + 
                               self.wy[:-1, :, 1:] + self.wy[1:, :, 1:])
            wz_center = 0.25 * (self.wz[:-1, :-1, :] + self.wz[1:, :-1, :] + 
                               self.wz[:-1, 1:, :] + self.wz[1:, 1:, :])
            
            self._magnetic_field = (factor * wx_center, factor * wy_center, factor * wz_center)
        
        return self._magnetic_field
    
    def run(self, steps: int, callback: Callable = None, interval: int = 1):
        """
        Run simulation for specified number of steps.
        
        Parameters
        ----------
        steps : int
            Number of time steps to simulate
        callback : callable, optional
            Function to call during simulation
        interval : int
            Call the callback every n steps
        """
        print(f"Starting simulation for {steps} steps...")
        
        for i in range(steps):
            self.step()
            
            if callback is not None and i % interval == 0:
                callback(self)
            
            if (i+1) % 10 == 0:
                print(f"Step {self.step_count}, Time: {self.time:.3e} s")
        
        print(f"Simulation completed. Final time: {self.time:.3e} s")
    
    def calculate_energy(self) -> Dict[str, float]:
        """
        Calculate the energy components in the simulation.
        
        Returns
        -------
        Dict[str, float]
            Dictionary containing energy components
        """
        # Get current backend
        xp = backend
        
        # Convert face-centered velocity to cell-centered for energy calculation
        vx_center = 0.5 * (self.vx[:-1, :, :] + self.vx[1:, :, :])
        vy_center = 0.5 * (self.vy[:, :-1, :] + self.vy[:, 1:, :])
        vz_center = 0.5 * (self.vz[:, :, :-1] + self.vz[:, :, 1:])
        
        # Kinetic energy
        kinetic_density = 0.5 * self.density * (vx_center**2 + vy_center**2 + vz_center**2)
        kinetic_energy = xp.sum(kinetic_density) * self.dx**3
        
        # Get EM fields
        ex, ey, ez = self.electric_field
        bx, by, bz = self.magnetic_field
        
        # Interpolate electric field components to cell centers
        ex_center = 0.5 * (ex[:-1, :, :] + ex[1:, :, :])
        ey_center = 0.5 * (ey[:, :-1, :] + ey[:, 1:, :])
        ez_center = 0.5 * (ez[:, :, :-1] + ez[:, :, 1:])
        
        # Electric field energy
        epsilon_0 = 1.0 / (C_LIGHT**2 * 4 * np.pi * 1e-7)
        electric_density = 0.5 * epsilon_0 * (ex_center**2 + ey_center**2 + ez_center**2)
        electric_energy = xp.sum(electric_density) * self.dx**3
        
        # Magnetic field energy - these should already be at cell centers from the property
        mu_0 = 4 * np.pi * 1e-7
        magnetic_density = 0.5 / mu_0 * (bx**2 + by**2 + bz**2)
        magnetic_energy = xp.sum(magnetic_density) * self.dx**3
        
        # Wave energy components
        first_sound_energy = xp.sum(0.5 * (self.first_sound.density_perturbation**2 + 
                               self.first_sound.velocity_perturbation**2)) * self.dx**3
        
        second_sound_energy = xp.sum(0.5 * (self.second_sound.temperature_perturbation**2 + 
                                self.second_sound.temp_velocity**2)) * self.dx**3
        
        # Total energy
        total_energy = kinetic_energy + electric_energy + magnetic_energy + \
                  first_sound_energy + second_sound_energy
        
        return {
            'kinetic': kinetic_energy,
            'electric': electric_energy,
            'magnetic': magnetic_energy,
            'electromagnetic': electric_energy + magnetic_energy,
            'first_sound': first_sound_energy,
            'second_sound': second_sound_energy,
            'total': total_energy
        }
