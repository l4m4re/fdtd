"""
Aether Grid: Space Time Potential Theory Simulation Framework
============================================================

This module implements a computational grid for simulating the kinetic substrate
according to Vector Laplacian Theory. It models both the action flux field (C)
and radiosity field (I) to capture electromagnetic and gravitational effects.

The simulation is based on the fundamental principle:
    a = -η/ρ·Δv

References:
- Theory/vector_laplacian/fields.md - Action flux and Radiosity fields
- Theory/equation_reference.tex - Complete equation reference
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any
from .operators import gradient, divergence, curl, vector_laplacian
from .wave_solvers import FirstSoundSolver, SecondSoundSolver
from .operators import helmholtz_decomposition, quantized_vortex

# Physical constants
VISCOSITY = 1/(4*np.pi*1e-7)  # η - dynamic viscosity [Pa·s]
C_LIGHT = 299792458.0  # Speed of light [m/s]
BACKGROUND_DENSITY = VISCOSITY/(C_LIGHT**2)  # ρ₀ = η/c² [kg/m³]

class AetherGrid:
    """
    Simulation grid for the kinetic substrate based on Space Time Potential Theory.
    
    This class implements the computational framework needed to simulate the 
    dynamics of the superfluid-like medium that constitutes space itself according
    to Vector Laplacian Theory. It tracks both first-order (velocity, vorticity) and
    second-order (acceleration, angular acceleration) fields.
    
    The simulation supports both incompressible (solenoidal) and compressible
    (irrotational) flow components, allowing it to model both electromagnetic
    and gravitational phenomena.
    """
    
    def __init__(self, 
                 nx: int, 
                 ny: int, 
                 nz: int,
                 dx: float = 1.0,
                 dt: float = 0.1,
                 viscosity: float = VISCOSITY,
                 background_density: float = BACKGROUND_DENSITY):
        """
        Initialize the aether grid simulation.
        
        Parameters
        ----------
        nx, ny, nz : int
            Grid dimensions
        dx : float
            Grid spacing (assumed uniform in all dimensions)
        dt : float
            Time step
        viscosity : float
            Dynamic viscosity η of the substrate [Pa·s]
        background_density : float
            Background density ρ₀ of the substrate [kg/m³]
        """
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx = dx
        self.dt = dt
        self.viscosity = viscosity
        self.background_density = background_density
        
        # First-order fields
        # PROMPT: Initialize these fields with appropriate initial conditions
        # based on the problem being studied (e.g., electromagnetic wave,
        # gravitational source, vortex structure)
        self.velocity = np.zeros((3, nx, ny, nz))  # Velocity field v [m/s]
        self.vorticity = np.zeros((3, nx, ny, nz))  # Vorticity ω = ∇×v [1/s]
        
        # Second-order fields
        self.acceleration = np.zeros((3, nx, ny, nz))  # Acceleration a [m/s²]
        self.angular_acceleration = np.zeros((3, nx, ny, nz))  # Angular acceleration α = ∇×a [1/s²]
        
        # Fundamental derived fields from theory
        # PROMPT: These fields are key to connecting with electromagnetic theory
        self.action_flux = np.zeros((3, nx, ny, nz))  # Action flux field C = η·v [kg/s², N/m]
        self.radiosity = np.zeros((3, nx, ny, nz))  # Radiosity field I = η·a [kg/s³, W/m²]
        
        # Potentials
        # PROMPT: The potentials connect to electromagnetic and gravitational potentials
        self.scalar_potential = np.zeros((nx, ny, nz))  # φ = ∇·(η·v) [Pa]
        self.vector_potential = np.zeros((3, nx, ny, nz))  # A = η/ρ·ω [m²/s]
        self.temperature_potential = np.zeros((nx, ny, nz))  # T = ∇·(η·a) [Pa/s]
        self.second_vector_potential = np.zeros((3, nx, ny, nz))  # χ = η/ρ·α [m²/s²]
        
        # Density variation for compressible dynamics
        self.density = np.ones((nx, ny, nz)) * background_density  # ρ [kg/m³]
        self.density_perturbation = np.zeros((nx, ny, nz))  # ρ' (for first sound)
        
        # Temperature variation for second sound effects
        self.temperature = np.zeros((nx, ny, nz))  # θ [K]
        self.temperature_perturbation = np.zeros((nx, ny, nz))  # θ' (for second sound)
        
        # Electromagnetic field equivalents
        # PROMPT: These are computed from the action flux field
        self.electric_field = np.zeros((3, nx, ny, nz))  # E = ∇·C/ε₀ [V/m, N/C]
        self.magnetic_field = np.zeros((3, nx, ny, nz))  # B = μ₀·∇×C [T, Wb/m²]
        
        # Add wave solvers
        self.first_sound_solver = FirstSoundSolver(
            nx, ny, nz, dx, dt, viscosity, background_density
        )
        self.second_sound_solver = SecondSoundSolver(
            nx, ny, nz, dx, dt, viscosity, background_density
        )
    
    def update_vorticity(self):
        """
        Update the vorticity field from the velocity field.
        
        Vorticity ω = ∇×v represents the local rotation of the substrate
        and is connected to the magnetic field in electromagnetic theory.
        """
        # PROMPT: Consider implementing vortex quantization constraints
        # for modeling particle-like structures
        
        # Extract velocity components
        vx = self.velocity[0]
        vy = self.velocity[1]
        vz = self.velocity[2]
        
        # Calculate curl
        wx, wy, wz = curl(vx, vy, vz, self.dx, self.dx, self.dx)
        
        # Update vorticity field
        self.vorticity[0] = wx
        self.vorticity[1] = wy
        self.vorticity[2] = wz
    
    def update_angular_acceleration(self):
        """
        Update the angular acceleration field from the acceleration field.
        
        Angular acceleration α = ∇×a is the time derivative of vorticity
        and relates to the time evolution of magnetic fields.
        """
        # PROMPT: This field is important for understanding dynamic
        # electromagnetic phenomena like induction
        
        # Extract acceleration components
        ax = self.acceleration[0]
        ay = self.acceleration[1]
        az = self.acceleration[2]
        
        # Calculate curl
        alphax, alphay, alphaz = curl(ax, ay, az, self.dx, self.dx, self.dx)
        
        # Update angular acceleration field
        self.angular_acceleration[0] = alphax
        self.angular_acceleration[1] = alphay
        self.angular_acceleration[2] = alphaz
    
    def update_action_flux(self):
        """
        Update the action flux field from velocity.
        
        The action flux field C = η·v is fundamental to the theory,
        corresponding to momentum flux density in the substrate.
        """
        # PROMPT: This field is directly related to electromagnetic
        # vector potential A and current density J
        
        self.action_flux[0] = self.viscosity * self.velocity[0]
        self.action_flux[1] = self.viscosity * self.velocity[1]
        self.action_flux[2] = self.viscosity * self.velocity[2]
    
    def update_radiosity(self):
        """
        Update the radiosity field from acceleration.
        
        The radiosity field I = η·a represents energy flux density
        and relates to Poynting vector in electromagnetics.
        """
        # PROMPT: This field is important for understanding energy
        # transfer and second-order effects like gravity
        
        self.radiosity[0] = self.viscosity * self.acceleration[0]
        self.radiosity[1] = self.viscosity * self.acceleration[1]
        self.radiosity[2] = self.viscosity * self.acceleration[2]
    
    def update_potentials(self):
        """
        Update the scalar and vector potentials from the fundamental fields.
        
        These potentials are key to connecting with standard electromagnetic
        and gravitational theories.
        """
        # PROMPT: The scalar potential directly relates to electric potential
        # and Newtonian gravitational potential
        
        # Action flux components
        cx = self.action_flux[0]
        cy = self.action_flux[1]
        cz = self.action_flux[2]
        
        # Scalar potential φ = ∇·C
        self.scalar_potential = divergence(cx, cy, cz, self.dx, self.dx, self.dx)
        
        # Vector potential A = η/ρ·ω
        factor = self.viscosity / self.background_density
        self.vector_potential[0] = factor * self.vorticity[0]
        self.vector_potential[1] = factor * self.vorticity[1]
        self.vector_potential[2] = factor * self.vorticity[2]
        
        # Radiosity components
        ix = self.radiosity[0]
        iy = self.radiosity[1]
        iz = self.radiosity[2]
        
        # Temperature potential T = ∇·I
        self.temperature_potential = divergence(ix, iy, iz, self.dx, self.dx, self.dx)
        
        # Second vector potential χ = η/ρ·α
        self.second_vector_potential[0] = factor * self.angular_acceleration[0]
        self.second_vector_potential[1] = factor * self.angular_acceleration[1]
        self.second_vector_potential[2] = factor * self.angular_acceleration[2]
    
    def update_electromagnetic_fields(self):
        """
        Update electromagnetic field equivalents from substrate fields.
        
        This method computes the electric and magnetic fields directly
        from operations on the action flux field C.
        """
        # PROMPT: These calculations demonstrate the direct connection
        # between substrate dynamics and Maxwell's equations
        
        # Action flux components
        cx = self.action_flux[0]
        cy = self.action_flux[1]
        cz = self.action_flux[2]
        
        # Electric field E = ∇·C/ε₀
        # Note: 1/ε₀ = c²μ₀ = c²/η
        div_C = divergence(cx, cy, cz, self.dx, self.dx, self.dx)
        factor_E = (C_LIGHT**2) / self.viscosity
        
        # Gradient of scalar potential
        ex, ey, ez = gradient(div_C, self.dx, self.dx, self.dx)
        
        self.electric_field[0] = factor_E * ex
        self.electric_field[1] = factor_E * ey
        self.electric_field[2] = factor_E * ez
        
        # Magnetic field B = μ₀·∇×C
        # Note: μ₀ = 1/η
        bx, by, bz = curl(cx, cy, cz, self.dx, self.dx, self.dx)
        factor_B = 1.0 / self.viscosity
        
        self.magnetic_field[0] = factor_B * bx
        self.magnetic_field[1] = factor_B * by
        self.magnetic_field[2] = factor_B * bz
    
    def update_acceleration(self):
        """
        Update the acceleration field using the vector Laplacian of velocity.
        
        According to the fundamental principle of the theory:
        a = -η/ρ·Δv
        """
        # PROMPT: This is the core equation of the theory, connecting
        # directly to Navier-Stokes with the characteristic sign reversal
        
        # Extract velocity components
        vx = self.velocity[0]
        vy = self.velocity[1]
        vz = self.velocity[2]
        
        # Calculate vector Laplacian
        lap_vx, lap_vy, lap_vz = vector_laplacian(vx, vy, vz, self.dx, self.dx, self.dx)
        
        # Update acceleration field
        factor = -self.viscosity / self.background_density
        self.acceleration[0] = factor * lap_vx
        self.acceleration[1] = factor * lap_vy
        self.acceleration[2] = factor * lap_vz
    
    def update(self):
        """
        Perform a complete update step for all fields.
        
        This integrates the velocity field using the current acceleration,
        then updates all derived fields.
        """
        # PROMPT: Consider splitting this into separate update methods
        # for first and second sound phenomena
        
        # Update velocity field using acceleration
        self.velocity[0] += self.acceleration[0] * self.dt
        self.velocity[1] += self.acceleration[1] * self.dt
        self.velocity[2] += self.acceleration[2] * self.dt
        
        # Update derived fields
        self.update_vorticity()
        self.update_action_flux()
        self.update_acceleration()
        self.update_angular_acceleration()
        self.update_radiosity()
        self.update_potentials()
        self.update_electromagnetic_fields()
        
        # Decompose velocity field
        self.decompose_velocity_field()
        
        # Update wave solvers
        self.update_wave_solvers()
        
        # Continue with remaining updates
    
    def add_charge(self, 
                   x: float, 
                   y: float, 
                   z: float, 
                   charge: float, 
                   radius: float = 1.0):
        """
        Add an electric charge at the specified position.
        
        In Space Time Potential Theory, charges are represented as
        divergences in the action flux field.
        
        Parameters
        ----------
        x, y, z : float
            Position of the charge
        charge : float
            Charge value [C]
        radius : float
            Characteristic radius of the charge distribution
        """
        # PROMPT: Implement proper charge representation consistent
        # with theory, possibly as a divergence in action flux
        
        # Convert grid coordinates to array indices
        ix = int(x / self.dx)
        iy = int(y / self.dx)
        iz = int(z / self.dx)
        
        # Ensure indices are within grid bounds
        if 0 <= ix < self.nx and 0 <= iy < self.ny and 0 <= iz < self.nz:
            # Simple point charge for now
            # PROMPT: Replace with proper Gaussian or analytical charge distribution
            self.scalar_potential[ix, iy, iz] += charge
    
    def add_vortex(self,
                   x: float,
                   y: float,
                   z: float,
                   circulation: float,
                   radius: float = 1.0,
                   direction: Tuple[float, float, float] = (0, 0, 1)):
        """
        Add a vortex at the specified position.
        
        In Space Time Potential Theory, particles can be represented as
        quantized vortices in the substrate.
        
        Parameters
        ----------
        x, y, z : float
            Position of the vortex center
        circulation : float
            Circulation strength [m²/s]
        radius : float
            Characteristic radius of the vortex
        direction : Tuple[float, float, float]
            Direction of the vortex axis
        """
        # PROMPT: Implement proper quantized vortex consistent with theory
        # Circulation should be quantized as h/m for particles
        
        # Convert grid coordinates to array indices
        ix = int(x / self.dx)
        iy = int(y / self.dx)
        iz = int(z / self.dx)
        
        # Normalize direction vector
        dx, dy, dz = direction
        norm = np.sqrt(dx**2 + dy**2 + dz**2)
        dx, dy, dz = dx/norm, dy/norm, dz/norm
        
        # Simple vortex model for now
        # PROMPT: Replace with proper analytical vortex solution
        
        # Calculate indices within radius of vortex center
        x_indices = np.arange(max(0, ix-radius), min(self.nx, ix+radius+1))
        y_indices = np.arange(max(0, iy-radius), min(self.ny, iy+radius+1))
        z_indices = np.arange(max(0, iz-radius), min(self.nz, iz+radius+1))
        
        # Set vorticity along the direction
        for i in x_indices:
            for j in y_indices:
                for k in z_indices:
                    dist = np.sqrt((i-ix)**2 + (j-iy)**2 + (k-iz)**2)
                    if dist <= radius:
                        intensity = circulation * (1 - dist/radius)
                        self.vorticity[0, i, j, k] = intensity * dx
                        self.vorticity[1, i, j, k] = intensity * dy
                        self.vorticity[2, i, j, k] = intensity * dz
    
    def decompose_velocity_field(self):
        """
        Decompose the velocity field into irrotational and solenoidal components.
        
        This is crucial for:
        - Understanding which components drive electromagnetic vs. gravitational effects
        - Separating first sound (density waves) from rotational dynamics
        """
        vx = self.velocity[0]
        vy = self.velocity[1]
        vz = self.velocity[2]
        
        vx_irrot, vy_irrot, vz_irrot, vx_sol, vy_sol, vz_sol = helmholtz_decomposition(
            vx, vy, vz, self.dx, self.dx, self.dx
        )
        
        # Store decomposed fields
        self.irrotational_velocity = np.array([vx_irrot, vy_irrot, vz_irrot])
        self.solenoidal_velocity = np.array([vx_sol, vy_sol, vz_sol])
        
        # Update density perturbation based on divergence of irrotational component
        div_irrot = divergence(vx_irrot, vy_irrot, vz_irrot, self.dx, self.dx, self.dx)
        self.first_sound_solver.density_perturbation = div_irrot
        
        # Update vorticity based on curl of solenoidal component
        wx, wy, wz = curl(vx_sol, vy_sol, vz_sol, self.dx, self.dx, self.dx)
        self.vorticity[0] = wx
        self.vorticity[1] = wy
        self.vorticity[2] = wz
    
    def update_wave_solvers(self):
        """
        Update both first and second sound wave solvers.
        """
        self.first_sound_solver.step()
        self.second_sound_solver.step()
        
        # Integrate wave solver results back into main simulation
        # This is crucial for capturing the interplay between different wave types
        
        # Add irrotational component from first sound solver
        self.velocity[0] += self.first_sound_solver.vx
        self.velocity[1] += self.first_sound_solver.vy
        self.velocity[2] += self.first_sound_solver.vz
        
        # Add acceleration component from second sound solver
        self.acceleration[0] += self.second_sound_solver.ax
        self.acceleration[1] += self.second_sound_solver.ay
        self.acceleration[2] += self.second_sound_solver.az
    
    def add_quantized_vortex(self, 
                             center: Tuple[float, float, float], 
                             circulation: float,
                             axis: Tuple[float, float, float] = (0, 0, 1),
                             core_radius: float = 3.0):
        """
        Add a quantized vortex to the velocity field.
        
        In Space Time Potential Theory, particles are represented as
        quantized vortices in the substrate.
        
        Parameters
        ----------
        center : Tuple[float, float, float]
            Position of the vortex center in grid coordinates
        circulation : float
            Strength of circulation [m²/s]
        axis : Tuple[float, float, float]
            Direction of the vortex axis
        core_radius : float
            Radius of the vortex core in grid cells
        """
        # Convert center to integer indices
        ix = int(center[0] / self.dx)
        iy = int(center[1] / self.dx)
        iz = int(center[2] / self.dx)
        
        # Generate quantized vortex velocity field
        vx_vortex, vy_vortex, vz_vortex = quantized_vortex(
            self.nx, self.ny, self.nz, (ix, iy, iz), 
            circulation, axis, int(core_radius)
        )
        
        # Add to velocity field
        self.velocity[0] += vx_vortex
        self.velocity[1] += vy_vortex
        self.velocity[2] += vz_vortex
        
        # Update derived fields
        self.update_vorticity()
        self.update_action_flux()
    
    def calculate_energy(self) -> Dict[str, float]:
        """
        Calculate the energy components in the simulation.
        
        Returns
        -------
        Dict[str, float]
            Dictionary with different energy components:
            - 'kinetic': Kinetic energy
            - 'electromagnetic': Electromagnetic energy
            - 'first_sound': Energy in density waves
            - 'second_sound': Energy in temperature waves
            - 'total': Total energy
        """
        # Calculate kinetic energy
        vx = self.velocity[0]
        vy = self.velocity[1]
        vz = self.velocity[2]
        kinetic_density = 0.5 * self.background_density * (vx**2 + vy**2 + vz**2)
        kinetic_energy = np.sum(kinetic_density) * self.dx**3
        
        # Calculate electromagnetic energy
        ex = self.electric_field[0]
        ey = self.electric_field[1]
        ez = self.electric_field[2]
        bx = self.magnetic_field[0]
        by = self.magnetic_field[1]
        bz = self.magnetic_field[2]
        
        # ε₀ = 1/(c²μ₀) = 1/(c²/η)
        epsilon_0 = 1.0 / (VISCOSITY * (C_LIGHT**2))
        
        e_energy_density = 0.5 * epsilon_0 * (ex**2 + ey**2 + ez**2)
        b_energy_density = 0.5 / VISCOSITY * (bx**2 + by**2 + bz**2)
        em_energy = np.sum(e_energy_density + b_energy_density) * self.dx**3
        
        # First sound energy
        rho_p = self.first_sound_solver.density_perturbation
        v_p = self.first_sound_solver.velocity_perturbation
        first_sound_energy = np.sum(0.5 * (rho_p**2 + v_p**2)) * self.dx**3
        
        # Second sound energy 
        theta_p = self.second_sound_solver.temperature_perturbation
        v_theta = self.second_sound_solver.temp_velocity
        second_sound_energy = np.sum(0.5 * (theta_p**2 + v_theta**2)) * self.dx**3
        
        # Total energy
        total_energy = kinetic_energy + em_energy + first_sound_energy + second_sound_energy
        
        return {
            'kinetic': kinetic_energy,
            'electromagnetic': em_energy,
            'first_sound': first_sound_energy,
            'second_sound': second_sound_energy,
            'total': total_energy
        }

# PROMPT: Additional classes to implement:
# 1. FirstSoundSolver - For simulating density perturbations
# 2. SecondSoundSolver - For simulating temperature waves
# 3. AetherParticle - For modeling quantized vortex structures in the substrate
