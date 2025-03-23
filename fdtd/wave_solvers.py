"""
Wave Solvers for Space Time Potential Theory
============================================

This module implements specialized solvers for first and second sound 
wave phenomena in the kinetic substrate. These solvers handle the distinct
wave equations that emerge from the vector Laplacian formulation.

References:
- Theory/vector_laplacian/wave_equations.md - Wave equations theory
- Theory/equation_reference.tex - Complete equation reference
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any
from .backend import backend
from .operators import scalar_laplacian, gradient, divergence

class FirstSoundSolver:
    """
    Solver for first sound (density waves) in the kinetic substrate.
    
    First sound represents oscillations in the substrate's density, analogous
    to pressure waves in classical fluids or longitudinal waves in elastic media.
    
    The wave equation is:
    ∂²ρ'/∂t² = (η/ρ₀)∇²ρ'
    
    With wave speed c₁ = √(η/ρ₀)
    """
    
    def __init__(self, nx: int, ny: int, nz: int, dx: float, dt: float, 
                 viscosity: float, background_density: float):
        """
        Initialize the first sound solver.
        
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
        # Add backend
        self.xp = backend  # Access backend directly
        
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dt = dt
        self.viscosity = viscosity
        self.background_density = background_density
        
        # Wave speed
        self.wave_speed = np.sqrt(viscosity / background_density)
        
        # CFL condition check
        cfl = self.wave_speed * dt / dx
        if cfl > 1.0:
            print(f"Warning: CFL condition not satisfied. CFL = {cfl} > 1.0")
            print(f"Consider reducing dt or increasing dx")
        
        # Fields - use backend's array creation
        self.density_perturbation = self.xp.zeros((nx, ny, nz))
        self.velocity_perturbation = self.xp.zeros((nx, ny, nz))  # ∂ρ'/∂t
        
        # Scalar potential from density perturbation
        self.scalar_potential = self.xp.zeros((nx, ny, nz))
        
        # Velocity field derived from density perturbation
        self.vx = self.xp.zeros((nx, ny, nz))
        self.vy = self.xp.zeros((nx, ny, nz))
        self.vz = self.xp.zeros((nx, ny, nz))
    
    def add_gaussian_perturbation(self, center: Tuple[int, int, int], 
                                 amplitude: float, width: float):
        """
        Add a Gaussian density perturbation to the field.
        
        Parameters
        ----------
        center : Tuple[int, int, int]
            Center of the Gaussian in grid coordinates
        amplitude : float
            Amplitude of the perturbation
        width : float
            Width (standard deviation) of the Gaussian in grid cells
        """
        x = np.arange(self.nx)
        y = np.arange(self.ny)
        z = np.arange(self.nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        cx, cy, cz = center
        r_squared = (X - cx)**2 + (Y - cy)**2 + (Z - cz)**2
        self.density_perturbation += amplitude * np.exp(-r_squared / (2 * width**2))
    
    def update_scalar_potential(self):
        """
        Update the scalar potential from density perturbation.
        
        The scalar potential is related to the density perturbation through
        the Poisson equation: ∇²φ = -ρ'
        """
        # Simple relaxation method for Poisson equation (could be optimized)
        self.scalar_potential = np.zeros_like(self.density_perturbation)
        for _ in range(20):  # Number of iterations
            lap_phi = scalar_laplacian(self.scalar_potential, self.dx, self.dx, self.dx)
            self.scalar_potential += 0.1 * (self.density_perturbation - lap_phi)  # Relaxation
    
    def update_velocity_field(self):
        """
        Update the velocity field from scalar potential.
        
        The velocity field is the gradient of the scalar potential:
        v = ∇φ
        """
        self.vx, self.vy, self.vz = gradient(
            self.scalar_potential, self.dx, self.dx, self.dx
        )
    
    def step(self):
        """
        Perform one time step update for the first sound wave equation.
        """
        # Calculate Laplacian of density perturbation
        lap_density = scalar_laplacian(
            self.density_perturbation, self.dx, self.dx, self.dx
        )
        
        # Update velocity perturbation (∂ρ'/∂t)
        wave_speed_squared = self.viscosity / self.background_density
        self.velocity_perturbation += self.dt * wave_speed_squared * lap_density
        
        # Update density perturbation
        self.density_perturbation += self.dt * self.velocity_perturbation
        
        # Update derived fields
        self.update_scalar_potential()
        self.update_velocity_field()

class SecondSoundSolver:
    """
    Solver for second sound (temperature waves) in the kinetic substrate.
    
    Second sound represents oscillations in the substrate's temperature or energy
    distribution, analogous to entropy waves in superfluids.
    
    The wave equation is:
    ∂²θ'/∂t² = -(η/ρ)∇²(∂θ'/∂t)
    
    This unusual wave equation, with a time derivative inside the Laplacian,
    is characteristic of second sound phenomena in superfluids.
    """
    
    def __init__(self, nx: int, ny: int, nz: int, dx: float, dt: float, 
                 viscosity: float, density: float):
        """
        Initialize the second sound solver.
        
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
        density : float
            Density ρ of the substrate [kg/m³]
        """
        # Add backend
        self.xp = backend  # Access backend directly
        
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dt = dt
        self.viscosity = viscosity
        self.density = density
        
        # Fields - use backend's array creation
        self.temperature_perturbation = self.xp.zeros((nx, ny, nz))
        self.temp_velocity = self.xp.zeros((nx, ny, nz))  # ∂θ'/∂t
        self.temp_acceleration = self.xp.zeros((nx, ny, nz))  # ∂²θ'/∂t²
        
        # Temperature potential (∇·(η·a))
        self.temperature_potential = self.xp.zeros((nx, ny, nz))
        
        # Derived acceleration field
        self.ax = self.xp.zeros((nx, ny, nz))
        self.ay = self.xp.zeros((nx, ny, nz))
        self.az = self.xp.zeros((nx, ny, nz))
    
    def add_gaussian_perturbation(self, center: Tuple[int, int, int], 
                                 amplitude: float, width: float):
        """
        Add a Gaussian temperature perturbation to the field.
        
        Parameters
        ----------
        center : Tuple[int, int, int]
            Center of the Gaussian in grid coordinates
        amplitude : float
            Amplitude of the perturbation
        width : float
            Width (standard deviation) of the Gaussian in grid cells
        """
        x = np.arange(self.nx)
        y = np.arange(self.ny)
        z = np.arange(self.nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        cx, cy, cz = center
        r_squared = (X - cx)**2 + (Y - cy)**2 + (Z - cz)**2
        self.temperature_perturbation += amplitude * np.exp(-r_squared / (2 * width**2))
    
    def update_temperature_potential(self):
        """
        Update the temperature potential from temperature velocity.
        
        The temperature potential T is related to the radiosity field:
        T = ∇·(η·a)
        """
        # For the temperature wave, the potential is related to the velocity
        # This is a simplified approximation
        self.temperature_potential = self.temp_velocity
    
    def update_acceleration_field(self):
        """
        Update the acceleration field from temperature potential.
        
        For second sound, the acceleration field is related to the gradient
        of the temperature potential.
        """
        self.ax, self.ay, self.az = gradient(
            self.temperature_potential, self.dx, self.dx, self.dx
        )
        
        # Scale by -1/(ρ) to get the actual acceleration
        factor = -1.0 / self.density
        self.ax *= factor
        self.ay *= factor
        self.az *= factor
    
    def step(self):
        """
        Perform one time step update for the second sound wave equation.
        
        This implements the unusual second sound wave equation which has
        a time derivative inside the Laplacian operator.
        """
        # Calculate Laplacian of temperature velocity
        lap_temp_velocity = scalar_laplacian(
            self.temp_velocity, self.dx, self.dx, self.dx
        )
        
        # Wave equation coefficient
        coefficient = -self.viscosity / self.density
        
        # Update temperature acceleration (∂²θ'/∂t²)
        self.temp_acceleration = coefficient * lap_temp_velocity
        
        # Update temperature velocity (∂θ'/∂t)
        self.temp_velocity += self.dt * self.temp_acceleration
        
        # Update temperature perturbation
        self.temperature_perturbation += self.dt * self.temp_velocity
        
        # Update derived fields
        self.update_temperature_potential()
        self.update_acceleration_field()
