"""
Vector Calculus Operators for Space Time Potential Theory Simulations
====================================================================

This module implements discrete versions of key vector calculus operators
required for simulating the kinetic substrate according to Vector Laplacian Theory.

The Space Time Potential Theory fundamentally relies on the vector Laplacian:
    Δv = ∇(∇·v) - ∇×(∇×v)

References:
- Theory/vector_laplacian/foundation.md - Mathematical foundations
- Theory/equation_reference.tex - Complete equation reference
"""

import numpy as np
from typing import Tuple

# PROMPT: Consider implementing advanced staggered grid techniques 
# for better representation of vector Laplacian components

def gradient(scalar_field: np.ndarray, dx: float = 1.0, dy: float = 1.0, dz: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the gradient of a scalar field.
    
    In Space Time Potential Theory, gradients are applied to:
    - Scalar potential φ = ∇·(η·v) to obtain force density
    - Temperature potential T = ∇·(η·a) to obtain yank density
    
    Parameters
    ----------
    scalar_field : np.ndarray
        3D scalar field
    dx, dy, dz : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the gradient vector field (gx, gy, gz)
    """
    # PROMPT: Evaluate boundary conditions for gradient operation
    # to maintain consistency with substrate physics
    
    gx = np.zeros_like(scalar_field)
    gy = np.zeros_like(scalar_field)
    gz = np.zeros_like(scalar_field)
    
    # x-component: central difference
    gx[1:-1, :, :] = (scalar_field[2:, :, :] - scalar_field[:-2, :, :]) / (2 * dx)
    # Forward/backward difference at boundaries
    gx[0, :, :] = (scalar_field[1, :, :] - scalar_field[0, :, :]) / dx
    gx[-1, :, :] = (scalar_field[-1, :, :] - scalar_field[-2, :, :]) / dx
    
    # y-component: central difference
    gy[:, 1:-1, :] = (scalar_field[:, 2:, :] - scalar_field[:, :-2, :]) / (2 * dy)
    # Forward/backward difference at boundaries
    gy[:, 0, :] = (scalar_field[:, 1, :] - scalar_field[:, 0, :]) / dy
    gy[:, -1, :] = (scalar_field[:, -1, :] - scalar_field[:, -2, :]) / dy
    
    # z-component: central difference
    gz[:, :, 1:-1] = (scalar_field[:, :, 2:] - scalar_field[:, :, :-2]) / (2 * dz)
    # Forward/backward difference at boundaries
    gz[:, :, 0] = (scalar_field[:, :, 1] - scalar_field[:, :, 0]) / dz
    gz[:, :, -1] = (scalar_field[:, :, -1] - scalar_field[:, :, -2]) / dz
    
    return gx, gy, gz

def divergence(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float = 1.0, dy: float = 1.0, dz: float = 1.0) -> np.ndarray:
    """
    Calculate the divergence of a vector field.
    
    In Space Time Potential Theory, divergence is applied to:
    - Action Flux field C = η·v to obtain scalar potential φ
    - Radiosity field I = η·a to obtain temperature potential T
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field
    dx, dy, dz : float
        Grid spacing in each dimension
        
    Returns
    -------
    np.ndarray
        Divergence scalar field
    """
    # PROMPT: Consider compressible vs incompressible treatment
    # for different substrate regimes
    
    div = np.zeros_like(vx)
    
    # Central difference for each component's derivative
    div[1:-1, :, :] += (vx[2:, :, :] - vx[:-2, :, :]) / (2 * dx)
    div[:, 1:-1, :] += (vy[:, 2:, :] - vy[:, :-2, :]) / (2 * dy)
    div[:, :, 1:-1] += (vz[:, :, 2:] - vz[:, :, :-2]) / (2 * dz)
    
    # Forward/backward difference at boundaries (x)
    div[0, :, :] += (vx[1, :, :] - vx[0, :, :]) / dx
    div[-1, :, :] += (vx[-1, :, :] - vx[-2, :, :]) / dx
    
    # Forward/backward difference at boundaries (y)
    div[:, 0, :] += (vy[:, 1, :] - vy[:, 0, :]) / dy
    div[:, -1, :] += (vy[:, -1, :] - vy[:, -2, :]) / dy
    
    # Forward/backward difference at boundaries (z)
    div[:, :, 0] += (vz[:, :, 1] - vz[:, :, 0]) / dz
    div[:, :, -1] += (vz[:, :, -1] - vz[:, :, -2]) / dz
    
    return div

def curl(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float = 1.0, dy: float = 1.0, dz: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the curl of a vector field.
    
    In Space Time Potential Theory, curl is applied to:
    - Velocity field v to obtain vorticity ω
    - Acceleration field a to obtain angular acceleration α
    - Action Flux field C = η·v to obtain magnetic field B
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field
    dx, dy, dz : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the curl vector field (curlx, curly, curlz)
    """
    # PROMPT: Implement quantized vorticity constraints
    # to model topological defects in substrate
    
    curlx = np.zeros_like(vx)
    curly = np.zeros_like(vy)
    curlz = np.zeros_like(vz)
    
    # x-component: dw/dy - dv/dz
    curly[:, 1:-1, :] += (vz[:, 2:, :] - vz[:, :-2, :]) / (2 * dy)
    curly[:, :, 1:-1] -= (vy[:, :, 2:] - vy[:, :, :-2]) / (2 * dz)
    
    # y-component: du/dz - dw/dx
    curly[:, :, 1:-1] += (vx[:, :, 2:] - vx[:, :, :-2]) / (2 * dz)
    curly[1:-1, :, :] -= (vz[2:, :, :] - vz[:-2, :, :]) / (2 * dx)
    
    # z-component: dv/dx - du/dy
    curlz[1:-1, :, :] += (vy[2:, :, :] - vy[:-2, :, :]) / (2 * dx)
    curlz[:, 1:-1, :] -= (vx[:, 2:, :] - vx[:, :-2, :]) / (2 * dy)
    
    # Boundary conditions (simplified version)
    # PROMPT: Implement proper boundary conditions for curl
    # consistent with substrate physics
    
    return curlx, curly, curlz

def vector_laplacian(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float = 1.0, dy: float = 1.0, dz: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the vector Laplacian of a vector field.
    
    The vector Laplacian is the fundamental operator in Space Time Potential Theory:
    Δv = ∇(∇·v) - ∇×(∇×v)
    
    Used to determine:
    - Acceleration: a = -η/ρ·Δv
    - Force density: f = -η·Δv
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field
    dx, dy, dz : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the vector Laplacian field
    """
    # PROMPT: Consider implementing both formulations:
    # 1. As gradient of divergence minus curl of curl
    # 2. As component-wise scalar Laplacian
    # and analyze differences in numerical stability
    
    # Method 1: Component-wise scalar Laplacian (more stable numerically)
    lap_vx = np.zeros_like(vx)
    lap_vy = np.zeros_like(vy)
    lap_vz = np.zeros_like(vz)
    
    # x-component
    lap_vx[1:-1, :, :] += (vx[2:, :, :] - 2*vx[1:-1, :, :] + vx[:-2, :, :]) / (dx**2)
    lap_vx[:, 1:-1, :] += (vx[:, 2:, :] - 2*vx[:, 1:-1, :] + vx[:, :-2, :]) / (dy**2)
    lap_vx[:, :, 1:-1] += (vx[:, :, 2:] - 2*vx[:, :, 1:-1] + vx[:, :, :-2]) / (dz**2)
    
    # y-component
    lap_vy[1:-1, :, :] += (vy[2:, :, :] - 2*vy[1:-1, :, :] + vy[:-2, :, :]) / (dx**2)
    lap_vy[:, 1:-1, :] += (vy[:, 2:, :] - 2*vy[:, 1:-1, :] + vy[:, :-2, :]) / (dy**2)
    lap_vy[:, :, 1:-1] += (vy[:, :, 2:] - 2*vy[:, :, 1:-1] + vy[:, :, :-2]) / (dz**2)
    
    # z-component
    lap_vz[1:-1, :, :] += (vz[2:, :, :] - 2*vz[1:-1, :, :] + vz[:-2, :, :]) / (dx**2)
    lap_vz[:, 1:-1, :] += (vz[:, 2:, :] - 2*vz[:, 1:-1, :] + vz[:, :-2, :]) / (dy**2)
    lap_vz[:, :, 1:-1] += (vz[:, :, 2:] - 2*vz[:, :, 1:-1] + vz[:, :, :-2]) / (dz**2)
    
    # PROMPT: Method 2 could be implemented as:
    # div = divergence(vx, vy, vz, dx, dy, dz)
    # grad_div_x, grad_div_y, grad_div_z = gradient(div, dx, dy, dz)
    # curl_x, curl_y, curl_z = curl(vx, vy, vz, dx, dy, dz)
    # curl_curl_x, curl_curl_y, curl_curl_z = curl(curl_x, curl_y, curl_z, dx, dy, dz)
    # lap_vx = grad_div_x - curl_curl_x
    # lap_vy = grad_div_y - curl_curl_y
    # lap_vz = grad_div_z - curl_curl_z
    
    return lap_vx, lap_vy, lap_vz

def scalar_laplacian(scalar_field: np.ndarray, dx: float = 1.0, dy: float = 1.0, dz: float = 1.0) -> np.ndarray:
    """
    Calculate the scalar Laplacian of a scalar field.
    
    In Space Time Potential Theory, this is used for:
    - Wave equations for density perturbations (first sound)
    - Wave equations for temperature perturbations (second sound)
    
    Parameters
    ----------
    scalar_field : np.ndarray
        3D scalar field
    dx, dy, dz : float
        Grid spacing in each dimension
        
    Returns
    -------
    np.ndarray
        Scalar Laplacian field
    """
    # PROMPT: Implement wave equation solver using this laplacian
    # for both first and second sound phenomena
    
    lap = np.zeros_like(scalar_field)
    
    lap[1:-1, :, :] += (scalar_field[2:, :, :] - 2*scalar_field[1:-1, :, :] + scalar_field[:-2, :, :]) / (dx**2)
    lap[:, 1:-1, :] += (scalar_field[:, 2:, :] - 2*scalar_field[:, 1:-1, :] + scalar_field[:, :-2, :]) / (dy**2)
    lap[:, :, 1:-1] += (scalar_field[:, :, 2:] - 2*scalar_field[:, :, 1:-1] + scalar_field[:, :, :-2]) / (dz**2)
    
    # Simplified boundary conditions
    # PROMPT: Improve boundary handling for scalar Laplacian
    
    return lap

def helmholtz_decomposition(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, 
                           dx: float = 1.0, dy: float = 1.0, dz: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Decompose a vector field into irrotational (curl-free) and solenoidal (divergence-free) components.
    
    In Space Time Potential Theory, this decomposition is fundamental for:
    - Separating compressible/irrotational flow (related to first sound and electric field)
    - Separating incompressible/solenoidal flow (related to vorticity and magnetic field)
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field
    dx, dy, dz : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Components of irrotational field (vx_irrot, vy_irrot, vz_irrot) and 
        solenoidal field (vx_sol, vy_sol, vz_sol)
    """
    # Calculate divergence
    div_v = divergence(vx, vy, vz, dx, dy, dz)
    
    # Solve Poisson equation for scalar potential
    # Using simple Jacobi iteration for demonstration (would use more efficient solver in practice)
    scalar_potential = np.zeros_like(div_v)
    for _ in range(100):  # Number of iterations
        lap_phi = scalar_laplacian(scalar_potential, dx, dy, dz)
        scalar_potential += 0.1 * (div_v - lap_phi)  # Relaxation factor 0.1
    
    # Calculate gradient of scalar potential to get irrotational component
    vx_irrot, vy_irrot, vz_irrot = gradient(scalar_potential, dx, dy, dz)
    
    # Solenoidal component is the difference
    vx_sol = vx - vx_irrot
    vy_sol = vy - vy_irrot
    vz_sol = vz - vz_irrot
    
    return vx_irrot, vy_irrot, vz_irrot, vx_sol, vy_sol, vz_sol

def first_sound_step(density_perturbation: np.ndarray, velocity_perturbation: np.ndarray, 
                    dt: float, dx: float, viscosity: float, background_density: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform one time step of the first sound wave equation.
    
    Implements the wave equation:
    ∂²ρ'/∂t² = (η/ρ₀)∇²ρ'
    
    Parameters
    ----------
    density_perturbation : np.ndarray
        Current density perturbation field
    velocity_perturbation : np.ndarray
        Current time derivative of density perturbation
    dt : float
        Time step
    dx : float
        Grid spacing
    viscosity : float
        Substrate viscosity η
    background_density : float
        Background density ρ₀
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Updated density perturbation and velocity perturbation
    """
    # Calculate Laplacian of density perturbation
    lap_density = scalar_laplacian(density_perturbation, dx, dx, dx)
    
    # Wave equation coefficient
    wave_speed_squared = viscosity / background_density
    
    # Update velocity perturbation (∂ρ'/∂t)
    velocity_perturbation += dt * wave_speed_squared * lap_density
    
    # Update density perturbation
    density_perturbation += dt * velocity_perturbation
    
    return density_perturbation, velocity_perturbation

def second_sound_step(temperature_perturbation: np.ndarray, temp_velocity: np.ndarray, temp_acceleration: np.ndarray,
                     dt: float, dx: float, viscosity: float, density: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform one time step of the second sound wave equation.
    
    Implements the wave equation:
    ∂²θ'/∂t² = -(η/ρ)∇²(∂θ'/∂t)
    
    Parameters
    ----------
    temperature_perturbation : np.ndarray
        Current temperature perturbation field
    temp_velocity : np.ndarray
        Current time derivative of temperature (∂θ'/∂t)
    temp_acceleration : np.ndarray
        Current second time derivative of temperature (∂²θ'/∂t²)
    dt : float
        Time step
    dx : float
        Grid spacing
    viscosity : float
        Substrate viscosity η
    density : float
        Substrate density ρ
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Updated temperature perturbation, velocity, and acceleration
    """
    # Calculate Laplacian of temperature velocity
    lap_temp_velocity = scalar_laplacian(temp_velocity, dx, dx, dx)
    
    # Wave equation coefficient
    coefficient = -viscosity / density
    
    # Update temperature acceleration (∂²θ'/∂t²)
    temp_acceleration = coefficient * lap_temp_velocity
    
    # Update temperature velocity (∂θ'/∂t)
    temp_velocity += dt * temp_acceleration
    
    # Update temperature perturbation
    temperature_perturbation += dt * temp_velocity
    
    return temperature_perturbation, temp_velocity, temp_acceleration

def quantized_vortex(nx: int, ny: int, nz: int, center: Tuple[int, int, int], 
                     circulation: float, axis: Tuple[float, float, float] = (0, 0, 1), 
                     core_radius: int = 3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate a quantized vortex velocity field.
    
    In Space Time Potential Theory, particles are represented as quantized vortices.
    This function creates an idealized vortex with quantized circulation.
    
    Parameters
    ----------
    nx, ny, nz : int
        Grid dimensions
    center : Tuple[int, int, int]
        Location of vortex center (ix, iy, iz)
    circulation : float
        Strength of circulation, typically quantized as h/m
    axis : Tuple[float, float, float]
        Direction of vortex axis, will be normalized
    core_radius : int
        Radius of vortex core in grid cells
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Velocity field components (vx, vy, vz) of the vortex
    """
    # Initialize velocity field
    vx = np.zeros((nx, ny, nz))
    vy = np.zeros((nx, ny, nz))
    vz = np.zeros((nx, ny, nz))
    
    # Normalize axis
    ax, ay, az = axis
    norm = np.sqrt(ax**2 + ay**2 + az**2)
    if norm > 0:
        ax, ay, az = ax/norm, ay/norm, az/norm
    else:
        ax, ay, az = 0, 0, 1  # Default to z-axis
    
    # Create coordinate grid
    x, y, z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')
    
    # Center coordinates
    cx, cy, cz = center
    x = x - cx
    y = y - cy
    z = z - cz
    
    # Calculate distance from vortex axis
    # For a z-axis aligned vortex, this is simply sqrt(x² + y²)
    # For arbitrary axis, we need to find the perpendicular distance to the axis
    
    # Project position vector onto axis to find parallel component
    parallel = (x*ax + y*ay + z*az)
    
    # Parallel component vector
    px = parallel * ax
    py = parallel * ay
    pz = parallel * az
    
    # Perpendicular component vector
    rx = x - px
    ry = y - py
    rz = z - pz
    
    # Distance from axis (perpendicular)
    r = np.sqrt(rx**2 + ry**2 + rz**2)
    
    # Avoid division by zero at center
    r = np.maximum(r, 0.001)
    
    # Azimuthal unit vector (perpendicular to both r and axis)
    theta_x = ay*rz - az*ry
    theta_y = az*rx - ax*rz
    theta_z = ax*ry - ay*rx
    
    # Normalize
    theta_norm = np.sqrt(theta_x**2 + theta_y**2 + theta_z**2)
    theta_norm = np.maximum(theta_norm, 0.001)  # Avoid division by zero
    theta_x /= theta_norm
    theta_y /= theta_norm
    theta_z /= theta_norm
    
    # Apply velocity profile: v = (circulation/(2π*r)) * θ̂ for r > core_radius
    # and v = (circulation/(2π*core_radius)) * (r/core_radius) * θ̂ for r <= core_radius
    velocity_magnitude = np.where(
        r > core_radius,
        circulation / (2 * np.pi * r),
        circulation * r / (2 * np.pi * core_radius**2)
    )
    
    # Set velocity components
    vx = velocity_magnitude * theta_x
    vy = velocity_magnitude * theta_y
    vz = velocity_magnitude * theta_z
    
    return vx, vy, vz

# PROMPT: Additional operators to consider:
# 1. Helmholtz decomposition function to separate vector fields into
#    curl-free and divergence-free components
# 2. Second-order vector Laplacian (Δ²) for jerk calculations
# 3. Wave equation integrators for both first and second sound
