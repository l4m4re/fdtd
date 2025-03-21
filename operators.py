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

# PROMPT: Additional operators to consider:
# 1. Helmholtz decomposition function to separate vector fields into
#    curl-free and divergence-free components
# 2. Second-order vector Laplacian (Δ²) for jerk calculations
# 3. Wave equation integrators for both first and second sound
