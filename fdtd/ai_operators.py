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
from typing import Tuple, Optional
from .backend import backend, get_backend_name

# PROMPT: Consider implementing advanced staggered grid techniques 
# for better representation of vector Laplacian components

def grad(field: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the gradient of a scalar field.
    
    In Space Time Potential Theory, gradients are applied to:
    - Scalar potential φ = ∇·(η·v) to obtain force density
    - Temperature potential T = ∇·(η·a) to obtain yank density
    
    Parameters
    ----------
    field : np.ndarray
        3D scalar field
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the gradient vector field (grad_x, grad_y, grad_z)
    """
    # PROMPT: Evaluate boundary conditions for gradient operation
    # to maintain consistency with substrate physics
    
    grad_x = np.zeros_like(field)
    grad_y = np.zeros_like(field)
    grad_z = np.zeros_like(field)
    
    # x-component: central difference
    grad_x[1:-1, :, :] = (field[2:, :, :] - field[:-2, :, :]) / (2 * dx)
    # Forward/backward difference at boundaries
    grad_x[0, :, :] = (field[1, :, :] - field[0, :, :]) / dx
    grad_x[-1, :, :] = (field[-1, :, :] - field[-2, :, :]) / dx
    
    # y-component: central difference
    grad_y[:, 1:-1, :] = (field[:, 2:, :] - field[:, :-2, :]) / (2 * dx)
    # Forward/backward difference at boundaries
    grad_y[:, 0, :] = (field[:, 1, :] - field[:, 0, :]) / dx
    grad_y[:, -1, :] = (field[:, -1, :] - field[:, -2, :]) / dx
    
    # z-component: central difference
    grad_z[:, :, 1:-1] = (field[:, :, 2:] - field[:, :, :-2]) / (2 * dx)
    # Forward/backward difference at boundaries
    grad_z[:, :, 0] = (field[:, :, 1] - field[:, :, 0]) / dx
    grad_z[:, :, -1] = (field[:, :, -1] - field[:, :, -2]) / dx
    
    return grad_x, grad_y, grad_z

def div(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float) -> np.ndarray:
    """
    Calculate the divergence of a vector field.
    
    In Space Time Potential Theory, divergence is applied to:
    - Action Flux field C = η·v to obtain scalar potential φ
    - Radiosity field I = η·a to obtain temperature potential T
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    np.ndarray
        Divergence scalar field
    """
    # PROMPT: Consider compressible vs incompressible treatment
    # for different substrate regimes
    
    div_v = np.zeros_like(vx)
    
    # Central difference for each component's derivative
    div_v[1:-1, :, :] += (vx[2:, :, :] - vx[:-2, :, :]) / (2 * dx)
    div_v[:, 1:-1, :] += (vy[:, 2:, :] - vy[:, :-2, :]) / (2 * dx)
    div_v[:, :, 1:-1] += (vz[:, :, 2:] - vz[:, :, :-2]) / (2 * dx)
    
    # Forward/backward difference at boundaries (x)
    div_v[0, :, :] += (vx[1, :, :] - vx[0, :, :]) / dx
    div_v[-1, :, :] += (vx[-1, :, :] - vx[-2, :, :]) / dx
    
    # Forward/backward difference at boundaries (y)
    div_v[:, 0, :] += (vy[:, 1, :] - vy[:, 0, :]) / dx
    div_v[:, -1, :] += (vy[:, -1, :] - vy[:, -2, :]) / dx
    
    # Forward/backward difference at boundaries (z)
    div_v[:, :, 0] += (vz[:, :, 1] - vz[:, :, 0]) / dx
    div_v[:, :, -1] += (vz[:, :, -1] - vz[:, :, -2]) / dx
    
    return div_v

def curl_face_to_edge(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, 
                     dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the curl vector field (curl_x, curl_y, curl_z)
    """
    # PROMPT: Implement quantized vorticity constraints
    # to model topological defects in substrate
    
    curl_x = np.zeros_like(vx)
    curl_y = np.zeros_like(vy)
    curl_z = np.zeros_like(vz)
    
    # x-component: dw/dy - dv/dz
    curl_y[:, 1:-1, :] += (vz[:, 2:, :] - vz[:, :-2, :]) / (2 * dx)
    curl_y[:, :, 1:-1] -= (vy[:, :, 2:] - vy[:, :, :-2]) / (2 * dx)
    
    # y-component: du/dz - dw/dx
    curl_y[:, :, 1:-1] += (vx[:, :, 2:] - vx[:, :, :-2]) / (2 * dx)
    curl_y[1:-1, :, :] -= (vz[2:, :, :] - vz[:-2, :, :]) / (2 * dx)
    
    # z-component: dv/dx - du/dy
    curl_z[1:-1, :, :] += (vy[2:, :, :] - vy[:-2, :, :]) / (2 * dx)
    curl_z[:, 1:-1, :] -= (vx[:, 2:, :] - vx[:, :-2, :]) / (2 * dx)
    
    # Boundary conditions (simplified version)
    # PROMPT: Implement proper boundary conditions for curl
    # consistent with substrate physics
    
    return curl_x, curl_y, curl_z

def curl_edge_to_face(wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, 
                     dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the curl of a vector field.
    
    In Space Time Potential Theory, curl is applied to:
    - Velocity field v to obtain vorticity ω
    - Acceleration field a to obtain angular acceleration α
    - Action Flux field C = η·v to obtain magnetic field B
    
    Parameters
    ----------
    wx, wy, wz : np.ndarray
        Components of the vector field
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the curl vector field (curl_x, curl_y, curl_z)
    """
    # PROMPT: Implement quantized vorticity constraints
    # to model topological defects in substrate
    
    curl_x = np.zeros_like(wx)
    curl_y = np.zeros_like(wy)
    curl_z = np.zeros_like(wz)
    
    # x-component: dw/dy - dv/dz
    curl_y[:, 1:-1, :] += (wz[:, 2:, :] - wz[:, :-2, :]) / (2 * dx)
    curl_y[:, :, 1:-1] -= (wy[:, :, 2:] - wy[:, :, :-2]) / (2 * dx)
    
    # y-component: du/dz - dw/dx
    curl_y[:, :, 1:-1] += (wx[:, :, 2:] - wx[:, :, :-2]) / (2 * dx)
    curl_y[1:-1, :, :] -= (wz[2:, :, :] - wz[:-2, :, :]) / (2 * dx)
    
    # z-component: dv/dx - du/dy
    curl_z[1:-1, :, :] += (wy[2:, :, :] - wy[:-2, :, :]) / (2 * dx)
    curl_z[:, 1:-1, :] -= (wx[:, 2:, :] - wx[:, :-2, :]) / (2 * dx)
    
    # Boundary conditions (simplified version)
    # PROMPT: Implement proper boundary conditions for curl
    # consistent with substrate physics
    
    return curl_x, curl_y, curl_z

def vector_laplacian(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray,
                    dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    dx : float
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
    lap_x = np.zeros_like(vx)
    lap_y = np.zeros_like(vy)
    lap_z = np.zeros_like(vz)
    
    # x-component
    lap_x[1:-1, :, :] += (vx[2:, :, :] - 2*vx[1:-1, :, :] + vx[:-2, :, :]) / (dx**2)
    lap_x[:, 1:-1, :] += (vx[:, 2:, :] - 2*vx[:, 1:-1, :] + vx[:, :-2, :]) / (dx**2)
    lap_x[:, :, 1:-1] += (vx[:, :, 2:] - 2*vx[:, :, 1:-1] + vx[:, :, :-2]) / (dx**2)
    
    # y-component
    lap_y[1:-1, :, :] += (vy[2:, :, :] - 2*vy[1:-1, :, :] + vy[:-2, :, :]) / (dx**2)
    lap_y[:, 1:-1, :] += (vy[:, 2:, :] - 2*vy[:, 1:-1, :] + vy[:, :-2, :]) / (dx**2)
    lap_y[:, :, 1:-1] += (vy[:, :, 2:] - 2*vy[:, :, 1:-1] + vy[:, :, :-2]) / (dx**2)
    
    # z-component
    lap_z[1:-1, :, :] += (vz[2:, :, :] - 2*vz[1:-1, :, :] + vz[:-2, :, :]) / (dx**2)
    lap_z[:, 1:-1, :] += (vz[:, 2:, :] - 2*vz[:, 1:-1, :] + vz[:, :-2, :]) / (dx**2)
    lap_z[:, :, 1:-1] += (vz[:, :, 2:] - 2*vz[:, :, 1:-1] + vz[:, :, :-2]) / (dx**2)
    
    # PROMPT: Method 2 could be implemented as:
    # div = div(vx, vy, vz, dx)
    # grad_div_x, grad_div_y, grad_div_z = grad(div, dx)
    # curl_x, curl_y, curl_z = curl_face_to_edge(vx, vy, vz, dx)
    # curl_curl_x, curl_curl_y, curl_curl_z = curl_edge_to_face(curl_x, curl_y, curl_z, dx)
    # lap_x = grad_div_x - curl_curl_x
    # lap_y = grad_div_y - curl_curl_y
    # lap_z = grad_div_z - curl_curl_z
    
    return lap_x, lap_y, lap_z

def quantize_circulation(field: np.ndarray, dx: float, 
                        kinematic_viscosity: float) -> np.ndarray:
    """
    Quantize the circulation of a vector field.
    
    In Space Time Potential Theory, particles are represented as quantized vortices.
    This function quantizes the circulation of the field.
    
    Parameters
    ----------
    field : np.ndarray
        Vector field
    dx : float
        Grid spacing
    kinematic_viscosity : float
        Kinematic viscosity of the substrate
        
    Returns
    -------
    np.ndarray
        Quantized field
    """
    # Initialize quantized field
    quantized_field = np.zeros_like(field)
    
    # Quantize circulation
    quantized_field = field / (2 * np.pi * kinematic_viscosity * dx)
    
    return quantized_field
