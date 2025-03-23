"""
Vector Calculus Operators for Space Time Potential Theory
========================================================

This module provides a collection of vector calculus operators implemented
on staggered grids for the Space Time Potential Theory simulations. These
operators maintain geometric clarity and ensure proper conservation properties.

Key operators include:
- Gradient (∇φ)
- Divergence (∇·v)
- Curl (∇×v)
- Vector Laplacian (Δv = ∇(∇·v) - ∇×(∇×v))
"""

import numpy as np
from typing import Tuple, Optional
from .backend import backend

def gradient(scalar_field: np.ndarray, dx: float = 1.0, dy: float = None, dz: float = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the gradient of a scalar field.
    
    In Space Time Potential Theory, gradients are applied to:
    - Scalar potential φ = ∇·(η·v) to obtain force density
    - Temperature potential T = ∇·(η·a) to obtain yank density
    
    Parameters
    ----------
    scalar_field : np.ndarray
        3D scalar field at cell centers
    dx : float
        Grid spacing in x dimension
    dy : float, optional
        Grid spacing in y dimension, defaults to dx if None
    dz : float, optional
        Grid spacing in z dimension, defaults to dx if None
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the gradient vector field (gx, gy, gz) at cell faces
    """
    # Get current backend
    xp = backend  # Changed from backend() to backend
    
    # Use dx for any unspecified dimensions
    if dy is None:
        dy = dx
    if dz is None:
        dz = dx
    
    # Initialize gradient components (on faces)
    gx = xp.zeros((scalar_field.shape[0]+1, scalar_field.shape[1], scalar_field.shape[2]))
    gy = xp.zeros((scalar_field.shape[0], scalar_field.shape[1]+1, scalar_field.shape[2]))
    gz = xp.zeros((scalar_field.shape[0], scalar_field.shape[1], scalar_field.shape[2]+1))
    
    # x-component (on x-faces)
    gx[1:-1, :, :] = (scalar_field[1:, :, :] - scalar_field[:-1, :, :]) / dx
    # Forward/backward difference at boundaries
    gx[0, :, :] = gx[1, :, :]  # Copy nearest neighbor
    gx[-1, :, :] = gx[-2, :, :]  # Copy nearest neighbor
    
    # y-component (on y-faces)
    gy[:, 1:-1, :] = (scalar_field[:, 1:, :] - scalar_field[:, :-1, :]) / dy
    # Forward/backward difference at boundaries
    gy[:, 0, :] = gy[:, 1, :]  # Copy nearest neighbor
    gy[:, -1, :] = gy[:, -2, :]  # Copy nearest neighbor
    
    # z-component (on z-faces)
    gz[:, :, 1:-1] = (scalar_field[:, :, 1:] - scalar_field[:, :, :-1]) / dz
    # Forward/backward difference at boundaries
    gz[:, :, 0] = gz[:, :, 1]  # Copy nearest neighbor
    gz[:, :, -1] = gz[:, :, -2]  # Copy nearest neighbor
    
    return gx, gy, gz

def divergence(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float = 1.0) -> np.ndarray:
    """
    Calculate the divergence of a vector field.
    
    In Space Time Potential Theory, divergence is applied to:
    - Action Flux field C = η·v to obtain scalar potential φ
    - Radiosity field I = η·a to obtain temperature potential T
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field on cell faces
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    np.ndarray
        Divergence scalar field at cell centers
    """
    # Get current backend
    xp = backend  # Changed from backend() to backend
    
    # Extract dimensions
    nx = vx.shape[0] - 1
    ny = vy.shape[1] - 1
    nz = vz.shape[2] - 1
    
    # Initialize divergence at cell centers
    div_v = xp.zeros((nx, ny, nz))
    
    # Calculate divergence
    div_v = (vx[1:, :, :] - vx[:-1, :, :] + 
             vy[:, 1:, :] - vy[:, :-1, :] + 
             vz[:, :, 1:] - vz[:, :, :-1]) / dx
    
    return div_v

def curl_face_to_edge(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, 
                     dx: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the curl of a face-centered vector field, producing an edge-centered result.
    
    In Space Time Potential Theory, curl is applied to:
    - Velocity field v to obtain vorticity ω
    - Acceleration field a to obtain angular acceleration α
    - Action Flux field C = η·v to obtain magnetic field B
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field on cell faces
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the curl vector field on cell edges
    """
    # Get current backend
    xp = backend  # Changed from backend() to backend
    
    # Extract dimensions
    nx = vx.shape[0] - 1
    ny = vy.shape[1] - 1
    nz = vz.shape[2] - 1
    
    # Initialize curl components (on edges)
    curl_x = xp.zeros((nx, ny+1, nz+1))
    curl_y = xp.zeros((nx+1, ny, nz+1))
    curl_z = xp.zeros((nx+1, ny+1, nz))
    
    # x-component (on x-edges)
    curl_x[:, 1:-1, 1:-1] = (vz[:, 1:, 1:-1] - vz[:, :-1, 1:-1]) / dx - \
                           (vy[:, 1:-1, 1:] - vy[:, 1:-1, :-1]) / dx
    
    # y-component (on y-edges)
    curl_y[1:-1, :, 1:-1] = (vx[1:-1, :, 1:] - vx[1:-1, :, :-1]) / dx - \
                           (vz[1:, :, 1:-1] - vz[:-1, :, 1:-1]) / dx
    
    # z-component (on z-edges)
    curl_z[1:-1, 1:-1, :] = (vy[1:, 1:-1, :] - vy[:-1, 1:-1, :]) / dx - \
                           (vx[1:-1, 1:, :] - vx[1:-1, :-1, :]) / dx
    
    # Apply boundary conditions (simplified)
    # x-component boundaries
    curl_x[:, 0, :] = curl_x[:, 1, :]
    curl_x[:, -1, :] = curl_x[:, -2, :]
    curl_x[:, :, 0] = curl_x[:, :, 1]
    curl_x[:, :, -1] = curl_x[:, :, -2]
    
    # y-component boundaries
    curl_y[0, :, :] = curl_y[1, :, :]
    curl_y[-1, :, :] = curl_y[-2, :, :]
    curl_y[:, :, 0] = curl_y[:, :, 1]
    curl_y[:, :, -1] = curl_y[:, :, -2]
    
    # z-component boundaries
    curl_z[0, :, :] = curl_z[1, :, :]
    curl_z[-1, :, :] = curl_z[-2, :, :]
    curl_z[:, 0, :] = curl_z[1, :]
    curl_z[:, -1, :] = curl_z[-2, :]
    curl_z[:, :, 0] = curl_z[:, :, 1]
    curl_z[:, :, -1] = curl_z[:, :, -2]
    
    return curl_x, curl_y, curl_z

def curl_edge_to_face(wx: np.ndarray, wy: np.ndarray, wz: np.ndarray, 
                     dx: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the curl of an edge-centered vector field, producing a face-centered result.
    
    This operation is the geometric dual of curl_face_to_edge.
    
    Parameters
    ----------
    wx, wy, wz : np.ndarray
        Components of the vector field on cell edges
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the curl vector field on cell faces
    """
    # Get current backend
    xp = backend
    
    # Get dimensions from input arrays
    nx_wx, ny_wx, nz_wx = wx.shape
    nx_wy, ny_wy, nz_wy = wy.shape
    nx_wz, ny_wz, nz_wz = wz.shape
    
    # Initialize curl components (on faces)
    curl_x = xp.zeros((nx_wx+1, ny_wx, nz_wx))
    curl_y = xp.zeros((nx_wy, ny_wy+1, nz_wy))
    curl_z = xp.zeros((nx_wz, ny_wz, nz_wz+1))
    
    # Compute each component of curl carefully with proper slicing
    # x-component (on x-faces)
    for i in range(1, nx_wx):
        for j in range(1, ny_wx-1):
            for k in range(1, nz_wx-1):
                # Compute ∂wz/∂y - ∂wy/∂z with proper indexing
                curl_x[i, j, k] = (wz[i-1, j, k] - wz[i-1, j-1, k]) / dx - \
                                 (wy[i-1, j-1, k] - wy[i-1, j-1, k-1]) / dx
    
    # y-component (on y-faces)
    for i in range(1, nx_wy-1):
        for j in range(1, ny_wy):
            for k in range(1, nz_wy-1):
                # Compute ∂wx/∂z - ∂wz/∂x with proper indexing
                curl_y[i, j, k] = (wx[i-1, j-1, k] - wx[i-1, j-1, k-1]) / dx - \
                                 (wz[i, j-1, k] - wz[i-1, j-1, k]) / dx
    
    # z-component (on z-faces)
    for i in range(1, nx_wz-1):
        for j in range(1, ny_wz-1):
            for k in range(1, nz_wz):
                # Compute ∂wy/∂x - ∂wx/∂y with proper indexing
                curl_z[i, j, k] = (wy[i, j-1, k-1] - wy[i-1, j-1, k-1]) / dx - \
                                 (wx[i-1, j, k-1] - wx[i-1, j-1, k-1]) / dx
    
    # Apply boundary conditions (simplified)
    # Fill boundary values with neighboring cells
    curl_x[0, :, :] = curl_x[1, :, :]
    curl_x[-1, :, :] = curl_x[-2, :, :]
    curl_x[:, 0, :] = curl_x[:, 1, :]
    curl_x[:, -1, :] = curl_x[:, -2, :]
    curl_x[:, :, 0] = curl_x[:, :, 1]
    curl_x[:, :, -1] = curl_x[:, :, -2]
    
    curl_y[0, :, :] = curl_y[1, :, :]
    curl_y[-1, :, :] = curl_y[-2, :, :]
    curl_y[:, 0, :] = curl_y[:, 1, :]
    curl_y[:, -1, :] = curl_y[:, -2, :]
    curl_y[:, :, 0] = curl_y[:, :, 1]
    curl_y[:, :, -1] = curl_y[:, :, -2]
    
    curl_z[0, :, :] = curl_z[1, :, :]
    curl_z[-1, :, :] = curl_z[-2, :, :]
    curl_z[:, 0, :] = curl_z[1, :]
    curl_z[:, -1, :] = curl_z[-2, :]
    curl_z[:, :, 0] = curl_z[:, :, 1]
    curl_z[:, :, -1] = curl_z[:, :, -2]
    
    return curl_x, curl_y, curl_z

def vector_laplacian(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray,
                    dx: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the vector Laplacian of a vector field.
    
    The vector Laplacian is the fundamental operator in Space Time Potential Theory:
    Δv = ∇(∇·v) - ∇×(∇×v)
    
    Implementation uses direct finite differencing to ensure dimension consistency.
    
    Used to determine:
    - Acceleration: a = -η/ρ·Δv
    - Force density: f = -η·Δv
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field on cell faces
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Components of the vector Laplacian field
    """
    # Get current backend
    xp = backend
    
    # Create output arrays with the same shape as inputs
    lap_vx = xp.zeros_like(vx)
    lap_vy = xp.zeros_like(vy)
    lap_vz = xp.zeros_like(vz)
    
    # Implement a direct finite difference approximation of the vector Laplacian
    # For x-component (on x-faces)
    lap_vx[1:-1, 1:-1, 1:-1] = (
        vx[2:, 1:-1, 1:-1] + vx[:-2, 1:-1, 1:-1] +  # x-direction
        vx[1:-1, 2:, 1:-1] + vx[1:-1, :-2, 1:-1] +  # y-direction
        vx[1:-1, 1:-1, 2:] + vx[1:-1, 1:-1, :-2] -  # z-direction
        6 * vx[1:-1, 1:-1, 1:-1]  # central point
    ) / (dx * dx)
    
    # For y-component (on y-faces)
    lap_vy[1:-1, 1:-1, 1:-1] = (
        vy[2:, 1:-1, 1:-1] + vy[:-2, 1:-1, 1:-1] +  # x-direction
        vy[1:-1, 2:, 1:-1] + vy[1:-1, :-2, 1:-1] +  # y-direction
        vy[1:-1, 1:-1, 2:] + vy[1:-1, 1:-1, :-2] -  # z-direction
        6 * vy[1:-1, 1:-1, 1:-1]  # central point
    ) / (dx * dx)
    
    # For z-component (on z-faces)
    lap_vz[1:-1, 1:-1, 1:-1] = (
        vz[2:, 1:-1, 1:-1] + vz[:-2, 1:-1, 1:-1] +  # x-direction
        vz[1:-1, 2:, 1:-1] + vz[1:-1, :-2, 1:-1] +  # y-direction
        vz[1:-1, 1:-1, 2:] + vz[1:-1, 1:-1, :-2] -  # z-direction
        6 * vz[1:-1, 1:-1, 1:-1]  # central point
    ) / (dx * dx)
    
    # Apply boundary conditions by copying nearest neighbor values
    # x-component boundaries
    lap_vx[0, :, :] = lap_vx[1, :, :]
    lap_vx[-1, :, :] = lap_vx[-2, :, :]
    lap_vx[:, 0, :] = lap_vx[:, 1, :]
    lap_vx[:, -1, :] = lap_vx[:, -2, :]
    lap_vx[:, :, 0] = lap_vx[:, :, 1]
    lap_vx[:, :, -1] = lap_vx[:, :, -2]
    
    # y-component boundaries
    lap_vy[0, :, :] = lap_vy[1, :, :]
    lap_vy[-1, :, :] = lap_vy[-2, :, :]
    lap_vy[:, 0, :] = lap_vy[:, 1, :]
    lap_vy[:, -1, :] = lap_vy[:, -2, :]
    lap_vy[:, :, 0] = lap_vy[:, :, 1]
    lap_vy[:, :, -1] = lap_vy[:, :, -2]
    
    # z-component boundaries
    lap_vz[0, :, :] = lap_vz[1, :, :]
    lap_vz[-1, :, :] = lap_vz[-2, :, :]
    lap_vz[:, 0, :] = lap_vz[:, 1, :]
    lap_vz[:, -1, :] = lap_vz[-2, :]
    lap_vz[:, :, 0] = lap_vz[:, :, 1]
    lap_vz[:, :, -1] = lap_vz[:, :, -2]
    
    return lap_vx, lap_vy, lap_vz

def scalar_laplacian(scalar_field: np.ndarray, dx: float = 1.0, dy: float = None, dz: float = None) -> np.ndarray:
    """
    Calculate the scalar Laplacian of a scalar field.
    
    Used for:
    - Wave equations for density perturbations (first sound)
    - Wave equations for temperature perturbations (second sound)
    
    Parameters
    ----------
    scalar_field : np.ndarray
        3D scalar field at cell centers
    dx : float
        Grid spacing in x dimension
    dy : float, optional
        Grid spacing in y dimension, defaults to dx if None
    dz : float, optional
        Grid spacing in z dimension, defaults to dx if None
        
    Returns
    -------
    np.ndarray
        Scalar Laplacian field at cell centers
    """
    # Get current backend
    xp = backend
    
    # Use dx for any unspecified dimensions
    if dy is None:
        dy = dx
    if dz is None:
        dz = dx
    
    # Initialize Laplacian field
    lap = xp.zeros_like(scalar_field)
    
    # Calculate Laplacian using central differences
    lap[1:-1, 1:-1, 1:-1] = (
        (scalar_field[2:, 1:-1, 1:-1] + scalar_field[:-2, 1:-1, 1:-1]) / (dx * dx) +
        (scalar_field[1:-1, 2:, 1:-1] + scalar_field[1:-1, :-2, 1:-1]) / (dy * dy) +
        (scalar_field[1:-1, 1:-1, 2:] + scalar_field[1:-1, 1:-1, :-2]) / (dz * dz) -
        2 * scalar_field[1:-1, 1:-1, 1:-1] * (1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz))
    )
    
    # Apply boundary conditions (simplified)
    lap[0, :, :] = lap[1, :, :]
    lap[-1, :, :] = lap[-2, :, :]
    lap[:, 0, :] = lap[:, 1, :]
    lap[:, -1, :] = lap[:, -2, :]
    lap[:, :, 0] = lap[:, :, 1]
    lap[:, :, -1] = lap[:, :, -2]
    
    return lap

def helmholtz_decomposition(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, 
                           dx: float = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                                   np.ndarray, np.ndarray, np.ndarray]:
    """
    Decompose a vector field into irrotational (curl-free) and solenoidal (divergence-free) components.
    
    The Helmholtz decomposition is fundamental for:
    - Separating compressible (irrotational) flow from incompressible (solenoidal) flow
    - Isolating electromagnetic (solenoidal) from gravitational (irrotational) effects
    
    Parameters
    ----------
    vx, vy, vz : np.ndarray
        Components of the vector field
    dx : float
        Grid spacing in each dimension
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Components of irrotational field (vx_irrot, vy_irrot, vz_irrot) and
        solenoidal field (vx_sol, vy_sol, vz_sol)
    """
    # Get current backend
    xp = backend  # Changed from backend() to backend
    
    # Calculate divergence
    div_v = divergence(vx, vy, vz, dx)
    
    # Solve Poisson equation for scalar potential
    # Using simple relaxation method (could be optimized)
    scalar_potential = xp.zeros_like(div_v)
    for _ in range(30):  # Number of iterations
        lap_phi = scalar_laplacian(scalar_potential, dx)
        scalar_potential += 0.1 * (div_v - lap_phi)  # Relaxation factor 0.1
    
    # Calculate gradient of scalar potential to get irrotational component
    vx_irrot, vy_irrot, vz_irrot = gradient(scalar_potential, dx)
    
    # Solenoidal component is the difference
    vx_sol = vx - vx_irrot
    vy_sol = vy - vy_irrot
    vz_sol = vz - vz_irrot
    
    return vx_irrot, vy_irrot, vz_irrot, vx_sol, vy_sol, vz_sol

def quantize_circulation(field: np.ndarray, dx: float, 
                       kinematic_viscosity: float) -> np.ndarray:
    """
    Quantize circulation in a vector field to multiples of h/m.
    
    In Space Time Potential Theory, particles are represented as quantized vortices
    with circulation quantized in units of h/m.
    
    Parameters
    ----------
    field : np.ndarray
        Vector field component
    dx : float
        Grid spacing
    kinematic_viscosity : float
        Kinematic viscosity (h/m) [m²/s]
        
    Returns
    -------
    np.ndarray
        Field with quantized circulation
    """
    # Get current backend
    xp = backend  # Changed from backend() to backend
    
    # Quantum of circulation
    quantum = kinematic_viscosity
    
    # Quantize to nearest multiple of the quantum
    quantized_field = xp.round(field / quantum) * quantum
    
    return quantized_field
