"""
Space Time Potential Theory Example
==================================

This example demonstrates the key capabilities of the Space Time Potential
Theory simulation, including:
1. Quantized vortex ring (representing a particle)
2. Electromagnetic wave propagation
3. Interaction between particle and field

This showcases the dual nature of matter and field in the unified framework.
"""

import os
import sys
# Add parent directory to Python path to allow imports when running from examples directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from fdtd import PotentialGrid, C_LIGHT, VISCOSITY

def visualize_fields(grid, z_slice=None, save_path=None):
    """
    Visualize the key fields in the simulation.
    
    Parameters
    ----------
    grid : PotentialGrid
        The simulation grid
    z_slice : int, optional
        Z-slice to visualize, defaults to middle of grid
    save_path : str, optional
        If provided, save figure to this path
    """
    if z_slice is None:
        z_slice = grid.nz // 2
        
    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # Get velocity components at z-slice
    vx_slice = 0.5 * (grid.vx[:-1, :, z_slice] + grid.vx[1:, :, z_slice])
    vy_slice = 0.5 * (grid.vy[:, :-1, z_slice] + grid.vy[:, 1:, z_slice])
    
    # Calculate velocity magnitude
    v_mag = np.sqrt(vx_slice**2 + vy_slice**2)
    
    # Plot velocity field
    im1 = axes[0, 0].pcolormesh(v_mag.T, cmap='viridis')
    plt.colorbar(im1, ax=axes[0, 0], label='Velocity magnitude [m/s]')
    axes[0, 0].set_title('Velocity Field')
    axes[0, 0].set_xlabel('X')
    axes[0, 0].set_ylabel('Y')
    
    # Add velocity vectors (downsampled)
    ds = 4  # Downsample factor
    X, Y = np.meshgrid(np.arange(0, grid.nx, ds), np.arange(0, grid.ny, ds))
    axes[0, 0].quiver(X, Y, vx_slice[::ds, ::ds].T, vy_slice[::ds, ::ds].T, 
                     color='white', scale=20)
    
    # Plot vorticity (z-component)
    # Interpolate from edges to cell centers
    wz_slice = 0.25 * (grid.wz[:-1, :-1, z_slice] + grid.wz[1:, :-1, z_slice] + 
                      grid.wz[:-1, 1:, z_slice] + grid.wz[1:, 1:, z_slice])
    
    im2 = axes[0, 1].pcolormesh(wz_slice.T, cmap='RdBu_r', 
                              vmin=-0.1*np.max(np.abs(wz_slice)+1e-10), 
                              vmax=0.1*np.max(np.abs(wz_slice)+1e-10))
    plt.colorbar(im2, ax=axes[0, 1], label='Vorticity (z) [1/s]')
    axes[0, 1].set_title('Vorticity Field')
    axes[0, 1].set_xlabel('X')
    axes[0, 1].set_ylabel('Y')
    
    # Plot scalar potential
    im3 = axes[0, 2].pcolormesh(grid.scalar_potential[:, :, z_slice].T, cmap='plasma')
    plt.colorbar(im3, ax=axes[0, 2], label='Scalar Potential [Pa]')
    axes[0, 2].set_title('Scalar Potential')
    axes[0, 2].set_xlabel('X')
    axes[0, 2].set_ylabel('Y')
    
    # Plot density field
    im4 = axes[1, 0].pcolormesh(grid.density[:, :, z_slice].T, cmap='viridis',
                              vmin=grid.background_density*0.9, 
                              vmax=grid.background_density*1.1)
    plt.colorbar(im4, ax=axes[1, 0], label='Density [kg/m³]')
    axes[1, 0].set_title('Density Field')
    axes[1, 0].set_xlabel('X')
    axes[1, 0].set_ylabel('Y')
    
    # Plot electric field
    # First get the components
    ex, ey, ez = grid.electric_field
    
    # Interpolate to cell centers for visualization
    ex_interp = 0.5 * (ex[:-1, :, z_slice] + ex[1:, :, z_slice])
    ey_interp = 0.5 * (ey[:, :-1, z_slice] + ey[:, 1:, z_slice])
    
    # Now calculate magnitude using the interpolated fields
    e_mag = np.sqrt(ex_interp**2 + ey_interp**2)
    
    im5 = axes[1, 1].pcolormesh(e_mag.T, cmap='inferno')
    plt.colorbar(im5, ax=axes[1, 1], label='Electric Field [V/m]')
    axes[1, 1].set_title('Electric Field Magnitude')
    axes[1, 1].set_xlabel('X')
    axes[1, 1].set_ylabel('Y')
    
    # Plot magnetic field
    bx, by, bz = grid.magnetic_field
    
    # Ensure correct dimensions for plotting
    if bz.shape[0] == grid.nx and bz.shape[1] == grid.ny:
        bz_plot = bz[:, :, z_slice]
    else:
        # Additional interpolation might be needed
        bz_plot = 0.25 * (bz[:-1, :-1, z_slice] + bz[1:, :-1, z_slice] +
                         bz[:-1, 1:, z_slice] + bz[1:, 1:, z_slice])
    
    im6 = axes[1, 2].pcolormesh(bz_plot.T, cmap='RdBu_r',
                              vmin=-np.max(np.abs(bz_plot)+1e-10), 
                              vmax=np.max(np.abs(bz_plot)+1e-10))
    plt.colorbar(im6, ax=axes[1, 2], label='Magnetic Field B_z [T]')
    axes[1, 2].set_title('Magnetic Field (z)')
    axes[1, 2].set_xlabel('X')
    axes[1, 2].set_ylabel('Y')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        plt.close(fig)
    else:
        plt.show()
    
    return fig

def main():
    # Create simulation grid
    nx, ny, nz = 64, 64, 16
    dx = 0.1  # Grid spacing in meters
    
    print("Initializing PotentialGrid simulation...")
    grid = PotentialGrid(nx, ny, nz, dx=dx)
    
    # Add a quantized vortex ring (representing a particle)
    grid.add_vortex_ring(
        center=(nx//2 * dx, ny//2 * dx, nz//2 * dx),
        radius=10 * dx,
        quantized=True
    )
    
    # Add a plane wave (electromagnetic wave)
    grid.add_plane_wave(
        amplitude=0.01,  # Small amplitude to avoid nonlinear effects
        frequency=1e9,   # 1 GHz
        direction=(1, 0, 0),
        polarization=(0, 1, 0)
    )
    
    # Visualize initial state
    print("Visualizing initial fields...")
    visualize_fields(grid, save_path="potential_theory_initial.png")
    
    # Run simulation
    print("Running simulation...")
    frames = []
    
    def save_frame(grid):
        fig = visualize_fields(grid)
        frames.append(fig)
    
    # Run for 100 steps
    grid.run(100, callback=save_frame, interval=10)
    
    # Visualize final state
    print("Visualizing final fields...")
    visualize_fields(grid, save_path="potential_theory_final.png")
    
    # Calculate and display energy components
    energy = grid.calculate_energy()
    print("\nEnergy components:")
    print(f"  Kinetic energy:       {energy['kinetic']:.5e} J")
    print(f"  Electric field energy: {energy['electric']:.5e} J")
    print(f"  Magnetic field energy: {energy['magnetic']:.5e} J")
    print(f"  First sound energy:   {energy['first_sound']:.5e} J")
    print(f"  Second sound energy:  {energy['second_sound']:.5e} J")
    print(f"  Total energy:         {energy['total']:.5e} J")
    
    print("\nSimulation completed successfully.")

if __name__ == "__main__":
    main()
