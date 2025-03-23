"""
Example of Space Time Potential Theory simulation using the fdtd framework.

This example demonstrates:
1. Setting up an AetherGrid simulation
2. Adding quantized vortices
3. Simulating both first and second sound phenomena
4. Analyzing the electromagnetic fields that emerge
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

# Add the parent directory to the path so we can import the fdtd module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Now we can import from fdtd
from fdtd.aethergrid import AetherGrid, VISCOSITY, C_LIGHT, BACKGROUND_DENSITY

# Grid parameters
nx, ny, nz = 50, 50, 50
dx = 1.0
dt = 0.1

# Create aether grid
grid = AetherGrid(nx, ny, nz, dx, dt, VISCOSITY, BACKGROUND_DENSITY)

# Add a quantized vortex (electron-like)
h = 6.62607015e-34  # Planck constant
m_e = 9.1093837015e-31  # Electron mass
circulation = h / m_e
grid.add_quantized_vortex(
    center=(nx/2, ny/2, nz/2), 
    circulation=circulation,
    axis=(0, 0, 1),
    core_radius=5.0
)

# Add density perturbation (first sound)
grid.first_sound_solver.add_gaussian_perturbation(
    center=(nx/4, ny/2, nz/2),
    amplitude=0.1,
    width=5.0
)

# Add temperature perturbation (second sound)
grid.second_sound_solver.add_gaussian_perturbation(
    center=(3*nx/4, ny/2, nz/2),
    amplitude=0.1,
    width=5.0
)

# Run simulation for 100 steps
energies = []
for step in range(100):
    grid.update()
    
    # Record energy
    energy = grid.calculate_energy()
    energies.append(energy)
    
    # Print progress every 10 steps
    if step % 10 == 0:
        print(f"Step {step}: Total energy = {energy['total']:.6e}")

# Plot energy evolution
steps = np.arange(len(energies))
energy_types = ['kinetic', 'electromagnetic', 'first_sound', 'second_sound', 'total']
energy_values = np.array([[e[key] for key in energy_types] for e in energies])

plt.figure(figsize=(10, 6))
for i, energy_type in enumerate(energy_types):
    plt.plot(steps, energy_values[:, i], label=energy_type)
plt.xlabel('Time Step')
plt.ylabel('Energy')
plt.title('Energy Evolution in Space Time Potential Theory Simulation')
plt.legend()
plt.grid(True)
plt.savefig('energy_evolution.png')

# Create a 2D slice visualization of electromagnetic fields
plt.figure(figsize=(12, 8))

# Electric field magnitude
plt.subplot(2, 2, 1)
E_mag = np.sqrt(
    grid.electric_field[0][:, ny//2, :]**2 + 
    grid.electric_field[1][:, ny//2, :]**2 + 
    grid.electric_field[2][:, ny//2, :]**2
)
plt.imshow(E_mag.T, origin='lower')
plt.colorbar(label='Electric Field Magnitude')
plt.title('Electric Field (y-slice)')

# Magnetic field magnitude
plt.subplot(2, 2, 2)
B_mag = np.sqrt(
    grid.magnetic_field[0][:, ny//2, :]**2 + 
    grid.magnetic_field[1][:, ny//2, :]**2 + 
    grid.magnetic_field[2][:, ny//2, :]**2
)
plt.imshow(B_mag.T, origin='lower')
plt.colorbar(label='Magnetic Field Magnitude')
plt.title('Magnetic Field (y-slice)')

# First sound (density perturbation)
plt.subplot(2, 2, 3)
plt.imshow(grid.first_sound_solver.density_perturbation[:, ny//2, :].T, origin='lower')
plt.colorbar(label='Density Perturbation')
plt.title('First Sound (y-slice)')

# Second sound (temperature perturbation)
plt.subplot(2, 2, 4)
plt.imshow(grid.second_sound_solver.temperature_perturbation[:, ny//2, :].T, origin='lower')
plt.colorbar(label='Temperature Perturbation')
plt.title('Second Sound (y-slice)')

plt.tight_layout()
plt.savefig('field_visualization.png')

print("Simulation complete. Results saved to 'energy_evolution.png' and 'field_visualization.png'")
