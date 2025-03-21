import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from fdtd.lgrid import LGrid    

# Initialize the grid
size = (100, 100, 100)
dt = 0.01
grid = LGrid(size, dt)

# Initialize the velocity field in a circular pattern
x, y, z = np.indices(size)
center = np.array(size) / 2
r = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
v_theta = 1 / r
grid.C = np.stack((np.zeros_like(r), v_theta, np.zeros_like(r)), axis=-1)

# Run the simulation
steps = 10
grid.run(steps)

# Visualize the results
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(x, y, z, grid.C[..., 0], grid.C[..., 1], grid.C[..., 2])
plt.show()