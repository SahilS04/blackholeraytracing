import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Step 1: Read the CSV file (assuming x, y, and intensity columns)
data = pd.read_csv('image.csv', header=None)  # Replace 'data.csv' with your file path

# Step 2: Extract x, y, and intensity values
x = data[0].values
y = data[1].values
intensity = data[2].values

# Step 3: Define the grid resolution (number of grid points in x and y direction)
grid_resolution_x = 1500  # Adjust as needed
grid_resolution_y = 1500  # Adjust as needed

# Step 4: Create a grid of points where you want to interpolate the data
x_grid = np.linspace(min(x), max(x), grid_resolution_x)
y_grid = np.linspace(min(y), max(y), grid_resolution_y)
x_grid, y_grid = np.meshgrid(x_grid, y_grid)

# Step 5: Interpolate the intensity data onto the grid
# Use griddata to interpolate non-integer (x, y) points onto a regular grid
grid_intensity = griddata((x, y), intensity, (x_grid, y_grid), method='linear')

# Step 6: Plot the intensity heatmap
plt.imshow(grid_intensity, extent=(min(x), max(x), min(y), max(y)), origin='lower', cmap='hot', aspect='auto')

# Add a color bar to show the intensity scale
#plt.colorbar(label='Intensity')

# Optionally, add labels and title
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Black Hole Image via Ray Tracing')

# Show the plot
plt.show()
