"""
Module Name: max_velocity_over_time.py

Description:
This module processes data from multiple NetCDF files to display the maximum
zonal velocity (u) within the latitude band 70 to 80 degrees over time.
Each file corresponds to a different ujet value, and the results are plotted together.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

# Constants
TIME_SCALE = 2707788
FONTSIZE = 15

# List of files corresponding to different ujet values
ujet_values = [10, 15, 20]
file_paths = [f"../../tests/output_ujet_{ujet}..nc" for ujet in ujet_values]

# Define the latitude range (degrees)
LAT_MIN = 70
LAT_MAX = 80

if not os.path.exists("../../output/velocities"):
    os.mkdir("../../output/velocities")

# Initialize the plot
plt.figure(figsize=(10, 6))

# Process each file
for ujet, file_path in zip(ujet_values, file_paths):
    nc = NetCDFFile(file_path)

    lats = nc.variables["theta"][:]
    u = nc.variables["u"][:]
    time = nc.variables["time"][:] * TIME_SCALE  # Time

    # Convert latitudes from radians to degrees
    if np.max(lats) < np.pi:
        lats = np.degrees(lats)
    # Apply the mask to select latitudes within the specified range
    lat_mask = (lats >= LAT_MIN) & (lats <= LAT_MAX)
    # Extract the corresponding u velocities
    selected_u = u[:, lat_mask, :]

    # Calculate the maximum u velocity within the band for each timestep
    max_u_per_timestep = np.max(selected_u, axis=(1, 2))

    # Plotting the maximum u velocity against time
    time_days = time / 86400  # Convert time to days
    plt.plot(
        time_days, max_u_per_timestep, label=f"ujet = {ujet} ($\\times 20 \\, \\text{{ms}}^{{-1}}$)"
    )

# Configure the plot
plt.xlabel("Time (days)", fontsize=FONTSIZE)
plt.ylabel(r"Maximum Zonal Velocity ($\times 20 \, \text{ms}^{-1}$)", fontsize=FONTSIZE)
plt.title("Maximum Zonal Velocity in the 70 to 80 Degree Band Over Time", fontsize=FONTSIZE)
plt.grid(True)
plt.legend(fontsize=FONTSIZE)
plt.xticks(fontsize=FONTSIZE)
plt.yticks(fontsize=FONTSIZE)
plt.savefig("../../output/velocities/max_velocity_over_time.pdf")
plt.show()
