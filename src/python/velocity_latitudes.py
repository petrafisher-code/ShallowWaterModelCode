"""
Module Name: velocity_latitudes.py

Description:
This module processes the data at the last time-step of a
simulation, finding and plotting the average zonal velocity
at each latitude.
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

# Load the NetCDF file
nc = NetCDFFile("../../tests/output.nc")

lats = nc.variables["theta"][:]
u = nc.variables["u"][:]

# Define the latitude range (65 to 86.5 degrees)
LAT_MIN = 65
LAT_MAX = round(lats[-1]*180/np.pi, 1)

if not os.path.exists("../../output/velocities"):
    os.mkdir("../../output/velocities")

# Check the range of latitudes available in the data
print(f"Latitude data range: {np.min(lats)*180/np.pi} to {np.max(lats)*180/np.pi} degrees")

# Convert latitudes from radians to degrees
if np.max(lats) < np.pi:
    lats = np.degrees(lats)

# Apply the mask again after converting to degrees
lat_mask = (lats >= LAT_MIN) & (lats <= LAT_MAX)
selected_lats = lats[lat_mask]

if selected_lats.size == 0:
    print("No latitudes found in the specified range.")
else:
    # Extract the corresponding u velocities
    selected_u = u[-1, lat_mask, :]

    # Calculate the average u velocity at each latitude
    average_u = np.mean(selected_u, axis=1)

    # Plotting the average u velocity against latitude
    plt.figure(figsize=(8, 6))
    plt.plot(average_u, selected_lats, color="purple", marker="o")
    plt.ylabel("Latitude (Degrees)", fontsize=FONTSIZE)
    plt.xlabel(r"Average Zonal Velocity ($\times 20 \, \text{ms}^{-1}$)", fontsize=FONTSIZE)
    plt.title(f"Average Zonal Velocities at Each Latitude (65 to ${LAT_MAX}$ degrees)",\
               fontsize=FONTSIZE)
    plt.grid(True)

    # Set tick parameters
    plt.xticks(fontsize=FONTSIZE)
    plt.yticks(fontsize=FONTSIZE)

    plt.savefig("../../output/velocities/velocity_latitudes.png")
    plt.show()
