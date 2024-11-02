"""
Module Name: end_v_val.py

Description:
This module processes meridional velocity data at a target latitude (75.8 degrees) from
a NetCDF file. It extracts min and max velocity values at the first, 5-day, and last
timesteps and plots velocity profiles across longitudes to analyze changes over time.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

# Use "DejaVu Serif" for LaTeX-like font
matplotlib.rcParams["font.family"] = "DejaVu Serif"
matplotlib.rcParams["font.serif"] = "DejaVu Serif"

# Set constants
FONTSIZE = 20
CASSINI_PERSPECTIVE = False
TIME_SCALE = 2707788

# Load the data from the NetCDF file
nc = NetCDFFile("../../tests/output_default_90days.nc")
lons = nc.variables["phi"][:]  # Longitudes
lats = nc.variables["theta"][:]  # Latitudes in radians
v = nc.variables["v"][:]  # v velocity component
time = nc.variables["time"][:] * TIME_SCALE

# Convert 75.9 degrees to radians
target_latitude_radians = np.radians(75.8)

# Find the index corresponding to latitude 75.8 degrees in radians
latitude_index = np.argmin(np.abs(lats - target_latitude_radians))

# Print the selected latitude in radians to verify
selected_latitude = lats[latitude_index]
print(
    f"Selected latitude is {selected_latitude:.4f} "
    f"(should be close to {target_latitude_radians:.4f} radians)"
)

# Calculate the timestep index for 5 days (5 days * 2 timesteps per day)
TIMESTEP_5_DAYS = 5 * 2

# Extract the v data for the first timestep, 5-day timestep, and last timestep
v_first_timestep = v[0, latitude_index, :]
v_5_day_timestep = v[TIMESTEP_5_DAYS, latitude_index, :]
v_last_timestep = v[-1, latitude_index, :]

# Print the min and max of v for each timestep to check for negative values
# First timestep
min_v_first, max_v_first = np.min(v_first_timestep), np.max(v_first_timestep)
print(f"First timestep (75.8°): Min v = {min_v_first}, Max v = {max_v_first}")

# 5-day timestep
min_v_5_day, max_v_5_day = np.min(v_5_day_timestep), np.max(v_5_day_timestep)
print(f"5-day timestep (75.8°): Min v = {min_v_5_day}, Max v = {max_v_5_day}")

# Last timestep
min_v_last, max_v_last = np.min(v_last_timestep), np.max(v_last_timestep)
print(f"Last timestep (75.8°): Min v = {min_v_last}, Max v = {max_v_last}")

# Define custom colors with a lighter green
colors = ["#66c2a5", "blue", "red"]  # '#66c2a5' is a pleasant light green

# Plotting the v data for all longitudes at the specified latitude for the selected timesteps
plt.figure(figsize=(12, 6))  # Increased width for better x-axis label fitting
plt.plot(lons, v_first_timestep, label="First Timestep (0 days)", color=colors[0])
plt.plot(lons, v_5_day_timestep, label="5-Day Timestep", color=colors[1])
plt.plot(
    lons, v_last_timestep, label=f"Last Timestep ({time[-1] / 86400:.1f} days)", color=colors[2]
)

plt.title(
    r"Meridional Velocity Component at Latitude 75.8° for Various Timesteps", fontsize=FONTSIZE
)
plt.xlabel("Longitude", fontsize=FONTSIZE)
plt.ylabel(r"$v \ (\times 20 \, \text{m/s})$", fontsize=FONTSIZE)
plt.grid(True)

# Set x-axis major ticks to display in degrees for better readability
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(np.pi / 6))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{np.degrees(x):.0f}°"))

# Set font size for tick labels
ax.tick_params(axis="both", which="major", labelsize=FONTSIZE)

plt.legend(fontsize=FONTSIZE)
plt.tight_layout()

# Save the plot as a PDF
OUTPUT_FILENAME = "../../output/velocities/v_velocity.pdf"
plt.savefig(OUTPUT_FILENAME)
print(f"Plot saved as {OUTPUT_FILENAME}")
