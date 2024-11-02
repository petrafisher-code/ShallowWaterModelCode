"""
Module Name: drift_rate.py

Description:
This module tracks a peak and plots the drift rate of the hexagon along the latitudes over time
from a single NetCDF file. One timestep equals 50 hours.

NOTE: This file is a DRAFT.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

# Constants
TIME_SCALE = 2707788
TARGET_LATITUDE = 75.9
HOURS_IN_A_DAY = 24  # Number of hours in a day
SIZE_MIN = 1e-5  # Minimum size to consider a peak significant
HOURS_PER_TIMESTEP = 2  # One timestep equals 50 hours


def load_data(file_name):
    """
    Loads velocity and time data from a NetCDF file.

    Parameters:
    file_name (str): The name of the NetCDF file.

    Returns:
    tuple: Tuple containing the velocity data at the target latitude, latitude array,
    and time array.
    """
    with NetCDFFile(file_name, "r") as nc:
        lats = nc.variables["theta"][:]
        v = nc.variables["v"][:]
        time = (
            nc.variables["time"][:] * TIME_SCALE / HOURS_PER_TIMESTEP
        )  # Convert time to timesteps

        # Find the index of the latitude closest to the target latitude
        lat_idx = np.abs(lats - TARGET_LATITUDE).argmin()

        # Extract the velocity data at the target latitude over time
        v_at_latitude = v[:, lat_idx, :]

    return v_at_latitude, lats, time


def track_peak_drift(file_name):
    """
    Tracks the peak of the hexagon and plots its drift rate over time.

    Parameters:
    file_name (str): The name of the NetCDF file to analyze.

    Returns:
    None
    """
    # Load the velocity data at the target latitude and the time array
    v_at_latitude, _, time = load_data(file_name)
    peak_positions = []

    # Track the position of the first peak at each timestep
    for v_longitude in v_at_latitude:
        # Find the peaks in the velocity data at this timestep
        peak_indices, _ = find_peaks(v_longitude, prominence=SIZE_MIN)

        # If a peak is found, save its position (longitude)
        if peak_indices.size > 0:
            peak_positions.append(peak_indices[0])

    # Convert the list of peak positions to a NumPy array
    peak_positions = np.array(peak_positions)

    # Calculate the movement of peaks between consecutive timesteps
    delta_positions = np.zeros_like(peak_positions[:-1])

    # Calculate the difference in peak positions between consecutive timesteps
    for i in range(1, len(peak_positions)):
        delta_positions[i - 1] = peak_positions[i] - peak_positions[i - 1]

    # Adjust delta positions to handle periodic boundary conditions in longitude (0 to 360 degrees)
    delta_positions[delta_positions > 180] -= 360
    delta_positions[delta_positions < -180] += 360

    # Calculate the drift rate as the change in position per timestep
    drift_rates = delta_positions / HOURS_PER_TIMESTEP

    # Adjust the time array to match the length of drift_rates and convert to days
    adjusted_time = time[1:] * HOURS_PER_TIMESTEP / HOURS_IN_A_DAY

    # Plot the drift rate over time in days
    plt.figure()
    plt.plot(adjusted_time[0:-1], drift_rates, label="Hexagon Drift Rate")
    plt.xlabel("Time (days)")
    plt.ylabel("Drift Rate (degrees/hour)")
    plt.title("Drift Rate of the Hexagon Along Latitudes Over Time")
    plt.legend()
    plt.show()

    # Print the final drift rate to 10 decimal places
    if drift_rates.size > 0:
        print(f"Final Drift Rate: {drift_rates[-1]:.10f} degrees/hour")
    else:
        print("No drift rates calculated.")


# File to analyze
FILE_NAME = "../../tests/output_default_90days.nc"

# Track and plot the drift rate of the hexagon
track_peak_drift(FILE_NAME)
