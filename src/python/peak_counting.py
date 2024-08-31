"""
Module Name: peak_counting.py

Description:
This module contains code to analyze the number
of polygon sides over time by counting peaks in
the velocity data at a specified latitude from multiple NetCDF files.
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
SECONDS_IN_A_DAY = 86400  # Number of seconds in a day
SIZE_MIN = 1e-5  # Minimum size to consider a peak significant


def load_data(file_name):
    """
    Loads velocity and time data from a NetCDF file.

    Parameters:
    file_name (str): The name of the NetCDF file.

    Returns:
    tuple: Tuple containing the velocity data and corresponding time array.
    """
    nc = NetCDFFile(file_name)
    lats = nc.variables["theta"][:]
    v = nc.variables["v"][:]
    time = nc.variables["time"][:] * TIME_SCALE / SECONDS_IN_A_DAY  # Convert time to days

    # Find the index of the latitude closest to the target latitude
    lat_idx = np.abs(lats - TARGET_LATITUDE).argmin()

    # Extract the velocity data at the target latitude over time
    v_at_latitude = v[:, lat_idx, :]

    return v_at_latitude, time


def count_peaks(v_longitude):
    """
    Counts the number of peaks in the given velocity data, considering periodic boundaries.

    Parameters:
    v_longitude (np.array): Velocity data along longitude for a specific time step.

    Returns:
    int: Number of peaks found in the velocity data.
    """
    peaks, _ = find_peaks(v_longitude, prominence=SIZE_MIN)

    # Handles periodic boundary in longitude data
    if v_longitude[-1] > v_longitude[-2] and v_longitude[-1] > v_longitude[0]:
        peaks = np.append(peaks, len(v_longitude) - 1)
    if v_longitude[0] > v_longitude[-1] and v_longitude[0] > v_longitude[1]:
        peaks = np.append(peaks, 0)
    peaks = np.unique(peaks)

    return len(peaks)


def analyze_peaks_over_time(data, time_array, label):
    """
    Analyzes the number of peaks over time and plots the results.

    Parameters:
    data (np.array): Velocity data at the target latitude over time.
    time_array (np.array): Time array converted to days.
    label (str): Label for the plot line.

    Returns:
    np.array: Array of the number of peaks over time.
    """
    peaks = []

    for t in range(data.shape[0]):
        v_longitude = data[t, :]
        num_peaks = count_peaks(v_longitude)
        peaks.append(num_peaks)
    peaks = np.array(peaks)

    # Plot the number of peaks over time
    plt.plot(time_array, peaks, label=label)


def plot_peaks_from_files(files):
    """
    Plots the number of peaks over time from multiple NetCDF files.

    Parameters:
    files (list): List of file names to process.
    """
    plt.figure()

    # Plot the horizontal dotted line behind the other lines
    plt.axhline(
        6, color="r", linestyle="--", label="Six sides", zorder=1
    )  # Lower zorder for the dotted line

    for file_name in files:
        v_at_latitude, time = load_data(file_name)
        if "ujet" in file_name:
            ujet_value = file_name.split("_ujet_")[1].split("..")[
                0
            ]  # Extract the ujet value from the filename
            label = f"ujet = {ujet_value}"
        else:
            label = "Polygon sides"
        analyze_peaks_over_time(v_at_latitude, time, label)

    plt.xlabel("Time (days)")
    plt.ylabel("Number of Peaks")
    plt.title("Number of sides in the polygon over time")
    plt.legend()
    plt.show()


# List of NetCDF files to analyze
# file_list = ["../../tests/output_ujet_10..nc", "../../tests/output_ujet_15..nc",
#               "../../tests/output_ujet_20..nc"]
file_list = ["../../tests/output.nc"]

# Analyze and plot peaks from the files
plot_peaks_from_files(file_list)
