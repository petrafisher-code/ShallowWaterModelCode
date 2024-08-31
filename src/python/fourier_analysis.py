"""
Module Name: fourier_analysis.py

Description:
This module contains code to analyse the number
of polygon sides over time with Fourier analysis.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftshift
from scipy.signal import find_peaks

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

# TODO: FIX AND TEST THIS MODULE pylint: disable=W0511

TIME_SCALE = 2707788

nc = NetCDFFile("../../tests/output.nc")
lons = nc.variables["phi"][:]
lats = nc.variables["theta"][:]
vort = nc.variables["vort"][:]
h = nc.variables["h"][:]
v = nc.variables["v"][:]
u = nc.variables["u"][:]
time = nc.variables["time"][:] * TIME_SCALE

# Function to compute the radial magnitude
def compute_radial_magnitude(magnitude_spectrum, r_bins, r):
    """
    Compute the radial magnitude by averaging the magnitude spectrum along each radius.

    Parameters:
        magnitude_spectrum : array
            The magnitude spectrum to analyze.
        r_bins : array
            The bins for radial distances.
        r : array
            The radial distance array.

    Returns:
        radial_magnitude : array
            The averaged magnitude along each radius.
    """
    radial_magnitude = np.zeros_like(r_bins)

    for i in range(len(r_bins) - 1):
        mask = (r >= r_bins[i]) & (r < r_bins[i + 1])
        radial_magnitude[i] = magnitude_spectrum[mask].mean()

    return radial_magnitude

# Function to compute the dominant mode in the data using Fourier analysis
def dominant_mode(data):
    """
    Compute the dominant mode (number of polygon sides) in the given data.

    Parameters:
        data : array
            The 2D data array to analyze.

    Returns:
        mode : int
            The number of polygon sides corresponding to the strongest mode.
    """
    data_fft = fft2(data)
    data_fft_shifted = fftshift(data_fft)
    magnitude_spectrum = np.abs(data_fft_shifted)

    height, width = magnitude_spectrum.shape
    center = (int(height / 2), int(width / 2))
    y, x = np.indices((height, width))
    r = np.sqrt((x - center[1]) ** 2 + (y - center[0]) ** 2)
    r_bins = np.linspace(0, np.max(r), height)

    radial_magnitude = compute_radial_magnitude(magnitude_spectrum, r_bins, r)

    peaks, _ = find_peaks(radial_magnitude)
    strongest_peak = peaks[np.argmax(radial_magnitude[peaks])]

    mode = strongest_peak * 2
    print(mode)

    return mode

# Function to analyze the strongest mode over time
def analyze_modes_over_time(data, time_array):
    """
    Analyze and plot the strongest mode in the data over time.

    Parameters:
        data : array
            The 3D data array (time, lat, lon) to analyze.
        time_array : array
            The array of time steps corresponding to the data.
    """
    modes = []

    for t in range(len(time_array)):
        mode = dominant_mode(data[t, :, :])
        modes.append(mode)

    plt.figure(figsize=(10, 5))
    plt.plot(time_array, modes, color="purple")
    plt.xlabel("Time (days)")
    plt.ylabel("Strongest Mode (Number of Polygon Sides)")
    plt.title("Strongest Mode Over Time")
    plt.grid(True)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.show()

# Example usage with velocity at all time steps
time_steps_in_days = time / 86400  # Convert time to days
analyze_modes_over_time(v, time_steps_in_days)
