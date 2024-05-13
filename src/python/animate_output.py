"""
Module Name: animate_output.py

Description:
This module contains code to generate frames and an animation
of the height, vorticity, and u and v velocities for the output.nc
netCDF file.
"""

import os
import getpass
import warnings

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

from create_map import create_map_func

# pylint: disable=E0611
# (pylint can't find Dataset in the netCDF4 package for some reason)
from netCDF4 import Dataset as NetCDFFile

matplotlib.use("Agg")

warnings.filterwarnings("ignore")

username = getpass.getuser()

nc = NetCDFFile("../../tests/output.nc")
lons = nc.variables["phi"][:]
lats = nc.variables["theta"][:]
vort = nc.variables["vort"][:]
h = nc.variables["h"][:]
v = nc.variables["v"][:]
u = nc.variables["u"][:]
time = nc.variables["time"][:]

if not os.path.exists("../../output/frames"):
    os.mkdir("../../output/frames")
if not os.path.exists("../../output/animations"):
    os.mkdir("../../output/animations")


def make_maps(data, ax, vmin_var, vmax_var, colourbar_label_var, title_var):  # pylint: disable=R0913
    """
    Generate an individual map plot with specified colorbar settings and labels.

    Parameters:
        data : array
            The data array to be plotted.
        ax : matplotlib.axes.Axes
            The axes to draw the map on.
        vmin_var : float
            The minimum value for the colorbar.
        vmax_var : float
            The maximum value for the colorbar.
        colourbar_label_var : str
            The label for the colorbar.
        title_var : str
            The title of the plot.

    Returns:
        ax : matplotlib.axes.Axes
            The modified axes containing the map plot.
    """
    # set up orthographic map projection
    x, y, basemap = create_map_func(lons, lats)
    # contour data over the basemap
    basemap.pcolor(x, y, data, cmap="jet", shading="auto", vmin=vmin_var, vmax=vmax_var)
    cbar = basemap.colorbar(location="bottom")
    cbar.ax.set_xticks(cbar.ax.get_xticks())
    cbar.ax.set_xticklabels(
        cbar.ax.get_xticklabels(), rotation="horizontal", fontsize=8
    )

    # Format colorbar tick labels using scientific notation
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, 3))
    cbar.ax.xaxis.set_major_formatter(formatter)
    power_text = cbar.ax.xaxis.get_offset_text()
    power_text.set_size(8)

    # Set label and title
    cbar.set_label(colourbar_label_var, fontsize=8)
    ax.set_title(title_var, fontsize=8)

    return ax


ITER = 0
for it1 in range(4, len(time) + 1, 4):
    it = it1 - 1
    f = plt.figure()

    ax1 = f.add_subplot(141)
    ax1 = make_maps(h[it, :, :], ax1, 60000, 64000, "h, m", "Height")

    ax2 = f.add_subplot(142)
    ax2 = make_maps(
        vort[it, :, :], ax2, -0.00005, 0.00005, "$\\zeta$, s$^{-1}$", "Vorticity"
    )

    ax3 = f.add_subplot(143)
    ax3 = make_maps(v[it, :, :], ax3, -7, 7, "v, m s$^{-1}$", "v")

    ax4 = f.add_subplot(144)
    ax4 = make_maps(u[it, :, :], ax4, -5, 60, "u, m s$^{-1}$", "u")

    ITER += 1
    plt.suptitle(f"t={time[it]/86400:.2f} days", fontsize=8, y=0.25)
    plt.savefig(f"../../output/frames/frame{ITER:03d}.png", format="png", dpi=300)

os.system(
    "ffmpeg -r 5 -f image2  -i ../../output/frames/frame%03d.png -vframes 34 -vcodec libx264\
          -crf 25 -pix_fmt yuv420p ../../output/animations/animation.mp4"
)
os.system(
    "ffmpeg -i ../../output/animations/animation.mp4  ../../output/animations/animation.gif"
)
