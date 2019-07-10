#!/usr/local/bin/python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# plots both potential and field from a fieldgen ev.dat file

def main():

    if (len(sys.argv) < 2):
        print("Usage: " + sys.argv[0] + " <filename>")
        exit()

    # make sure that matplotlib does not have other images somehow already open
    plt.close("all")

    data = np.loadtxt(sys.argv[1])
    # figure out the size of grid the X-Y points
    x = set(data[:,0])  # sets only have one copy of anything, so if there are repeats, they are removed
    y = set(data[:,1])
    bounds = (min(x), max(x), min(y), max(y))

    # first plot the potential in subplot 2
    z_index = 4
    z = data[:,z_index]
    # now reshape the zvals array into the appropriate shape, and find the boundaries
    zvals = z.reshape(len(x), len(y))
    # mask 0 values to prevent from plotting
    zvals[zvals==0] = np.nan
    # imshow plots columns and rows opposite to how you'd expect; so transpose them
    zvals = zvals.T

    # get a figure object in subplot 1
    fig = plt.subplot(1, 2, 1)
    ip = plt.imshow(zvals,
                    extent=bounds,   # set the boundaries of the edges of the 'image' data
                    origin="lower",  # tell matplotlib where [0,0] is in the bottom
                    cmap='jet')      # use the 'jet color map scheme, there are a bunch of options
                                     # see: https://matplotlib.org/examples/color/colormaps_reference.html
    # label axes
    #plt.xlabel("Radial position [mm]", size=15)
    #plt.ylabel("Axial position [mm]",  size=15, labelpad=10)
    plt.xlabel("Radial position [mm]")
    plt.ylabel("Axial position [mm]", labelpad=10)

    # make the color legend
    cbar = plt.colorbar(ip, fraction=0.046, pad=0.04)
    label_text = ["Potential", "Field", "Radial Field", "Longitudinal Field"]
    units_text = [" [V]", " [V/cm]", " [V/cm]", " [V/cm]"]
    #cbar.set_label(label_text[z_index-2] + units_text[z_index-2], size=16, labelpad=5)
    #plt.title(label_text[z_index-2] + " vs. Position\n", fontsize=18, linespacing=0.5)
    cbar.set_label(label_text[z_index-2] + units_text[z_index-2], labelpad=5)
    plt.title(label_text[z_index-2] + " vs. Position\n", linespacing=0.5)

    # now plot the field in subplot 2
    z_index = 5
    z = data[:,z_index]
    zvals = z.reshape(len(x), len(y))
    zvals[zvals==0] = np.nan
    zvals = zvals.T
    fig = plt.subplot(1, 2, 2)
    ip = plt.imshow(zvals,
                    extent=bounds,   # set the boundaries of the edges of the 'image' data
                    origin="lower",  # tell matplotlib where [0,0] is in the bottom
                    cmap='jet')      # use the 'jet color map scheme, there are a bunch of options
                                     # see: https://matplotlib.org/examples/color/colormaps_reference.html
    # label axes
    plt.xlabel("Radial position [mm]")
    plt.ylabel("Axial position [mm]", labelpad=10)

    # make the color legend
    cbar = plt.colorbar(ip, fraction=0.046, pad=0.04)
    label_text = ["Potential", "Field", "Radial Field", "Longitudinal Field"]
    units_text = [" [V]", " [V/cm]", " [V/cm]", " [V/cm]"]
    cbar.set_label(label_text[z_index-2] + units_text[z_index-2], labelpad=5)
    plt.title(label_text[z_index-2] + " vs. Position\n", linespacing=0.5)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
