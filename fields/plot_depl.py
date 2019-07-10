#!/usr/local/bin/python3
import sys
import numpy as np
import matplotlib.pyplot as plt


def main():

    double = 1  # set to 0/1 to skip/add plotting of negative r values if they are not in the data file

    if (len(sys.argv) < 2):
        print("Usage: " + sys.argv[0] + " <filename>")
        exit()
    z_index = 2
    max_z = 2

    data = np.loadtxt(sys.argv[1])
    # figure out the size of grid the X-Y points
    x = set(data[:,0])  # sets only have one copy of anything, so if there are repeats, they are removed
    y = set(data[:,1])
    z = data[:,z_index]
    if ("abs" in sys.argv[2:]):
        z = np.absolute(z)  # take absolute value

    # now reshape the zvals array into the appropriate shape, and find the boundaries
    zvals = z.reshape(len(x), len(y))
    # mask 0 values to prevent from plotting
    zvals[zvals==0] = np.nan
    if (max_z != None): zvals[zvals>max_z] = np.nan

    # imshow plots columns and rows opposite to how you'd expect; so transpose them
    zvals = zvals.T

    if (double == 1 and min(x) == 0):
        # stack so we can plot the data from one half of the detector (positive r-values only)
        zvals_neg = np.fliplr(zvals)
        zvals_full = np.hstack((zvals_neg,zvals))
        bounds = (-1*max(x), max(x), min(y), max(y))
    else:
        zvals_full = zvals
        bounds = (min(x), max(x), min(y), max(y))

    # make sure that matplotlib does not have other images somehow already open
    # plt.close("all")
    # get a figure object and  the axes object associated with it
    fig = plt.figure() # (figsize=(10,5))
    # axes = plt.gca()
    # plot the image
    ip = plt.imshow(zvals_full,
                    vmax=max_z,
                    extent=bounds,   # set the boundaries of the edges of the 'image' data
                    origin="lower",  # tell matplotlib where [0,0] is in the bottom
                    cmap='jet')      # use the 'jet color map scheme, there are a bunch of options
                                     # see: https://matplotlib.org/examples/color/colormaps_reference.html
    # tell the axes what their limits are so it does not try to autoscale them
    # axes.set_xlim(bounds[:2])
    # axes.set_ylim(bounds[2:])
    # label axes
    plt.xlabel("r [mm]", size=13)
    plt.ylabel("z [mm]", labelpad=8,  size=13)
    fname = sys.argv[1][-9:-4]
    while (fname[0] == '0'):
        fname = fname[1:]
    plt.title(fname + "\n", fontsize=14, linespacing=0.4)

    plt.tight_layout()
    #plt.show()
    fname = sys.argv[1]+".png"
    print('Saving frame', fname)
    plt.savefig(fname)
    #files.append(fname)


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False


if __name__ == "__main__":
    main()
