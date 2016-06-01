#!/usr/bin/env python -W ignore::DeprecationWarning
# Will Yingling
# April 2016


from mpl_toolkits.mplot3d import Axes3D  # for 3D
import numpy as np
import matplotlib.pyplot as plt  # for 2D
from matplotlib.colors import LogNorm  # for colors and scale
import densitysmoothing as ds
import matplotlib.animation as animation
import warnings
import matplotlib.mlab as mlab


# this reads in the COM files we made with Assigment 4
def read_file(filename):
    # skip header is only 0 because
    # it's the file I made with only 1 header line
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=0)

    diff = data['diff']


    indexID = np.where(diff < 1000.)

    diff = data['diff'][indexID]
    # store y position of particles
    fronts = data['f'][indexID]
    # store z position of particles
    backs = data['b'][indexID]



    return diff, fronts, backs



#################
# main function #
#################
def main():

    print("Starting...")

    Histo_file = [ "Stdev_Histos.txt"]

    diff, front_part, back_part = read_file(Histo_file[0])

    snap_time = np.arange(0., 700., 1)



    mean_stdev = np.mean(diff)

    mean_temp = np.ones(len(snap_time))

    mean_stdev_arr = mean_temp * mean_stdev

    mean_label_string = "Mean Standard Deviation = " + str(round(mean_stdev,4))

    line_stdev, = plt.plot(snap_time, diff, label='Standard Deviation Difference')
    line_mean, = plt.plot(snap_time, mean_stdev_arr, '--', color='red',  label=mean_label_string)
    plt.title("Difference in Sigma Across All Snaps")
    plt.xlabel("Snapshot Number")
    plt.ylabel("Difference between Stdev of front and back particles")
    plt.legend(handles=[line_stdev, line_mean], loc=2)
    plt.savefig("Diff_Stdev.png")
    plt.show()


    print('end')



if __name__ == '__main__':
    main()
