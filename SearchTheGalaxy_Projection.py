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


with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)


# setting the global variable for gravitational constant
def set_big_g():
    global BIG_G
    BIG_G = 4.498768e-6  # in kpc^3/M_sun/Gyr

# setting the global variable for gravitational constant
def set_RADIUS():
    global RADIUS
    RADIUS = 10.  # in kpc

def set_image_data(all_data):
    global COUNT
    COUNT = 0

    global ImData
    ImData = all_data

"""
def find_image_data():
    current_plot_data = ImData[COUNT]
    old_count = COUNT
    global COUNT
    COUNT = old_count + 1
    return current_plot_data
"""

# this reads in the COM files we made with Assigment 4
def read_COM(filename):
    # skip header is only 0 because
    # it's the file I made with only 1 header line
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=0)

    time = data['time']

    new_x = data['x']
    # store y position of particles
    new_y = data['y']
    # store z position of particles
    new_z = data['z']
    # store x position of particles
    new_vx = data['vx']
    # store y position of particles
    new_vy = data['vy']
    # store z position of particles
    new_vz = data['vz']

    return time, new_x, new_y, new_z, new_vx, new_vy, new_vz


# this reads in the COM files we made with Assigment 4
def read_VHighRes(filename, COMfile):
    # skip header is only 0 because
    # it's the file I made with only 1 header line

    PType = 1.

    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)

    time_COM, x_COM, y_COM, z_COM, vx_COM, vy_COM, vz_COM = read_COM(COMfile)

    x_COM_M33 = x_COM[700]
    y_COM_M33 = y_COM[700]
    z_COM_M33 = z_COM[700]

    indexID = np.where(data['type'] == PType)

    Px = data['x'][indexID]
    # store y position of particles
    Py = data['y'][indexID]
    # store z position of particles
    Pz = data['z'][indexID]
    # store x position of particles
    Pvx = data['vx'][indexID]
    # store y position of particles
    Pvy = data['vy'][indexID]
    # store z position of particles
    Pvz = data['vz'][indexID]

    #TODO make this a for loop over all time.
    new_x = Px #- x_COM[700]
    new_y = Py #- y_COM[700]
    new_z = Pz #- z_COM[700]
    new_vx = Pvx #- vx_COM[700]
    new_vy = Pvy #- vy_COM[700]
    new_vz = Pvz #- vz_COM[700]


    """
    if PType == 1.:
        M31_file = "Halo_of_M31.txt"
    elif PType == 2.:
        M31_file = "Disk_of_M31.txt"
    elif PType == 3.:
        M31_file = "Bulge_of_M31.txt"

    open_M31_Bulge = open(M31_file, '+w')

    open_M31_Bulge.write("x y z vx vy vz xCOM yCOM zCOM")
    open_M31_Bulge.write("\n")

    for i in range(len(new_x)):
        # write info with COM at snap 000.
        COMP = [new_x[i], new_y[i], new_z[i], new_vx[i], new_vy[i], new_vz[i], x_COM[0], y_COM[0], z_COM[0]]
        # converting COMP to a string
        COMP_to_str = " ".join(str(j) for j in COMP)
        open_M31_Bulge.write(COMP_to_str)
        open_M31_Bulge.write("\n")

    open_M31_Bulge.close()
    """

    return new_x, new_y, new_z, new_vx, new_vy, new_vz  # , x_COM, y_COM, z_COM


# calc magnitude
# xx, yy, zz are arrays of the same size
def find_wake(xx, yy, zz, vxx, vyy, vzz, x_COM, y_COM, z_COM, radius):
    local_x = []
    local_y = []
    local_z = []

    local_vxx = []
    local_vyy = []
    local_vzz = []

    upper_bound_x = x_COM + radius
    lower_bound_x = x_COM - radius

    upper_bound_y = y_COM + radius
    lower_bound_y = y_COM - radius

    upper_bound_z = z_COM + radius
    lower_bound_z = z_COM - radius


    for i in range(len(xx)):

        m = np.sqrt((x_COM-xx[i])**2 + (y_COM-yy[i])**2 + (z_COM-zz[i])**2)

        #if xx[i] <= upper_bound_x and xx[i] >= lower_bound_x:
        if m <= radius:
            x = xx[i] - x_COM
            local_x.append(x)
            local_vxx.append(vxx[i])

        # if yy[i] <= upper_bound_y and yy[i] >= lower_bound_y:
            y = yy[i] - y_COM
            local_y.append(y)
            local_vyy.append(vyy[i])

        # if zz[i] <= upper_bound_z and zz[i] >= lower_bound_z:
            z = zz[i] - z_COM
            local_z.append(z)
            local_vzz.append(vzz[i])

    return local_x, local_y, local_z, local_vxx, local_vyy, local_vzz


# return and array of the hernquist profile
# a is buffer parameter
def hernquist(mass, radius, a):
    # initialize array
    hern_mass = []

    # loop through each radius step
    for i in range(len(radius)):
        # derived H mass profile
        profile = mass * radius[i]**2 / (radius[i] + a)**2
        hern_mass.append(profile)

    return hern_mass


# calculating the circular speed with an array of masses over each radii step
def circular_speed_hern(mass, scale_length, radius):
    v_circ = []

    for i in range(len(radius)):
        v_circ_at_r = np.sqrt(BIG_G * mass / (radius[i] + scale_length)**2 * radius[i])
        v_circ.append(v_circ_at_r)

    return v_circ


def circular_speed_MN(mass, scale_radius, radius):
    v_circ = []
    scale_height = scale_radius / 5.

    for i in range(len(radius)):
        v_circ_at_r = np.sqrt(BIG_G * mass / (radius[i]**2 + (scale_radius + scale_height)**2)**1.5 * radius[i]**2)
        v_circ.append(v_circ_at_r)

    return v_circ


# calc the combined rotation curve
# with arrays of each component and radii
def rotation_curve(v_circ_halo, v_circ_disk, v_circ_bulge):
    rot_curve = []

    for i in range(len(v_circ_halo)):
        rot_curve_at_r = np.sqrt(v_circ_halo[i]**2 + v_circ_disk[i]**2 + v_circ_bulge[i]**2)
        rot_curve.append(rot_curve_at_r)

    return rot_curve


def MagnitudeVector(a, b, c, d, e, f):

    mag_vec = np.sqrt((a-d)**2 + (b-e)**2 + (c-f)**2)

    return mag_vec


def find_path(snap, x_COM, y_COM, z_COM):

    #initialize with the starting position
    x_track = [0]
    y_track = [0]
    z_track = [0]

    for i in range(7):
        # stepping back in time

        # find dist between current pos and prev ten places
        x_prev = -x_COM[snap] + x_COM[snap-i]
        y_prev = -y_COM[snap] + y_COM[snap-i]
        z_prev = -z_COM[snap] + z_COM[snap-i]

        x_track.append(x_prev)
        y_track.append(y_prev)
        z_track.append(z_prev)

    return x_track, y_track, z_track


def find_M31(snap, x_COM_M31, z_COM_M31, y_COM_M31, x_COM_M33, y_COM_M33, z_COM_M33):

    x_current_M31 = x_COM_M31[snap]
    y_current_M31 = y_COM_M31[snap]
    z_current_M31 = z_COM_M31[snap]
    x_current_M33 = x_COM_M33[snap]
    y_current_M33 = y_COM_M33[snap]
    z_current_M33 = z_COM_M33[snap]

    mag = MagnitudeVector(x_current_M31,  y_current_M31, z_current_M31, x_current_M33, y_current_M33, z_current_M33,)

    elongate = 1.

    x_normal = -x_current_M33 * elongate / mag
    y_normal = -y_current_M33 * elongate / mag
    z_normal = -z_current_M33 * elongate / mag

    return x_normal, y_normal, z_normal

# setting the global variable for gravitational constant
def set_RADIUS():
    global RADIUS
    RADIUS = 45.  # in kpc


#################
# main function #
#################
def main():

    # initialize the global variables
    set_big_g()
    set_RADIUS()
    print("Starting...")

    #TODO find velocity vector
    #TODO plot just the density of M33
    # look near half mass radius
    # cosider ohase diagram
    # find how people ave looked for it.

    # list with the necessary COM files generated from Assignment 4
    # if these haven't been made, find A4 and make them now
    COM_file = [ "COM_M33.txt" , "COM_M31.txt", "COM_MW.txt"]

    # define where to start
    snap = 215
    dir_path = "/Volumes/GALAXYM33/"

    images = []
    images_data = []
    diff_stdev = []

    while snap < 700:
        if len(str(snap)) == 1:
            MW_file = dir_path + "MW_00" + str(snap) + ".txt"
            M31_file = dir_path + "M31_00" + str(snap) + ".txt"
            M33_file = dir_path + "M33_00" + str(snap) + ".txt"
        elif len(str(snap)) == 2:
            MW_file = dir_path + "MW_0" + str(snap) + ".txt"
            M31_file = dir_path + "M31_0" + str(snap) + ".txt"
            M33_file = dir_path + "M33_0" + str(snap) + ".txt"
        else:
            MW_file = dir_path + "MW_" + str(snap) + ".txt"
            M31_file = dir_path + "M31_" + str(snap) + ".txt"
            M33_file = dir_path + "M33_" + str(snap) + ".txt"


        time, x_COM_M33, y_COM_M33, z_COM_M33, vx_COM_M33, vy_COM_M33, vz_COM_M33 = read_COM("COM_M33.txt")
        time31, x_COM_M31, y_COM_M31, z_COM_M31, vx_COM_M31, vy_COM_M31, vz_COM_M31 = read_COM("COM_M31.txt")
        time33, x_COM_M33, y_COM_M33, z_COM_M33, vx_COM_M33, vy_COM_M33, vz_COM_M33 = read_COM("COM_M33.txt")

        x_M31, y_M31, z_M31, vx_M31, vy_M31, vz_M31 = read_VHighRes(M31_file, "COM_M33.txt")
        x_MW, y_MW, z_MW, vx_MW, vy_MW, vz_MW = read_VHighRes(MW_file, "COM_M33.txt")
        x_M33, y_M33, z_M33, vx_M33, vy_M33, vz_M33 = read_VHighRes(M33_file, "COM_M33.txt")

        x_local_M31, y_local_M31, z_local_M31, vx_local_M31, vy_local_M31, vz_local_M31 = find_wake(x_M31, y_M31, z_M31, vx_M31, vy_M31, vz_M31, x_COM_M33[snap], y_COM_M33[snap], z_COM_M33[snap], RADIUS)
        x_local_MW, y_local_MW, z_local_MW, vx_local_MW, vy_local_MW, vz_local_MW = find_wake(x_MW, y_MW, z_MW, vx_MW, vy_MW, vz_MW, x_COM_M33[snap], y_COM_M33[snap], z_COM_M33[snap], RADIUS)
        x_local_M33, y_local_M33, z_local_M33, vx_local_M33, vy_local_M33, vz_local_M33 = find_wake(x_M33, y_M33, z_M33, vx_M33, vy_M33, vz_M33, x_COM_M33[snap], y_COM_M33[snap], z_COM_M33[snap], RADIUS)

        x_local = x_local_M31 + x_local_MW
        y_local = y_local_M31 + y_local_MW
        z_local = z_local_M31 + z_local_MW
        vx_local = vx_local_M31 + vx_local_MW
        vy_local = vy_local_M31 + vy_local_MW
        vz_local = vz_local_M31 + vz_local_MW

        if snap > 10:
            x_prev_steps, y_prev_steps, z_prev_steps = find_path(snap, x_COM_M33, y_COM_M33, z_COM_M33)


        x_to_M31, y_to_M31, z_to_M31 = find_M31(snap, x_COM_M31, y_COM_M31, z_COM_M31, x_COM_M33, y_COM_M33, z_COM_M33)

        mag_vec = MagnitudeVector(vx_COM_M33[snap], vy_COM_M33[snap], vz_COM_M33[snap], vx_COM_M31[snap], vy_COM_M31[snap], vz_COM_M31[snap])

        x_proj_particles = []
        y_proj_particles = []
        z_proj_particles = []
        vx_proj_particles = []
        vy_proj_particles = []
        vz_proj_particles = []



        for i in range(len(x_local)):
            #dot_prod_vel = np.dot([vx_COM_M33[snap], vy_COM_M33[snap], vz_COM_M33[snap]],[x_local[i], y_local[i], z_local[i]])

            dot_prod_M31 = np.dot([x_to_M31, y_to_M31, z_to_M31] ,[x_local[i], y_local[i], z_local[i]])

            #print(dot_prod_M31)

            if dot_prod_M31 > -5. and dot_prod_M31 < 5.:
                x_proj_particles.append(x_local[i])
                vx_proj_particles.append(vx_local[i])
                y_proj_particles.append(y_local[i])
                vy_proj_particles.append(vy_local[i])
                z_proj_particles.append(z_local[i])
                vz_proj_particles.append(vz_local[i])


        #vel_mag_front = np.sqrt(np.square(vx_front_particles) + np.square(vy_front_particles) + np.square(vz_front_particles))
        #vel_mag_back = np.sqrt(np.square(vx_back_particles) + np.square(vy_back_particles) + np.square(vz_back_particles))

        #normalized_back = vel_mag_back / np.sum(vel_mag_back) * 100
        #normalized_front = vel_mag_front / np.sum(vel_mag_front) * 100
        """
        mu_front = np.mean(vel_mag_front)
        mu_back = np.mean(vel_mag_back)

        stdev_front = round(np.std(vel_mag_front), 2)
        stdev_back = round(np.std(vel_mag_back), 2)

        diff_stdev.append(stdev_front - stdev_back)


        stdev_string = r'Front $\sigma=$' + str(stdev_front) + '\nBack $\sigma=$' + str(stdev_back)

        n, bins, patches = plt.hist([vel_mag_front, vel_mag_back], 25, normed=True, color=['blue', 'red'], label=['Front Particles', 'Back Particles'])
        plt.xlabel("Velocity (km/s)")
        plt.axis([0, 500, 0, 0.01 ])
        plt.ylabel("Count")
        plt.title("Histogram of Particle Velocities at Snap " + str(snap))
        front_line = mlab.normpdf(bins, mu_front, stdev_front)
        back_line = mlab.normpdf(bins, mu_back, stdev_back)
        #plt.legend(handles=[], loc=1)

        fline, = plt.plot(bins, front_line, '--', c='b', label="Front Gaussain Fit")
        bline, = plt.plot(bins, back_line, '--', c='red', label="Back Gaussian Fit")

        plt.text(15, 0.0085, stdev_string)
        plt.legend(handles=[fline, bline], loc=1)
        #lt.legend(stdev_string, loc=2)

        if snap < 10.:
            snap = "00" + str(snap)
        elif snap < 100. and int(snap) > 9.:
            snap = "0" + str(snap)

        file_name = "Velocities_Histo_Snap_" + str(snap) + ".png"
        plt.savefig(file_name)
        plt.show()
        """


        print("num of local particles enclosed ", len(x_local))
        print("num of particles enclosed ", len(x_proj_particles))

        rho_test = ds.grid(x_proj_particles, y_proj_particles, z_proj_particles, 100)

        #images_data.append(np.log10(rho_test.T))
        #print('hewe5re')
        # add mass component
        #width, height
        plt.figure(figsize=(25,10))
        plt.subplot(131)
        #f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
        plt.imshow(np.log10(rho_test.T), origin='lower',extent=[min(x_local), max(x_local), min(y_local), max(y_local)], cmap='spectral', animated=True, vmin=-4., vmax=0.)
        plt.plot(0., 0., '*', c='k')
        plt.plot(vx_COM_M33[snap]*2./mag_vec, vy_COM_M33[snap]*2./mag_vec, '+', c='b')
        plt.plot([0., vx_COM_M33[snap]*2./mag_vec], [0., vy_COM_M33[snap]*2./mag_vec], c='b')
        if snap > 10:
            plt.plot(x_prev_steps, y_prev_steps, c='k')
        #plt.plot([0., x_to_M31], [0., y_to_M31], c='c')
        plt.xlabel("X (kpc)")
        plt.axis([-40, 40, -40, 40])
        plt.ylabel("Y (kpc)")
        #plt.colorbar()  # logarithmic
        #plt.scatter(y_local_M33, z_local_M33)

        rho_test = ds.grid(y_proj_particles, z_proj_particles, x_proj_particles, 100)
        #images_data.append(np.log10(rho_test.T))
        #print('hewe5re')
        # add mass component
        plt.subplot(132)
        #f, axs = plt.subplots(figsize=(15,15))
        plt.imshow(np.log10(rho_test.T), origin='lower',extent=[min(x_local), max(x_local), min(y_local), max(y_local)], cmap='spectral', animated=True, vmin=-4, vmax=0)
        plt.plot(0., 0., '*', c='k')
        plt.plot(vy_COM_M33[snap]*2./mag_vec, vz_COM_M33[snap]*2./mag_vec, '+', c='b')
        plt.plot([0., vy_COM_M33[snap]*2./mag_vec], [0., vz_COM_M33[snap]*2./mag_vec], c='b')
        if snap > 10:
            plt.plot(y_prev_steps, z_prev_steps, c='k')
        #plt.plot([0., x_to_M31], [0., y_to_M31], c='c')
        plt.xlabel("Y (kpc)")
        plt.axis([-40, 40, -40, 40])
        plt.ylabel("Z (kpc)")
        plt.title("Density of Halo Particles at Snap" + str(snap) + "  :" + str(len(x_proj_particles)) + " Particles")


        rho_test = ds.grid(z_proj_particles, x_proj_particles, y_proj_particles, 100)
        #images_data.append(np.log10(rho_test.T))
        #print('hewe5re')
        # add mass component
        plt.subplot(133)
        #f, axs = plt.subplots(figsize=(15,15))
        plt.imshow(np.log10(rho_test.T), origin='lower',extent=[min(x_local), max(x_local), min(y_local), max(y_local)], cmap='spectral', animated=True, vmin=-4, vmax=0)
        plt.plot(0., 0., '*', c='k')
        plt.plot(vz_COM_M33[snap]*2./mag_vec, vx_COM_M33[snap]*2./mag_vec, '+', c='b')
        plt.plot([0., vz_COM_M33[snap]*2./mag_vec], [0., vx_COM_M33[snap]*2./mag_vec], c='b')
        if snap > 10:
            plt.plot(z_prev_steps, x_prev_steps, c='k')
        #plt.plot([0., x_to_M31], [0., y_to_M31], c='c')
        plt.xlabel("Z (kpc)")
        plt.axis([-40, 40, -40, 40])
        plt.ylabel("X (kpc)")
        plt.colorbar()

        #Contours for the plot
        if snap < 10.:
            snap = "00" + str(snap)
        elif snap < 100. and int(snap) > 9.:
            snap = "0" + str(snap)

        file_name = "Wake_Projected_Snap_" + str(snap) + ".png"

        print(file_name)
        plt.tight_layout()
        plt.savefig(file_name, dpi=600)

        # plt.show() has to come AFTER savefig()
        plt.show()
        plt.close()
        print()


        #images.append([figXY])
        #print(images)

        snap = int(snap) + 1

    """
    snap_time = np.arange(0., 10., len(diff_stdev)/10.)

    line_stdev, = plt.plot(snap_time, diff_stdev, label='Standard Deviation Difference')
    plt.title("Difference in Sigma Across All Snaps")
    plt.xlabel("Snapshot Number")
    plt.ylabel("Difference between standard deviation of front and back particles")
    plt.legend(handles=[line_stdev])
    plt.savefig("Diff_Stdev.png")
    plt.show()
    """

    print('end')



if __name__ == '__main__':
    main()
