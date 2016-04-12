# Will Yingling
# April 2016


import numpy as np

# for 2D graphs
import matplotlib.pyplot as plt

# setting the global variable for gravitational constant
def set_big_g():
    global BIG_G
    BIG_G = 4.498768e-6  # in kpc^3/M_sun/Gyr


# Hernquist definition of acceleration
# mass of the component we are looking at. scalar
# some characeristic scale_radius
# x, y, z, are scalars
# dummy _var is just to ID which axis we are looking along
def HernquistAccel(Mass, scale_rad, x, y, z, dummy_var):

    radius = np.sqrt(x**2 + y**2 + z**2)

    # if x
    if dummy_var == 0:
        coord = x
    # if y
    if dummy_var == 1:
        coord = y
    # if z
    if dummy_var == 2:
        coord = z

    # calc the hernquist accel
    HA = -1.*BIG_G*Mass / (radius * (scale_rad + radius)**2) * coord

    return HA


# Miyamoto definition of acceleration
# mass of the component we are looking at. scalar
# some characeristic scale_radius
# x, y, z, are scalars
# dummy _var is just to ID which axis we are looking along
def MiyamotoNagaiAccel(Mass, rd, x, y, z, dummy_var):
    #scale height
    zd = rd / 5.

    #buffer
    B = rd + np.sqrt(z*z + zd*zd)

    R = np.sqrt(x**2 + y**2)

    # if x
    if dummy_var == 0:
        coord = x
    # if y
    if dummy_var == 1:
        coord = y
    # if z we need to adjust the accel eqn because of disk
    if dummy_var == 2:
        equation_for_disk = B / np.sqrt(z**2 + zd**2)
        coord = z * equation_for_disk

    # calc MN acceleration
    MNA = -1.*BIG_G*Mass / (R*R + B*B)**(1.5) * coord

    return MNA


# Calculating Frictional force felt
# x, y, z, vx, vy, vz are scalars of COM
# Msat is mass of M33
def DynamicalFriction(Msat, x, y, z, vx, vy, vz, dummy):

    r = np.sqrt(x**2 + y**2 + z**2)
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    # v_circ = HernquistCircularSpeed(Msat, 62., r)
    v_circ = 240.
    b_max = r
    b_min = G * Msat / v_circ**2  # 250 is circ speed

    #TODO change the v_circ to changing Hern profile


    # in x dir
    if dummy == 0:
        component = vx
    # in y dir
    if dummy == 1:
        component = vy
    # in z dir
    if dummy == 2:
        component = vz

    coloumb_log = np.log(b_max / b_min)

    DF_accel = -0.428 * G * Msat * coloumb_log / (r**2 * v) * component


    fudge = 0.75  # grumble grumble
    DF_accel = DF_accel * fudge

    return DF_accel


# funct to bring all accels together.
def M31Accel(x, y, z, vx, vy, vz, dummy_var):
    # disk scale radius and mass
    rd = 5.
    mass_disk = 1.2e11

    # bulge mass and radius
    mass_bulge = 1.9e10
    radius_bulge = 1.

    # halo scale radius and mass
    mass_halo = 1.9e12
    radius_halo = 62.

    sat_mass = 1.96e11

    halo_accel = HernquistAccel(mass_halo, radius_halo, x, y, z, dummy_var)
    bulge_accel = HernquistAccel(mass_bulge, radius_bulge, x, y, z, dummy_var)
    disk_accel = MiyamotoNagaiAccel(mass_disk, rd, x, y, z, dummy_var)
    fric_accel = DynamicalFriction(sat_mass, x, y, z, vx, vy, vz, dummy_var)

    total_acceleration = halo_accel + bulge_accel + disk_accel + fric_accel

    return total_acceleration


# Looks ahead half a step and predict the orbit
# in one full step
def LeapFrog(dt, x, y, z, vx, vy, vz):

    x_half = x + vx * dt / 2.
    y_half = y + vy * dt / 2.
    z_half = z + vz * dt / 2.

    # print("Halfs",  x_half, y_half, z_half)
    vx_new = vx + M31Accel(x_half, y_half, z_half, vx, vy, vz, 0) * dt
    vy_new = vy + M31Accel(x_half, y_half, z_half, vx, vy, vz, 1) * dt
    vz_new = vz + M31Accel(x_half, y_half, z_half, vx, vy, vz, 2) * dt

    # print("Vels",  vx_new, vy_new, vz_new)
    x_new = x + (vx + vx_new) * dt / 2.
    y_new = y + (vy + vy_new) * dt / 2.
    z_new = z + (vz + vz_new) * dt / 2.
    # print("Step",  x_new, y_new, z_new)

    return x_new, y_new, z_new, vx_new, vy_new, vz_new


# Integrates the orbit from 0 to 10Gyr
def OrbitIntegrator(t_init, dt, t_final, file):
    time = t_init

    # diff set ups for ideal and real cases
    # also writes to a file with the pos and vel and time of the integration
    if file == 0:

        # Setup1
        x = 30.
        y = 0.
        z = 0.
        vx = 0.
        vy = 215.13
        vz = 0.
        open_future_M33 = open("Future_M33_Opt1.txt", '+w')

    if file == 1:

     # SetUp2
        x = -476. + 378.
        y = 491. - 611.
        z = -412. + 285.
        vx = 43. - 74.
        vy = 102. + 72.
        vz = 142. - 49.
        open_future_M33 = open("Future_M33_Opt2.txt", '+w')

    open_future_M33.write("time x y z vx vy vz")
    open_future_M33.write("\n")
    # COMP = [x, y, z, vx, vy, vz]
    COMP = [time, x, y, z, vx, vy, vz]

    # converting COMP to a string
    COMP_to_str = " ".join(str(j) for j in COMP)
    open_future_M33.write(COMP_to_str)
    open_future_M33.write("\n")

    while time < t_final-dt:
        x_new, y_new, z_new, vx_new, vy_new, vz_new = LeapFrog(dt, x, y, z, vx, vy, vz)
        time += dt
        # print("Integerator",  x_new, y_new, z_new)
        COMP = [time, x_new, y_new, z_new, vx_new, vy_new, vz_new]

        # converting COMP to a string
        COMP_to_str = " ".join(str(j) for j in COMP)
        open_future_M33.write(COMP_to_str)
        open_future_M33.write("\n")

        x = x_new
        y = y_new
        z = z_new
        vx = vx_new
        vy = vy_new
        vz = vz_new

    open_future_M33.close()

    return


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


# file: txt file name
# COM : lists of x, y, z COMs
# snap : snap number from which to eval M31 and M33
def mass_enclosed(file, COM, PType, snap, radius):

    d = open(file, 'r')
    data = np.genfromtxt(file, dtype=None, names=True, skip_header=3)

    # array of disk particle indexes
    indexID = np.where(data['type'] == PType)

    # pull out mass array and individual mass
    mass = data['m'][indexID]
    mass_per_P = mass[0]
    # assign x, y, z, values of the unique particles into arrays
    Px = data['x'][indexID]
    Py = data['y'][indexID]
    Pz = data['z'][indexID]

    # look in the COM file, and extract the x, y, z position at desired snap
    x_COM_var = COM[0][snap]
    y_COM_var = COM[1][snap]
    z_COM_var = COM[2][snap]

    # make an array containing the same COM for every element
    x_COM = [x_COM_var]*len(Px)
    y_COM = [y_COM_var]*len(Py)
    z_COM = [z_COM_var]*len(Pz)

    # find the relative distances between particles and COM
    sep_x = Px - x_COM
    sep_y = Py - y_COM
    sep_z = Pz - z_COM

    # calc the magnitude of the separation
    # an array is returned for the mag of sep for every particle
    mag = magnitude(sep_x, sep_y, sep_z)

    # initializing mass enclosed
    mass_enc = 0.

    # search through every particle
    # if the sep of the investigated particle is <= current radius,
    # then add an individual mass
    # to the ever growing enclosed mass
    for j in range(len(mag)):
        if mag[j] <= radius:
            mass_enc += mass_per_P

    tot_mass = np.sum(mass)

    return tot_mass, mass_enc


# calc magnitude
# xx, yy, zz are arrays of the same size
def magnitude(xx, yy, zz):
    mag = []

    for i in range(len(xx)):
        m = np.sqrt(xx[i]**2 + yy[i]**2 + zz[i]**2)
        mag.append(m)

    return mag


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

    # returns sep_in_x, sep_in_y, sep_in_z
    return mag_vec


#################
# main function #
#################
def main():
    set_big_g()

    # list with the necessary COM files generated from Assignment 4
    # if these haven't been made, find A4 and make them now
    COM_file = ["COM_M31.txt", "COM_M33.txt"]  #, "COM_M33.txt"]

    # gal_file = ["/Volumes/GALAXYM33/VLow_Res/M31_000.txt", "/Volumes/GALAXYM33/VLow_Res/M33_000.txt"]

    M31_halo_mass = 1.9e12
    M31_disk_mass = 1.2e11
    M31_bulge_mass = 1.9e10

    R = np.arange(0., 30., 0.5)

    time, x_COM, y_COM, z_COM, vx_COM, vy_COM, vz_COM = read_COM(COM_file[0])

    M31_halo_scale_radius = 62.
    M31_bulge_scale_radius = 1.0
    M31_disk_scale_radius = 5.0


    whole_hog = [  [ [[],[],[]], [[],[],[]], [[],[],[]] ], [ [[],[],[]], [[],[],[]] ]   ]

    M33_M31_r = []
    M33_M31_v = []
    M33_M31_r2 = []
    M33_M31_v2 = []


    # loop over the galaxies
    for galaxy in range(len(COM_file)):

        # print("in", COM_file[galaxy])
        # read in the COM files we made earlier and assign variables
        time, x_COM, y_COM, z_COM, vx_COM, vy_COM, vz_COM = read_COM(COM_file[galaxy])

        # define the com with arrays of x, y, z
        gal_COM_position = [x_COM, y_COM, z_COM]
        gal_COM_velocity = [vx_COM, vy_COM, vz_COM]

        if galaxy == 0:
            M31_position = gal_COM_position
            M31_velocity = gal_COM_velocity
        if galaxy == 1:
            M33_position = gal_COM_position
            M33_velocity = gal_COM_velocity

        """ this is for if you want a rot curve (takes a few minutes to run)
        # loop throught the Ptypes
        for PType in range(3):
            # This processes takes a while to run
            # so this print is just a visual confirmation that the script is running okay
            print("Ptype", PType)

            # pass into here for everything except for when the loop gets to M33's Bulge
            if galaxy == 0 or galaxy == 1 and PType != 2:
                # starting our calc at a reasonable value
                kpc = 1.
                # loop out to radius RADII
                RADII = 31.

                # loop from intial kpc to almost RADII
                while kpc < RADII:

                    radius = kpc

                    # this is where the magic happens. or at least gets executed
                    # this returns an array tot_mass(kpc) and mass_enclosed(kpc)
                    # mass_enclosed(galaxy file, galaxy COM list, PType, snap, radius array)
                    Pmass, mass_enc = mass_enclosed(gal_file[galaxy], gal_COM_position, float(PType)+1., 0, radius)

                    # set the calc'd values to their corresponding locations in whole_hog
                    whole_hog[galaxy][PType][0].append(radius)
                    whole_hog[galaxy][PType][1].append(mass_enc)
                    whole_hog[galaxy][PType][2].append(Pmass)

                    # kpc is our counter so add one to prev value
                    kpc += 1.
            else:
                pass

    # define the corresponding values for M31
    M31_DM_radius = whole_hog[0][0][0]
    M31_DM_mass = whole_hog[0][0][1]
    M31_DM_tot_mass = whole_hog[0][0][2]
    M31_Disk_radius = whole_hog[0][1][0]
    M31_Disk_mass = whole_hog[0][1][1]
    M31_Disk_tot_mass = whole_hog[0][1][2]
    M31_Bulge_radius = whole_hog[0][2][0]
    M31_Bulge_mass = whole_hog[0][2][1]
    M31_Bulge_tot_mass = whole_hog[0][2][2]

    # calc the herquist profile of the galaxies
    # hernquist(summed mass, radius array, buffer parameter a)
    #    M33_hern = hernquist(M33_DM_tot_mass[0], M33_DM_radius, 24.5)
    M31_hern = hernquist(M31_DM_tot_mass[0], M31_DM_radius, M31_halo_scale_radius)
    #    MW_hern = hernquist(MW_DM_tot_mass[0], MW_DM_radius, 62.5)

    # calculate the circular velocities of each PType for M31
    # circular_speed(mass_enclosed array, radius array)
    M31_DM_v_circ = circular_speed_hern(M31_halo_mass, M31_halo_scale_radius, M31_DM_radius)
    M31_Disk_v_circ = circular_speed_hern(M31_disk_mass, M31_disk_scale_radius, M31_Disk_radius)
    M31_Bulge_v_circ = circular_speed_hern(M31_bulge_mass, M31_bulge_scale_radius, M31_Bulge_radius)

    # Calculate the rotation curves
    # rotation_curve(DM ENCLOSED ARRAY, DISK ENCLOSED ARRAY, BULGE ARRAY, radius array)
    #    M33_rot_curve = rotation_curve(M33_DM_mass, M33_Disk_mass, M33_Bulge_mass, M33_DM_radius)
    M31_rot_curve = rotation_curve(M31_DM_v_circ, M31_Disk_v_circ, M31_Bulge_v_circ)
    #    MW_rot_curve = rotation_curve(MW_DM_mass, MW_Disk_mass, MW_Bulge_mass, MW_DM_radius)

    # plot the rotation curves of each component of M31
    # DM, Disk, Bulge, Combined rotation
    line15, = plt.plot(M31_DM_radius, M31_DM_v_circ, label="DM Hern")
    line16, = plt.plot(M31_Disk_radius, M31_Disk_v_circ, label="Disk Hern")
    line17, = plt.plot(M31_Bulge_radius, M31_Bulge_v_circ, label="Bulge Hern")
    line18, = plt.plot(M31_DM_radius, M31_rot_curve, c='k', label="Rotation Curve Disk MN")
    legend = plt.legend(handles=[line15, line16, line17, line18], loc=5)
    plt.title("M31 Rotation curve")
    plt.xlabel("M31 radius kpc")
    plt.ylabel("M31 Vel km/s")
    plt.ylim(40,250)
    plt.show()
    """

    # for the two set ups
    # Start at 0., dt=0.01, 10 Gyr, which set up
    OrbitIntegrator(0., 0.01, 10., 0)
    OrbitIntegrator(0., 0.01, 10., 1)

    integration_file = "Future_M33_Opt1.txt"

    time, xx, yy, zz, vxx, vyy, vzz = read_COM(integration_file)

    integration_file2 = "Future_M33_Opt2.txt"

    time2, xx2, yy2, zz2, vxx2, vyy2, vzz2 = read_COM(integration_file2)

    total_position_circ = magnitude(xx, yy, zz)
    total_velocity_circ = magnitude(vxx, vyy, vzz)
    total_position = magnitude(xx2, yy2, zz2)
    total_velocity = magnitude(vxx2, vyy2, vzz2)

    fake_time = np.arange(0, 10, 0.014268)

    r_sim = MagnitudeVector(M31_position[0], M31_position[1], M31_position[2], M33_position[0], M33_position[1], M33_position[2])
    v_sim = MagnitudeVector(M31_velocity[0], M31_velocity[1], M31_velocity[2], M33_velocity[0], M33_velocity[1], M33_velocity[2])

    line25, = plt.plot(fake_time, r_sim, label="M33 Belsa Sim")
    line35, = plt.plot(time2, total_position, label="M33 Integrator")
    line45, = plt.plot(time2, total_position_circ, label="M33 Position Circ")
    legend = plt.legend(handles=[line25, line35,  line45], loc=5)
    plt.title("M33 Pos vs Time")
    plt.xlabel("time in Myr")
    plt.ylabel("x in kpc")
    # plt.xlim(0,10000)
    plt.ylim(0,500)
    plt.show()

    line24, = plt.plot(fake_time, v_sim, label="M33 Besla Sim")
    line34, = plt.plot(time2, total_velocity, label="M33 Integrator")
    line44, = plt.plot(time2, total_velocity_circ, label="M33 Vel Circ")
    legend = plt.legend(handles=[line24, line34, line44], loc="upper right")
    plt.title("M33 Vel vs time")
    plt.xlabel("time in Myr")
    plt.ylabel("velocity in km/s")
    # plt.xlim(0,10000)
    plt.ylim(0,1000)
    plt.show()

if __name__ == '__main__':
    main()
