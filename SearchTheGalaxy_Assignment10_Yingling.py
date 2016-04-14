# Will Yingling
# April 2016

from mpl_toolkits.mplot3d import Axes3D  # for 3D
import numpy as np
import matplotlib.pyplot as plt  # for 2D
from matplotlib.colors import LogNorm  # for colors and scale


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

    # calc the radius
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
    # scale height
    zd = rd / 5.

    # buffer parameter
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
    b_min = BIG_G * Msat / v_circ**2

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

    # calc the coloumb logarithm
    coloumb_log = np.log(b_max / b_min)

    # calc the accel felt due to dyn fric
    DF_accel = -0.428 * BIG_G * Msat * coloumb_log / (r**2 * v) * component

    fudge = 0.75  # "fudge factor" to account for tidal forces instead of a point mass
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

    # find all the forces of friction
    halo_accel = HernquistAccel(mass_halo, radius_halo, x, y, z, dummy_var)
    bulge_accel = HernquistAccel(mass_bulge, radius_bulge, x, y, z, dummy_var)
    disk_accel = MiyamotoNagaiAccel(mass_disk, rd, x, y, z, dummy_var)
    fric_accel = DynamicalFriction(sat_mass, x, y, z, vx, vy, vz, dummy_var)

    # sum the accels. DF is neg
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
    COMP = [time, x, y, z, vx, vy, vz]

    # converting COMP to a string
    COMP_to_str = " ".join(str(j) for j in COMP)
    open_future_M33.write(COMP_to_str)
    open_future_M33.write("\n")

    while time < t_final-dt:
        # Integrate, yo
        x_new, y_new, z_new, vx_new, vy_new, vz_new = LeapFrog(dt, x, y, z, vx, vy, vz)
        time += dt

        COMP = [time, x_new, y_new, z_new, vx_new, vy_new, vz_new]

        # converting COMP to a string then writing to file
        COMP_to_str = " ".join(str(j) for j in COMP)
        open_future_M33.write(COMP_to_str)
        open_future_M33.write("\n")

        # define calc'd values as current values
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


# this reads in the COM files we made with Assigment 4
def read_VHighRes(filename, COMfile):
    # skip header is only 0 because
    # it's the file I made with only 1 header line

    PType = 2.

    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)

    time_COM, x_COM, y_COM, z_COM, vx_COM, vy_COM, vz_COM = read_COM(COMfile)

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

    new_x = Px - x_COM[0]
    new_y = Py - y_COM[0]
    new_z = Pz - z_COM[0]
    new_vx = Pvx - vx_COM[0]
    new_vy = Pvy - vy_COM[0]
    new_vz = Pvz - vz_COM[0]

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

    return new_x, new_y, new_z, new_vx, new_vy, new_vz, x_COM, y_COM, z_COM


def read_bulge(filename):
    # skip header is only 0 because
    # it's the file I made with only 1 header line
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=0)

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

    x_COM = data['xCOM']
    # store y position of particles
    y_COM = data['yCOM']
    # store z position of particles
    z_COM = data['zCOM']

    return new_x, new_y, new_z, new_vx, new_vy, new_vz, x_COM, y_COM, z_COM


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


# to orient the viewing angle to something else
# x and y are arrays
# k is the amount of rotation in radians
def tilt(x, y, k):
    # initialize arrays
    x_tilt = np.zeros(len(x))
    y_tilt = np.zeros(len(y))

    # rotate the viewing angle of the particles
    for i in range(len(x)):
        x_tilt[i] = x[i]*np.cos(k) - y[i]*np.sin(k)
        y_tilt[i] = y[i]*np.cos(k) + x[i]*np.sin(k)

    return x_tilt, y_tilt

# calc the mass to half light ratio
# std is standard deviation
def mass_to_half_light(std):

    r_half_light = 2.  # kpc
    eta = 1.

    mass_light = 6. * eta * r_half_light * std**2 / BIG_G

    return mass_light


#################
# main function #
#################
def main():
    set_big_g()

    # list with the necessary COM files generated from Assignment 4
    # if these haven't been made, find A4 and make them now
    COM_file = ["COM_M31.txt", "COM_M33.txt"]  # , "COM_M33.txt"]

    gal_file = "M31_000.txt"

    # After the code is run once, it creates a file with the variable info
    # so things go quicker. comment out the line below and uncomment the line with the corresponding file
    xx, yy, zz, vxx, vyy, vzz, x_COM, y_COM, z_COM = read_VHighRes(gal_file, COM_file[0])
    # xx, yy, zz, vxx, vyy, vzz, x_COM, y_COM, z_COM = read_bulge("Halo_of_M31.txt")
    # xx, yy, zz, vxx, vyy, vzz, x_COM, y_COM, z_COM = read_bulge("Disk_of_M31.txt")
    # xx, yy, zz, vxx, vyy, vzz, x_COM, y_COM, z_COM = read_bulge("Bulge_of_M31.txt")

    # define the amount of stretching in each dir for isophote contours
    a_axisXY = [.6, 1.2, 1.8, 2.4, 3.6, 4.8]
    b_axisXY = [0.5, 1., 1.5, 2., 3., 4.]

    a_axisYZ = [0.5, 1., 1.5, 2., 3., 4.]
    b_axisYZ = [0.55, 1.1, 1.65, 2.2, 3.3, 4.4]

    # q = minor_axis / major_axis
    qXY = b_axisXY[0] / a_axisXY[0]
    qYZ = a_axisYZ[0] / b_axisYZ[0]

    # execute the tilting of the view
    # last number is amount or rotation in radians
    x_tiltXY, y_tiltXY = tilt(xx, yy, 0.05)
    y_tiltYZ, z_tiltYZ = tilt(yy, zz, -0.15)

    # for the contours
    angle = np.linspace(0, 2*np.pi, 100)

    # plot a density hist for XY plane with colors and contours
    plt.hist2d(x_tiltXY, y_tiltXY, bins=10000, norm=LogNorm())
    plt.plot(0.8, 0.1, '*', c='k')
    plt.title("Plotting Bulge Particles Around M31 COM XY Plane q=" + str(qXY))
    plt.xlabel("X (kpc)")
    plt.ylabel("Y (kpc)")
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.plot(a_axisXY[0]*np.cos(angle)+0.8, b_axisXY[0]*np.sin(angle)+0.1, c='r')
    plt.plot(a_axisXY[1]*np.cos(angle)+0.8, b_axisXY[1]*np.sin(angle)+0.1, c='r')
    plt.plot(a_axisXY[2]*np.cos(angle)+0.8, b_axisXY[2]*np.sin(angle)+0.1, c='r')
    plt.plot(a_axisXY[3]*np.cos(angle)+0.8, b_axisXY[3]*np.sin(angle)+0.1, c='r')
    plt.plot(a_axisXY[4]*np.cos(angle)+0.8, b_axisXY[4]*np.sin(angle)+0.1, c='r')
    plt.plot(a_axisXY[5]*np.cos(angle)+0.8, b_axisXY[5]*np.sin(angle)+0.1, c='r')
    plt.colorbar()
    plt.show()

    # plot a density hist for YZ plane with colors and contours
    plt.hist2d(y_tiltYZ, z_tiltYZ, bins=10000, norm=LogNorm())
    plt.plot(0., .5, '*', c='k')
    plt.title("Plotting Bulge Particles Around M31 COM YZ Plane q=" + str(qYZ))
    plt.xlabel("Y (kpc)")
    plt.ylabel("Z (kpc)")
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.plot(a_axisYZ[0]*np.cos(angle), b_axisYZ[0]*np.sin(angle)+0.5, c='r')
    plt.plot(a_axisYZ[1]*np.cos(angle), b_axisYZ[1]*np.sin(angle)+0.5, c='r')
    plt.plot(a_axisYZ[2]*np.cos(angle), b_axisYZ[2]*np.sin(angle)+0.5, c='r')
    plt.plot(a_axisYZ[3]*np.cos(angle), b_axisYZ[3]*np.sin(angle)+0.5, c='r')
    plt.plot(a_axisYZ[4]*np.cos(angle), b_axisYZ[4]*np.sin(angle)+0.5, c='r')
    plt.plot(a_axisYZ[5]*np.cos(angle), b_axisYZ[5]*np.sin(angle)+0.5, c='r')
    plt.colorbar()
    plt.show()

    # for next ~ 30 lines are 3D plot code copied from
    # http://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib

    # to make a 3D plot of the axial ratios
    fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
    ax = fig.add_subplot(111, projection='3d')

    coefs = (1.2, 1., 1.1)  # are our unique axial ratios
    # Radii corresponding to the coefficients:
    rx, ry, rz = 1. / np.sqrt(coefs)

    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))

    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')

    # Adjustment of the axes, so that they all have the same span:
    max_radius = max(rx, ry, rz)
    for axis in 'xyz':
        getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

    plt.show()
    # end copy

    # print ellipticity
    print("ellipticity of XY plane is e = ", round(1.-qXY, 2))
    print("ellipticity of YZ plane is e = ", round(1.-qYZ, 2), "\n")

    # variable args are arrays
    radius_to_particle = np.sqrt(xx**2 + yy**2 + zz**2)
    # find all particles within a given radius
    indexID = np.where(radius_to_particle <= 10.)

    # redefine vels with only particle within a given kpc of COM
    vx_within_radius = vxx[indexID]
    vy_within_radius = vyy[indexID]
    vz_within_radius = vzz[indexID]

    # calc avg of each vel component
    vx_avg = np.average(vx_within_radius)
    vy_avg = np.average(vy_within_radius)
    vz_avg = np.average(vz_within_radius)

    print("Average Velocity in x-dir = ", round(vx_avg, 2), "km/s")
    print("Average Velocity in y-dir = ", round(vy_avg, 2), "km/s")
    print("Average Velocity in z-dir = ", round(vz_avg, 2), "km/s", "\n")

    # calc the standard deviation
    vx_std = np.std(vx_within_radius)
    vy_std = np.std(vy_within_radius)
    vz_std = np.std(vz_within_radius)

    print("Standard deviation in x-dir = ", round(vx_std, 2))
    print("Standard deviation in y-dir = ", round(vy_std, 2))
    print("Standard deviation in z-dir = ", round(vz_std, 2), "\n")

    # calculate the mass to half light ratio
    x_mass_light = mass_to_half_light(vx_std)
    y_mass_light = mass_to_half_light(vy_std)
    z_mass_light = mass_to_half_light(vz_std)

    print("Mass to half light in x-dir = %e" % round(x_mass_light, 2), "M_sun")
    print("Mass to half light in y-dir = %e" % round(y_mass_light, 2), "M_sun")
    print("Mass to half light in z-dir = %e" % round(z_mass_light, 2), "M_sun\n")

    # plot phase diagrams to show if rotation occurs
    plt.scatter(xx, vyy, marker=".", label="Phase for x vs vy")
    plt.title("Phase Diagram for all Disk particles x vx vy")
    plt.xlabel("x (kpc)")
    plt.ylabel("vy (km/s)")
    plt.show()

    plt.scatter(yy, vxx, marker=".", label="Phase for y vs vx")
    plt.title("Phase Diagram for all Disk particles in y vx vx")
    plt.xlabel("y (kpc)")
    plt.ylabel("vx (km/s)")
    plt.show()

if __name__ == '__main__':
    main()
