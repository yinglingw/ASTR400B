# Will Yingling
# Apr 2016
# code put together by Nico Garavito Camargo
# http://nbviewer.jupyter.org/github/jngaravitoc/shapes/blob/master/Density.ipynb

from sklearn.neighbors import KDTree
import warnings
warnings.filterwarnings('ignore')
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


#Finding the 32 neighboors
def nearest_neighboors(x, y, z, r):

    D = np.array([x, y, z])
    # transpose the array
    D = D.T
    # scikit. machine learning. Classification. nearest neighbor
    tree = KDTree(D, leaf_size=2500)
    # k refers to number of nearby particles to look at
    dist, ind = tree.query(r, k=33)

    return dist[0], ind[0]


# Evaluating the kernel
def kernel(r, h):

    if r<h:
        W = 1. - 3./2. * (r/h)**2.0 + 3./4. * (r/h)**3.0
    elif((r>h) & (r<2.0*h)):
        W = 1./4. * (2. - r/h)**3.
    else:
        W = 0.0

    return W/(np.pi*h**3.0)


# Computing the local density
def density(x, y, z, mass, r):

    dn, idn = nearest_neighboors(x, y, z, r)
    h = np.max(dn)/2.0
    rho = np.zeros(33.0)
    m = mass[idn]

    for i in range(len(dn)):
        W = kernel(dn[i], h)
        rho[i] = m[i]*W

    return np.sum(rho)


#Making a grid
# This is the part gets passed to first with
# arrays of particle locations
# res is resolution (100-500 is good for wake)
# the time to execute increases with res
def grid(X, Y, Z, res):

    mass = np.ones(len(X))
    rho = np.zeros((res, res))

    rx = np.linspace(min(X)+min(X)*0.2, max(X)+max(X)*0.2, res)
    ry = np.linspace(min(Y)+min(Y)*0.2, max(Y)+max(Y)*0.2, res)

    for i in range(res):
        for j in range(res):
            rho[i][j] = density(X, Y, Z, mass, [rx[i], ry[j], 0])

    return rho

