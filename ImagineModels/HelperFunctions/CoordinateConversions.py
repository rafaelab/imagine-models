import numpy as np


def cyl2cart(coordinate):
    """
    Converts cylindrical coordinates to Cartesian coordinates.
    INPUT:
    -----
     - coordinate: (3,npts)-array with cylindrical coordinates
                   [rho, z, phi] is assumed
    OUTPUT:
    -------
     - ouput : (3,npts)-array of cartesian coordinates
    Created on Oct 21 2016
    @author: V.Pelgrims
    """
    x = coordinate[0] * np.cos(coordinate[2])
    y = coordinate[0] * np.sin(coordinate[2])
    z = coordinate[1]

    return np.array([x,y,z])