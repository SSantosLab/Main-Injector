from typing import Tuple
import numpy as np
import matplotlib.path
import os


def get_footprint_RaDec(file: str = "") -> Tuple[float, float]:
    if file.startswith('desi'):
        file = 'desi_footprint_1.txt'
    else:
        file = 'des_footprint.txt'
    if file == "":
        file='des_footprint.txt'
    gw_data_dir = os.environ["DESGW_DATA_DIR"]
    print(gw_data_dir)
    print(file)
    footFile = gw_data_dir+'/des_footprint.txt'
    print(footFile)
    ra_and_dec = np.genfromtxt(footFile, unpack=False, comments="#")
    ra = ra_and_dec[:, 0]
    dec = ra_and_dec[:,1]

    return ra, dec


def inside_footprint(ra, dec):
    ix = ra > 180
    ra[ix] = ra[ix]-360.
    footprint = get_footprint()
    ix = footprint.contains_points(list(zip(ra, dec)))
    return ix


def get_footprint():
    ra, dec = get_footprint_RaDec()
    footprint = des_path(ra, dec)
    return footprint


def des_path(ra_des, dec_des):
    footprint = matplotlib.path.Path(list(zip(ra_des, dec_des)))
    #footprint = list(zip(ra_des, dec_des))
    return footprint
