from typing import Tuple
import numpy as np
import matplotlib.path
import os


def get_footprint_RaDec(file: str) -> Tuple[float, float]:

    gw_data_dir = os.environ["DESGW_DATA_DIR"]
    footFile = gw_data_dir + file
    ra, dec = np.genfromtxt(footFile, unpack=True, comments="#")
    return ra, dec


def inside_footprint(ra, dec):
    ix = ra > 180
    ra[ix] = ra[ix]-360.
    footprint = get_footprint()
    ix = footprint.contains_points(zip(ra, dec))
    return ix


def get_footprint():
    ra, dec = get_footprint_RaDec()
    footprint = des_path(ra, dec)
    return footprint


def des_path(ra_des, dec_des):
    footprint = matplotlib.path.Path(zip(ra_des, dec_des))
    return footprint
