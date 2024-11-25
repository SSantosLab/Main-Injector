import healpy as hp
import numpy as np
import shutil

from os.path import basename
from subprocess import run

license="""
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast
   & light style.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

#
# from a healpix map, 
# return  ra, dec, val
#   fluxConservation = False => averaging maps when changing resolution
#   fluxConservation = True => sums maps when changing resolution
#
def flatten(skymap: str, nside: int = 512, backup : bool = False) -> None:
    """
    Flatten an multi-order skymap given and skymap filepath name and nside.
    The output is an flatten skymap with the same name as the input.
    If backup=True, savy a copy of the moc skymap (default is False).
    Arguments:
    ----------
        skymap: str
            multi-order skymap fullpath name.
        nside : int
          nside for flattened skymap None.
        backup: bool (default: False)
           if set to True, save a copy of the moc skymap.
    Returns:
    --------
        None.
    """

    filename = basename(skymap)
    filename = skymap.split('.')[0]
    out = f'{skymap}_flatten.fits.gz'
    command = f'ligo-skymap-flatten '
    command +=f'--nside {nside} {skymap} {out}'
    run(command, shell=True)

    backup_skymap(skymap, backup, filename, out)

def backup_skymap(skymap, backup, filename, out):
    if backup:
        print(f'Saving a copy of {basename(skymap)}')
        shutil.copy(skymap, f'{filename}_moc.fits.gz')

    shutil.move(out, skymap)

def hp2np (hp_map_file, nan=True, degrade=False, fluxConservation=True, field=0):
    hm = hp.read_map(hp_map_file, field=field)
    ra,dec,vals = map2np(hm, resolution=degrade, fluxConservation=fluxConservation)
    return ra,dec,vals

#
# healpix maps have implicit ra,dec-
# one is supposed to know from position along
# vector what the ra,dec is.
# Often this is inconvenient. 
#
# Give me, ra, dec, val.
def map2np(hp_map, resolution=False, fluxConservation=True) :
    nside = hp.npix2nside(len(hp_map))
    if resolution :
        if fluxConservation :
            hp_map = hp.ud_grade(hp_map, resolution, power=-2)
        else :
            hp_map = hp.ud_grade(hp_map, resolution)
        nside = hp.npix2nside(len(hp_map))
    ix = range(0,hp_map.size)
    # pix2and wants the indicies numbers of the pixel to get the coord of
    theta,phi = hp.pix2ang(nside,ix)
    theta = theta*360./2./np.pi;
    phi = phi*360./2./np.pi
    ra = phi
    dec = 90-theta
    ix=np.nonzero(ra > 180); ra[ix]=ra[ix]-360.
    return ra, dec, hp_map

#
# convert a non-hp map into a hp map
#
def convert(template_ra, template_dec, ra, dec, vals) :
    import scipy.spatial
    tree = scipy.spatial.KDTree( list(zip(ra,dec)) )
    newmap = []
    for i in range(0,template_ra.size) :
        data = tree.query(np.array([template_ra[i], template_dec[i]]))
        index = data[1]
        newmap.append(vals[index])
    newmap = np.array(newmap)
    return newmap

def radec2hp(ra, dec, nsides=1024) :
# coordinate conversion
    phi = ra*2*np.pi/360.;
    theta = (90-dec)*2*np.pi/360.
    pix = hp.ang2pix(nsides,theta,phi)
# build an empty healpix map
    npix = hp.nside2npix(nsides)
    vector = np.empty(npix)
    vector[:] = np.NAN
# fill it
    unique_pixels  = np.unique(pix)
    for up in unique_pixels:
        ix = np.nonzero(pix==up)
        count = ix[0].size
        vector[up] = count
    return vector
