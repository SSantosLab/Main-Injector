import numpy as np
from sky.atmosphere import extinctionModel

license = """
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


def scattering_airmass(zenith_distance):
    """
    Return airmass condition given zenith distance.

    Parameter:
    ----------
    zenith_distance: list (?)
        current zenith distance

    Returns:
        airmass (list)
    """

    rad_per_deg = (2*np.pi/360.)
    sin_zd = np.sin(zenith_distance)
    airmass = 1/(np.sqrt(1 - 0.96*sin_zd**2))
    # the function is symmetrical about zd=90, not what we want
    if airmass.size > 1:
        ix = np.nonzero(zenith_distance > 90*rad_per_deg)
        if ix[0].size > 0:
            airmass[ix] = airmass.max()
    return airmass


def moon_mag(phase, filter):
    """
    Return moon's magnitude given phase and filter.


    Parameters:
    -----------
        phase: float
            TBD
        filter: str
            current observation filter.

    Returns:
        mag ( float / list ?)

    """
    vmag = 3.84 + 0.026*phase + 4.0e-9 * phase**4
    vmag = 10**(-0.4*vmag)
    vmag = 26.5754 - 1.0958*np.log(vmag)  

    delgV = 0.34  # ~ from sdss.org filter transforms page
    delgr = 0.44  # ~ from sdss.org filter transforms page
    delri = 0.11  # ~ from sdss.org filter transforms page
    deliz = 0.03  # ~ from sdss.org filter transforms page
    delzy = 0.00  # assume
    if filter == "V":
        mag = 0
    elif filter == "g":
        mag = delgV
    elif filter == "r":
        mag = delgV - delgr
    elif filter == "i":
        mag = delgV - delgr - delri
    elif filter == "z":
        mag = delgV - delgr - delri - deliz
    elif filter == "y":
        mag = delgV - delgr - delri - deliz - delzy
    else:
        raise Exception("no such filter {}".format(filter))
    mag = mag+vmag
    return mag

def sky_brightness_model(filter, zenith_distance):

    extinction = extinctionModel(filter)

    if filter == "g":
        sky = 22.0
        slope = 0.0114
    elif filter == "r":
        sky = 20.9
        slope = 0.0131
    elif filter == "i":
        sky = 19.96
        slope = 0.0081
    elif filter == "z":
        sky = 18.63
        slope = 0.0081
    elif filter == "y":
        sky = 18.0
        slope = 0.0150
    elif filter == "Y":
        sky = 18.0
        slope = 0.0150
    else:
        raise Exception("no such filter {}".format(filter))

    # those values are defined such that delta skybrightness = 0 at 1.3 airmass
    # so remove that
    sky = sky + slope*(zenith_distance-(40.*2*np.pi/360.))
    # and add it tback in
    sky = ks_sky_brightness(sky, extinction, zenith_distance)
    return sky

def scattering_function(moon_sep):
    """

    Parameters:
    -----------
    moon_sep

    Returns:
    --------
    """

    rad_per_deg = (2*np.pi/360.)
    term_1a = 10**5.36
    term_1b = 1.06 + np.cos(moon_sep)**2
    term_2 = 10**(6.15-moon_sep/(40.*rad_per_deg))
    scattering = -2.5*np.log10(term_1a * term_1b + term_2)
    return scattering

def ks_sky_brightness(sky, k, zenith_distance):
    """
    sky:
    k:
    zenith_distance:
    """
    
    airmass = scattering_airmass(zenith_distance)
    sky_zd = -2.5*np.log10(airmass) + k*(airmass - 1)
    mag = sky + sky_zd
    return mag

def moon_brightness_model(filter, extinction, phase,
                        moon_sep, moon_zd, obj_zd):
    """
    """

    test = 0
    lunar_airmass = scattering_airmass(moon_zd)
    obj_airmass = scattering_airmass(obj_zd)
    if test:
        print("\t airmass of moon & object : ", lunar_airmass, obj_airmass)
    lunar_airmass = extinction*lunar_airmass
    obj_airmass = -2.5*np.log10(1-10**(-0.4*extinction*obj_airmass))

    scattering = scattering_function(moon_sep)
    if moon_zd > 93:
        phase = 0.0
    lunar_mag = moon_mag(phase, filter)
    sky = scattering + lunar_mag + lunar_airmass + obj_airmass
    return sky

def skyFiducial(filter: str) -> float:
    """
    """
    zenith_distance = 40.  # degrees = 1.3 airmasses
    sky = sky_brightness_model(filter, zenith_distance)
    return sky



def sky_brightness_at_time(filter, zenith_distance, moon_zd,
                            moon_sep, moon_phase):
    """
    Returns the sky brightness os the sky.

    Parameters:
    -----------
    filter: str
        Observation filter.
    zenith_distance: float / list?
        Zenith distance.
    moon_zd: flot
        Zenith distance of the moon
    moon_sep: float

    """

    extinction = extinctionModel(filter)

    sky_mag = total_sky_brightness(filter, zenith_distance, extinction,
                                   moon_phase, moon_sep, moon_zd)
    return sky_mag

def total_sky_brightness(filter, zenith_distance, extinction,
                         phase, moon_sep, moon_zd):
    """
    """
    sky = sky_brightness_model(filter, zenith_distance)

    moon = moon_brightness_model(filter,
                                 extinction,
                                 phase,
                                 moon_sep,
                                 moon_zd,
                                 zenith_distance)
    # these are additive
    linear_sky = 10**(-0.4*sky + 9)
    linear_moon = 10**(-0.4*moon + 9)
    total = sky + -2.5*np.log10((linear_sky+linear_moon)/linear_sky)
    return total
