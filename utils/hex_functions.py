import hex_object
import numpy as np
import hp2np 
import hexalate
import decam2hp
import jsonMaker
from os import getenv
import os
import pandas as pd
import sys
import  matplotlib.pyplot as plt
import mags
from pyslalib import slalib
import datetime

#get a readable date/time from mjd
def mjd_to_date_time(mjd, with_date = False):
    mjd_integer = int(mjd)
    mjd_decimal = mjd - mjd_integer

    # Convert MJD integer to date
    reference_date = datetime.datetime(1858, 11, 17, 12, 0, 0)  # Reference MJD date
    target_date = reference_date + datetime.timedelta(days=mjd_integer)

    # Convert MJD decimal to time
    total_seconds = mjd_decimal * 24 * 60 * 60
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    time_str = f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"
    
    if with_date:
        return target_date.strftime("%Y-%m-%d") + " " + time_str
    else:
        return time_str
    
def calculate_zenith_position(mjd, latitude):
    jd = mjd + 2400000.5
    t = (jd - 2451545.0) / 36525.0
    lst = 280.46061837 + 360.98564736629 * (jd - 2451545) + t * t * (0.000387933 - t / 38710000)
    lst = lst % 360
    ra_zenith = lst
    dec_zenith = latitude
    
    return ra_zenith, dec_zenith

def get_hexinfo(ra, dec, prob, expTime, filt, mjd, get_sunrise_sunset = False, camera="decam"):
    """
    Gets information about given hexes, primarily rise and set times, then initializes
    a hex object for each (from hex_object.py) and stores them in a list. Rise and set 
    time calculations, along with sunrise/sunset calculations, based on Jim's 
    simplicity.calc code and draw heavily from his mags.py code. 
    

    Parameters:
    -----------
    ra: np.vector
        vector of hex ra
    dec: np.vector
        vector of hex dec
    prob: np.vector
        vector of hex summed prob from LIGO map
    expTime: np.vector
        vector of exposure times
    filt: str
        telescope filt.
    mjd: float
        modified julian date.
    get_sunrise_sunset
        True/False if find and return sunrise/sunset for given mjd night
    camera: str
        telescope camera name to take in account.
        
    Outputs:
    --------
    list_hexes: list
        list of hex objects
    sunrise: np.float64
        mjd decimal of given night's sunrise, given if get_sunrise_sunset is true
    sunset: np.float64
        mjd decimal of given night's sunset, given if get_sunrise_sunset is true
    """

    # Find night statistics
    night, sunset, sunrise = mags.findNightDuration(mjd, camera)
    night = np.float64(night)*24.
    sunset = np.float64(sunset)
    sunrise = np.float64(sunrise)
#     print(sunset)

    # Work every 30 seconds from sunset to sunrise, commented out is working every min
    #mjd_list = np.arange(sunset, sunrise+1./1440., 1./1440.)
    mjd_list = np.arange(sunset, sunrise + (30. / 86400.), (30. / 86400.))

    # Calculate the limiting magnitudes at each mjd
    limit_mag = []
    moon_sep = []
    moon_phase = []
    obs = mags.observed(ra,dec,prob, sunset, doMaps=False, verbose=False)
    for mjd in mjd_list :
        obs.resetTime(mjd)
        obs.limitMag(filt,exposure=expTime)
        limit_mag.append(obs.maglim)
    limit_mag = np.vstack(limit_mag)
    
    # Set any limiting magnitudes less than 0 to 0
    ix = limit_mag < 0; limit_mag[ix] = 0

    # Initialize lists to store rise and set times 
    set_list = [ ]
    rise_list = [ ]

    # Find rise and set times. to find rise time, for each object it looks for the 
    # index where the previous limiting magnitude was 0 but the current index 
    # is nonzero and saves the corresponding mjd, or if this never occurs, assume
    # it rises at sunset. to find set time, for each object look for index where
    # previous limiting mag was nonzero but current limiting mag is 0 and saves the
    # corresponding mjd. if this doesn't occur assume it sets at sunrise. 
    for i in range(len(limit_mag)):
        for j in range(len(limit_mag[i])):
            if limit_mag[i][j] !=0 and limit_mag[i-1][j] == 0:
                rise_list.append([mjd_list[i], ra[j], dec[j]])
            elif limit_mag[i-1][j] !=0 and limit_mag[i][j] == 0:
                set_list.append([mjd_list[i], ra[j], dec[j]])
                #print(f'setting j={j}, i={i}, {ra[j]}, {mjd_list[i]}')
    
    set_times = []
    rise_times = []
    
    #find correct rise and set time for each hex
    for i in range(len(ra)):
        sets =[]
        rises=[]
        
        #find rise and set times corresponding to each ra and dec
        for s in set_list:
            if s[1]==ra[i] and s[2]==dec[i]:
                sets.append(s[0])

        #find rise time for each ra and dec set (in other words, each hex)
        for r in rise_list:
            if r[1]==ra[i] and r[2]==dec[i]:
                rises.append(r[0])

        if len(sets) != 0:
            set_times.append(np.max(sets))
        elif len(sets) == 0:
            set_times.append(mjd_list[-1])
        if len(rises) !=0:
            rise_times.append(np.min(rises))
        elif len(rises) == 0:
            rise_times.append(mjd_list[0])
    
    list_hexes = []
    for i in range(len(ra)):
        thishex = hex_object.HexObject(ra[i], dec[i], prob[i], rise_times[i], set_times[i], expTime, index = i)
        list_hexes.append(thishex)
        
    if get_sunrise_sunset:
        return list_hexes, sunrise, sunset
    else:
        return list_hexes
    