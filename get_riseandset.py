import numpy as np
import mags
import obsSlots
import hp2np 
import hexalate

import decam2hp
import jsonMaker
from os import getenv
import pandas as pd
import sys
import  matplotlib.pyplot as plt

def get_riseset(ra, dec, prob, expTime, filt, mjd, best=False, camera="decam"):
    """
    Gets information about given hexes. Based on Jim's simplicity.calc
    code.

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
    best: bool
        True/False if use best or not
    camera: str
        telescope camera name to take in account.
        
    Outputs:
    --------
    hexinfo: pd.dataframe
        dataframe with each row representing a hex. contains ra, dec,
        prob, and rise and set times in MJD (minutes), and is sorted by 
        probability.
    """

    # find the highest prob hex
    best_ix = np.argmax(prob)

    # find night statistics
    night, sunset, sunrise = mags.findNightDuration(mjd, camera)
    night = np.float64(night)*24.
    sunset = np.float64(sunset)
    sunrise = np.float64(sunrise)

    # work every minute from sunset to sunrise
    mjd_list = np.arange(sunset, sunrise+1./1440., 1./1440.)
    
    # calculate the limiting magnitudes
    limit_mag = []
    best_limit_mag = []
    moon_sep = []
    moon_phase = []
    obs = mags.observed(ra,dec,prob, sunset, doMaps=False, verbose=False)
    
    #for every mjd, find limiting magnitude for every hex
    for mjd in mjd_list :
        obs.resetTime(mjd)
        obs.limitMag(filt,exposure=expTime)
        limit_mag.append(obs.maglim)
        best_limit_mag.append(obs.maglim[best_ix])
        moon_sep.append(obs.moonSep[0]*360./2/np.pi)
        moon_phase.append(obs.moonPhase)
    limit_mag = np.vstack(limit_mag)
    
    #set negative limiting magnitudes (meaning magnitudes at which hex isn't visible) to 0
    ix = limit_mag < 0; limit_mag[ix] = 0
    
    #copied from Jim's code, left it here in case someone wants moon information
    moon_sep = np.array(moon_sep)
    moon_phase = np.array(moon_phase)
    # moon phase: where 0 = full, 90 equals half, and  180 = new
    # convert to %full
    moon_phase = ((180.-moon_phase)/180.)*100.
    bins = np.arange(21,25.5,0.5)
    
    set_list = [ ]
    rise_list = [ ]
    
    #find rise and set times -- if the value before is 0 but the value itself isn't 0, it's rise time
    #if value before is not 0 and value itself is 0, it's set time
    for i, row in enumerate(limit_mag):
        for j, value in ((j, row[j]) for j in range(len(row))):
            if value != 0 and limit_mag[i-1][j] == 0:
                rise_list.append([mjd_list[i], ra[j], dec[j]])
        #     #if hex has already risen, count first observable time as rise time. this is not implemented
        #     rise_list.append([mjd_list[i], ra[0], dec[0]])
            if limit_mag[i-1][j] != 0 and value == 0:
                set_list.append([mjd_list[i-1], ra[j], dec[j]])

    
    set_times = []
    rise_times = []
    #print(rise_list)
    
    
    for i in range(len(ra)):
        sets =[]
        rises=[]
        
        #find rise and set times corresponding to each ra and dec
        for s in set_list:
            #print(s[1])
            if s[1]==ra[i] and s[2]==dec[i]:
                #print(s[0], s[1], ra[i], dec[i])
                sets.append(s[0])
                #print(s[1], s[2], ra[i], dec[i], s[0])
        #find rise time for each ra and dec
        for r in rise_list:
            if r[1]==ra[i] and r[2]==dec[i]:
#                print(f'For {ra[i]}, have {r[0]}')
                rises.append(r[0])

    #print(rise_times, set_times)
    
    #create the dataframe
    d = {'RA (deg)': ra, 'DEC (deg)': dec, 'PROBABILITY': prob, 'RISES (MJD)': rise_times, 'SETS (MJD)': set_times}
    info = pd.DataFrame.from_dict(d)
    
    #sort by probability 
    hexinfo = info.sort_values(by='PROBABILITY', ascending=False)
#    print(hexinfo)
    
    return hexinfo
