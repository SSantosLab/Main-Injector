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
import healpy as hp
from matplotlib.path import Path
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches 
from astropy.io import fits
from astropy import units as u
import ligo.skymap.plot
from matplotlib import pyplot as plt
import healpy as hp
import numpy as np
import ligo.skymap
import datetime
import os

def mjd_decimal_to_minutes_since_sunset(mjd_decimal_list, sunset_mjd_decimal):
    # Constants
    MINUTES_IN_DAY = 24 * 60
    
    def calculate_minutes(mjd_decimal):
        # Convert decimal portion to minutes
        decimal_part = mjd_decimal - int(mjd_decimal)
        minutes_from_decimal = decimal_part * MINUTES_IN_DAY
        
        # Convert whole days to minutes
        whole_days = int(mjd_decimal)
        minutes_from_days = whole_days * MINUTES_IN_DAY
        
        # Total minutes
        total_minutes = minutes_from_days + minutes_from_decimal
        return total_minutes
    
    sunset_minutes = calculate_minutes(sunset_mjd_decimal)
    
    minutes_since_sunset = []
    for mjd_decimal in mjd_decimal_list:
        total_minutes = calculate_minutes(mjd_decimal)
        
        # Calculate minutes since sunset
        time_since_sunset = total_minutes - sunset_minutes
        minutes_since_sunset.append(time_since_sunset)
    
    return minutes_since_sunset

def get_ccd_corners(ccd):
    # Open the file
    with open('zeroed_corners.list', 'r') as file:
        # Read the header line
        header = file.readline().strip().split(',')

        # Find the column index for the specified CCD
        ccd_index = header.index('CCD')
        ccd = f'{ccd}.0'
        
        # Find the column indices for the corner coordinates
        corner1_ra_index = header.index('Corner 1 Ra')
        corner1_dec_index = header.index('Corner 1 Dec')
        corner2_ra_index = header.index('Corner 2 Ra')
        corner2_dec_index = header.index('Corner 2 Dec')
        corner3_ra_index = header.index('Corner 3 Ra')
        corner3_dec_index = header.index('Corner 3 Dec')
        corner4_ra_index = header.index('Corner 4 Ra')
        corner4_dec_index = header.index('Corner 4 Dec')

        # Iterate through the file lines
        for line in file:
            # Split the line into columns
            columns = line.strip().split(',')

            # Check if the CCD matches the specified value
            #print(columns[ccd_index], ccd)
            if columns[ccd_index] == ccd:
                # Extract the corner coordinates
                value = 1
                corner1 = [value*float(columns[corner1_ra_index]), value*float(columns[corner1_dec_index])]
                corner2 = [value*float(columns[corner2_ra_index]), value*float(columns[corner2_dec_index])]
                corner3 = [value*float(columns[corner3_ra_index]), value*float(columns[corner3_dec_index])]
                corner4 = [value*float(columns[corner4_ra_index]), value*float(columns[corner4_dec_index])]
                # Return the list of corner coordinates
                return [corner1, corner2, corner3, corner4]

    # If the CCD is not found, return None
    return None

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
    
    if len(ra)==0:
    	print('No hexes in region')
    	return []

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
        obs.limitMag(filt[0],exposure=expTime[0])
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
        thishex = hex_object.HexObject(ra[i], dec[i], prob[i], rise_times[i], set_times[i], expTime, filt, index = i)
        list_hexes.append(thishex)
        
    if get_sunrise_sunset:
        return list_hexes, sunrise, sunset
    else:
        return list_hexes
    
def calculate_zenith_position(mjd, latitude):
    jd = mjd + 2400000.5
    t = (jd - 2451545.0) / 36525.0
    lst = 280.46061837 + 360.98564736629 * (jd - 2451545) + t * t * (0.000387933 - t / 38710000)
    lst = lst % 360
    ra_zenith = lst
    dec_zenith = latitude
    
    return ra_zenith, dec_zenith

def order_hexes(hex_list):
    new_list = sorted(hex_list, key=lambda x: x.awesomeness_factor, reverse=True)
    return new_list

#FOR TESTING
def plot_hexes_on_skymap(ra_cent, dec_cent, all_cor_ras, all_cor_decs, skymap, eventname="N/A", eventcoords=None):
    total_ras = []
    total_decs = []

    for this_ra, this_dec, these_ras, these_decs in zip(ra_cent, dec_cent, all_cor_ras, all_cor_decs):
        i_possible = decam2hp.isCatalogInHex(int(this_ra), int(this_dec), np.array(these_ras), np.array(these_decs))
        full_shape = (these_ras[i_possible], these_decs[i_possible])

        corners = (np.min(full_shape[0]), np.max(full_shape[0]), np.min(full_shape[1]), np.max(full_shape[1]))

        total_ras.append([corners[0], corners[0], corners[1], corners[1], corners[0]])
        total_decs.append([corners[2], corners[3], corners[3], corners[2], corners[2]])
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518".'''
    name = eventname
    url = skymap
    date = str(datetime.date.today())
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)

    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame
    # print(maxprob_ra)

    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.box(False)
    plt.xticks([])
    plt.yticks([])

    ax = plt.axes(
        projection='astro zoom',
        center=center,
        radius=10*u.deg)

    for key in ['ra', 'dec']:
        ax.coords[key].set_ticklabel_visible(True)
        ax.coords[key].set_ticks_visible(True)


    ax.scalebar((0.1, 0.1), 5 * u.deg).label()

    ct = ax.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])

    ax.imshow_hpx(url, cmap='cylon')
    ax.plot(
        maxprob_ra, maxprob_dec,
        transform=ax.get_transform('world'),
        marker=ligo.skymap.plot.reticle(),
        markersize=30,
        markeredgewidth=3)
    for i in range(len(total_ras)):
        ax.plot(total_ras[i], total_decs[i], transform=ax.get_transform('world'))
    if eventcoords != None:
        ax.scatter(eventcoords, transform=ax.get_transform('world'), color='pink')


def plot_testhex(hex_test, sunset):
    fig, axes = plt.subplots(2,3, figsize=(15,10))
    # axes[1][2].set_visible(False)
    for ax_row in axes:
        for ax in ax_row:
            ax.set_xlabel('Minutes since sunset')


    minutes_since_sunset = mjd_decimal_to_minutes_since_sunset(hex_test.mjd_list, sunset)
    # axes[1][0].set_position([0.24,0.125,0.228,0.343])
    # axes[1][1].set_position([0.55,0.125,0.228,0.343])
    axes[1][0].ticklabel_format(axis="x", style="plain")

    axes[0][0].plot(minutes_since_sunset, hex_test.airmass_Factor)
    axes[0][0].set_title('Airmass Factor')
    axes[0][1].plot(minutes_since_sunset, hex_test.hex_VisibilityFactor)
    axes[0][1].set_title('Hex Visibility Factor')
    axes[0][2].plot(minutes_since_sunset, hex_test.slewTime_Factor)
    axes[0][2].set_title('Slew Time Factor')

    axes[1][0].plot(minutes_since_sunset, hex_test.lunarSeparation_Factor)
    axes[1][0].set_title('Lunar Separation Factor')
    axes[1][1].plot(minutes_since_sunset, hex_test.awesomeness_factor_list)
    axes[1][1].set_title('Awesomeness Factor')

    axes[1][2].plot(minutes_since_sunset, hex_test.coverage_list)
    axes[1][2].set_title('Sky Coverage Factor')

    axes[0][1].axvline(mjd_decimal_to_minutes_since_sunset([hex_test.rise_mjd], sunset), color='red', linestyle='--')
    axes[0][1].axvline(mjd_decimal_to_minutes_since_sunset([hex_test.set_mjd], sunset), color='red', linestyle='--')
    axes[0][0].axvline(mjd_decimal_to_minutes_since_sunset([hex_test.rise_mjd], sunset), color='red', linestyle='--')
    axes[0][0].axvline(mjd_decimal_to_minutes_since_sunset([hex_test.set_mjd], sunset), color='red', linestyle='--')

    
def extract_coordinates(hexlist):
    ras = [hexy.ra for hexy in hexlist]
    decs = [hexy.dec for hexy in hexlist]
    return ras, decs

def generate_coordinates(ras, decs):
    radius = 1.4
    n=100
    xs = np.linspace(-radius, radius, n)
    ys = np.linspace(-radius, radius, n)
    all_x = []
    all_y = []
    for n in range(len(ras)):
        xlist = []
        ylist = []
        for x in xs:
            for y in ys:
                xlist.append(x + ras[n])
                ylist.append(y + decs[n])
        all_x.append(np.array(xlist))
        all_y.append(np.array(ylist))
    return all_x, all_y


def run_or(
    skymap,
    probArea_outer,
    probArea_inner,
    flt,
    expTime_inner,
    expTime_outer,
    mjd,
    resolution=64,
    hexFile=getenv('DATA_DIR')+"/all-sky-hexCenters-decam.txt", 
    trigger_id="LIGO/Virgo", 
    trigger_type="bright", 
    propid='propid', 
    jsonFilename="des-gw.json",
    test=False,
    plot_numbers=False,
    huh=False
):
    
    camera = "decam"
    if (probArea_inner > 1) or (probArea_outer > 1) : raise Exception("probArea_outer or inner is > 1, impossible")

    # get ligo map
    ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)
    # imprint decam hexes onto sky map
    raCenter,decCenter,junk = hexalate.getHexCenters(hexFile)
    nHexes= raCenter.size
    ix, = np.where(ligo > 10e-6) ## AG: 12-13-22 consider a different mode to the hard cut we are using now.
    ra, dec, ligo = ra[ix], dec[ix], ligo[ix]
    sum_prob_in_hex = np.zeros(nHexes)

    treedata = decam2hp.buildtree(ra,dec,nsides=resolution,recompute=True)
    tree = treedata[2]
    sum_prob_in_hex = decam2hp.hexalateMapWithoutOverlap(ra,dec,ligo,tree, raCenter, decCenter, "decam", verbose=False)

    # now we have our data!
    # sort on probability so that max is first on new vector
    index_on_max_prob = np.argsort(sum_prob_in_hex)[::-1]

    # lets find the n'th and 90th percentile list of hexes
    cumsum = np.cumsum(sum_prob_in_hex[index_on_max_prob]) 
    outer_percentile_index = np.argmax( cumsum >= probArea_outer )
    inner_percentile_index = np.argmax( cumsum >= probArea_inner )
    print("N hexes in inner= {}  and in outer= {}, total number of hexes={}".format(
        inner_percentile_index, outer_percentile_index, inner_percentile_index + outer_percentile_index ))

    # now lets get the list of outer and inner hexes
    inner_ra = []; outer_ra=[]
    inner_dec = []; outer_dec=[]
    inner_prob = []; outer_prob=[]
    for i in range(0,inner_percentile_index) :
        # the raCenters sorted on index_on_max_prob are in decending order of prob
        # so we just have to grab them one by on till we reach inner_percentile_index
        inner_ra.append(raCenter[index_on_max_prob][i] )
        inner_dec.append(decCenter[index_on_max_prob][i] )
        inner_prob.append(sum_prob_in_hex[index_on_max_prob][i] )
    for i in range(inner_percentile_index, outer_percentile_index) :
        # and here, too, from inner_percentile_index to outer_percentile_index
        outer_ra.append(raCenter[index_on_max_prob][i] )
        outer_dec.append(decCenter[index_on_max_prob][i] )
        outer_prob.append(sum_prob_in_hex[index_on_max_prob][i] )
    inner_ra = np.array(inner_ra); inner_dec = np.array(inner_dec); 
    outer_ra = np.array(outer_ra); outer_dec = np.array(outer_dec); 
    inner_prob = np.array(inner_prob); outer_prob = np.array(outer_prob); 
    if huh == True:
        return ra, dec, inner_ra, outer_ra, inner_dec, outer_dec, outer_prob, inner_prob,expTime_inner, expTime_outer
    return inner_ra, outer_ra, inner_dec, outer_dec, outer_prob, inner_prob,expTime_inner, expTime_outer

def data_for_onering(strategy_csv):
    #input strategy csv location
    df = pd.read_csv(strategy_csv)
    print(df)
    df.sort_values(by='Deprob1', ascending=False, inplace=True)
    optimal_strategy = df.iloc[0]
    print(optimal_strategy)
    outer = optimal_strategy['Region Coverage']
    inner = optimal_strategy['Region Coverage_deep']
    filt = optimal_strategy['Filter_comb'][0]
    exposure_outer = optimal_strategy['Exposure01']
    exposure_inner = optimal_strategy['Exposure01_deep']
    print(exposure_outer, exposure_inner)
    return outer, inner, filt, exposure_outer, exposure_inner

##FROM NORA FOR PLOTTING SKYMAPS
def make_alert_skymap(map_path):
    skymap = hp.read_map(map_path, field=range(3))
    prob, distmu, distsigma = skymap
    
    npix = len(prob)
    nside = hp.npix2nside(npix)


    maxprobcoord_tup = hp.pix2ang(nside, np.argmax(prob))
    maxprobcoord = [0, 0]
    maxprobcoord[1] = np.rad2deg(0.5*np.pi-maxprobcoord_tup[0])
    maxprobcoord[0] = np.rad2deg(maxprobcoord_tup[1])
    

    def find_area(prob_array, contours, nside = nside):
        areas = []
        for contour in contours:
            sortedprob = np.sort(prob_array)
            probsum = 0
            probcutoff = 1
            npix = 0
            while probsum<contour:
                probsum = probsum+sortedprob[-1]
                probcutoff = sortedprob[-1]
                sortedprob = sortedprob[:-1]
                npix = npix+1

            area = npix * hp.nside2pixarea(nside, degrees=True)
            areas.append(area)
        return ([int(round(i,0)) for i in areas])
    
    
    def ci(level, sorted_prob=np.flip(np.sort(prob))):
        csum = 0
        c = 0
        index = 0
        while csum < level:
            csum += sorted_prob[index]
            c = sorted_prob[index]
            index += 1
        return c

    c50 = ci(.5)
    c90 = ci(.9)
    levels = [c90, c50]
        
    area50, area90 = find_area(prob, contours = [.5, .9])
    
    maxprob_ra = round(maxprobcoord[0],2)
    maxprob_dec = round(maxprobcoord[1],2)
    
    maxprob_dist = int(round(distmu[np.argmax(prob)],0))
    maxprob_distsigma = int(round(distsigma[np.argmax(prob)],0))
    
    return (area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob)

    
def make_skymap_plot(url, name):
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518".'''
    
    
    date = str(datetime.date.today())
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)
    
    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame


    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.annotate('Event Name: {}'.format(os.path.basename(name)) + '\n'
                 + r'50% Area: {} deg$^2$'.format(area50) + '\n'
                 + r'90% Area: {} deg$^2$'.format(area90) + '\n'
                 + r'Max Prob Coordinates (degrees): ({},{})'.format(maxprob_ra, maxprob_dec) + '\n'
                 + r'Max Prob Distance: {}$\pm${} Mpc'.format(maxprob_dist, maxprob_distsigma)
                 ,(0.9,0.8))
    plt.box(False)
    plt.xticks([])
    plt.yticks([])

    ax = plt.axes(projection='astro hours mollweide')

    ax_inset = plt.axes(
        [0.9, 0.2, 0.2, 0.2],
        projection='astro zoom',
        center=center,
        radius=10*u.deg)

    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(False)
        ax_inset.coords[key].set_ticks_visible(True)
   

    ax.mark_inset_axes(ax_inset)
    ax.connect_inset_axes(ax_inset, 'upper right')
    ax.connect_inset_axes(ax_inset, 'lower left')
    ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()

    ax.imshow_hpx(url, cmap='cylon')
    cs = ax.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])
    ct = ax_inset.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])

    ax_inset.imshow_hpx(url, cmap='cylon')
    ax_inset.plot(
        maxprob_ra, maxprob_dec,
        transform=ax_inset.get_transform('world'),
        marker=ligo.skymap.plot.reticle(),
        markersize=30,
        markeredgewidth=3)
    ax.plot(
        maxprob_ra, maxprob_dec,
        transform=ax.get_transform('world'),
        marker=ligo.skymap.plot.reticle(inner=0),
        markersize=10,
        markeredgewidth=3, label = "Max Prob Coord")
    ax.legend(loc = (0.1,1))
