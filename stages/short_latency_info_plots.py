"""
Done by Thomas Ruch
"""

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import ligo.skymap.plot
from matplotlib import pyplot as plt
import healpy as hp
import numpy as np
import ligo.skymap
import pandas as pd
from astroplan import Observer, FixedTarget
from astroplan import (AltitudeConstraint, AirmassConstraint,AtNightConstraint)
from astroplan import is_observable
from astropy.time import Time
from astroplan.plots import plot_airmass
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun
from astropy.coordinates import get_body
import datetime
import ephem
import json
import os

### Function for reading + parsing the skymap 
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


def moon_airmass(event_name, todays_date, target_coords):
    date = datetime.date.today()
    m = ephem.Moon(todays_date)
    phase = round(m.moon_phase, 2)
    
    plt.style.use(astropy_mpl_style)
    quantity_support()
    
    
    
    CTIO = EarthLocation(lat=-30.17*u.deg, lon=-70.80*u.deg, height=3000*u.m)
    utcoffset = -4*u.hour  # Eastern Daylight Time

    mytime = todays_date
    midnight = Time(mytime) - utcoffset
    
    
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_tonight = midnight + delta_midnight
    frame = AltAz(obstime=times_tonight, location=CTIO)
    sunaltazs = get_sun(times_tonight).transform_to(frame)
    
    max_prob_coord = SkyCoord(target_coords[0], target_coords[1], unit="deg")
    max_prob_coord_altazs = max_prob_coord.transform_to(frame)
    
    
    moon = get_body("moon", times_tonight)
    moonaltazs = moon.transform_to(frame)
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    t = max_prob_coord_altazs.secz
    condition = [i for i in range(len(t)-1) if abs(t[i] - t[1+1]) > 10]
    t[condition] = np.nan
    
    
    l1 = ax1.plot(delta_midnight, moonaltazs.alt, color='blue', ls='--', label='Moon')
    l2 = ax2.plot(delta_midnight, t, color = 'red', ls='-', label='Max Prob Coord', alpha = 0.5)
    label = ['Moon (Alt)', 'Max Prob Coord (Airmass)']
    fig.legend(labels=label,
           loc= (0.1, 0.82), fontsize = 9)
    
  
    ax2.grid(visible = False)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs.alt < -0*u.deg, color='0.9', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs.alt < -6*u.deg, color='0.8', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs.alt < -12*u.deg, color='0.7', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs.alt < -18*u.deg, color='0.6', zorder=0)
    
    ax1.text(-3, 80, 'Moon phase = {}'.format(phase),bbox=dict(facecolor='blue', alpha = 0.3))
    ax1.set_title(todays_date)
    ax1.set_xlim(-12*u.hour, 12*u.hour)
    ax1.set_xticks((np.arange(13)*2-12)*u.hour)
    ax1.set_ylim(0*u.deg, 90*u.deg)
    ax1.set_xlabel('Hours from CTIO Local Midnight (UTC-4)')
    ax1.set_ylabel('Altitude [deg]')
    ax2.set_ylabel('Airmass')
    ax2.set_ylim(4,1)
    #plt.savefig(event_name+'_Moon.png',dpi=300, bbox_inches = "tight")
    
def make_plots_initial(url, name):
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518". 
    
    This function should return 2 plots: One of the skymap, and another showing the moon phase + altitude and airmass of the max probability coordinate from the skymap. Currently, these are saved in the working directory.'''
    
    date = str(datetime.date.today())
    
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)
    
    moon_airmass(name, date, [maxprob_ra, maxprob_dec])
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
   
    #ax.grid()
    #ax_inset.grid()
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

    
    #plt.savefig(name+'_initial_skymap.png',dpi=300, bbox_inches = "tight")
    
####################### Functions needed post strategy-code ##################################
def ra_dec2theta_phi(ra,dec):
    theta = 0.5 * np.pi - np.pi*dec/180
    phi = np.deg2rad(ra)
    return theta, phi

def get_prob_from_observing_json(NSIDE, json_data, prob_array):
    hex_number = []
    total_prob = []
    for i in range(len(json_data)):
        ra = json_data[i]['RA']
        dec = json_data[i]['dec']
        theta_hex, phi_hex = ra_dec2theta_phi(ra,dec)
        vec = hp.ang2vec(theta_hex, phi_hex)
        decam_hex_disc = hp.query_disc(NSIDE, vec, radius=np.radians(0.9772))
        #decam_hex_disc = hp.query_disc(nside, vec, radius=np.radians(0.75))
        prob_covered = np.sum(prob_array[decam_hex_disc])
        hex_number.append(i)
        total_prob.append(prob_covered)
    
    probcum = [sum(total_prob[:i+1]) for i in range(len(total_prob))]
   
    prob_percent = [i*100 for i in probcum]
    return hex_number, prob_percent

def airmass(event_name,target_coords):
    ''' Target coords is a list of tuples containing ra, dec, and the name of the target'''

    CTIO = Observer(longitude=-70.80*u.deg, latitude=-30.17*u.deg,
                  elevation=3000*u.m, name="CTIO",timezone='America/Santiago')
    
    HSC = Observer(longitude=-155.4761*u.deg, latitude=19.825*u.deg,
                  elevation=4139*u.m, name="Subaru", timezone="US/Hawaii")
    
    tscope = CTIO ; tscope_str = 'CTIO'

    if len(target_coords) > 1:
        targets=[FixedTarget(coord=SkyCoord(ra=coords[0]*u.deg,dec=coords[1]*u.deg),name=coords[2]) for coords in target_coords]
    else:
        targets=FixedTarget(coord=SkyCoord(ra = target_coords[0][0]*u.deg,dec=target_coords[0][1]*u.deg),name='Max Prob Coord')

    constraints = [AltitudeConstraint(10*u.deg, 90*u.deg), AtNightConstraint.twilight_civil()]

    time_range=(Time.now(),Time.now()+1*u.day)
    #time_range=(Time.now()*u.day,Time.now()+1*u.day)
    ever_observable = is_observable(constraints, tscope, targets, time_range=time_range)
    
    

    fig,ax=plt.subplots(figsize=(12,8))
    ax.set_title('DECam Hex Airmasses')
    plot_airmass(targets,ax=ax,
            observer=tscope,
            time=tscope.midnight(Time.now()).to_datetime(timezone=tscope.timezone),
            use_local_tz=True,
            brightness_shading=True,
            max_region=3,min_region=1.5)
    plt.legend(loc='best')
    plt.savefig(event_name+'_Airmass_Hexes',dpi=300, bbox_inches = "tight")
