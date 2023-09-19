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
import sys
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun
from astropy.coordinates import get_body
import datetime
import ephem
import json

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

    c90 = ci(.5)
    c50 = ci(.9)
    levels = [c50, c90]
        
    area50, area90 = find_area(prob, contours = [.5, .9])
    
    maxprob_ra = round(maxprobcoord[0],2)
    maxprob_dec = round(maxprobcoord[1],2)
    
    maxprob_dist = int(round(distmu[np.argmax(prob)],0))
    maxprob_distsigma = int(round(distsigma[np.argmax(prob)],0))
    
    return (area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob)


def moon_airmass(event_name, todays_date, target_coords):
    date = datetime.date.today()
    #print(date)
    m = ephem.Moon(todays_date)
    #m = ephem.Moon(date)
    phase = round(m.moon_phase, 2)
    
    plt.style.use(astropy_mpl_style)
    quantity_support()
    
    
    
    CTIO = EarthLocation(lat=-30.17*u.deg, lon=-70.80*u.deg, height=3000*u.m)
    HSC = EarthLocation(lat=19.8255*u.deg, lon=-155.476*u.deg, height=4139*u.m)
    utcoffset = -4*u.hour  # Eastern Daylight Time
    
    
    
    #midnight = Time('2012-7-13 00:00:00') - utcoffset
    # mytime = todays_date + ' 00:00:00'
    mytime = todays_date
    midnight = Time(mytime) - utcoffset
    
    
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=CTIO)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
    
    m33 = SkyCoord(target_coords[0], target_coords[1], unit="deg")
    #m33 = SkyCoord.from_name('M33')
    
    m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)
    
    moon_July12_to_13 = get_body("moon", times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    slope = m33altazs_July12_to_13.secz / delta_midnight
    slope = np.array([abs(i.value) for i in slope])
    condition = [i for i in range(len(slope)) if slope[i] > 10]
    n1 = [i == 0 if i in condition else m33altazs_July12_to_13.secz[i] for i in range(len(m33altazs_July12_to_13.secz)) ]
    t = m33altazs_July12_to_13.secz
    t[condition] = np.nan
    
    l1 = ax1.plot(delta_midnight, moonaltazs_July12_to_13.alt, color='blue', ls='--', label='Moon')
    #ax2.plot(delta_midnight, m33altazs_July12_to_13.secz, color='red', ls='--', label='M33')
    l2 = ax2.plot(delta_midnight, t, color = 'red', ls='-', label='Max Prob Coord')
    label = ['Moon (Alt)', 'Max Prob Coord (Airmass)']
    fig.legend([l1, l2], labels=label,
           loc= (0.1, 0.82), fontsize = 9)
    
    ax2.grid(visible = False)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -0*u.deg, color='0.9', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -6*u.deg, color='0.8', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -12*u.deg, color='0.7', zorder=0)
    ax1.fill_between(delta_midnight, 0*u.deg, 90*u.deg, 
                     sunaltazs_July12_to_13.alt < -18*u.deg, color='0.6', zorder=0)
    
    ax1.text(-3, 80, 'Moon phase = {}'.format(phase),bbox=dict(facecolor='blue', alpha=0.3))
    ax1.set_title(todays_date)
    #ax1.legend(loc='upper left')
    #ax2.legend(loc='upper left')
    ax1.set_xlim(-12*u.hour, 12*u.hour)
    ax1.set_xticks((np.arange(13)*2-12)*u.hour)
    ax1.set_ylim(0*u.deg, 90*u.deg)
    ax1.set_xlabel('Hours from CTIO Local Midnight (UTC-4)')
    ax1.set_ylabel('Altitude [deg]')
    ax2.set_ylabel('Airmass')
    ax2.set_ylim(4,1)
    plt.savefig(event_name+'_Moon',dpi=300, bbox_inches = "tight")
    
def make_plots_initial(url, name):
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518". 
    
    This function should return 2 plots: One of the skymap, and another showing the moon phase + altitude and airmass of the max probability coordinate from the skymap. Currently, these are saved in the working directory.'''
    
    date = str(datetime.date.today())
    
    print('URL: ',url)
    print('Name: ',name)
    print('Date: ',date)
    
    
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)
    
    moon_airmass(name, date, [maxprob_ra, maxprob_dec])
    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame

    #airmass(name, [(maxprob_ra, maxprob_dec)])

    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.annotate('Event Name: {}'.format(name) + '\n'
                 + r'50% Area: {} deg$^2$'.format(area50) + '\n'
                 + r'90% Area: {} deg$^2$'.format(area90) + '\n'
                 + r'Max Prob Coordinates: ({},{})'.format(maxprob_ra, maxprob_dec) + '\n'
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
    
    plt.savefig(name+'_initial_skymap',dpi=300, bbox_inches = "tight")
    
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
#######################
    
def make_plots_post_strat(url, name, jsonloc):
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518", and jsonloc is the path to the json file produced by the strategy code. 
    This function should return 3 plots: One of the skymap with hex tilings, another showing the airmass of the coordinates of all hexes, and finally one showing how much probability each hex covers.'''


    print('URL: ',url)
    print('Name: ',name)
    print('jsonloc: ',jsonloc)
    
    f = open(jsonloc)
    data = json.load(f)

    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)

    airmass_input = [(data[i]['RA'], data[i]['dec'], 'Hex_'+str(i)) for i in range(len(data))]
    airmass(name, airmass_input)

    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame


    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.annotate('Event Name: {}'.format(name) + '\n'
                 + r'50% Area: {} deg$^2$'.format(area50) + '\n'
                 + r'90% Area: {} deg$^2$'.format(area90) + '\n'
                 + r'Max Prob Coordinates: ({},{})'.format(maxprob_ra, maxprob_dec) + '\n'
                 + r'Max Prob Distance: {}$\pm${} Mpc'.format(maxprob_dist, maxprob_distsigma)
                 ,(0.9,0.8))
    plt.box(False)
    plt.xticks([])
    plt.yticks([])

    hex_centers = []
    for i in range(len(data)):
        ra = data[i]['RA']
        dec = data[i]['dec']
        centers = SkyCoord(ra, dec, unit="deg")
        hex_centers.append(centers)
        
        
    ax = plt.axes(
        #[0.05, 0.05, 0.9, 0.9],
        projection='astro hours mollweide')
    
    ax_inset = plt.axes(
        [0.9, 0.2, 0.2, 0.2],
        projection='astro zoom',
        center=hex_centers[0],
        radius=25*u.deg)
    
    ax_hexes = []
    for i in range(len(data)):
        axis = plt.axes(projection='astro degrees zoom')
        ra = data[i]['RA']
        dec = data[i]['dec']
        centers = SkyCoord(ra, dec, unit="deg")
        ax_hexes.append(axis)
        
        
        ax.mark_inset_circle(axis, centers, '.9772 deg', edgecolor = 'green')
        ax_inset.mark_inset_circle(axis, centers, '.9772 deg', edgecolor = 'green')
        #ax.mark_inset_circle(axis, centers, '.75 deg', edgecolor = 'green')
        #ax_inset.mark_inset_circle(axis, centers, '.75 deg', edgecolor = 'green')
        axis.remove()
   
    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(False)
        ax_inset.coords[key].set_ticks_visible(True)
        
        
    ax.grid()
    ax_inset.grid()
    ax.mark_inset_axes(ax_inset)
    ax.connect_inset_axes(ax_inset, 'upper left')
    ax.connect_inset_axes(ax_inset, 'lower left')
    #ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()
    #ax_inset.compass(0.9, 0.1, 0.2)

    ax.imshow_hpx(url, cmap='cylon')
    cs = ax.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])
    ct = ax_inset.contour_hpx(url, levels = levels, colors = ['black'], linewidths = [1,0.5])
    ax_inset.imshow_hpx(url, cmap='cylon')
    
    '''ax_inset.plot(
        maxprob_ra, maxprob_dec,
        transform=ax_inset.get_transform('world'),
        marker=ligo.skymap.plot.reticle(),
        markersize=10,
        markeredgewidth=3)'''
    ax.plot(maxprob_ra, maxprob_dec,
        transform=ax.get_transform('world'),
        marker=ligo.skymap.plot.reticle(inner=0),
        markersize=10,
        markeredgewidth=3, label = "Max Prob Coord")
    ax.legend(loc = (0.1,1))
    plt.savefig(name + '_w_Hexes',dpi=300, bbox_inches = "tight")
    
    plt.figure()
    plt.bar(get_prob_from_observing_json(nside, data, prob)[0], get_prob_from_observing_json(nside, data, prob)[1], width = 1, edgecolor = 'k')
    plt.xlabel('Hex #')
    plt.ylabel(r'Cumulative $\%$ of Probability Covered')
    plt.title(name)
    plt.savefig(name+'_Cumulative_Hex_Prob',dpi=300, bbox_inches = "tight")

def parser():
    pass

if __name__ == '__main__':
    pass
