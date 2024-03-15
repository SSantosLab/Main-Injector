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

    moon_plot = event_name+'/Moon.png'
    plt.savefig(moon_plot, dpi=300, bbox_inches = "tight")
    
    return moon_plot
    
def make_plots_initial(url, name):
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518". 
    
    This function should return 2 plots: One of the skymap, and another showing the moon phase + altitude and airmass of the max probability coordinate from the skymap. Currently, these are saved in the working directory.'''
    
    date = str(datetime.date.today())
    
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)
    
    moon_plot = moon_airmass(name, date, [maxprob_ra, maxprob_dec])
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
    
    initial_skymap_plot = name+'/initial_skymap.png'
    plt.savefig(initial_skymap_plot,dpi=300, bbox_inches = "tight")

    return initial_skymap_plot, moon_plot