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
from astropy.coordinates import AltAz, EarthLocation, SkyCoord,get_sun,get_moon
from astropy.coordinates import get_sun
from astropy.coordinates import get_body
import datetime
import ephem
import json
import os
from astropy.coordinates import ICRS
import pytz
import subprocess

def getSunTimes(date,location,elevations=[-14*u.deg,-18*u.deg]):
    """
    A function to give the sunset, sunrise, and time of sun elevations as a function of the night of observation
    """
    sunset = getSunset(date,location)
    sunrise = getSunrise(date,location)
    sunElTimes = {}
    for el in elevations:
        sunElTimes[el] = getSunElevations(date,location,el)
    return sunset,sunrise,sunElTimes

def getSunset(date,location):
    """
    Get the sunset time (UTC) for a given date and location.
    
    Parameters:
    - date (np.datetime64): The date of interest (only the date part is used).
    - location (EarthLocation): Observer's location on Earth.
    
    Returns:
    - np.datetime64: Sunset time in UTC.
    """
    # Convert date to Astropy Time object (at noon to avoid timezone edge cases)
    astropy_time = Time(date.astype('datetime64[D]') + np.timedelta64(12, 'h'), scale='utc') # this might need adjustment

    # AltAz frame for that location and time
    altaz_frame = AltAz(obstime=astropy_time, location=location)

    # Find the sunset by stepping in small intervals
    times = astropy_time + np.linspace(-12, 12, 1000) * u.hour
    sun_altazs = get_sun(times).transform_to(AltAz(obstime=times, location=location))
    
    # Sunset is where altitude crosses from positive to negative
    altitudes = sun_altazs.alt.deg
    crossing_indices = np.where((altitudes[:-1] > 0) & (altitudes[1:] <= 0))[0]

    if len(crossing_indices) == 0:
        raise ValueError("Sunset not found for the given date and location.")
        return None
    
    sunset_time_astropy = times[crossing_indices[0]]

    # Convert back to np.datetime64
    sunset_time_np = np.datetime64(sunset_time_astropy.utc.iso, 's')
    return sunset_time_np

def getSunrise(date,location):
    """
    Get the sunrise time (UTC) for a given date and location.

    Parameters:
    - date (np.datetime64): The date of interest.
    - location (EarthLocation): Observer's location on Earth.

    Returns:
    - np.datetime64: Sunrise time in UTC.
    """
    # Convert input date to an astropy Time object at midnight UTC
    midnight = Time(date.astype('datetime64[D]'), scale='utc') # This might need some adjustment...

    # Create a range of times around the date
    times = midnight + np.linspace(0, 24, 1000) * u.hour  # full day coverage

    # Calculate the Sun's altitude at these times
    frame = AltAz(obstime=times, location=location)
    altitudes = get_sun(times).transform_to(frame).alt.deg

    # Find where the sun crosses from below 0° to above 0° (sunrise)
    crossing_indices = np.where((altitudes[:-1] <= 0) & (altitudes[1:] > 0))[0]

    if len(crossing_indices) == 0:
        raise ValueError("Sunrise not found for the given date and location.")
        return None

    # Sunrise time is at the first crossing
    sunrise_time_astropy = times[crossing_indices[0]]

    # Convert back to np.datetime64
    sunrise_time_np = np.datetime64(sunrise_time_astropy.utc.iso, 's')
    return sunrise_time_np

def getSunElevations(date,location,target_altitude):
    """
    Find the UTC times when the Sun rises and sets relative to a specified altitude 
    for a given date and location.

    Parameters:
    - date (np.datetime64): The date of interest (UTC).
    - location (EarthLocation): Observer's location.
    - target_altitude (Quantity): Desired Sun altitude (e.g., 0*u.deg for sunrise/sunset).

    Returns:
    - [np.datetime64, np.datetime64]: Tuple of (rising_time, setting_time) in UTC.
    """
    if not target_altitude.unit.is_equivalent(u.deg):
        raise ValueError("target_altitude must have angular units (e.g., degrees)")

    # Start from midnight UTC of the given date
    midnight = Time(date.astype('datetime64[D]'), scale='utc')

    # Create a full day grid of times
    times = midnight + np.linspace(0, 24, 2000) * u.hour

    # Calculate Sun altitudes
    frame = AltAz(obstime=times, location=location)
    altitudes = get_sun(times).transform_to(frame).alt

    # Find crossings: where sign of (altitude - target) changes
    alt_diff = altitudes - target_altitude
    crossing_indices = np.where(alt_diff[:-1] * alt_diff[1:] <= 0)[0]

    if len(crossing_indices) < 2:
        raise ValueError("Could not find both rising and setting times for this altitude on this date.")
        return None

    # Rising: first crossing; Setting: second crossing
    rising_idx = crossing_indices[0]
    setting_idx = crossing_indices[1]

    # Linear interpolation for more precise crossing times
    def interpolate_time(idx):
        t1, t2 = times[idx], times[idx+1]
        a1, a2 = altitudes[idx].value, altitudes[idx+1].value
        frac = (target_altitude.value - a1) / (a2 - a1)
        return t1 + frac * (t2 - t1)

    rising_time = interpolate_time(rising_idx)
    setting_time = interpolate_time(setting_idx)

    # Convert to np.datetime64
    rising_time_np = np.datetime64(rising_time.utc.iso, 's')
    setting_time_np = np.datetime64(setting_time.utc.iso, 's')

    return [rising_time_np, setting_time_np]

def getMoonTimes(maxProbCoord,date,location,crossing=30*u.deg):
    """
    A function to give the moonrise, moonset, and moon closest separation measurements. Also, crossing of 30 deg separation if it is present on the night of observation
    """
    moonrise = getMoonrise(date,location)
    moonset = getMoonset(date,location)
    moonClosest = getMoonClosest(date,location,maxProbCoord)
    moonCrossing = getMoonCrossing(date,location,maxProbCoord,crossing)
    return moonrise,moonset,moonClosest,moonCrossing

def getMoonrise(date,location):
    """
    Find the UTC time of moonrise for a given date and location.

    Parameters:
    - date (np.datetime64): The date of interest (UTC).
    - location (EarthLocation): Observer's location.

    Returns:
    - np.datetime64: Moonrise time in UTC.
    """
    # Start from midnight UTC of the given date
    midnight = Time(date.astype('datetime64[D]'), scale='utc')

    # Create a full day grid of times
    times = midnight + np.linspace(0, 24, 2000) * u.hour

    # Calculate Moon altitudes
    frame = AltAz(obstime=times, location=location)
    altitudes = get_moon(times).transform_to(frame).alt

    # Find crossings: Moon altitude crosses 0 degrees upwards
    alt_diff = altitudes - 0*u.deg
    crossing_indices = np.where((alt_diff[:-1] < 0) & (alt_diff[1:] >= 0))[0]

    if len(crossing_indices) == 0:
        raise ValueError("No moonrise found on this date at this location.")
        return None

    # First moonrise
    rising_idx = crossing_indices[0]

    # Linear interpolation for more precise crossing time
    def interpolate_time(idx):
        t1, t2 = times[idx], times[idx+1]
        a1, a2 = altitudes[idx].value, altitudes[idx+1].value
        frac = (0 - a1) / (a2 - a1)
        return t1 + frac * (t2 - t1)

    moonrise_time = interpolate_time(rising_idx)

    # Convert to np.datetime64
    moonrise_time_np = np.datetime64(moonrise_time.utc.iso, 's')

    return moonrise_time_np

def getMoonset(date,location):
    """
    Find the UTC time of moonset for a given date and location.

    Parameters:
    - date (np.datetime64): The date of interest (UTC).
    - location (EarthLocation): Observer's location.

    Returns:
    - np.datetime64: Moonset time in UTC.
    """
    # Start from midnight UTC of the given date
    midnight = Time(date.astype('datetime64[D]'), scale='utc')

    # Create a full day grid of times
    times = midnight + np.linspace(0, 24, 2000) * u.hour

    # Calculate Moon altitudes
    frame = AltAz(obstime=times, location=location)
    altitudes = get_moon(times).transform_to(frame).alt

    # Find crossings: Moon altitude crosses 0 degrees downward
    alt_diff = altitudes - 0*u.deg
    crossing_indices = np.where((alt_diff[:-1] > 0) & (alt_diff[1:] <= 0))[0]

    if len(crossing_indices) == 0:
        raise ValueError("No moonset found on this date at this location.")
        return None

    # First moonset
    setting_idx = crossing_indices[0]

    # Linear interpolation for more precise crossing time
    def interpolate_time(idx):
        t1, t2 = times[idx], times[idx+1]
        a1, a2 = altitudes[idx].value, altitudes[idx+1].value
        frac = (0 - a1) / (a2 - a1)
        return t1 + frac * (t2 - t1)

    moonset_time = interpolate_time(setting_idx)

    # Convert to np.datetime64
    moonset_time_np = np.datetime64(moonset_time.utc.iso, 's')

    return moonset_time_np

def getMoonClosest(date,location,target_coord):
    """
    Find the UTC time when the Moon is closest to a given sky coordinate.

    Parameters:
    - date (np.datetime64): The date of interest (UTC).
    - location (EarthLocation): Observer's location.
    - target_coord (SkyCoord): Target sky coordinate.

    Returns:
    - np.datetime64: Time of minimum angular separation (UTC).
    """
    # Start from midnight UTC of the given date
    midnight = Time(date.astype('datetime64[D]'), scale='utc')

    # Sample the whole day finely
    times = midnight + np.linspace(0, 24, 2000) * u.hour

    # Get Moon positions at these times
    moon_coords = get_moon(times, location=location)

    # Compute angular separation from target
    separations = moon_coords.separation(target_coord)

    # Find minimum separation
    min_idx = np.argmin(separations)
    min_sep = np.min(separations)

    # Optionally refine the minimum by fitting around it
    # For now, we just interpolate linearly between neighboring points
    if 0 < min_idx < len(times)-1:
        t1, t2 = times[min_idx-1], times[min_idx+1]
        s1, s2 = separations[min_idx-1].arcsec, separations[min_idx+1].arcsec
        sm = separations[min_idx].arcsec

        # Fit a parabola and find the vertex
        denom = (s1 - 2*sm + s2)
        if denom != 0:
            frac = 0.5 * (s1 - s2) / denom
            best_time = times[min_idx] + frac * (t2 - t1)/2
        else:
            best_time = times[min_idx]
    else:
        best_time = times[min_idx]

    # Convert to np.datetime64
    best_time_np = np.datetime64(best_time.utc.iso, 's')

    return best_time_np,min_sep


def getMoonCrossing(date,location,target_coord,crossing):
    """
    Find the UTC time when the Moon's angular separation from a target
    drops below a given threshold.

    Parameters:
    - date (np.datetime64): The date of interest (UTC).
    - location (EarthLocation): Observer's location.
    - target_coord (SkyCoord): Target sky coordinate.
    - threshold (Quantity): Threshold separation (e.g., 5*u.deg).

    Returns:
    - np.datetime64: Time of crossing below threshold (UTC).
    """
    # Start from midnight UTC of the given date
    midnight = Time(date.astype('datetime64[D]'), scale='utc')

    # Create time grid for the full day
    times = midnight + np.linspace(0, 24, 2000) * u.hour

    # Get Moon coordinates at these times
    moon_coords = get_moon(times, location=location)

    # Calculate angular separations
    separations = moon_coords.separation(target_coord)

    # Find where separation crosses below threshold
    crossing_indices = np.where((separations[:-1] > threshold) & (separations[1:] <= threshold))[0]

    if len(crossing_indices) == 0:
        raise ValueError("No crossing below the threshold found on this date at this location.")

    # Take the first crossing
    idx = crossing_indices[0]

    # Interpolate time for more precise crossing
    t1, t2 = times[idx], times[idx+1]
    s1, s2 = separations[idx].deg, separations[idx+1].deg
    frac = (threshold.value - s1) / (s2 - s1)
    crossing_time = t1 + frac * (t2 - t1)

    # Convert to np.datetime64
    crossing_time_np = np.datetime64(crossing_time.utc.iso, 's')

    return crossing_time_np

def make_agn_plot(plots_path,trigger_id,mass_chirp,dist,distsigma):   
 
    print(plots_path,type(plots_path))
    print(trigger_id,type(trigger_id))
    print(mass_chirp,type(mass_chirp))
    print(dist,type(dist))
    print(distsigma,type(distsigma))
    plotString =    'python ' +\
                    f'./handlers/PhysicalObservable.py ' +\
                    f'--triggerid {trigger_id} ' +\
                    f'--chirpmass {mass_chirp} ' +\
                    f'--distance {dist} ' +\
                    f'--sigma_dist {distsigma} ' +\
                    f'--save_path {plots_path}'

    subprocess.run(plotString,shell=True)

    return str(plots_path)+'/RealEventsLum_'+trigger_id+'.png'

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


def moon_airmass(event_name, todays_date, target_coords,return_many=False):
    date = datetime.date.today()
    m = ephem.Moon(todays_date)
    phase = round(m.moon_phase, 2)
    
    plt.style.use(astropy_mpl_style)
    quantity_support()
    
    
    
    CTIO = EarthLocation.of_site('Cerro Tololo Interamerican Observatory')
    chile_now = datetime.datetime.now(pytz.timezone('Chile/Continental'))
    utcoffset =  u.hour * int(chile_now.utcoffset().total_seconds()/60/60)
    # utcoffset = 0
    # print("UTC offset:",utcoffset)

    mytime = todays_date
    midnight = Time(mytime) - utcoffset # - -> plus?
    
    
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_tonight = midnight + delta_midnight

    print(todays_date,times_tonight)

    frame = AltAz(obstime=times_tonight, location=CTIO)
    sunaltazs = get_sun(times_tonight).transform_to(frame)
    
    max_prob_coord = SkyCoord(target_coords[0], target_coords[1], unit="deg",frame='icrs')
    max_prob_coord_altazs = max_prob_coord.transform_to(frame)
    
    moon = get_body("moon", times_tonight)
    moonaltazs = moon.transform_to(frame)

    moon_separation = max_prob_coord_altazs.separation(moonaltazs)
    
    # print("Max moon sep: {}".format(np.min(moon_separation.deg)))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    t = max_prob_coord_altazs.secz
    condition = [i for i in range(len(t)-1) if abs(t[i] - t[i+1]) > 10]
    t[condition] = np.nan
    
    
    l1 = ax1.plot(delta_midnight, moon_separation.deg, color='blue', ls='--', label='Moon separation')
    l2 = ax2.plot(delta_midnight, t, color = 'red', ls='-', label='Max Prob Coord', alpha = 0.5)
    label = ['Moon separation (deg)', 'Max Prob Coord (Airmass)']
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
    ax1.set_ylim(60*u.deg, 90*u.deg)
    ax1.set_xlabel("Hours from CTIO Local Midnight (UTC{})".format(utcoffset))
    ax1.set_ylabel('Altitude [deg]')
    ax2.set_ylabel('Airmass')
    ax2.set_ylim(2,1)

    print("Peak airmass: {}\n Peak airmass time (local): {}".format(np.nanmin(t),delta_midnight[np.nanargmin(t)]))

    moon_plot = event_name+'/Moon.png'
    moon_plot = f'/data/des70.a/data/desgw/O4/Main-Injector-O4b/utils/Moon_{todays_date}.jpg' #uncomment this line if you are using the moonplot figure in utils
    plt.savefig(moon_plot, dpi=300, bbox_inches = "tight")
    os.chmod(moon_plot, 0o0777)
    
    # Clear the current axes.
    plt.cla() 
    # Clear the current figure.
    plt.clf() 
    # Closes all the figure windows.
    plt.close('all')


    if return_many:
        return moon_plot,sunaltazs,delta_midnight,moon_separation,t
    else:   
        return moon_plot
 
def __format_sepMessage(dateOfInterest,sunset,sunrise,sunElTimes,moonrise,moonset,moonClosest,moonCrossing):
    msg = "**Observability statistics for the night of {}**\n\n".format(dateOfInterest) +\
           "Sunset time: {}\n".format(sunset) +\
           "Sunrise time: {}\n".format(sunrise)
    for el,times in zip(sunElTimes.keys(),sunElTimes.values()):
        msg += "Sun crossing of {} deg at {} and {}\n".format(el,times[0],times[1])
    msg += "Moonrise time: {}\n".format(moonrise) +\
           "Moonset time: {}\n".format(moonset) +\
           "Moon closest angular separation to max probability coordinate at {}, {:.2f} degrees separation\n".format(moonClosest[0],moonClosest[1]) +\
           "Moon crosses 30 deg separation at {}".format(moonCrossing)
    return msg
 
 
def make_plots_initial(url, name,chirpEstimate,dist,distsigma,slackBot=None):
    '''url is either the skymap url or the local path to the skymap, name is something like "S230518". 
    
    This function should return 2 plots: One of the skymap, and another showing the moon phase + altitude and airmass of the max probability coordinate from the skymap. Currently, these are saved in the working directory.'''
    
    date = str(datetime.date.today())
    
    
    area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(url)
    
    moon_plot = moon_airmass(name, date, [maxprob_ra, maxprob_dec])
    center = SkyCoord(maxprob_ra, maxprob_dec, unit="deg")  # defaults to ICRS frame

     for dateAdjust in np.arange(3):
        dateOfInterest = np.datetime64(datetime.date.today())+np.timedelta64(dateAdjust,"D")
        sunset,sunrise,sunElTimes = getSunTime(dateOfInterest,CTIO)
        moonrise,moonset,moonClosest,moonCrossing = getMoonTimes(maxProbCoord,dateOfInterest,CTIO)
        msg = __format_sepMessage(dateOfInterest,sunset,sunrise,sunElTimes,moonrise,moonset,moonClosest,moonCrossing)
        if slackBot!=None:
            slackBot.post_message("",msg)# Post the message
        else:
            print(msg)
        

    event_id = os.path.basename(os.path.abspath(os.path.join(name ,"../..")))
    
    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.annotate('Event Name: {}'.format(event_id) + '\n'
                 + r'50% Area: {} deg$^2$'.format(area50) + '\n'
                 + r'90% Area: {} deg$^2$'.format(area90) + '\n'
                 + r'Max Prob Coordinates (degrees): ({},{})'.format(maxprob_ra, maxprob_dec) + '\n'
                 + r'Weighted average distance: {}$\pm${} Mpc'.format(dist, distsigma) + '\n'
                 + 'Chirp mass estimate: {:.2f}$M$'.format(float(chirpEstimate)),(0.9,0.8))
    plt.box(False)
    plt.xticks([])
    plt.yticks([])

    ax = plt.axes(projection='astro hours mollweide')

    # print("Limits:",ax.get_ylim(),ax.get_xlim())

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


    ### Add galactic plane and +- 15 deg to skymap plot 
    
    seanLimit = 15 # The upper and lower limit on the galactic latitude range - typically, this is 15 degrees
    galacticLatitude =np.append(np.arange(121,474,step=1), []) 

    galacticCenterline = np.full(np.shape(galacticLatitude),0)
    galacticLowerLimit = np.full(np.shape(galacticLatitude),-seanLimit)
    galacticUpperLimit = np.full(np.shape(galacticLatitude),seanLimit)

    galacticCenterlineCoords = SkyCoord(l=galacticLatitude*u.degree,b=galacticCenterline*u.degree,frame='galactic')
    galacticLowerLimitCoords = SkyCoord(l=galacticLatitude*u.degree,b= galacticLowerLimit*u.degree,frame='galactic')
    galacticUpperLimitCoords = SkyCoord(l=galacticLatitude*u.degree,b= galacticUpperLimit*u.degree,frame='galactic')

    # plot it

    galaxyKwargs = {"Center": {'ls':"--",'color':'black','label':"Galactic centerline"},
                    "Upper limit": {'ls':"--","color":'black','alpha':0.5,'label':"Galactic latitude limit +/- {} deg".format(seanLimit)},
                    "Lower limit": {'ls':"--","color":'black','alpha':0.5}}

    for coord,label,galkey in zip([galacticCenterlineCoords,galacticLowerLimitCoords,galacticUpperLimitCoords],["Galactic center","Galactic lower limit","Galactic upper limit"],galaxyKwargs.keys()):
        # coord = 
        # wrapLoc = 267*u.deg
        
        # print("WrapLocation:", wrapLoc)
        galRa = coord.icrs.ra
        galDec = coord.icrs.dec
        ax.plot(galRa,galDec,transform=ax.get_transform('world'),**galaxyKwargs[galkey])

    # print("Galactic coords after transformation:",galacticCenterlineCoords.ra,galacticCenterlineCoords.dec)
    # print("Galactic coords after transformation:",galacticCenterlineCoords.l,galacticCenterlineCoords.b) 

    ### End galactic plane plot

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
    os.chmod(initial_skymap_plot, 0o0777)

    # Clear the current axes.
    plt.cla() 
    # Clear the current figure.
    plt.clf() 
    # Closes all the figure windows.
    plt.close('all')

    return initial_skymap_plot, moon_plot
    # return  initial_skymap_plot, moon_plot, galacticCenterlineCoords,maxprob_ra, maxprob_dec 
