from typing import Tuple
import os
import sys
import numpy as np
import healpy as hp
from pyslalib import slalib

import hp2np
import atmosphere
import telescope
import dust
import seeingModel
import sky_model
import mcbryde


license = """
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style. (Note the rate at which people who stand on the
   summit of K2 successfully make it down.)

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""


class observed(object):
    """
    Constuct the observed map.
    This produces the limiting magnitude map for a 30 sec i-band
    exposure using DECam/Blanco over the whole sky given a time in MJD.

    Internally all coordinates, time and space, are in radians, except MJD.

    Attributes:
    -----------
    ra: float
        Right Accesion to observed class instance.
    dec: float
        Declination of observed class instance.
    values: float
    mjd: float
        Modified Julian Date of observed class instance.
    alpha: float (Default: 0.)
    doMaps: bool (Default: True)
        If True, make skymap of observed class instance.
    camera: str (Default: "decam")
        Camera type for which observed class instance.
    verbose: bool (Default: True)
        If True, givens output to terminel.
    """

    def __init__(self,
                 ra,
                 dec,
                 values,
                 mjd,
                 alpha=0.,
                 degradeRes=True,  # change map resolution
                 doMaps=True,  # don't think of ra,dec,vals as a healpy map; don't do map work
                 camera="decam",
                 verbose=True,   # be very loud or no
                 single_obj = False # be even less loud or no, added by Elise
                 ):

        self.verbose = verbose
        self.single_obj = single_obj
        print(self.single_obj)
        data_dir = os.environ["DATA_DIR"]
        ctio_lat = -30.16527778
        ctio_lon = -70.8125
        ctio_height = 2215.
        kpno_lat = 31.9600784
        kpno_lon = -111.598169
        kpno_height = 2067.
        maunakea_lat = 19.826667
        maunakea_lon = -155.471667
        maunakea_height = 4215.

        if camera == "decam" or camera == "des":
            obs_lat = ctio_lat
            obs_lon = ctio_lon
            obs_height = ctio_height
        elif camera == "desi":
            obs_lat = kpno_lat
            obs_lon = kpno_lon
            obs_height = kpno_height
        elif camera == "hsc":
            obs_lat = maunakea_lat
            obs_lon = maunakea_lon
            obs_height = maunakea_height
        self.observatory = camera

        self.date = slalib.sla_djcl(mjd)

        self.degToRad = 0.0174532925199
        self.lat = obs_lat*self.degToRad
        self.lon = obs_lon*self.degToRad
        self.height = obs_height

        self.mjd = mjd
        self.ra = ra*self.degToRad
        self.dec = dec*self.degToRad
        self.map = values
        if doMaps:
            self.nside = hp.npix2nside(values.size)
        else:
            self.nside = 0

        # observational astronomy
        self.lst = self.mjdToLST(self.mjd, self.lon)
        self.ha, self.zd = self.equatorialToObservational(
            self.ra, self.dec, self.lst, self.lat)
        self.airmass = self.airmassModel(self.zd)
        self.sunData = self.getSolarPosition(
            self.mjd, self.lon, self.lat, self.lst)
        self.sunZD = self.sunData[2]
        self.moonData = self.getLunarPosition(
            self.mjd, self.lon, self.lat, self.lst)
        self.moonZD = self.moonData[2]
        self.moonSep = self.getLunarSeparation(
            self.ra, self.dec, self.moonData)
        self.moonPhase = self.getLunarPhase()
        self.moonRa = self.moonData[0]
        self.moonDec = self.moonData[1]

        # telescope limits
        self.limits = 0.0

        # galactic dust
        dust_dir = data_dir
        dust_file = "plank-ebv-HFI_CompMap_ThermaDustModel.fits"
        ra, dec, ebv = dust.loadDust(dust_dir, dust_file)
        if degradeRes:
            ra, dec, ebv = hp2np.map2np(
                ebv, resolution=self.nside, fluxConservation=False)
        if not doMaps:
            ra = self.ra/self.degToRad
            dec = self.dec/self.degToRad
            ebv = hp.get_interp_val(ebv, ra, dec, lonlat=True)
            ra = self.ra
            dec = self.dec
        self.dust_ra = ra
        self.dust_dec = dec
        self.ebv = ebv

        # stellar density in the form of a probability of recognizing object map
        if not self.single_obj:
            print(
                "loading inverse stellar density = probability of recognition map"
            )

        ra, dec, precog = hp2np.hp2np(os.path.join(data_dir, "precognize.fits"))

        if degradeRes:
            ra, dec, precog = hp2np.map2np(
                precog,
                resolution=self.nside,
                fluxConservation=False
            )
        if not doMaps:
            ra = self.ra/self.degToRad
            dec = self.dec/self.degToRad
            precog = hp.get_interp_val(precog, ra, dec, lonlat=True)
            ra = self.ra
            dec = self.dec
        self.pra = ra
        self.pdec = dec
        self.precog = precog

        # limiting mag
        self.maglim = ""
        self.maglimall = ""

        # where to put the central meridan in the  mcbryde-thomas projection
        # - the ra which one wants on the meridian
        self.mcbryde_alpha = alpha
        # set the figure axes so that the aspect ratio is set by data
        # for the right Mcbryde plot, which projects onto a flat x,y page
        # plt.axes().set_aspect('equal')

    # change to a new time
    def resetTime(self, mjd):
        self.mjd = mjd
        self.lst = self.mjdToLST(self.mjd, self.lon)
        self.ha, self.zd = self.equatorialToObservational(
            self.ra, self.dec, self.lst, self.lat)
        self.airmass = self.airmassModel(self.zd)
        self.sunData = self.getSolarPosition(
            self.mjd, self.lon, self.lat, self.lst)
        self.sunZD = self.sunData[2]
        self.moonData = self.getLunarPosition(
            self.mjd, self.lon, self.lat, self.lst)
        self.moonZD = self.moonData[2]
        self.moonSep = self.getLunarSeparation(
            self.ra, self.dec, self.moonData)
        self.moonRa = self.moonData[0]
        self.moonDec = self.moonData[1]

    def limitMag(self, filt, exposure=30, verbose=True):
        # print "\t JTA limiting magnitude", self.mjd, self.mjdToLST(self.mjd, self.lon)
        # print self.ha[int(self.ha.size/2.)]
        mjd = self.mjd
        zd = self.zd
        sun_zd = self.sunZD
        moon_zd = self.moonZD
        moon_phase = self.moonPhase
        moon_sep = self.moonSep
        airmass = self.airmass
        ebv = self.ebv
        alpha = self.mcbryde_alpha
        sunIsUp = self.sunBrightnessModel(sun_zd)
        if sunIsUp:
            m = np.copy(ebv)*0.1+-10.0
            mglobal = np.copy(ebv)*0.1+-10.0
            if not self.single_obj:
                print("\t ... the sun is up")
        else:
            dust = self.dustTransmission(filt, ebv)
            atmo = self.atmosphereTransmission(zd, airmass, filt, moon_sep)
            seeing = self.seeing(airmass, filt, seeingAtZenith=0.9)
            sky = self.skyBrightness(zd, moon_zd, moon_sep, moon_phase, filt)
            skyFid = self.skyBrightnessFid(filt)

            telescope = self.telescopeLimits(self.ha, self.dec)
            self.limits = telescope

            if filt == "g":
                m_zp = 23.3
            elif filt == "r":
                m_zp = 23.4
            elif filt == "i":
                m_zp = 22.9
            elif filt == "z":
                m_zp = 22.5
            elif filt == "y":
                m_zp = 20.6
            else:
                raise Exception("no such filt")

            SN = telescope*dust*atmo*(1./seeing)*np.sqrt(exposure/30.)
            ix = np.nonzero(SN <= 0)
            SN[ix] = 1.0e-12
            m = m_zp + 2.5*np.log10(SN) + 0.5*(sky - skyFid)
            # arbitarily limit limiting mag to that of moon
            ix = np.nonzero(m < 0)
            m[ix] = 0
            ix, = np.where(telescope == 0)
            m[ix] = -11.0

            SNglobal = dust*atmo*(1./seeing)*np.sqrt(exposure/30.)
            ix = np.nonzero(SNglobal <= 0)
            SNglobal[ix] = 1.0e-12
            mglobal = m_zp + 2.5*np.log10(SNglobal) + 0.5*(sky - skyFid)
            # arbitarily limit limiting mag to that of moon
            ix = np.nonzero(mglobal < 0)
            mglobal[ix] = 0

        # project into equal area map
        x, y = mcbryde.mcbryde(self.ra/self.degToRad,
                               self.dec/self.degToRad, alpha=alpha)
        self.x = x
        self.y = y
        hx, hy = mcbryde.mcbryde(
            self.ha/self.degToRad, self.dec/self.degToRad, alpha=alpha)
        self.hx = hx
        self.hy = hy
        self.maglim = m
        self.maglimall = mglobal

    def dustTransmission(self, filt, ebv):
        if self.verbose:
            print("\t ... dust")
        dustTransmission = dust.dustTransmission(filt, ebv)
        return dustTransmission

    def atmosphereTransmission(self, zd, airmass,
                               filt, moon_sep, refAirmass=1.3):
        if self.verbose:
            print("\t ... atmosphere")
        atransmission = atmosphere.transmission(airmass, filt, refAirmass)
        if self.verbose:
            print("\t ... earth")
        dtransmission = atmosphere.dirtTransmission(zd)
        if self.verbose:
            print("\t ... moon")
        mtransmission = atmosphere.lunarDirtTransmission(moon_sep)
        transmission = atransmission*dtransmission*mtransmission
        return transmission

    def telescopeLimits(self, ha, dec):
        if self.verbose:
            print("\t ... telescope limits")
        if self.observatory == "decam" or self.observatory == "des":
            #limits = telescope.blancoLimits(ha, dec)
            # print "HACK HACK HACK"
            limits = np.zeros(ha.size)
            ix = self.zd*360./2/np.pi <= 67.5
            limits[ix] = 1.0
        elif self.observatory == "desi" or self.observatory == "hsc":
            limits = np.zeros(ha.size)
            ix = self.zd*360./2/np.pi <= 80.
            limits[ix] = 1.0
        return limits

    def seeing(self, airmass, wavelength=775., seeingAtZenith=1.0):
        if self.verbose:
            print("\t ... seeing")
        seeing = seeingModel.seeingWithAirmassAndLambda(
            airmass, wavelength, seeingAtZenith)
        return seeing

    def skyBrightness(self, zd, moon_zd, moon_sep, moon_phase, filt):
        if self.verbose:
            print("\t ... sky brightness")
        sky = sky_model.sky_brightness_at_time(
            filt, zd, moon_zd, moon_sep, moon_phase)
        return sky

    def skyBrightnessFid(self, filt):
        skyfiducial = sky_model.skyFiducial(filt)
        return skyfiducial

    #
    # First, let's turn RA Dec into local sideral time, the start of all observations
    #
    # equation of Equinoxes is an ~1 second effect
    # "east longitude", where positive numbers increase as one moves east
    #   recall that lat, long is in radians
    def mjdToLST(self, mjd, eastLongitude):
        if self.verbose:
            print("\t MJD to LST")
        gmst = slalib.sla_gmst(mjd)
        eqEquinoxes = slalib.sla_eqeqx(mjd)
        lst = gmst + eqEquinoxes + eastLongitude
        return lst

    #
    # Now we turn LST and dec into Hour Angle and zenith distance
    #
    # Technically this is apparent ha,dec; i.e. we ignore diurnal abberation
    #   recall that lat, long is in radians
    # The best way to view this output is as hexbin(ra,dec,zd) or ha
    def equatorialToObservational(self, ra, dec, lst, latitude):
        if self.verbose:
            print("\t LST to HA,zenith distance")
    #Elise's addition to help handling single hexes
        try:
            if len(self.ra) > 1:
                ha = lst - ra
                zd = self.zenithDistance(ha, dec, latitude)
                ix = np.nonzero(ha > 180.*2*np.pi/360.)
                ha[ix] = ha[ix] - 360*2*np.pi/360.
                return ha, zd
        
        except:
            ha = lst - self.ra
            zd = self.zenithDistance(ha, dec, latitude)
            return ha, zd

    #
    # Calculating zenith distance
    #
    # sin(HC) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(LHA)
    def zenithDistance(self, ha, dec, latitude):
        degToRad = 2.*np.pi/360.
        sinAltRad = np.sin(latitude)*np.sin(dec) + \
            np.cos(latitude)*np.cos(dec)*np.cos(ha)
        altRad = np.arcsin(sinAltRad)
        zenithDist = 90*degToRad - altRad
        return zenithDist

    #
    # Calculating airmass
    #
    # from Young 1994 "Air mass and refraction", Applied Optics 33:1108-1110
    # max error = 0.0037 airmasses at zd=90 degrees. (!)
    def airmassModel(self, zd):
        if self.verbose:
            print("\t zenith distance to airmass")
        airmass = zd*0+45  # a hack to set the values of airmass under the horizon
        ix = np.nonzero(zd <= 89.9999*2*np.pi/360.)
        coszd = np.cos(zd)
        coszdsq = coszd*coszd
        numerator = 1.002432*coszdsq + 0.148386*coszd + 0.0096467
        denominator = coszdsq*coszd + 0.149864*coszdsq + 0.0102963*coszd \
                      + 0.000303978
    # Another Elise addition so it can handle a single hex 
        try:
            airmass[ix] = numerator[ix]/denominator[ix]
        except:
            airmass = numerator/denominator
        return airmass

    #
    # check lunar position
    #

    def getLunarPosition(self, mjd, eastLongitude, latitude, lst):
        ra, dec, diam = slalib.sla_rdplan(mjd, 3, eastLongitude, latitude)
        ha = lst - ra
        zd = self.zenithDistance(ha, dec, latitude)
        return ra, dec, zd

    #
    # returns moon phase in degrees
        # K&S want moon phase as angle in degrees,
        #   where 0 = full, 90 equals half, and  180 = new
    #
    def getLunarPhase(self):
        moon_ra, moon_dec = self.moonData[0], self.moonData[1]
        sun_ra, sun_dec = self.sunData[0], self.sunData[1]
        # moon is full when elongation = 180, new when elongation = 0
        moon_elongation = self.gc_separation(
            sun_ra, sun_dec, moon_ra, moon_dec)
        # K&S want moon phase as angle in degrees,
        #   where 0 = full, 90 equals half, and  180 = new
        phase = (180*(2*np.pi/360.) - moon_elongation)*360./2/np.pi
        return phase

    def getLunarSeparation(self, ra, dec, moonData):
        moon_ra, moon_dec = self.moonData[0], self.moonData[1]
        moon_sep = self.gc_separation(ra, dec, moon_ra, moon_dec)
        return moon_sep

    #
    # check solar position
    #
    def getSolarPosition(self, mjd, eastLongitude, latitude, lst):
        ra, dec, diam = slalib.sla_rdplan(mjd, 0, eastLongitude, latitude)
        ha = lst - ra
        zd = self.zenithDistance(ha, dec, latitude)
        return ra, dec, zd
    #
    # check to see if sun near the horizen
    #

    def sunBrightnessModel(self, sunZD):
        twilight = 100.*2*np.pi/360.
        if sunZD <= twilight:
            bright = True
        else:
            bright = False
        return bright

    #
    # great circle separations
    #
    def gc_separation(self, ra1, dec1, ra2, dec2):
        delDec = dec1-dec2
        delRa = ra1-ra2
        dhav = self.haversine(delDec)
        rhav = self.haversine(delRa)
        hav = dhav + np.cos(dec1)*np.cos(dec2)*rhav
        gc_distance = self.ahaversine(hav)
        return gc_distance

    def haversine(self, theta):
        hav = np.sin(theta/2.)**2
        return hav

    def ahaversine(self, x):
        ahav = 2*np.arcsin(np.sqrt(x))
        return ahav


def findNightDuration(mjd: float,
                      camera="decam") -> Tuple[float, float, float]:
    """
    return
    """
    ctio_lat = -30.16527778
    ctio_lon = -70.8125
    ctio_height = 2215.
    kpno_lat = 31.9600784
    kpno_lon = -111.598169
    kpno_height = 2067.
    maunakea_lat = 19.826667
    maunakea_lon = -155.471667
    maunakea_height = 4215.
    if camera == "decam" or camera == "des":
        obs_lat = ctio_lat
        obs_lon = ctio_lon
        obs_height = ctio_height
    elif camera == "desi":
        obs_lat = kpno_lat
        obs_lon = kpno_lon
        obs_height = kpno_height
    elif camera == "hsc":
        obs_lat = maunakea_lat
        obs_lon = maunakea_lon
        obs_height = maunakea_height

    degToRad = 2.*np.pi/360.
    lat = obs_lat*degToRad
    lon = obs_lon*degToRad
    height = obs_height

    imjd = np.int64(mjd)
    start_mjd = imjd - 6./24.  # before sunset at CTIO

    sunset = ""
    sunrise = ""
    
    print(f'mjd is: {mjd}')
    
    #check every minute
    for i in np.arange(0, 1., 1./(24.*60)):
        mjd = start_mjd + i
        gmst = slalib.sla_gmst(mjd)
        eqEquinoxes = slalib.sla_eqeqx(mjd)
        lst = gmst + eqEquinoxes + lon

        sunra, sundec, diam = slalib.sla_rdplan(mjd, 0, lon, lat)
        sunha = lst - sunra
        sinAltRad = np.sin(lat)*np.sin(sundec) + \
            np.cos(lat)*np.cos(sundec)*np.cos(sunha)
        altRad = np.arcsin(sinAltRad)
        zenithDist = 90*degToRad - altRad

        twilight = 100.*2*np.pi/360.
        if zenithDist <= twilight:
            bright = True
        else:
            bright = False
        if sunset == "" and bright == False:
            sunset = mjd
        if sunset != "" and sunrise == "" and bright == True:
            sunrise = mjd
    print(sunrise, sunset)
    duration = float(sunrise)-float(sunset)
    return duration, sunset, sunrise
