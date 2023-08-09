import numpy as np
import hp2np 
import hexalate
import decam2hp
import jsonMaker
from os import getenv
import pandas as pd
import sys
import  matplotlib.pyplot as plt
import mags
from pyslalib import slalib
from enum import Enum

CTIO_LAT = -30.16527778
CTIO_LON = -70.8125
CTIO_HEIGHT = 2215.
KPNO_LAT = 31.9600784
KPNO_LON = -111.598169
KPNO_HEIGHT = 2067.
MAUNAKEA_LAT = 19.826667
MAUNAKEA_LON = -155.471667
MAUNAKEA_HEIGHT = 4215.

class CameraType(Enum):
    DECam = "decam"
    DESI = "desi"
    HSC = "hsc"
    
class hex_object:
    """
    Adapted heavily from Jim's mags.py code (at least for calculations of lunar separation). 

    Internally all coordinates, time and space, are in radians, except MJD, but outputs are in degrees.

    Attributes:
    -----------
    ra: float
        Right Accesion of hex. 
    dec: float
        Declination of hex.
    prob: float
        Given hex's probability.
    rise_est: string
        EST and date of hex's rise time.
    set_est: string
        EST and date of hex's set time. 
    rise_mjd: float
        Modified Julian Date of hex rise time.
    set_mjd: float
        Modified Julian Date of hex set time.
    camera: str (Default: "decam")
        Camera type for observations. 
    """
    def __init__(self, ra, dec, prob, rise_est, set_est, rise_mjd, set_mjd, camera="decam"):
        self.prob = prob
        self.rise_est = rise_est
        self.set_est = set_est
        self.rise_mjd = rise_mjd
        self.set_mjd = set_mjd

        #get conversion factors
        self.degToRad = 2. * np.pi / 360.
        self.radToDeg = 360 / 2 / np.pi

        #determine whether using decam, desi, or hsc (not sure if necessary, jim's mags just does this 
        #so I put it in 
        self.observatory = CameraType(camera)

        #calculations for moon separation/distance from moon
        self.lat = self.obs_lat * self.degToRad
        self.lon = self.obs_lon * self.degToRad
        self.height = self.obs_height

        self.ra = ra * self.degToRad
        self.dec = dec * self.degToRad

        #get moon separation at rise time of hex and set time of hex, saved in
        #self.rise_moon_sep and self.set_moon_sep, respectively 
        self.calculate_lunar_separations()

#         print(self.rise_moon_sep, self.set_moon_sep)


###WILL WANT TO FILL IN CODE HERE####
    def airmass_factor(self, time):
        #airmass factor as a function of time and position 
    def distance_from_hex(idk):
        #get distance from hex to another hex
    def slew_time(distance_from_hex i think):
        #get slew time
    def time_to_set(self, current_time):
        return current_time - self.set_mjd
    def get_awesomeness_factor(self):
        #get the awesomeness factor filling in previous factors here
        
#####REST OF THIS CODE THEORETICALLY KINDA DONE#########

    def calculate_lunar_separations(self):
            self.resetTime(self.rise_mjd)
            rise_moon_sep = self.radToDeg * self.getLunarSeparation(self.ra, self.dec, self.moonData)

            self.resetTime(self.set_mjd)
            set_moon_sep = self.radToDeg * self.getLunarSeparation(self.ra, self.dec, self.moonData)

            self.rise_moon_sep = rise_moon_sep
            self.set_moon_sep = set_moon_sep

    def getLunarPosition(self, mjd, eastLongitude, latitude, lst):
        ra, dec, diam = slalib.sla_rdplan(mjd, 3, eastLongitude, latitude)
        ha = lst - ra
        zd = self.zenithDistance(ha, dec, latitude)
        return ra, dec, zd
    
    def getLunarSeparation(self, ra, dec, moonData):
        moon_ra, moon_dec = self.moonData[0], self.moonData[1]
        moon_sep = self.gc_separation(ra, dec, moon_ra, moon_dec)
        return moon_sep

    def obs_lat(self):
        if self.observatory == CameraType.DECam or self.observatory == CameraType.DESI:
            return CTIO_LAT
        elif self.observatory == CameraType.HSC:
            return MAUNAKEA_LAT

    def obs_lon(self):
        if self.observatory == CameraType.DECam or self.observatory == CameraType.DESI:
            return CTIO_LON
        elif self.observatory == CameraType.HSC:
            return MAUNAKEA_LON
        
    def obs_height(self):
        if self.observatory == CameraType.DECam or self.observatory == CameraType.DESI:
            return CTIO_HEIGHT
        elif self.observatory == CameraType.HSC:
            return MAUNAKEA_HEIGHT
        
    def mjdToLST(self, mjd, eastLongitude):
        gmst = slalib.sla_gmst(mjd)
        eqEquinoxes = slalib.sla_eqeqx(mjd)
        lst = gmst + eqEquinoxes + eastLongitude
        return lst
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
    def equatorialToObservational(self, ra, dec, lst, latitude):
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
        
    def zenithDistance(self, ha, dec, latitude):
        degToRad = 2.*np.pi/360.
        sinAltRad = np.sin(latitude)*np.sin(dec) + \
            np.cos(latitude)*np.cos(dec)*np.cos(ha)
        altRad = np.arcsin(sinAltRad)
        zenithDist = 90*degToRad - altRad
        return zenithDist
    
    def resetTime(self, mjd):
        self.mjd = mjd
        self.lst = self.mjdToLST(self.mjd, self.lon)
        self.ha, self.zd = self.equatorialToObservational(
            self.ra, self.dec, self.lst, self.lat)
        self.moonData = self.getLunarPosition(
            self.mjd, self.lon, self.lat, self.lst)
        self.moonZD = self.moonData[2]
        self.moonSep = self.getLunarSeparation(
            self.ra, self.dec, self.moonData)
        self.moonRa = self.moonData[0]
        self.moonDec = self.moonData[1]