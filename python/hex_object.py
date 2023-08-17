import datetime
from pyslalib import slalib
import numpy as np

def zenith_distance(ha, dec, latitude):
    # Input and Output in degrees
    sinAltRad = np.sin(np.deg2rad(latitude))*np.sin(np.deg2rad(dec)) + np.cos(np.deg2rad(latitude))*np.cos(np.deg2rad(dec))*np.cos(np.deg2rad(ha))
    altRad = np.arcsin(sinAltRad)
    zenithDist = np.pi/2 - altRad
    return np.rad2deg(zenithDist)

def mjd_to_LST(mjd, eastLongitude):
    # Input and Output in degrees
    gmst = slalib.sla_gmst(mjd)
    eqEquinoxes = slalib.sla_eqeqx(mjd)
    lst = np.rad2deg(gmst) + np.rad2deg(eqEquinoxes) + eastLongitude
    return lst

def get_moon_pos(mjd, camera='decam'):
    # Coordinates in degrees
    if camera == "decam":
        lat = -30.16527778
        lon = -70.8125
    elif camera == "desi":
        lat = 31.9600784
        lon = -111.598169
    elif camera == "hsc":
        lat = 19.826667
        lon = -155.471667
    else:
        print("Camera not supported")
    
    # These should all be in degrees
    moon_ra_rad, moon_dec_rad, moon_dia_rad = slalib.sla_rdplan(mjd, 3, np.deg2rad(lon), np.deg2rad(lat))
    moon_ha = mjd_to_LST(mjd, lon) - np.rad2deg(moon_ra_rad)
    moon_zd = zenith_distance(moon_ha, np.rad2deg(moon_dec_rad), lat)
    return np.rad2deg(moon_ra_rad), np.rad2deg(moon_dec_rad), moon_ha, moon_zd

class HexObject:
    """
    Adapted heavily from Jim's mags.py code (at least for calculations of lunar separation). 

    Internally all coordinates, time and space, are in radians, except MJD, but outputs are in degrees.

    Attributes:
    -----------
    ra: float
        Right Accesion of hex (deg)
    dec: float
        Declination of hex (deg)
    prob: float
        Given hex's probability.
    rise_mjd: float
        Modified Julian Date of hex rise time.
    set_mjd: float
        Modified Julian Date of hex set time.
    """
    def __init__(self, ra, dec, prob, rise_mjd, set_mjd, expTime, index = None):
        self.prob = prob
        self.rise_mjd = rise_mjd
        self.set_mjd = set_mjd
        self.ra = ra
        self.dec = dec
        self.index = index
        self.expTime = expTime
        self.awesomeness_factor = None
        
#for testing
        self.awesomeness_factor_list = []
        self.airmass_Factor = []
        self.hex_VisibilityFactor = []
        self.slewTime_Factor = []
        self.lunarSeparation_Factor = []
        self.mjd_list = []

    def getAwesomenessFactor(self, mjd, last_ra, last_dec, moon_ra, moon_dec, moon_zd, detailed = False):
        # Ranks hex based on attributes. Inputs are current time (MJD), current scope coordinates (deg), and moon coordinates (deg).
        if detailed:
            print(f'Probability: {self.prob} Airmass factor: {self.airmassFactor(mjd)} hexVisibility factor: {self.hexVisibilityFactor(mjd)}, slew time: {self.slewTimeFactor(last_ra, last_dec)}, lunar sep: {self.lunarSeparationFactor(moon_ra, moon_dec, moon_zd)}')
        self.airmass_Factor.append(self.airmassFactor(mjd))
        self.hex_VisibilityFactor.append(self.hexVisibilityFactor(mjd))
        self.slewTime_Factor.append(self.slewTimeFactor(last_ra, last_dec))
        self.lunarSeparation_Factor.append(self.lunarSeparationFactor(moon_ra, moon_dec, moon_zd))
        self.mjd_list.append(mjd)
        
        AwesomenessFactor = self.prob * self.airmassFactor(mjd) * self.hexVisibilityFactor(mjd) * \
                self.slewTimeFactor(last_ra, last_dec) * self.lunarSeparationFactor(moon_ra, moon_dec, moon_zd)
        self.awesomeness_factor = AwesomenessFactor
        self.awesomeness_factor_list.append(AwesomenessFactor)
        
    def airmassFactor(self, mjd, camera='decam'):
        # Coordinates in degrees
	if (mjd<self.rise_mjd) or (mjd>self.set_mjd):
	    # hex not in sky
	    return 0.
	
        if camera == "decam":
            lat = -30.16527778
            lon = -70.8125
            alt = 2215
        elif camera == "desi":
            lat = 31.9600784
            lon = -111.598169
            alt = 2067
        elif camera == "hsc":
            lat = 19.826667
            lon = -155.471667
            alt = 4215
        else:
            print("Camera not supported")

        # Calculate the hour angle for the hex
        hex_ha = mjd_to_LST(mjd, lon) - self.ra

        # Calculate the zenith distance for the hex
        hex_zd = zenith_distance(hex_ha, self.dec, lat)

        # Calculate the airmass X using the formula from Rozenberg (1966), with an airmass value of 40 at the horizon
        X = (np.cos(np.deg2rad(hex_zd)) + 0.025*np.exp(-11*np.cos(np.deg2rad(hex_zd))))**-1

        # Calculate the airmass factor. Airmass = 1 (source at zenith) is given a score of 1. A source with a ~50 degree
        # zenith angle (airmass = 1.6) is given a score of 0.7. Sources at zenith angles greater than this are given a score
        # which decays exponentially, with a 56 degree zenith angle (airmass ~ 1.8: limit of camera) given a score of ~0.21
        # and all airmasses above 2.3 getting a default score of 0.01
        if 0.9<=X<=1.6:
            return 1 - 0.5*(X-1)
        elif 1.6<X<=2.3:
            return 0.7*np.exp(-6*(X-1.6))
        elif 2.3<X<=40: # Airmass of 2.3 corresponds to a source 26 degs from the horizon
            return 0.01
        else:
            # Not a valid airmass
            return 0.
 
    def hexVisibilityFactor(self, mjd):
        # Favors hexes that are close to setting and deweights hexes that are not visible. Input is current time (MJD).
        if (mjd < self.rise_mjd) or (mjd > self.set_mjd):
            return 0. # Hex not available if not in sky
        time_to_set = self.set_mjd - mjd
        return self.visFactorFunc(time_to_set)

    def visFactorFunc(self, x, a=1.83, xmin=1./24, pmin=0.05):
        # Tuned so that the median is at 2 hours to set. Plateau value defaults to 1 hour to set. Input is time to set in days.
        c = xmin**a
        xmax = xmin * pmin**(-1./a)
        if x <= xmin:
            return 1.
        if x >= xmax:
            return pmin
        return c * x**(-a)
        
    def slewTimeFactor(self, last_ra, last_dec):
        # Find the slew time from current scope position to this hex. Inputs must be current scope position in degrees.
        hex_sep = self.getAngSep(self.ra, self.dec, last_ra, last_dec)
        t_slew = hex_sep / 1. # DECam slews at 1 degree per second. This line is for clarity.
        return self.slewFactorFunc(t_slew)
        
    def slewFactorFunc(self, x, a=0.188, xmin=0.5):
        # Numbers tuned such that the factor at 20s slew time is 0.5. Input is slew time in seconds.
        if x <= xmin:
            return 1.
        c = xmin**a
        return c * x**(-a)
        
    def lunarSeparationFactor(self, moon_ra, moon_dec, moon_zd):
        # Inputs must be moon RA, Dec, and Zenith Distance in degrees.
        moon_sep = self.getAngSep(self.ra, self.dec, moon_ra, moon_dec)
        if moon_zd > 90:
            return 1. # If the moon is set, then it shouldn't be a problem. Might need to tweak horizon degree value.
        elif moon_sep > 15:
            return 1.
        else:
            return 0.

    def getAngSep(self, ra1, dec1, ra2, dec2):
        # Combined the haversine functions to be more coherent. Inputs must be in degrees.
        d2r = np.pi/180
        a1 = ra1 * d2r
        d1 = dec1 * d2r
        a2 = ra2 * d2r
        d2 = dec2 * d2r
        return 2 * np.arcsin(np.sqrt(np.sin((d2-d1)/2)**2 + np.cos(d1)*np.cos(d2)*np.sin((a2-a1)/2)**2)) / d2r

    def getRiseTimeReadable(self, with_date=False):
        mjd_integer = int(self.rise_mjd)
        mjd_decimal = self.rise_mjd - mjd_integer

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

    def getSetTimeReadable(self, with_date=False):
        mjd_integer = int(self.set_mjd)
        mjd_decimal = self.set_mjd - mjd_integer

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
      
    def display_hexinfo(self):
        print(f'RA: {self.ra:.3f}, DEC: {self.dec:0.3}, PROB: {self.prob:0.3}, AWESOMENESS: {self.awesomeness_factor:0.3}')
