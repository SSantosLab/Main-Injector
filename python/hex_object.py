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
    def __init__(self, ra, dec, prob, rise_mjd, set_mjd, expTime, filt, detP, coverage = 0, camera='decam', index=None):
        self.prob = prob
        self.rise_mjd = rise_mjd
        self.set_mjd = set_mjd
        self.ra = ra
        self.dec = dec
        self.detP = detP # tuple of (detP_1,detP_2), where detP_1 is first pass and detP_2 is second pass
        
        if camera == "decam":
            self.lat = -30.16527778
            self.lon = -70.8125
            self.alt = 2215
        elif camera == "desi":
            self.lat = 31.9600784
            self.lon = -111.598169
            self.alt = 2067
        elif camera == "hsc":
            self.lat = 19.826667
            self.lon = -155.471667
            self.alt = 4215
        else:
            print("Camera not supported")
        self.limit_mjd = self.getUnobservableTime() # Time when hex reaches telescope airmass limit
        self.index = index
        self.awesomeness_factor = None
        self.coverage = 0
        
        base_hex_coverage = 0.8522 #percent of a single hex covered without a dither (ie at 0, 0)
        self.coverage_factor = base_hex_coverage #percent of hex that would be covered initially
        
        self.dither = [0.00, 0.00] #ra and dec of dither on this current observation; updates if observed again at a new dither

#for testing
        self.awesomeness_factor_list = []
        self.airmass_Factor = []
        self.hex_VisibilityFactor = []
        self.slewTime_Factor = []
        self.lunarSeparation_Factor = []
        self.mjd_list = []
        self.coverage_list = []

        self.expTime = expTime
        self.filt = filt
        self.slewTime = 0
        
        
    
    def getAwesomenessFactor(self, mjd, last_ra, last_dec, moon_ra, moon_dec, moon_zd, dither = [0, 0], detailed = False):
        # Ranks hex based on attributes. Inputs are current time (MJD), current scope coordinates (deg), and moon coordinates (deg).
        if detailed:
            print(f'Probability: {self.prob} Airmass factor: {self.airmassFactor(mjd)} hexVisibility factor: {self.hexVisibilityFactor(mjd)}, slew time: {self.slewTimeFactor(last_ra, last_dec)}, lunar sep: {self.lunarSeparationFactor(moon_ra, moon_dec, moon_zd)}')
            
        self.airmass_Factor.append(self.airmassFactor(mjd))
        self.hex_VisibilityFactor.append(self.hexVisibilityFactor(mjd))
        self.slewTime_Factor.append(self.slewTimeFactor(last_ra, last_dec))
        self.lunarSeparation_Factor.append(self.lunarSeparationFactor(moon_ra, moon_dec, moon_zd))
        self.mjd_list.append(mjd)
        self.updateCoverageFactor()
        self.coverage_list.append(self.coverage_factor)
        
        AwesomenessFactor = 2*self.prob * self.airmassFactor(mjd) * self.hexVisibilityFactor(mjd) * \
                self.slewTimeFactor(last_ra, last_dec) * self.lunarSeparationFactor(moon_ra, moon_dec, moon_zd) * self.coverage_factor
        
        self.awesomeness_factor = AwesomenessFactor
        self.awesomeness_factor_list.append(AwesomenessFactor)
        
    def airmassFactor(self, mjd):
        # Rates the hex based on the airmass of the current position.
        if (mjd<self.rise_mjd) or (mjd>self.set_mjd):
            # Hex not in sky
            return 0.
        
        airmass = self.getAirmass(mjd)
        
        if airmass<=1.6 and airmass >= 1: # Airmass of 1.6 is where we wanted decline to become more steep
            return 2 - airmass
        elif airmass <= 1.8: # Airmass of 1.8 is DECam limit
            return 3.6 - 2*airmass
        elif airmass<1:
            return 1.
        else:
            return 0.
    
    def getAirmass(self, mjd):
        # Calculate the airmass X using the formula from Rozenberg (1966), with an airmass value of 40 at the horizon
        hex_ha = mjd_to_LST(mjd, self.lon) - self.ra # Calculate the hour angle for the hex
        hex_zd = zenith_distance(hex_ha, self.dec, self.lat) # Calculate the zenith distance for the hex
        X = (np.cos(np.deg2rad(hex_zd)) + 0.025*np.exp(-11*np.cos(np.deg2rad(hex_zd))))**-1
        return X
    
    def getUnobservableTime(self):
        # Gets the MJD at which the airmass for this hex reaches 1.8. If the hex is always below above 1.8, then this returns None.
        t = np.linspace(self.rise_mjd+(self.set_mjd-self.rise_mjd)/2, self.set_mjd, 100)
        for time in t:
            if self.getAirmass(time)>=1.8:
                return time
        return self.set_mjd

    def hexVisibilityFactor(self, mjd):
        # Favors hexes that are close to airmass limit and deweights hexes that are not visible. Input is current time (MJD).
        if (mjd < self.rise_mjd) or (mjd > self.set_mjd):
            return 0. # Hex not available if not in sky
        if self.limit_mjd==None:
            return 0.
        else:
            return self.visFactorFunc(self.limit_mjd - mjd)

    def visFactorFunc(self, time, a=1.1, tmin=1/48., pmin=0.1):
        # Input is time until an event. Increases until event-tmin, at which point this function starts to return 1.
        c = tmin**a
        tmax = tmin * pmin**(-1./a)
        if time <= tmin:
            return 1.
        return c * time**(-a)
    
    def observe_hex(self):
        #update coverage if hex has been observed
        if self.dither == [0.00, 0.00]:
            #update dither to be next most coverage
            self.dither = [0.06389, 0.287436]
        elif self.dither == [0.06389, 0.287436]:
            self.dither = [0.0484, -0.6725]
            self.coverage = self.coverage_factor
        elif self.dither == [0.0484, -0.6725]:
            self.dither = [-0.2395875, -0.135264]
            self.coverage = self.coverage_factor
        elif self.dither == [-0.2395875, -0.135264]:
            self.dither = [0.9423775, -0.405792]
            self.coverage = self.coverage_factor
        elif self.dither == [0.9423775, 0.405792]:
            self.dither = [-0.76668, 0.473424]
            self.coverage = self.coverage_factor
        elif self.dither == [-0.76668, 0.473424]:
            self.dither = [-0.543065, -0.828492]
            self.coverage = self.coverage_factor
        elif self.dither == [-0.543065, -0.828492]:
            self.dither = [-0.0479175, 0.388884]
            self.coverage = self.coverage_factor
    
    def updateCoverageFactor(self):
        #determine sky coverage factor for a given dither. this number is additional percent of sky covered, which will be multiplied by hex probability in awesomeness factor calculations
        if self.dither == [0.00, 0.00]:
            #if hex has not been covered, use sky coverage from no dither
            pass
        elif self.dither == [0.06389, 0.287436]:
            added_coverage = 0.111371 #coverage from dither [-0.2395875, -0.135264], which adds the most probability space
            
            self.coverage = self.coverage_factor
            self.coverage_factor = added_coverage
        elif self.dither == [0.0484, -0.6725]:
            added_coverage = 0.022066
            self.coverage = self.coverage_factor
            self.coverage_factor = added_coverage
        elif self.dither == [-0.2395875, -0.135264]:
            added_coverage = 0.007529
            self.coverage = self.coverage_factor
            self.coverage_factor = added_coverage
        elif self.dither == [0.9423775, -0.405792]:
            added_coverage = 0.004673
            self.coverage = self.coverage_factor
            self.coverage_factor = added_coverage  
        elif self.dither == [-0.543065, -0.828492]:
            added_coverage = 0.001298
            self.coverage = self.coverage_factor
            self.coverage_factor = added_coverage   
        elif self.dither == [-0.76668, 0.473424]:
            added_coverage = 0.000519
            self.coverage = self.coverage_factor
            self.coverage_factor = added_coverage  
        else:
            added_coverage = 0.0 #coverage from dither [-0.2395875, -0.135264], which adds the most probability space
#         else:
#             added_coverage = 0.5 #make it bigger to prioritize dithers more
            self.coverage_factor = added_coverage
    
            
        

    def slewTimeFactor(self, last_ra, last_dec):
        # Find the slew time from current scope position to this hex. Inputs must be current scope position in degrees.
        hex_sep = self.getAngSep(self.ra, self.dec, last_ra, last_dec)
        t_slew = hex_sep / 1. # DECam slews at 1 degree per second. This line is for clarity.
        self.slewTime = t_slew
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
