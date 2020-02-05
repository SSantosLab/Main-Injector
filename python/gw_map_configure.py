"""
Form the classes for handling the many different
inputs to the map making code.
"""
import fitsio
import numpy as np
import hp2np
import healpy as hp

class trigger(object):
    """
    """
    def __init__(self, skymap, trigger_id, trigger_type, resolution, days_since_burst=0.) :
        """
        """
        self.skymap = skymap
        self.trigger_id = trigger_id
        self.trigger_type = trigger_type
        self.days_since_burst = days_since_burst

        hdr = fitsio.read_header(skymap,1)
        burst_mjd = np.float(hdr["mjd-obs"])
        try:
            distance = hdr["distmean"]
        except:
            print('distmean was not in payload... setting distance to 60mpc')
            distance = 60
        self.distance  = distance
        self.burst_mjd = burst_mjd

        print "\ngw_map_trigger: read map {}, {} at {:.2f} Mpc\n".format(
            skymap, trigger_type, distance)

        ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)
        ligo_dist, ligo_dist_sig, ligo_dist_norm  = \
            distance*np.ones(ra.size), np.zeros(ra.size), np.zeros(ra.size)
        try :
            ligo_dist      = hp.read_map(skymap, field=1, verbose=False)
            ligo_dist_sig  = hp.read_map(skymap, field=2, verbose=False)
            ligo_dist_norm = hp.read_map(skymap, field=3, verbose=False)
            # there are infinite distance pixels in the map. deal with these
            ligo_dist[np.isinf(ligo_dist)] = 10000.

            # change resolution
            ligo_dist      = hp.ud_grade(ligo_dist, resolution)
            ligo_dist_sig  = hp.ud_grade(ligo_dist, resolution)
            ligo_dist_norm = hp.ud_grade(ligo_dist, resolution)
        except:
            print "\n\t !!!!!!!! ------- no distance information in skymap ------ !!!!!!!!\n"
    
        # JTA Hack
        # GW190425
        #print "HACK HACK HACK"
        #ix = ((ra > 150) | (( ra > -180) & (ra < -120)) ) & (dec > -10) & (dec < 40)
        #ix = np.invert(ix)
        #ligo[ix] = 0.0
        #print "HACK HACK HACK"
    
        # GW170217 hack JTA
        #ix = (ra > 0) & ( ra < 180) & (dec >= -30)
        #ix = np.invert(ix)
        #ligo[ix] = 0.0
        # GW170225 hack JTA
        #ix = (dec >= 2)
        #ligo[ix] = 0.0
        # GW170814 hack JTA
        #ix = (ra > -10) & ( ra < 60) & (dec < -20)
        #ix = np.invert(ix)
        #ligo[ix] = 0.0
        self.ligo_ra = ra
        self.ligo_dec = dec
        self.ligo = ligo
        self.ligo_dist = ligo_dist
        self.ligo_dist_sig = ligo_dist_sig
        self.ligo_norm = ligo_dist_norm



class strategy(object) :
    """
    """
    def __init__(self, camera, exposure_list, filter_list, maxHexesPerSlot, hoursAvailable, propid):
        """
        """
        self.camera = camera
        self.exposure_list = exposure_list
        self.filter_list = filter_list
        self.maxHexesPerSlot = maxHexesPerSlot
        self.hoursAvailable = hoursAvailable
        self.propid = propid

        if camera == "decam" :
            self.overhead =  30. # seconds
            self.area_per_hex = 3.0 # sq-degrees
        elif camera == "hsc" :
            self.overhead = 20.
            self.area_per_hex = 1.5
        else: raise Exception("camera {} not handled".format(camera))

class control(object):
    """
    """
    def __init__(self, resolution, data_dir, debug=False, allSky=False,
            snarf_mi_maps=False, mi_map_dir="/data/des41.a/data/desgw/O3FULL/Main-Injector/OUTPUT/O3REAL/") :
        """
        """
        self.resolution = resolution
        self.debug = debug
        self.this_tiling = []
        self.reject_hexes= []
        self.allSky = allSky
        self.datadir = data_dir
        # find hexes starting with slot start_slot. Carries on to jsons
        self.start_slot = -1
        # find hexes for the do_nslots starting with slot start_slot. Carries on to jsons
        self.do_nslots = -1
        self.just_sort_by_ra = False # if true, take the calculated hexes and just sort
                                    # the time to be observed as increasing ra
        # should i bypass the mapmaking by copying over from the MI directories?
        self.snarf_mi_maps = snarf_mi_maps
        # what directory to holds the existing maps? Copy over to output_dir
        self.mi_map_dir = mi_map_dir 
 

class results(object):
    """
    """
    def __init__(self) :
        """
        """
        self.probability_per_slot = False
        self.time_of_slot = False
        self.isdark = False
        self.slotDuration = False
        self.hoursPerNight = False
        self.n_slots = False
        self.first_slot = False
        self.best_slot = False

