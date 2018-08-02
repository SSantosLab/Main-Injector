import numpy as np
import matplotlib.path

# the intent is that these two main routines,
#   hexesOnMap & hexalateMap
# are widely useful ways to deal with DECam and DES
# data with healpix maps

# nside=32 ~ camera   3.36 sq-degree
# nside=256 ~ ccd     0.052 sq-degrees
# nside=512 resolves ccds, (47 sq-arcmin vs 165 sq-arcmin for a ccd)
# nsides = 1024 = 3.4x3.4 arcmin vs  9x18 arcmin

# keeps count of times a hex, any hex, overlies a map pixel
#   camera outline for the hexes given
def hexesOnMap(ra,dec, raHexen, decHexen, camera) :
    print "    checking {} hexes : ".format(raHexen.size),
    count = np.zeros(ra.size)
    for i in range(0,raHexen.size) :
        ix = radecInHex( raHexen[i], decHexen[i], ra, dec, camera)
        count[ix] += 1
    print " "
    return count


# I think doing this using a tree requires a projection we will
# use the Sanson-Flamsteed projection, aka sinusoidal projection (x=ra*cos(dec), y=dec)
def hexalateMap(ra, dec, vals, tree, raHexen, decHexen,  camera, verbose=1) :
    if verbose : print "\t hexalateMap \t nhex = {},".format(raHexen.size),
    if verbose: print " npix = {}".format(ra.size)
    hexVal = np.zeros(raHexen.size)
    counter = 0
    for i in range(0,raHexen.size) :
        ix = radecInHex( raHexen[i], decHexen[i], ra, dec, tree, camera)
        if not ix.size: continue
        try  :
            hexVal[i] = vals[ix].sum()
        except Exception: 
            print "why are there exceptions in hexalateMap?"
            hexVal[i] = 0
    return hexVal


# I think doing this using a tree requires a projection we will
# use the Sanson-Flamsteed projection, aka sinusoidal projection (x=ra*cos(dec), y=dec)
#
# return an index that is true/exists if ra,dec is inside
# the camera outline for the hex 

def radecInHex ( raCenter, decCenter, ra,dec,tree, camera,mapCenter= 0.) :
    if camera == 'decam':
        camera_radius = 1.1
    elif camera == 'hsc':
        camera_radius = 0.75
    else: 
        raise Exception('No such camera')
    radius = camera_radius + 0.1
    sf_scale = 2.6 *np.abs(decCenter)/90. ;# the distortion at high dec means we need wider "circle"
    radius += sf_scale
    data = tree.query_ball_point(((raCenter-mapCenter)*np.cos(decCenter*2*np.pi/360), decCenter), radius)
    near_ix = np.array(data)
    if not near_ix.size: return near_ix

    # find the answer
    if camera == 'decam':
        hex = hexPath(raCenter, decCenter, camera)
        inside_ix = radecInHexPath( hex, ra[near_ix], dec[near_ix])
    elif camera == 'hsc':
        inside_ix = radecInCircle(raCenter, decCenter, ra[near_ix], dec[near_ix], radius=camera_radius)
    near_ix = near_ix[inside_ix]

    return near_ix



#
# Use a circle of area pi sq-degrees as a reasonable
# approximation to the camera outline.
#
def radecInCircle (raCenter, decCenter, ra, dec, radius=1.0) :
    raDist  = (ra  - raCenter )
    decDist = (dec - decCenter)
    raDist  = raDist * np.cos(dec * 2 * np.pi/360.)
    distance = np.sqrt(raDist**2 + decDist**2)

    ix = distance <= radius

    return ix

def buildtree(ra,dec,nsides=1024,recompute=False, \
        wrapForMetricEval=False, dumpTree = False) :
    # I think doing this using a tree requires a projection we will
    # use a Sanon-Flamsteed projection (aka sinusoidal, x=ra*cos(dec), y=dec)
    import scipy.spatial
    import os.path
    import cPickle

    # zero is center of map (i.e, 180 is at singularity)
    ix = ra > 180; ra[ix]=ra[ix]-360.
    if wrapForMetricEval :
        # buffer each side of the singularity so to avoid issues at boundry
        ix = ra >= 0; left_ra=ra[ix]-360.; left_dec=dec[ix]
        ix = ra < 0; right_ra=ra[ix]+360.; right_dec=dec[ix]
        ra = np.append(left_ra, ra)
        ra = np.append(ra,right_ra)
        dec = np.append(left_dec, dec)
        dec = np.append(dec,right_dec)
    
    #file = "/data/des30.a/data/annis/test/healpix_{}_ra_dec_tree_wrap_180.pickle".format(nsides)
    file = "./healpix_{}_ra_dec_tree_wrap_180.pickle".format(nsides)
    if not recompute and os.path.exists(file):
        ra,dec,tree = cPickle.load(open(file,"rb"))
    else :
        tree = scipy.spatial.cKDTree(zip(ra*np.cos(dec*2*np.pi/360.),dec))
        if dumpTree: cPickle.dump([ra,dec,tree],open(file,"wb"))

    # 180 is center of map (i.e, 0 is at singularity)
    ra0 = np.copy(ra); dec0 = dec
    ix = ra < 0; ra0[ix]=ra0[ix]+360.
    if wrapForMetricEval :
        # buffer each side of the singularity so to avoid issues at boundry
        ix = ra0 >= 180; left_ra=ra0[ix]-360.; left_dec=dec0[ix]
        ix = ra0 < 180; right_ra=ra0[ix]+360.; right_dec=dec0[ix]
        ra0 = np.append(left_ra, ra0)
        ra0 = np.append(ra0,right_ra)
        dec0 = np.append(left_dec, dec0)
        dec0 = np.append(dec0,right_dec)
    #file = "/data/des30.a/data/annis/test/healpix_{}_ra_dec_tree_wrap_0.pickle".format(nsides)
    file = "./healpix_{}_ra_dec_tree_wrap_0.pickle".format(nsides)
    if not recompute and os.path.exists(file):
        ra0,dec0,tree0 = cPickle.load(open(file,"rb"))
    else :
        tree0 = scipy.spatial.cKDTree(zip((ra0-180.)*np.cos(dec0*2*np.pi/360.),dec0))
        if dumpTree: cPickle.dump([ra0,dec0,tree0],open(file,"wb"))
    return ra,dec,tree, ra0, dec0, tree0
    
def radecInHexPath ( hex_path, ra, dec) :
    ix = hex_path.contains_points( zip(ra,dec) )
    ix = np.array(ix)
    return ix

def hexPath (raCenter, decCenter, camera) :
    ra,dec = cameraOutline(raCenter, decCenter, camera) 
    hex = matplotlib.path.Path(zip(ra,dec))
    return hex

def cameraOutline (raCenter, decCenter, camera) :
    if camera == 'decam':
        ra, dec = cameraOutlineAtZero(camera)
    elif camera == 'hsc':
        ra, dec = cameraOutlineAtZero(camera)
    else:
        raise Exception('No such Camera')
    dec = decCenter + dec
    ra = raCenter + ra/np.cos(dec*2*np.pi/360.)
    return ra,dec

# two dead chips
def cameraOutlineAtZero (camera) :
    if camera == 'decam':
        ra  = [-0.47256, 
            -0.47256, -0.15752,  -0.15752,  0.15752,  0.15752, 0.47256, 
            0.47256, 0.63008, 0.63008,
            0.7876, 0.7876, 0.94512, 0.94512, 0.94512, 1.10264, 1.10264,
            1.10264, 0.94512, 0.94512, 0.94512, 0.7876, 0.7876, 0.63008,
            0.63008, 0.47256, 
            0.47256, 0.15752,  0.15752,  -0.15752,  -0.15752, -0.47256, 
            -0.47256, -0.63008, -0.63008,
            -0.7876, -0.7876, -0.94512, -0.94512, -0.94512, -1.10264, -1.10264,
            -1.10264, -0.94512, -0.94512, -0.94512, -0.7876, -0.7876, -0.63008,
            -0.63008, -0.47256]
        dec = [ 0.824266666667, 
            0.98912, 0.98912,  0.824266666667, 0.824266666667, 0.98912, 0.98912, 
            0.824266666667, 0.824266666667, 0.659413333333,
            0.659413333333, 0.49456, 0.49456, 0.329706666667, 
            0.164853333333, 0.164853333333, 0.0,
            -0.164853333333, -0.164853333333, -0.329706666667, 
            -0.49456, -0.49456, -0.659413333333, -0.659413333333,
            -0.824266666667, -0.824266666667, 
            -0.98912, -0.98912,  -0.824266666667, -0.824266666667, -0.98912, -0.98912, 
            -0.824266666667, -0.824266666667, -0.659413333333,
            -0.659413333333, -0.49456, -0.49456, -0.329706666667, 
            -0.164853333333, -0.164853333333, 0.0,
            0.164853333333, 0.164853333333, 0.329706666667, 0.49456, 
            0.49456, 0.659413333333, 0.659413333333,
            0.824266666667, 0.824266666667]
    if camera == 'hsc':
        ra = [ 0.750000, 0.748490, 0.743966, 0.736447, 0.725962, 0.712553, 0.696276, 0.677195, 0.655387, 0.630940, 0.603953, 0.574533, 0.542801, 0.508882, 0.472915, 0.435043, 0.395419, 0.354203, 0.311561, 0.267665, 0.222690, 0.176819, 0.130236, 0.083129, 0.035686, -0.011899, -0.059437, -0.106736, -0.153605, -0.199855, -0.245301, -0.289759, -0.333050, -0.375000, -0.415440, -0.454207, -0.491146, -0.526106, -0.558948, -0.589540, -0.617757, -0.643488, -0.666627, -0.687081, -0.704769, -0.719620, -0.731572, -0.740579, -0.746604, -0.749622, -0.749622, -0.746604, -0.740579, -0.731572, -0.719620, -0.704769, -0.687081, -0.666627, -0.643488, -0.617757, -0.589540, -0.558948, -0.526106, -0.491146, -0.454207, -0.415440, -0.375000, -0.333050, -0.289759, -0.245301, -0.199855, -0.153605, -0.106736, -0.059437, -0.011899, 0.035686, 0.083129, 0.130236, 0.176819, 0.222690, 0.267665, 0.311561, 0.354203, 0.395419, 0.435043, 0.472915, 0.508882, 0.542801, 0.574533, 0.603953, 0.630940, 0.655387, 0.677195, 0.696276, 0.712553, 0.725962, 0.736447, 0.743966, 0.748490, 0.750000]
        dec = [ 0.000000, 0.047568, 0.094944, 0.141938, 0.188361, 0.234025, 0.278747, 0.322346, 0.364648, 0.405481, 0.444681, 0.482091, 0.517559, 0.550944, 0.582110, 0.610932, 0.637294, 0.661090, 0.682224, 0.700611, 0.716177, 0.728859, 0.738606, 0.745379, 0.749151, 0.749906, 0.747641, 0.742366, 0.734102, 0.722882, 0.708751, 0.691766, 0.671995, 0.649519, 0.624427, 0.596821, 0.566812, 0.534521, 0.500077, 0.463619, 0.425295, 0.385258, 0.343670, 0.300698, 0.256515, 0.211299, 0.165233, 0.118501, 0.071292, 0.023796, -0.023796, -0.071292, -0.118501, -0.165233, -0.211299, -0.256515, -0.300698, -0.343670, -0.385258, -0.425295, -0.463619, -0.500077, -0.534521, -0.566812, -0.596821, -0.624427, -0.649519, -0.671995, -0.691766, -0.708751, -0.722882, -0.734102, -0.742366, -0.747641, -0.749906, -0.749151, -0.745379, -0.738606, -0.728859, -0.716177, -0.700611, -0.682224, -0.661090, -0.637294, -0.610932, -0.582110, -0.550944, -0.517559, -0.482091, -0.444681, -0.405481, -0.364648, -0.322346, -0.278747, -0.234025, -0.188361, -0.141938, -0.094944, -0.047568, -0.000000]
    ra = np.array(ra)
    dec=np.array(dec)
    return ra,dec

# full camera
def cameraOutlineAtZeroFull () :
    ra = [-0.47256, -0.47256, 0, 0.47256, 0.47256, 0.63008, 0.63008,
        0.7876, 0.7876, 0.94512, 0.94512, 0.94512, 1.10264, 1.10264,
        1.10264, 0.94512, 0.94512, 0.94512, 0.7876, 0.7876, 0.63008,
        0.63008, 0.47256, 0.47256, 0, -0.47256, -0.47256, -0.63008, -0.63008,
        -0.7876, -0.7876, -0.94512, -0.94512, -0.94512, -1.10264, -1.10264,
        -1.10264, -0.94512, -0.94512, -0.94512, -0.7876, -0.7876, -0.63008,
        -0.63008, -0.47256]
    dec = [ 0.824266666667, 0.98912, 0.98912, 0.98912, 
        0.824266666667, 0.824266666667, 0.659413333333,
        0.659413333333, 0.49456, 0.49456, 0.329706666667, 
        0.164853333333, 0.164853333333, 0.0,
        -0.164853333333, -0.164853333333, -0.329706666667, 
        -0.49456, -0.49456, -0.659413333333, -0.659413333333,
        -0.824266666667, -0.824266666667, -0.98912, -0.98912, 
        -0.98912, -0.824266666667, -0.824266666667, -0.659413333333,
        -0.659413333333, -0.49456, -0.49456, -0.329706666667, 
        -0.164853333333, -0.164853333333, 0.0,
        0.164853333333, 0.164853333333, 0.329706666667, 0.49456, 
        0.49456, 0.659413333333, 0.659413333333,
        0.824266666667, 0.824266666667]

    ra = np.array(ra)
    dec=np.array(dec)
    return ra,dec

