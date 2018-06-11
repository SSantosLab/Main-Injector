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
def hexesOnMapTested(ra,dec, raHexen, decHexen) :
    print "    checking {} hexes : ".format(raHexen.size),
    count = np.zeros(ra.size)
    for i in range(0,raHexen.size) :
        ix = radecInHexTested( raHexen[i], decHexen[i], ra, dec)
        count[ix] += 1
    print " "
    return count
# keeps count of times a hex, any hex, overlies a map pixel
#   camera outline for the hexes given
def hexesOnMap(ra,dec, raHexen, decHexen) :
    print "    checking {} hexes : ".format(raHexen.size),
    count = np.zeros(ra.size)
    for i in range(0,raHexen.size) :
        ix = radecInHex( raHexen[i], decHexen[i], ra, dec)
        count[ix] += 1
    print " "
    return count

# if one wants to know the sum of the ligo probability in
# a set of observed hexes, hexalateMap is the routine to use.
#
# return the sum of the map vals inside the hexMap hexes
#
def hexalateMapTested(ra,dec, vals, raHexen, decHexen, verbose=1, useCircle=False, radius=1.0) :
    if verbose : print "\t hexalateMap \t nhex = {},".format(raHexen.size),
    if verbose: print " npix = {}".format(ra.size)
    hexVal = np.zeros(raHexen.size)
    for i in range(0,raHexen.size) :
        if useCircle :
            ix = radecInCircle( raHexen[i], decHexen[i], ra, dec, radius=radius)
        else :
            ix = radecInHexTested( raHexen[i], decHexen[i], ra, dec)
        try  :
            if not ix.size: continue
            hexVal[i] = vals[ix].sum()
        except Exception: 
            print "why are there exceptions in hexalateMap?"
            hexVal[i] = 0
    return hexVal

# I think doing this using a tree requires a projection we will
# use the Sanson-Flamsteed projection, aka sinusoidal projection (x=ra*cos(dec), y=dec)
def hexalateMap(ra,dec, vals, tree, raHexen, decHexen, verbose=1) :
    if verbose : print "\t hexalateMap \t nhex = {},".format(raHexen.size),
    if verbose: print " npix = {}".format(ra.size)
    hexVal = np.zeros(raHexen.size)
    counter = 0
    for i in range(0,raHexen.size) :
        ix = radecInHex( raHexen[i], decHexen[i], ra, dec, tree)
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
def radecInHex ( raCenter, decCenter, ra,dec,tree, mapCenter=0.) :
    camera_radius = 1.2
    sf_scale = 2.6 *np.abs(decCenter)/90. ;# the distortion at high dec means we need wider "circle"
    camera_radius += sf_scale
    data = tree.query_ball_point(
        ((raCenter-mapCenter)*np.cos(decCenter*2*np.pi/360), decCenter), camera_radius)
    near_ix = np.array(data)
    if not near_ix.size: return near_ix

    # find the answer
    hex = hexPath(raCenter, decCenter)
    inside_ix = radecInHexPath( hex, ra[near_ix], dec[near_ix])
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


# return an index that is true/exists if ra,dec is inside
# the camera outline for the hex 
def radecInHexTested ( raCenter, decCenter, ra, dec) :
    # cut away hopless areas in ra, dec
    cosDec = np.cos(decCenter*2*np.pi/360.)
    near_ix = (abs(ra-raCenter)*cosDec <= 1.1) & (dec-decCenter <= 1.1)
    ## just doing the dec sped it up by x2
    #near_ix = (dec-decCenter < 1.5)

    # if there is nothing near, no reason to check further
    if np.all(~near_ix) :
        return near_ix

    # find the answer
    hex = hexPath(raCenter, decCenter)
    inside_ix = radecInHexPath( hex, ra[near_ix], dec[near_ix])
    
    # the issue here is that near_ix is of size ~1,000,000
    # although only ~3000 are True, so inside_ix is of size ~3,000
    # and about ~200 are True. One wants to adjust near_ix using
    # the truth values of inside_ix.
    near_ix[near_ix] = inside_ix
    # simple, but not obvious

    return near_ix


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
    
    file = "/data/des30.a/data/annis/test/healpix_{}_ra_dec_tree_wrap_180.pickle".format(nsides)
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
    file = "/data/des30.a/data/annis/test/healpix_{}_ra_dec_tree_wrap_0.pickle".format(nsides)
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

def hexPath (raCenter, decCenter) :
    ra,dec = cameraOutline(raCenter, decCenter) 
    hex = matplotlib.path.Path(zip(ra,dec))
    return hex

def cameraOutline (raCenter, decCenter ) :
    ra, dec = cameraOutlineAtZero()
    dec = decCenter + dec
    ra = raCenter + ra/np.cos(dec*2*np.pi/360.)
    return ra,dec

# two dead chips
def cameraOutlineAtZero () :
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
    ra = np.array(ra)
    dec=np.array(dec)
    return ra,dec

# full camera
def cameraOutlineAtZeroFull () :
    ra  = [-0.47256, -0.47256, 0, 0.47256, 0.47256, 0.63008, 0.63008,
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


