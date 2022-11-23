import numpy as np
import hp2np 
import hexalate
import decam2hp
import jsonMaker
import simplicity


# Code to demonstrate a very simple way to go from Strategy (from the strategy paper)
# to Json file (mocked up here by printing hexes, filters, exposure times to screen).
#
# This is fast.
#
# Next steps:
#   done 1) get json file writing into place
#   prototyped 2) check that the hexes are above the horizon, probably using code in simplicity.py
#   3) consider more sophisticated sorting of hexes
#       what is in place now is: go to highest hex first, then do nearest neighbors in inner region,
#       then go to highest probability in outer region, then do nearest neighbors in outer region,
#   4) can we use galaxy catalog, hexelation, and rank our hexes?
#
#   Jim Annis and Nora Sherman, 
#   Nov 7, 2022
#
# PS, there are less then 100 lines of code in this file!!!! We declare success!
#
# OneRing.py "it is corrupting us, if it does nothing else!"

# first pass only!
def nike ( skymap, probArea_outer, probArea_inner, filter, expTime_inner, expTime_outer,
            mjd,
            hexFile="all-sky-hexCenters-decam.txt", 
            trigger_id="LIGO/Virgo", 
            trigger_type="bright", 
            propid='propid', 
            jsonFilename="des-gw.json") :

    camera = "decam"
    if (probArea_outer > 1) or ( probArea_inner) : raise Exception("probArea_outer,inner is > 1, impossible")

    # get ligo map
    ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)

    # imprint decam hexes onto sky map
    raCenter,decCenter,junk = hexalate.getHexCenters (hexFile)
    nHexes= raCenter.size
    sum_prob_in_hex = np.zeros(nHexes)
    for i in nHexes:
        hexPath = decam2hp.hexPath(raCenter[i], decCenter[i], camera)
        inside_ix = decam2hp.radecInHexPath( hexPath, ra, dec)
        sum_prob_in_hex[i] = ligo[inside_ix].sum()

    # now we have our data!
    # sort on probability so that max is first on new vector
    index_on_max_prob = np.argsort(sum_prob_in_hex)[::-1]

    # lets find the n'th and 90th percentile list of hexes
    cumsum = np.cumsum(sum_prob_in_hex[index_on_max_prob]) 
    outer_percentile_index = np.argmax( cumsum >= probArea_outer )
    inner_percentile_index = np.argmax( cumsum >= probArea_inner )
    print("N hexes in inner= {}  and in outer= {}".format(
        inner_percentile_index, outer_percentile_index ))

    # now lets get the list of outer and inner hexes
    inner_ra = []; outer_ra=[]
    inner_dec = []; outer_dec=[]
    inner_prob = []; outer_prob=[]
    for i in range(0,inner_percentile_index) :
        # the raCenters sorted on index_on_max_prob are in decending order of prob
        # so we just have to grab them one by on till we reach inner_percentile_index
        inner_ra.append(raCenter[index_on_max_prob][i] )
        inner_dec.append(decCenter[index_on_max_prob][i] )
        inner_prob.append(sum_prob_in_hex[index_on_max_prob][i] )
    for i in range(inner_percentile_index, outer_percentile_index) :
        # and here, too, from inner_percentile_index to outer_percentile_index
        outer_ra.append(raCenter[index_on_max_prob][i] )
        outer_dec.append(decCenter[index_on_max_prob][i] )
        outer_prob.append(sum_prob_in_hex[index_on_max_prob][i] )
    inner_ra = np.array(inner_ra); inner_dec = np.array(inner_dec); 
    outer_ra = np.array(outer_ra); outer_dec = np.array(outer_dec); 
    inner_prob = np.array(inner_prob); outer_prob = np.array(outer_prob); 

    # Sort in the Alyssa approved way!
    # done separately in inner and outer areas.
    new_inner_ra, new_inner_dec, new_inner_prob = \
        sort_nearest_neighbor (inner_ra, inner_dec, inner_prob) 
    new_outer_ra, new_outer_dec, new_outer_prob = \
        sort_nearest_neighbor (outer_ra, outer_dec, outer_prob) 

    # having resorted on optimal slew time, combine and prep for JSON writing
    inner_ra = new_inner_ra; inner_dec = new_inner_dec
    outer_ra = new_outer_ra; outer_dec = new_outer_dec
    inner_prob = new_inner_prob; outer_prob = new_outer_prob
    ra = np.concatenate(inner_ra,outer_ra)
    dec = np.concatenate(inner_dec,outer_dec)
    prob = np.concatenate(inner_prob,outer_prob)

    expTime = np.concatenate( np.ones(ra.size)*expTime_inner , np.ones(ra.size)*expTime_outer)

    # Now we have the hexes:  can we observe them tonight?
    # this prints details to the screen, a prototype for some action on them
    obs_object = simplicity.calc(ra, dec, prob, expTime, filter, mjd)

    # Now we have the hexes we're going to observe tonight!
    # JSON writing
    jsonMaker.writeJson(ra,dec,expTime,filter, 
            trigger_id, trigger_type, propid, skymap, jsonFilename ) 


#
# one of many possible sorting metrics
#
def sort_nearest_neighbor (ra, dec, prob) :
    # Sort in the Alyssa approved way!
    # This particular implementation just chooses nearest neighbor for minimal slew time
    # done separately in inner and outer areas.
    nHexes_inner = ra.size 
    new_ra = np.copy(ra)
    new_dec = np.copy(dec)
    new_prob = np.copy(prob)
    search_index = np.array(range(1,nHexes_inner))
    # having put the max probability onto the stack, find the nearest hex to it
    for i in range(1,nHexes_inner) :
        distances = np.sqrt( 
            (ra[search_index]-new_ra[i-1])**2 + 
            (dec[search_index]-new_dec[i-1])**2   )

        ix = np.argsort(distances)[0] 

        new_ra[i] = ra[search_index][ix]
        new_dec[i] = dec[search_index][ix]
        new_prob[i] = prob[search_index][ix]
        
        search_index = np.delete(search_index, ix)
    return new_ra, new_dec, new_prob

