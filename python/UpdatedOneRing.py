import numpy as np
import hp2np 
import hex_object
import hexalate
import decam2hp
import jsonMaker
from os import getenv
import pandas as pd
import sys
import copy
import matplotlib.pyplot as plt
import hex_functions
import awesomeness_functions as af
from hex_functions import get_hexinfo 

# Code to demonstrate a very simple way to go from Strategy (from the strategy paper)
# to Json file.
#
# This is fast.
#
#   Jim Annis and Nora Sherman, 
#   Nov 7, 2022
#
#   Major Updates by Elise Kesler, January 2024. 
#   All credit to the Main Injectors (Nora, Elise, Andre, Thomas, and Isaac) for the 
#   Awesomeness Factor (which is housed in hex_object)
#
# OneRing.py "it is corrupting us, if it does nothing else!"

# first pass only!
def run_or(
    skymap,
    probArea_outer,
    probArea_inner,
    flt,
    expTime_inner,
    expTime_outer,
    mjd,
    resolution=64,
    hexFile=getenv('DATA_DIR')+"/all-sky-hexCenters-decam.txt", 
    trigger_id="LIGO/Virgo", 
    trigger_type="bright", 
    propid='propid', 
    jsonFilename="des-gw.json",
    test=False,
    plot_numbers=False
):
    
    camera = "decam"
    inner_exptime = expTime_inner
    outer_exptime = expTime_outer
    filt=flt
    if (probArea_inner > 1) or (probArea_outer > 1) : raise Exception("probArea_outer or inner is > 1, impossible")

    # get ligo map
    ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)

    # imprint decam hexes onto sky map
    raCenter,decCenter,junk = hexalate.getHexCenters(hexFile)
    nHexes= raCenter.size
    ix, = np.where(ligo > 10e-6) ## AG: 12-13-22 consider a different mode to the hard cut we are using now.
    ra, dec, ligo = ra[ix], dec[ix], ligo[ix]
    sum_prob_in_hex = np.zeros(nHexes)

    treedata = decam2hp.buildtree(ra,dec,nsides=resolution,recompute=True)
    tree = treedata[2]
    sum_prob_in_hex = decam2hp.hexalateMapWithoutOverlap(ra,dec,ligo,tree, raCenter, decCenter, "decam", verbose=False)

    # now we have our data!
    # sort on probability so that max is first on new vector
    index_on_max_prob = np.argsort(sum_prob_in_hex)[::-1]

    # lets find the n'th and 90th percentile list of hexes
    cumsum = np.cumsum(sum_prob_in_hex[index_on_max_prob]) 
    outer_percentile_index = np.argmax( cumsum >= probArea_outer )
    inner_percentile_index = np.argmax( cumsum >= probArea_inner )
    print("N hexes in inner= {}  and in outer= {}, total number of hexes={}".format(
        inner_percentile_index, outer_percentile_index, inner_percentile_index + outer_percentile_index ))

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

    
    # From here on is implementation of The Main Injectors (consisting of Elise, Nora, Isaac, Thomas, and Andre) new hex sorting code, with sorting based on awesomeness factor
    og_inner_hexlist, sunrise, sunset = af.get_hexinfo(inner_ra, inner_dec, inner_prob, inner_exptime, filt, mjd, True)
    og_outer_hexlist = af.get_hexinfo(outer_ra, outer_dec, outer_prob, outer_exptime, filt, mjd)
    
    # Create copies of original inner and outer lists of hexes in order to remove hexes as they're sorted
    inner_hexlist = og_inner_hexlist.copy()
    outer_hexlist = og_outer_hexlist.copy()
    
    #find starting time for observing by finding the earliest time a hex is observable
    all_risetimes = []
    for thishex in inner_hexlist:
        all_risetimes.append(thishex.rise_mjd)
    for thishex in outer_hexlist:
        all_risetimes.append(thishex.rise_mjd)
    
    starttime = np.min(all_risetimes)
        
    # initialize ra and dec for slew time factor and list of observing mjds and corresponding hexes
    last_ra = None
    last_dec = None
    observe_mjds = []
    obs_order = []
    current_mjd = sunset
    
    
    def order_hexes(hex_list):
        """
        Takes list of hex objects and sorts the list by awesomeness factor (highest 
        awesomeness factor first)
        """
        new_list = sorted(hex_list, key=lambda x: x.awesomeness_factor, reverse=True)
        return new_list
    
    #do hex sorting
    last_ra = None
    last_dec = None
    observe_mjds = []
    obs_order = []

    def sort_hexes(hex_list, current_mjd, last_ra, last_dec, dithering = False):
        """
        Takes list of hex objects, sorts it by awesomeness factor at the current mjd, 
        and observes the best hex. 
        ---
        Inputs
        ---
        hex_list (list): List of hex objects to sort.
        current_mjd (float): The MJD we are observing at.
        last_ra (float): The last ra observed at. Used to calculate slew time. 
        last_dec (float): The last dec observed at. Used to calculate slew time. 
        dithering (boolean): True if we care about observing the hex again with another dither later.
        Outputs
        ---
        hex_list (list): New list of hex objects to sort. If a hex was observed, it's
        gone from the list.
        new_mjd (float): Updated MJD used for the next observation.
        last_ra (float): Updated last ra observed at. 
        last_dec (float): Updated last dec observed at. 
        found_hex (boolean): True if a good hex was found and observed, false if not. Basically
        just tells you if it found an observable hex. 
        """
        
        #get moon position
        moon_ra, moon_dec, moon_ha, moon_zd = hex_object.get_moon_pos(current_mjd)
        
        #determine the awesomeness factor for each hex in the list
        for i in range(len(hex_list)):
            if len(obs_order) == 0:
                last_ra, last_dec = hex_list[i].ra, hex_list[i].dec
            hex_list[i].getAwesomenessFactor(current_mjd, last_ra, last_dec, moon_ra, moon_dec, moon_zd)
            
        #figures out if any of the awesomeness factors are nonzero (ie determines if any observable)
        has_nonzero_awesomeness = any(hex_obj.awesomeness_factor != 0 for hex_obj in hex_list)

        #if there are observable hexes, sort them by awesomeness factor by calling order_hexes
        if has_nonzero_awesomeness:
            hexes_list = order_hexes(hex_list)
            # Get the best hex. has to be a copy for dither purposes
            best_hex = copy.deepcopy(hexes_list[0])
            # Add it to the observation plan
            obs_order.append(hexes_list[0])
            observe_mjds.append(current_mjd)

            # Remove elements at the specified index from all lists (remove observed hex)
            new_mjd = current_mjd + (hexes_list[0].expTime / 86400)
            # update last observed ra and dec based on hex just observed
            last_ra, last_dec = hexes_list[0].ra, hexes_list[0].dec
            
            #if you're observing something again, ie seeing it with a dither, say we dithered
            if hexes_list[0].dither != [0.00, 0.00]:
                print(f'Observed hex {hexes_list[0].ra}, {hexes_list[0].dec} multiple times')
            
            # update the observed hex object so it's possible to see it with dither
            best_hex.observe_hex()
            #Add the observed hex to hexes_list_outer so that it may be observed again with new dither
            if dithering:
                outer_hexlist.append(best_hex)
            index_to_remove = hex_list.index(hexes_list[0])
            #remove the observed hex from hexlist
            hex_list.pop(index_to_remove)

            all_mjds.append(new_mjd)
    #         print(f'after observing a hex, have {af.mjd_to_date_time(new_mjd)}')
            return hex_list, new_mjd, last_ra, last_dec, True

        else: 
    #         print(f'found no hex')
            #if didn't find any hex, just skip forward 30 seconds so the code can try again
        #this is actually done outside this function, which is why it is commented out
            new_mjd = current_mjd #+ (30 / 86400)
    #         all_mjds.append(new_mjd)
            return hex_list, new_mjd, last_ra, last_dec, False

    current_mjd = np.float64(starttime)


    while current_mjd < sunrise:
        #feel bad there's a while loop, but this works. this goes through each mjd, finds best
        #hex, and makes the json. 
        
        #prioritize inner hexes that have not been observed yet -- if there are any, sort those
        #and observe best one. 
        if len(inner_hexlist) != 0:
            found_any_hex = True
            inner_hexlist, current_mjd, last_ra, last_dec, found_hex = sort_hexes(inner_hexlist, current_mjd, last_ra, last_dec, dithering=True)
            found_any_hex = found_hex
    #         if found_any_hex:
    #                 print('yay i found one')
            #if didn't find any good inner hexes, look for non-priority hexes.
            if len(inner_hexlist) != 0 and not found_hex:
                outer_hexlist, current_mjd, last_ra, last_dec, found_outer_hex = sort_hexes(outer_hexlist, current_mjd, last_ra, last_dec, dithering=False)
                found_any_hex = found_outer_hex
                
                #if didn't find any hexes at all, just skip forward and try again
                if not found_any_hex:
                    current_mjd += (30 / 86400)
                    #print(f'first guy bud at {af.mjd_to_date_time(current_mjd, with_date=True)}')
        #if there's no inner hexes left, observe nonpriority hexes (outer hexes and inner hexes 
        #with a dither)
        elif len(outer_hexlist) != 0:
                found_any_hex = True
                outer_hexlist, current_mjd, last_ra, last_dec, found_outer_hex = sort_hexes(outer_hexlist, current_mjd, last_ra, last_dec)
                found_any_hex = found_outer_hex
                if not found_any_hex:
                    current_mjd += (30 / 86400)
                #if found_outer_hex:
                    #print(f'after finding hex outer hexlist has {len(outer_hexlist)}')

        else:
            print(len(outer_hexlist), len(inner_hexlist))
            current_mjd = sunrise
            break

    # Now we have the hexes we're going to observe tonight!
    ra = []
    dec = []
    expTime = []

    for i in range(len(obs_order)):
        print(f'At {hex_functions.mjd_to_date_time(observe_mjds[i])} will observe this hex:')
        obs_order[i].display_hexinfo()
        ra.append(obs_order[i].ra)
        dec.append(obs_order[i].dec)
        expTime.append(obs_order[i].expTime)

    if plot_numbers:
        plt.clf()
        mymarker = ['$'+str(int(x))+'$' for x in df['order'].values]
        myra = df['ra'].values
        mydec = df['dec'].values
        for i in range(len(myra)):
            plt.scatter(myra[i], mydec[i], marker=mymarker[i], s=100)
            plt.savefig('plotnumbers.jpeg')
    ra = np.array(ra)
    dec = np.array(dec)
    expTime = np.array(expTime)
    
    print(f'Writing observing plan to {jsonFilename}')
    # JSON writing
    jsonMaker.writeJson(ra,dec,expTime,flt, 
            trigger_id, trigger_type, propid, skymap, jsonFilename ) 

    return


##FROM HERE ON IS OLD ONERING ##
#
# one of many possible sorting metrics
#
def sort_nearest_neighbor (ra, dec, prob) :
    # Sort in the Alyssa approved way!
    # This particular implementation just chooses nearest neighbor for minimal slew time
    # done separately in inner and outer areas.
    nHexes = ra.size 
    new_ra = np.copy(ra)
    new_dec = np.copy(dec)
    new_prob = np.copy(prob)
    search_index = np.array(range(1,nHexes))
    # having put the max probability onto the stack, find the nearest hex to it
    for i in range(1,nHexes) :
        distances = np.sqrt( 
            (ra[search_index]-new_ra[i-1])**2 + 
            (dec[search_index]-new_dec[i-1])**2   )
    
        ix = np.argsort(distances)[0]

        new_ra[i] = ra[search_index][ix]
        new_dec[i] = dec[search_index][ix]
        new_prob[i] = prob[search_index][ix]
        
        search_index = np.delete(search_index, ix)
    return new_ra, new_dec, new_prob


def sort_nearest_highest(ra, dec, prob, startorder=0., radThreshold=3.):
    # Sort by nearest neighbor, find all neighbors within some radius (radThreshold), and sort by probability.
    
    order = np.ones(ra.size)*-999
    df = pd.DataFrame({'ra':ra, 'dec':dec, 'probs':prob, 'order':order})

    i = startorder ## order counter
    idxmaxprob = df['probs'].idxmax()
    while any(df['order'] == -999):
        # calculate distance from max prob point, and discretize by threshold radius
        mydist = np.sqrt((df.iloc[idxmaxprob]['ra'] - df['ra'])**2 + (df.iloc[idxmaxprob]['dec'] - df['dec'])**2)
        df['distance'] = mydist // radThreshold

        # sort  by  distance with prob breaking ties
        df.sort_values(['distance', 'probs'], ascending=[True, False])
        
        groupdf = df.groupby(['distance']).get_group(0) # get all rows where quantized distance = 0 
        groupdf = groupdf[groupdf.order == -999] # make sure we are only ordering points that haven't been observed 

        # this is a check to ensure we don't get stuck on banana or lobe
        # eg everything in the closest radius has already been observed, check for groups of points further away
        j = 0
        while groupdf.empty:
            try: # Try b/c disconnected blobs make empty groups
                groupdf = df.groupby(['distance']).get_group(j)
                groupdf = groupdf[groupdf['order'] == -999]
            except(KeyError): pass
            j += 1

        neworder = np.arange(i, i + len(groupdf)) 
        groupdf['order'] = neworder
        i += len(groupdf)

        df.loc[groupdf.index.values ,'order'] = groupdf['order']
        idmaxprob = groupdf['probs'].idxmin() # start next iteration with the last value in group

    df['distance'] = mydist
    return df

#main(args)

def sort_spiral(inner_ra, outer_ra, inner_dec, outer_dec, inner_prob, outer_prob, radThreshold=3.):
    # Use descrete sets of inner and outer hexes, within each of those sets, sort by distance then by prob.
    # JA: Cartoon code. 
    ra, dec, prob = inner_ra, inner_dec, inner_prob
    order = np.ones(ra.size)*-999
    df = pd.DataFrame({'ra':ra, 'dec':dec, 'probs':prob, 'order':order})

    i = 0 ## order counter                                                                             
    idxmaxprob = df['probs'].idxmax()

    while any(df['order'] == -999):
        # calculate distance from max prob point, and discretize by threshold radius                   
        mydist = np.sqrt((df.iloc[idxmaxprob]['ra'] - df['ra'])**2 + (df.iloc[idxmaxprob]['dec'] - df['dec'])**2)
        df['distance'] = mydist // radThreshold

        # sort  by  distance with prob breaking ties                                                   
        df.sort_values(['distance', 'probs'], ascending=[True, False])

        groupdf = df.groupby(['distance']).get_group(0) # get all rows where quantized distance = 0    
        groupdf = groupdf[groupdf.order == -999] # make sure we are only ordering points that haven't been observed                                                                                          

        # this is a check to ensure we don't get stuck on banana or lobe                               
        # eg everything in the closest radius has already been observed, check for groups of points further away                                                                                             
        j = 0
        while len(groupdf) == 0:
            try:
                groupdf = df.groupby(['distance']).get_group(j)
                groupdf = groupdf[groupdf['order'] == -999]
            except(KeyError): pass
            j += 1

        neworder = np.arange(i, i + len(groupdf))
        groupdf['order'] = neworder
        i += len(groupdf)

        df.loc[groupdf.index.values ,'order'] = groupdf['order']
        idmaxprob = groupdf['probs'].idxmin() # start next iteration with the last value in group       
    df['distance'] = mydist
    df['ra'] = df['ra'] + 360
    innerdf = df



#### outer
    ra, dec, prob = outer_ra, outer_dec, outer_prob
    order = np.ones(ra.size)*-999
    df = pd.DataFrame({'ra':ra, 'dec':dec, 'probs':prob, 'order':order})

    i = len(innerdf) ## order counter                                                                                  
    idxmaxprob = df['probs'].idxmax()

    while any(df['order'] == -999):
        # calculate distance from max prob point, and discretize by threshold radius                       
        mydist = np.sqrt((df.iloc[idxmaxprob]['ra'] - df['ra'])**2 + (df.iloc[idxmaxprob]['dec'] - df['dec'])**2)
        df['distance'] = mydist // radThreshold
        # sort  by  distance with prob breaking ties                                                       
        df.sort_values(['distance', 'probs'], ascending=[True, False])
        groupdf = df.groupby(['distance']).get_group(0) # get all rows where quantized distance = 0        
        groupdf = groupdf[groupdf.order == -999] # make sure we are only ordering points that haven't been observed                                                                                                    

        # this is a check to ensure we don't get stuck on banana or lobe                                  
        # eg everything in the closest radius has already been observed, check for groups of points furthe away                                                                                                       
        j = 0
        while len(groupdf) == 0:
            try:
                groupdf = df.groupby(['distance']).get_group(j)
                groupdf = groupdf[groupdf['order'] == -999]
            except(KeyError): pass
            j += 1
        neworder = np.arange(i, i + len(groupdf))
        groupdf['order'] = neworder
        i += len(groupdf)
        df.loc[groupdf.index.values ,'order'] = groupdf['order']
        idmaxprob = groupdf['probs'].idxmin() # start next iteration with the last value in group      
    df['distance'] = mydist

    # concat inner and outer df, but make sure the inner is first  
    df = pd.concat([innerdf, df], axis=0)
    
    return df


def sort_zone(inner_ra, outer_ra, inner_dec, outer_dec, inner_prob, outer_prob, radThreshold=3.):
    
    if len(inner_ra) == 0:
        inner_df = pd.DataFrame()
    else:
        inner_df = sort_nearest_highest(inner_ra, inner_dec, inner_prob, radThreshold=radThreshold)
    if len(outer_ra) == 0:
        outer_df = pd.DataFrame()
    else:
        outer_df = sort_nearest_highest(outer_ra, outer_dec, outer_prob, startorder=len(inner_df), radThreshold=radThreshold)
    return pd.concat([inner_df, outer_df], axis=0)

def order_hexes(hex_list):
    new_list = sorted(hex_list, key=lambda x: x.awesomeness_factor, reverse=True)
    return new_list

last_ra = None
last_dec = None
observe_mjds = []
obs_order = []
#for testing purposes
all_mjds = []
all_awesomeness_factors = []

def sort_hexes(hex_list, current_mjd, last_ra, last_dec):
    moon_ra, moon_dec, moon_ha, moon_zd = hex_object.get_moon_pos(current_mjd)
    for i in range(len(hex_list)):
        if len(obs_order) == 0:
            last_ra, last_dec = hex_list[i].ra, hex_list[i].dec
        hex_list[i].getAwesomenessFactor(current_mjd, last_ra, last_dec, moon_ra, moon_dec, moon_zd)
    has_nonzero_awesomeness = any(hex_obj.awesomeness_factor != 0 for hex_obj in hex_list)

    if has_nonzero_awesomeness:
        hexes_list = order_hexes(hex_list)
#         print(hexes_list[0].expTime)
        obs_order.append(hexes_list[0])
        observe_mjds.append(current_mjd)

        # Remove elements at the specified index from all lists
        new_mjd = current_mjd + (hexes_list[0].expTime / 86400)
        last_ra, last_dec = hexes_list[0].ra, hexes_list[0].dec
        hex_list.remove(hexes_list[0])
        all_mjds.append(new_mjd)
        return new_mjd, last_ra, last_dec, True
            
    else: 
#         print(f'found no hex')
        new_mjd = current_mjd #+ (30 / 86400)
#         all_mjds.append(new_mjd)
        return new_mjd, last_ra, last_dec, False
    


