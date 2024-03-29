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
import time

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
    detP,
    resolution=64,
    hexFile=getenv('DATA_DIR')+"/all-sky-hexCenters-decam.txt", 
    trigger_id="LIGO/Virgo", 
    trigger_type="bright", 
    propid='propid', 
    jsonFilename="des-gw.json",
    test=False,
    plot_numbers=False
):
    
    timer_start = time.perf_counter()
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
    og_inner_hexlist, sunrise, sunset = af.get_hexinfo(inner_ra, inner_dec, inner_prob, inner_exptime, filt, mjd, detP,True)
    og_outer_hexlist = af.get_hexinfo(outer_ra, outer_dec, outer_prob, outer_exptime, filt, detP,mjd)

    # def find_total_sky(self):
        # self.get_hexlists()
        # total_sky = len(og_inner_hexlist)+len(og_outer_hexlist) # literally just the number of hexes
    total_prob = np.sum([hexy.prob for hexy in og_inner_hexlist]) + np.sum([hexy.prob for hexy in og_outer_hexlist]) # Total prob coverage *possible* from strategy code, not total prob possible for the skymap
    
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


    def sort_hexes(hex_list, current_mjd, last_ra, last_dec, obs_order, observe_mjds, filt_list, exp_list, dithering=False, secondpass=False):
        moon_ra, moon_dec, moon_ha, moon_zd = hex_object.get_moon_pos(current_mjd)

        # Update awesomeness factor for each hex
        for hex_obj in hex_list:
            # If obs_order (the observing plan) is empty, make last ra and dec just current one bc no slew time
            if not obs_order:
                last_ra, last_dec = hex_obj.ra, hex_obj.dec
            hex_obj.getAwesomenessFactor(current_mjd, last_ra, last_dec, moon_ra, moon_dec, moon_zd)

        # Figure out if any hexes have nonzero awesomeness factor
        has_nonzero_awesomeness = any(hex_obj.awesomeness_factor != 0 for hex_obj in hex_list)

        if not has_nonzero_awesomeness:
            # If no hex has a nonzero awesomeness, return early with minimal updates (and saying didn't find any)
            return hex_list, current_mjd, last_ra, last_dec, obs_order, observe_mjds, False, filt_list, exp_list

        hexes_list = order_hexes(hex_list)
        #deep copy is necessary to preserve dithers, even though computationally heavy
        best_hex = copy.deepcopy(hexes_list[0])
        obs_order.append(hexes_list[0])
        observe_mjds.append(current_mjd)

        # Decide on exposure time and filter based on pass
        exptime, filt = (hexes_list[0].expTime[1], hexes_list[0].filt[1]) if secondpass else (hexes_list[0].expTime[0], hexes_list[0].filt[0])
        exp_list.append(exptime)
        filt_list.append(filt)
        
        #get slew time for overhead time calculation
        slewtime = hexes_list[0].slewTime

        # Calculate new MJD with exposure time and overhead
        new_mjd = current_mjd + exptime / 86400 + 28.6 / 86400 + slewtime / 86400 #add overhead
        #28.6 for overhead. 20.6s is the time between readouts, and 8s is the max for filter changes/hexapod movement. additional time is from slew time.
        last_ra, last_dec = best_hex.ra, best_hex.dec
        best_hex.observe_hex()
        print(best_hex.dither)

        # Handle dithering, if this is the first pass send it to be observed with second pass, else put it in the nonpriority 
        #hexlist for further dithers
        if dithering:
            (outer_hexlist if secondpass else dithers_hexlist).append(best_hex)

        # Remove observed hex from queue 
        hex_list.pop(hex_list.index(hexes_list[0]))

        return hex_list, new_mjd, last_ra, last_dec, obs_order, observe_mjds, True, filt_list, exp_list


    current_mjd = np.float64(starttime)



    def plan_observations(starttime, sunrise, inner_hexlist, outer_hexlist, last_ra, last_dec, obs_order, dithers_hexlist):
        current_mjd = np.float64(starttime)
        observed_twice = []
        obs_order = []
        observe_mjds = []
        filt_list = []
        exp_list = []
        prob_list = []
        current_seconds = 0
        while current_mjd < sunrise:
            found_any_hex = False
            #if it has been 2-3 min, take second pass of hexes already covered (in other words, take first dither with 2nd pass filter)
            if current_seconds >= 120:
                    while len(dithers_hexlist) != 0:
                        dithers_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, found_secondpass_hex, filt_list, exp_list = sort_hexes(dithers_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, filt_list, exp_list, dithering=True, secondpass=True)
                        if found_secondpass_hex:
                            prob_list.append(obs_order[-1].detP[1])
                        elif not found_secondpass_hex:
                            print(f'Could not do second pass on some inner hexes.')
                            break
                    current_seconds = 0
            elif len(inner_hexlist) != 0 and current_seconds < 120:
                    old_mjd = current_mjd
                    inner_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, found_hex, filt_list, exp_list = sort_hexes(inner_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, filt_list, exp_list, dithering=True)
                    
                    found_any_hex = found_hex
                    if found_any_hex:
                        current_seconds += (current_mjd - old_mjd)*86400
                        prob_list.append(obs_order[-1].detP[0])
    #                     print(current_seconds)
                    if len(inner_hexlist) != 0 and not found_hex:
                        outer_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, found_outer_hex, filt_list, exp_list = sort_hexes(outer_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, filt_list, exp_list, dithering=False)
                        found_any_hex = found_outer_hex
                        if found_outer_hex:
                            prob_list.append(obs_order[-1].detP[0])
                        if not found_any_hex:
                            current_mjd += (30 / 86400)
            elif len(outer_hexlist) != 0:
                outer_hexlist, current_mjd, last_ra, last_dec,  obs_order, observe_mjds, found_outer_hex, filt_list, exp_list = sort_hexes(outer_hexlist, current_mjd, last_ra, last_dec, obs_order, observe_mjds, filt_list, exp_list, dithering=False)
                found_any_hex = found_outer_hex
                if found_outer_hex:
                    prob_list.append(obs_order[-1].detP[0])
                if not found_any_hex:
                    current_mjd += (30 / 86400)
            else:
                print(len(outer_hexlist), len(inner_hexlist))
                break

        sky_covered = 0.
        prob_covered = 0.
        disc_prob_covered = 0.
        coverage_list = []
        localization_prob_list = [] # discovery prob based on strategy
        disc_prob_list = [] # total discovery prob based on strategy + model

        for i in range(len(obs_order)):
            print(f'At {af.mjd_to_date_time(observe_mjds[i])} will observe this hex:')
            obs_order[i].display_hexinfo()
            print(f'With filter {filt_list[i]}, exptime {exp_list[i]}, and dither {obs_order[i].dither}')
            sky_covered += obs_order[i].coverage_factor
            prob_covered += obs_order[i].coverage_factor * obs_order[i].prob
            disc_prob_covered += obs_order[i].coverage_factor * obs_order[i].prob * prob_list[i]
            coverage_list.append(sky_covered)
            localization_prob_list.append(prob_covered)
            disc_prob_list.append(disc_prob_covered)
            
        
        final_local_prob = prob_covered
        final_disc_prob = disc_prob_covered

        return obs_order, observe_mjds, coverage_list, prob_list, filt_list, exp_list,find_local_prob,final_disc_prob



    # Now we have the hexes we're going to observe tonight!
    obs_order = []
    dithers_hexlist = []
    # inner_hexlist = og_inner_hexlist.copy()
    # outer_hexlist = og_outer_hexlist.copy()
    obs_order, obs_mjds, coverage_list, prob_list, filt_list, exp_list = plan_observations(starttime, sunrise, inner_hexlist, outer_hexlist, None, None, obs_order, dithers_hexlist)
    ra = []
    dec = []
    for i in range(len(obs_order)):
        # print(f'At {hex_functions.mjd_to_date_time(observe_mjds[i])} will observe this hex:')
        # obs_order[i].display_hexinfo()
        ra.append(obs_order[i].ra)
        dec.append(obs_order[i].dec)

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
#    expTime = np.array(expTime)
    print("Printing values for posterity in order: ra,dec,exp_list,filt_list, trigger_id, trigger_type, propid, skymap, jsonFilename",ra,dec,exp_list,filt_list,trigger_id, trigger_type, propid, skymap, jsonFilename,sep="\n\n")
    print(f'Writing observing plan to {jsonFilename}')
    # JSON writing
    jsonMaker.writeJson(ra,dec,exp_list,filt_list, 
            trigger_id, trigger_type, propid, skymap, jsonFilename ) 

    timer_end = time.perf_counter()
    print(f"Finished OneRing in {timer_end - timer_start:0.4f} seconds")

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

# def sort_hexes(hex_list, current_mjd, last_ra, last_dec):
#     moon_ra, moon_dec, moon_ha, moon_zd = hex_object.get_moon_pos(current_mjd)
#     for i in range(len(hex_list)):
#         if len(obs_order) == 0:
#             last_ra, last_dec = hex_list[i].ra, hex_list[i].dec
#         hex_list[i].getAwesomenessFactor(current_mjd, last_ra, last_dec, moon_ra, moon_dec, moon_zd)
#     has_nonzero_awesomeness = any(hex_obj.awesomeness_factor != 0 for hex_obj in hex_list)

#     if has_nonzero_awesomeness:
#         hexes_list = order_hexes(hex_list)
# #         print(hexes_list[0].expTime)
#         obs_order.append(hexes_list[0])
#         observe_mjds.append(current_mjd)

#         # Remove elements at the specified index from all lists
#         new_mjd = current_mjd + (hexes_list[0].expTime / 86400)
#         last_ra, last_dec = hexes_list[0].ra, hexes_list[0].dec
#         hex_list.remove(hexes_list[0])
#         all_mjds.append(new_mjd)
#         return new_mjd, last_ra, last_dec, True
            
#     else: 
# #         print(f'found no hex')
#         new_mjd = current_mjd #+ (30 / 86400)
# #         all_mjds.append(new_mjd)
#         return new_mjd, last_ra, last_dec, False
    

