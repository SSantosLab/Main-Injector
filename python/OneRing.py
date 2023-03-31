import numpy as np
import hp2np 
import hexalate
import decam2hp
import jsonMaker
import simplicity
from os import getenv
import pandas as pd
import sys
import  matplotlib.pyplot as plt


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
def run_or ( skymap, probArea_outer, probArea_inner, flt, expTime_inner, expTime_outer,
            mjd,
            resolution=64,
             hexFile=getenv('DATA_DIR')+"/all-sky-hexCenters-decam.txt", 
            trigger_id="LIGO/Virgo", 
            trigger_type="bright", 
            propid='propid', 
            jsonFilename="des-gw.json",
             test=False,
             plot_numbers=False) :
    
    camera = "decam"
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
        inner_percentile_index, outer_percentile_index, outer_percentile_index ))

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
    # Nora, Elise, Andre (hereafter NEA) update sorting to *only* use sort_zone (the bestest sorting function in OneRing as of 31.03.2023)
   # new_inner_ra, new_inner_dec, new_inner_prob = \
   #     sort_nearest_neighbor (inner_ra, inner_dec, inner_prob) 
   # new_outer_ra, new_outer_dec, new_outer_prob = \
   #     sort_nearest_neighbor (outer_ra, outer_dec, outer_prob) 

    # having resorted on optimal slew time, combine and prep for JSON writing
    #inner_ra = new_inner_ra; inner_dec = new_inner_dec
    #outer_ra = new_outer_ra; outer_dec = new_outer_dec
    #inner_prob = new_inner_prob; outer_prob = new_outer_prob
    #ra = np.concatenate([inner_ra,outer_ra])
    #dec = np.concatenate([inner_dec,outer_dec])
    #prob = np.concatenate([inner_prob,outer_prob])

    #expTime = np.concatenate([np.ones(ra.size)*expTime_inner , np.ones(ra.size)*expTime_outer])
    #### test new nearest_highest:                                     
#    df = sort_nearest_highest(ra, dec, prob)
#    print('nearest highest')
#    print(df)
    
    ### test sort by nearst only
#    df = sort_nearest_neighbor(ra, dec, prob)
#    print('nearest neighbor')
#    print(df)

    ## test spiral
    df = sort_zone(inner_ra, outer_ra, inner_dec, outer_dec, inner_prob, outer_prob, radThreshold=4.)

    #NEA 31.03.2023
    ra = df['ra'].values
    dec = df['dec'].values
    prob = df['probs'].values
    expTime = np.concatenate([np.ones(ra.size)*expTime_inner , np.ones(ra.size)*expTime_outer])
#    print('spiral')
   # print(df)
    
    if plot_numbers:
        plt.clf()
        mymarker = ['$'+str(int(x))+'$' for x in df['order'].values]
        myra = df['ra'].values
        mydec = df['dec'].values
        for i in range(len(myra)):
            plt.scatter(myra[i], mydec[i], marker=mymarker[i], s=100)
            plt.savefig('plotnumbers.jpeg')

    # Now we have the hexes:  can we observe them tonight?
    # this prints details to the screen, a prototype for some action on them
    #obs_object = [simplicity.calc(ra, dec, prob, exp, flt, mjd) for exp in expTime]
    #obs_object = simplicity.calc(ra[0:1], dec[0:1], prob[0:1], expTime[0:1], flt, mjd)

    # Now we have the hexes we're going to observe tonight!
    # JSON writing
    jsonMaker.writeJson(ra,dec,expTime,flt, 
            trigger_id, trigger_type, propid, skymap, jsonFilename ) 

    return
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

    


