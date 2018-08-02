import numpy as np
import healpy as hp
import os

import sourceProb
import modelRead

#
# mapsAtTimeT
#
# This file contains routines to compute probability maps
#   for a given burst time and given observation time 
#   in days since the burst.
#
# probabilityMaps
#   Use the mags and sourceProb objects to calculate
#   the limiting magnitude and probability maps at a given time
#   in a given filter and exposure time
#   The models needed as input are found via:
#       models = modelRead.getModels()
#   and then models are ready to roll
#
# totalProbability
#   Return the total probability of in the map at a given time
#   in a given filter and exposure. 
#       This is a thin wrapper about probabilityMaps
#
# manyDaysOfTotalProbability 
#   Return a list of times and total probabilites
#   at deltaTime starting at startOfDays and ending at endOfDays
#   using totalProbabilityAtTimeT
#       The deltaTime parameter sets the duration of a "slot",
#       in the langauge of getHexObservations
#
# oneDayOfTotalProbability 
#   A thin wrapper about manyDaysOfTotalProbability to do one 24 hour period
#
# probabilityMapSaver 
#   Given a list of times and total probabilites,
#   use probabilityMapsAtTimeT (if the probability is high enough)
#   to re-make the maps and save them to disk
#

# Over 11 days
#   calculate sourceProb.map and given a model, sm.calculateProb()
#   sum the probability in the map, append that number to an array,
#   return
#
#   On the subject of deltaTime:
#    deltaTime = 1./24.  ;# once per hour
#    deltaTime = 0.0223  ;# 32 minutes, = 1 hex (i,z,z,i) 180s+30s = 8 minutes 
#       so 4 hexes per 32 minute slot. 
#       More convenient for survey strategy planning (as 8 doesnt go into 60)
#   deltaTime = 0.223*2 ;# twice as long slots, (8 hexes/slot)


def oneDayOfTotalProbability (obs, mjd, spatial, distance, distance_sig, 
        models, deltaTime=0.0223, start_mjd = 0, 
        probTimeFile="probTime.txt", trigger_type="NS",
        halfNight = False, firstHalf= True) :

    # the work.
    start_of_days=0
    end_of_days=1
    #print "JTA debugging  ====="
    #end_of_days=0.1
    totalProbs,times = manyDaysOfTotalProbability(
        obs, mjd, spatial, distance, distance_sig, models, 
        start_mjd = start_mjd, 
        startOfDays=start_of_days, endOfDays=end_of_days,
        deltaTime=deltaTime, probTimeFile=probTimeFile,
        trigger_type=trigger_type, 
        halfNight = halfNight, firstHalf= firstHalf) 

    return totalProbs,times

def manyDaysOfTotalProbability (
        obs, burst_mjd, spatial, distance, distance_sig, 
        models, start_mjd=0,
        startOfDays=0, endOfDays=11, deltaTime=0.0223, 
        probTimeFile="probTime.txt", trigger_type="NS",
        halfNight = False, firstHalf= True) :
    times = []
    totalProbs = []

    # one might want to delay the start of computations
    # till some time after the burst
    if start_mjd == 0: start_mjd = burst_mjd
    delayTime = start_mjd - burst_mjd

    # in the language of getHexObservations:
    #   each slot is 32 minutes, each slot can hold 4 hexes
    dark = False
    isDark = []
    for time in np.arange(startOfDays,endOfDays,deltaTime) :
        time = time + delayTime
        if time < (1.5/24.) : continue
        print "================================== ",
        print "hours since Time Zero: {:.1f}".format(time*24.),
        totalProb, sunIsUp = totalProbability(obs, burst_mjd, time, \
            spatial, distance, distance_sig, models, \
            trigger_type=trigger_type)
        if sunIsUp: print "\t ... the sun is up"
        times.append(time)
        totalProbs.append(totalProb)
        if not dark and not sunIsUp:
            dark = True ;# model sunset
        if dark and sunIsUp:
            dark = False ;# model sunrise
        isDark.append(dark)
    totalProbs =np.array(totalProbs)
    times = np.array(times)
    isDark = np.array(isDark).astype(bool)
    # JTA successful half nights
    # an attempt to deal with half nights
    # if successful, these two booleans should go in yaml
    # somehow and be passed down to here
    #halfNight = True; firstHalf= True
    darkCount = np.nonzero(isDark)[0]
    if halfNight:
        darkHalfSize = np.int(np.round(darkCount.size/2.))
        if firstHalf:
            darkCount = darkCount[:darkHalfSize]
        else :
            darkCount = darkCount[darkHalfSize:]
        tp = totalProbs
        totalProbs = totalProbs*0.0
        totalProbs[darkCount] = tp[darkCount]
    # informational
    print "total all-sky summed probability of detection (list1) and daysSinceBurst (list2)"
    print totalProbs,"\n",times,"\n",isDark.astype(int)
    print "dark slots=",darkCount.size
    
#    print "===== times with total prob > 10**-2"
#    ix = totalProbs > 10**-2; 
#    if np.nonzero(ix)[0].size == 0 :
#        totalProbs = np.array([0,])
#        times = np.array([times[0],])
#    else :
#        totalProbs = totalProbs[ix]
#        times = times[ix]
#    print "total all-sky summed probability of detection (list1) and daysSinceBurst (list2)"
#    print totalProbs,"\n",times

    data = np.array([totalProbs, times, isDark.astype(int)]).T
    np.savetxt(probTimeFile, data, "%f %f %d")
    return totalProbs,times

#==============================================================
#
# core computation
#
def totalProbability(obs, mjdOfBurst, daysSinceBurst, \
        spatial, distance, distance_sig, models,
        filter="i", exposure=180, trigger_type="NS") :
    obs,sm,sunIsUp = probabilityMaps(obs, mjdOfBurst, daysSinceBurst, \
        spatial, distance, distance_sig,
        models, filter, exposure, trigger_type=trigger_type)
    if sunIsUp:
        totalProb = 0.0
    else :
        totalProb = (obs.map * sm.probMap).sum()
    return totalProb, sunIsUp

# drive the probability map calculations. In the end, distance only is used here
def probabilityMaps(obs, mjdOfBurst, daysSinceBurst, \
        spatial, distance, distance_sig, models,
        filter="i", exposure=180, trigger_type="NS") :
    obs.resetTime(mjdOfBurst+daysSinceBurst)
    sm=sourceProb.map(obs, type=trigger_type);  

    sunIsUp = obs.sunBrightnessModel(obs.sunZD)
    if sunIsUp: return obs, sm, sunIsUp

    obs.limitMag(filter, exposure=exposure)
    if trigger_type == "NS" :
        # we may need to rescale the light curve in the models
        models_at_t = modelRead.modelsAtTimeT (models, daysSinceBurst)
        model_scale, model_scale_time = models["scale"]
        AbsMag = sm.modelAbsoluteMagnitude
        new_scale = AbsMag - model_scale
        abs_mag = models_at_t[0] + new_scale
        sm.absMagMean = abs_mag
    sm.searchDistance = np.array([distance,])
    result = sm.calculateProb(spatial, distance, distance_sig)
    if not result:
        sunIsUp = 1
    return obs,sm, sunIsUp

#==============================================================
#
# very  expensive hexelation of maps
#   thus it doublesa a routine that saves 
#   maps and hexelated maps to disk.
#
# for each time in the times,
# calculate the probabilities, 
#   and hexalate
#   then save the maps
#       times,probabilities are the output of manyDaysOfTotalProbability
#
# a onlyHexesAlreadyDone can point at a file that contains
# ra,decs that are already done, and these will replace the
# all sky hexes
# 
def probabilityMapSaver (obs, sim, mjd, ligo, distance, distance_sig,
        models, times, probabilities, data_dir, debug, camera,
        onlyHexesAlreadyDone="", reject_hexes="",
        performHexalatationCalculation=True, trigger_type="NS") :
    import decam2hp
    import hexalate
    import os
    # one reads the tiling 9 hex centers as that is our default position
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    hexFile = gw_data_dir + "all-sky-hexCenters-"+camera+".txt"
    keep_flag = performHexalatationCalculation
    
    prob_slots = np.percentile(probabilities, 95)
    print("95th percentile",prob_slots)
    counter = -1
    for time,prob  in zip(times, probabilities) :
        counter += 1
        performHexalatationCalculation = keep_flag
        if prob <= 0 : 
            performHexalatationCalculation = False
# terrible hack
        if debug:
            if prob <= prob_slots : 
                performHexalatationCalculation = False
        #print "probabilityMapSaver: counter, time= ", counter, time
        if time < 0.06: time = 0.06 ;# if less than 1.5 hours, set to 1.5 hours
        print "================== map save =====>>>>>>>>===== ",
        print "slot {} | hours since Time Zero: {:.1f}".format(counter, time*24.),
        if prob <= 0 : 
            print "\t total probability is zero"
            continue
        obs,sm, isDark = \
            probabilityMaps( obs, mjd, time, ligo, distance, distance_sig,
            models, trigger_type=trigger_type)
        if obs.sunBrightnessModel (obs.sunZD ): 
            #print "\t\tThe sun is up. Continue"
            continue

        # obs.ra, obs.dec, obs.map  = ligo map
        # obs.hx, obs.hy   = mcbryde projection of houar angle, dec
        # obs.maglim = limiting mag
        # sm.prob = limiting mag convolve abmag convolve volume
        # sm.probMap = total prob map
        # hexRa,hexDec,hexVals
        nameStem = os.path.join(data_dir, str(sim) + "-{}".format(str(counter)))
        print "\t Writing files as {}".format(nameStem)

        name = nameStem + "-ra.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.ra)
        name = nameStem + "-dec.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.dec)
        name = nameStem + "-ha.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.ha)
        name = nameStem + "-hx.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.hx)
        name = nameStem + "-hy.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.hy)
        name = nameStem + "-x.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.x)
        name = nameStem + "-y.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.y)
        name = nameStem + "-map.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.map)
        name = nameStem + "-maglim.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.maglim)
        name = nameStem + "-maglim-global.hp"
        if os.path.exists(name): os.remove(name)
        hp.write_map(name, obs.maglimall)
# do we need this?
        #name = nameStem + "-prob.hp"
        #if os.path.exists(name): os.remove(name)
        #hp.write_map(name, sm.prob)
        name = nameStem + "-probMap.hp"
        if os.path.exists(name): os.remove(name)
        #sm.probMap = sm.probMap*gal
        hp.write_map(name, sm.probMap)

        treedata = decam2hp.buildtree(obs.ra*360./2/np.pi,obs.dec*360./2/np.pi,\
            nsides=hp.get_nside(obs.ra), recompute=True)
        tree = treedata[2]
        if performHexalatationCalculation :
            raHexen, decHexen, idHexen = hexalate.getHexCenters(hexFile)
            if len(onlyHexesAlreadyDone) > 0 :
                do_these = np.in1d(idHexen, onlyHexesAlreadyDone)
                do_these = np.nonzero(do_these)
                raHexen, decHexen, idHexen = \
                    raHexen[do_these], decHexen[do_these], idHexen[do_these]
            if len(reject_hexes) > 0 :
                dont_do_these = np.in1d(idHexen, reject_hexes)
                dont_do_these = np.nonzero(np.invert(dont_do_these))
                raHexen, decHexen, idHexen = \
                    raHexen[do_these], decHexen[do_these], idHexen[do_these]

            raHexen, decHexen, idHexen, hexVals, rank = \
                hexalate.cutAndHexalateOnRaDec ( obs, sm, raHexen, decHexen, idHexen, tree, camera)

            # where rank is to be understood as the indicies of the
            # ranked hexes in order; i.e., they have nothing to do with
            # raHexen, decHexen, hexVals except as a sorting key
            name = nameStem + "-hexVals.txt"
            if os.path.exists(name): os.remove(name)
            #print mjd,time
            
            f = open(name,'w')
            for j in range(0,raHexen.size) :
                f.write("{:.6f}, {:.5f}, {:s}, {:.4e}, {:d}, {:.4f}\n".format(
                    raHexen[j],decHexen[j],idHexen[j],hexVals[j],rank[j],(np.asfarray(rank*0.)+(mjd+time))[j]))
            f.close()
            #np.savetxt(name,data.T,fmt="%.6f, %.5f, %s, %.4e, %d, %.4f")
    

# Get the saved maps for each day and hour.
def readMaps (data_dir, simNumber, slot) :
    name = os.path.join(data_dir, str(simNumber) + "-{}".format(str(slot)))

    ra=hp.read_map(name+"-ra.hp");
    dec=hp.read_map(name+"-dec.hp");
    ha=hp.read_map(name+"-ha.hp");
    map=hp.read_map(name+"-map.hp");
    maglim=hp.read_map(name+"-maglim.hp");
    prob=hp.read_map(name+"-prob.hp");
    probMap=hp.read_map(name+"-probMap.hp");
    hx=hp.read_map(name+"-hx.hp");
    hy=hp.read_map(name+"-hy.hp");
    x=hp.read_map(name+"-x.hp");
    y=hp.read_map(name+"-y.hp");
    return ra, dec, map, maglim, prob, probMap, x,y, hx,hy

