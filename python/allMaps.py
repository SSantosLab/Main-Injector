import numpy as np
import healpy as hp
import hp2np
import sourceProb
import mags
import modelRead
import mapsAtTimeT


# sims, mjds, distances, models = allMaps.veni(); 
# allMaps.vidi(sims, mjds, distances, models)

# sims, mjds, distances, models = allMaps.veni(False); 
# allMaps.vidi(sims, mjds, distances, models, False)

# ra, dec, map, maglim, maglimg, prob, probMap, hx,hy = allMaps.readem( 789064, 0)
# sim,mjd,distance,snr, hlv, cra, cdec, p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2015.txt", unpack=True);
#
#=================================================
#
# RUN the LIGO Simulations!!!
#
#=================================================

def selectSingleSim (
        simNumber, data_dir="/data/des30.a/data/annis/des-gw/ligo/sims/") :
    sims, mjds, distances, models = veni()
    ix = sims==simNumber
    sims, mjds, distances = sims[ix], mjds[ix], distances[ix] 
    print "found ",sims[0]
    simfile="bayestar-{:d}.fits.gz".format(sims[0])
    ligoMapFile = data_dir+simfile
    return sims, mjds, distances, ligoMapFile

#
#==== Big routine # 1: collect metadata
#
# Get the list of Ligo Maps
# Get the NS models
#    Definitively a prep into memory routine
#
def veni( ) :
    type = "BBH"
    if type == "2015" :
        dir = "/data/des30.a/data/annis/des-gw/ligo/"
        simsFile = "2015_inj.txt"
        sims, mjds, distances = np.genfromtxt(dir+simsFile, unpack=True, skiprows=40, usecols=(0,2,8))
        sims = sims.astype("int")
    elif type == "2016" :
        dir = "/data/des30.a/data/annis/des-gw/ligo/"
        simsFile = "2016_inj.txt"
        sims, mjds, distances = np.genfromtxt(dir+simsFile, unpack=True, skiprows=40, usecols=(0,2,8))
        sims = sims.astype("int")
    elif type == "BBH" :
        import glob
        dir = "/data/des30.a/data/annis/des-gw/ligo/burst/2015/"
        file_list = glob.glob(dir+"BBH_LIB/*LIB_C.fits.gz")
        gps = np.genfromtxt(dir+"all_BF2Y-2015-BBH_index.txt", unpack=True, skip_header=1, usecols=0)
        sims  = np.genfromtxt(dir+"all_BF2Y-2015-BBH_index.txt", unpack=True, skip_header=1, usecols=1, dtype="str")
        distances= np.ones(sims.size)*20.0
        mjds = gpsToMjd(gps)
    else :
        raise Exception("no sim package known of type {}".format(type))
    models = modelRead.getModels()
    return sims, mjds, distances, models

def gpsToMjd (gps):
    """ 
The BBH sims are given a time in gps seconds.
This is seconds since 0:0:0 6-Jan-1980 and as of 2016 is 17 seconds ahead of UTC
We take this zeropoint time to be MJD= 44244.0.
    This routine corrects for the 17 seconds, which means it becomes progressively
    less accuate as the date moves away from 2016.
    """
    days = (gps-17.0)/3600./24.
    mjd = days + 44244.0
    return mjd

#==== Big routine # 2: find the probabilities over 10 days
#
#  For each sim build the  mags.observed object
#
#   sims, mjds, distances, models=allMaps.veni()
#   allMaps.vidi(sims, mjds, distances, models)
#           test: ix=np.nonzero((sims==10934)|(sims==1087))
#           allMaps.vidi(sims[ix], mjds[ix], distances[ix], models)
#
def vidi(sims, mjds, distances, models, quick=False) :
    import os.path
    type = "BBH"
    if type == "2015" :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out/"
        #odata_dir = "/data/des30.a/data/annis/des-gw/ligo/nsims-2015-out/"
        odata_dir = "/data/des30.a/data/annis/des-gw/nsims-2015-out/"
        file = "bayestar-{:d}.fits.gz"
        trigger_type = "NS"
    elif type == "2016" :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out/"
        #odata_dir = "/data/des30.a/data/annis/des-gw/ligo/nsims-2016-out/"
        odata_dir = "/data/des30.a/data/annis/des-gw/nsims-2016-out/"
        file = "bayestar-{:d}.fits.gz"
        trigger_type = "NS"
    elif type == "BBH" :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/burst/2015/"
        odata_dir = "/data/des30.a/data/annis/des-gw/ligo/burst/out2015/"
        file = "BBH_LIB/{:s}-LIB_C.fits.gz"
        trigger_type = "BH"
    else :
        raise Exception("no sim package known of type {}".format(type))
    outfile = odata_dir + "mainInjector-sim-mjd-dist-bslot-nslot-probcovered-econ_area-need_area-quality.txt"
    fd = open(outfile,"w"); fd.close()
    counter = 0
    for sim, mjd, distance in zip(sims,mjds,distances) :
        simfile = file.format(sim)
        print "sim, distance: ", sim, distance
        simfile = data_dir + simfile
        outdir = odata_dir + str(sim) + "/"
        name = outdir + str(sim)+"-probabilityPlot.png"
        if not quick  and os.path.exists(name) : continue
        if not os.path.exists(simfile): 
            print ".... skipping as no such file"; 
            continue

        best_slot, n_slots, first_slot, \
            econ_prob, econ_area, area_need, quality = \
            mainInjector (sim, simfile, mjd, distance, trigger_type, outdir, 
            recycler_mjd=mjd+(0.5/24.), quick=quick)
        fd=open(outfile,"a")
        fd.write("{} {} {} {} {} {} {} {} {}\n".format(sim, mjd, distance, \
            best_slot, n_slots,  econ_prob, econ_area, area_need,\
            quality))
        fd.close()
        counter += 1
        #if counter >= 21 : raise Exception("im done here")
        
#
#  so what hexes did we do, given the json files?
def gethexIDfromJson (dir = "/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_night1_json/",
        verbose=True) :
    import json; import glob
    import os
    import hexalate
    cwd = os.getcwd(); print cwd
    os.chdir(dir)
    files = glob.glob("*.json")
    nfiles=len(files);
    ra = np.array([])
    dec = np.array([])
    slot = np.array([])
    for n in range(nfiles) :
        print files[n]
        fd=open(files[n],"r"); data=json.load(fd);fd.close()
        for i in range(len(data)): 
            ra = np.append(ra, data[i]["RA"])
            dec = np.append(dec, data[i]["dec"])
            slot = np.append(slot, files[n][9:11])
            if verbose: print data[i]["RA"],data[i]["dec"],files[n][9:11]

    os.chdir(cwd)
    id = hexalate.getHexId(ra,dec)
    return ra,dec,id, slot
def gethexIDfromDB() :
    import hexalate
    file="/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_exp_table_night1.txt"
    file="/data/des30.a/data/annis/des-gw/Jan4-2017-event/GW170104_exp_table_night1_2.txt" 
    ra,dec = np.genfromtxt(file, unpack=True,usecols=(1,2)); 
    ii=ra>180; ra[ii] = ra[ii]-360
    id = hexalate.getHexId(ra,dec)
    id = np.unique(id)
    return id
#
# shall we run the hexes from a ra-dec-id-prob-mjd-slot.tx file?
def hexID (file = "../Jan4-2017-event/GW170104-ra-dec-id-prob-mjd-slot.txt") :
    ra,dec = np.genfromtxt(file, unpack=True,usecols=(0,1) )
    id = np.genfromtxt(file, unpack=True,usecols=(2), dtype = "str" )
    # for GW170104; doing this cut does the obs right, not doing makes the gif cover better area
    ix = ra < 100; ra = ra[ix]; dec = dec[ix]; id = id[ix]
    return ra,dec,id
#
# trigger_type = "NS" or "BH"
#
def mainInjector (trigger_id, skymap, trigger_type,\
        outputDir, recycler_days_since_burst=0.3, 
        hexID_to_do=[], hexID_to_reject=[], 
        hours_used_by_NS=0, hours_used_by_BH=0,
        halfNight = False, firstHalf= True,
        quick=False, debug=True, camera='decam') :

    import os
    import yaml
    import getHexObservations
    skipEcon = False
    skipEcon = True

    with open("maininjector.yaml","r") as f:
        config = yaml.safe_load(f); 
    
    # debug
    debug = config["debug"]    

    # camera
    camera   = config["camera"]

    #resolution
    resolution = config["resolution"]

    # season parameters
    start_of_season   = config["start_of_season"]
    end_of_season     = config["end_of_season"]

    # static observation details
    overhead          = config["overhead"]
    area_per_hex      = config["area_per_hex"]

    # strategy
    exposure_length_ns= config["exposure_length_NS"]
    filter_list_ns    = config["exposure_filter_NS"]
    maxHexesPerSlot_ns= config["maxHexesPerSlot_NS"]
    exposure_length_bh= config["exposure_length_BH"]
    filter_list_bh    = config["exposure_filter_BH"]
    maxHexesPerSlot_bh= config["maxHexesPerSlot_BH"]

    # economics analysis for NS and for BH
    hoursAvailable_ns = config["time_budget_for_NS"]
    hoursAvailable_bh = config["time_budget_for_BH"]
    lostToWeather_ns  = config["hours_lost_to_weather_for_NS"]
    lostToWeather_bh  = config["hours_lost_to_weather_for_BH"]
    rate_bh           = config["rate_of_bh_in_O2"];# events/year
    rate_ns           = config["rate_of_ns_in_O2"];# events/year


    # configure strategy for the event type
    if trigger_type == "NS" :
        hoursAvailable       = hoursAvailable_ns - lostToWeather_ns - hours_used_by_NS
        rate                 = rate_ns
        exposure_length      = exposure_length_ns
        filter_list          = filter_list_ns
        maxHexesPerSlot      = maxHexesPerSlot_ns
        nepochs           = config["nepochs_NS"]
        epochs = np.array([])
        for i in range (1,nepochs+1) :
            epochs        = np.append(epochs, config["epoch{}_NS".format(i)])
        end_date          = config["enddate_NS"]
        distance          = -999
    elif trigger_type == "BH" :
        hoursAvailable       = hoursAvailable_bh - lostToWeather_bh - hours_used_by_BH
        rate                 = rate_bh
        exposure_length      = exposure_length_bh
        filter_list          = filter_list_bh 
        maxHexesPerSlot      = maxHexesPerSlot_bh
        nepochs           = config["nepochs_BH"]
        epochs = np.array([])
        for i in range (1,nepochs+1) :
            epochs        = np.append(epochs, config["epoch{}_BH".format(i)])
        end_date          = config["enddate_BH"]
        distance          = 1.0
    else :
        raise Exception(
            "trigger_type={}  ! Can only compute BH or NS".format(trigger_type))
    exposure_length   = np.array(exposure_length)
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    if camera == "hsc":
        overhead = 20.
        area_per_hex = 1.5

    # make the maps
    if quick :
        skipAll = True
    else :
        skipAll = False
    #resolution = 256 ;# default, resolution element on order of ccd area size
    #resolution = 128 ;# roughly 4 ccds
    #resolution = 64 ;# very fast, debuging, roughly 1/4 of the camera size
    #if debug :
       # return getHexObservations.prepare(
       # skymap, trigger_id, outputDir, outputDir, distance=distance,
       # exposure_list=exposure_length, filter_list=filter_list,
       # overhead=overhead, maxHexesPerSlot=maxHexesPerSlot,
       # start_days_since_burst = recycler_days_since_burst, 
       # skipHexelate=False, skipAll = skipAll,
       # this_tiling = hexID_to_do, reject_hexes = hexID_to_reject, 
       # resolution=resolution, trigger_type=trigger_type, 
       # halfNight = halfNight, firstHalf= firstHalf, debug=debug,camera=camera)
    probs,times,slotDuration,hoursPerNight = getHexObservations.prepare(
        skymap, trigger_id, outputDir, outputDir, distance=distance,
        exposure_list=exposure_length, filter_list=filter_list,
        overhead=overhead, maxHexesPerSlot=maxHexesPerSlot,
        start_days_since_burst = recycler_days_since_burst, 
        skipHexelate=False, skipAll = skipAll,
        this_tiling = hexID_to_do, reject_hexes = hexID_to_reject,
        resolution=resolution, trigger_type=trigger_type,
        halfNight = halfNight, firstHalf = firstHalf,debug=debug, camera=camera)

#    print skymap, trigger_id, outputDir, outputDir, "distance=",distance, \
#        "exposure_list=",exposure_length, "filter_list=",filter_list, \
#        "overhead=",overhead, "maxHexesPerSlot=",maxHexesPerSlot, \
#        "start_days_since_burst = ",recycler_days_since_burst,  \
#        "skipHexelate=",False, "skipAll =", skipAll, \
#        "this_tiling =", hexID_to_do, "reject_hexes =", hexID_to_reject,  \
#        "resolution=",resolution, "trigger_type=",trigger_type, \
#        "halfNight =", halfNight, "firstHalf=", firstHalf
#
#    print probs,times,slotDuration,hoursPerNight  
                
    # figure out how to divide the night
    n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
        probs, times, hoursPerNight=hoursPerNight,
        hoursAvailable=hoursAvailable)

    if quick :
        skipJson = True; 
    else :
        skipJson = False; 
    # compute the best observations
    best_slot = getHexObservations.now( 
        n_slots, mapDirectory=outputDir, simNumber=trigger_id, 
        maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot, 
        exposure_list=exposure_length, filter_list=filter_list, 
        trigger_type = trigger_type, skipJson =skipJson)


    if n_slots > 0 :
#   area_left is the number of hexes we have left to observe this season
#   T_left is the number of days left in the season
#   rate is the effective rate of triggers
#
        time_cost_per_hex = nepochs * np.sum(overhead + exposure_length) #sec 
        area_left =  area_per_hex * (hoursAvailable * 3600)/(time_cost_per_hex)
        time_left = end_of_season - start_of_season

        if skipEcon :
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!     SKIPPING ECON ANALYSIS!  !!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            econ_prob, econ_area, quality = 0,0, 1.0
            need_area = 11734.0 
        # we may want to evaluate each event independently
        if not skipEcon: 
            # do Hsun-yu Chen's 
            print "======================================>>>>>>>>>>>>>>>>>>"
            print " economics "
            print "getHexObservations.economics (", trigger_id, ",",\
                best_slot, ", mapDirectory= \"",outputDir, "\" ,",\
                "area_left=",area_left, ", days_left=",time_left, ",rate=",rate,") "
            print "======================================>>>>>>>>>>>>>>>>>>"
            econ_prob, econ_area, need_area, quality = \
                getHexObservations.economics (trigger_id, 
                    best_slot, mapDirectory=outputDir, 
                    area_left=area_left, days_left=time_left, rate=rate) 
    
            if econ_area > 0.0 :
                hoursOnTarget = (econ_area/area_per_hex ) * (time_cost_per_hex/3600.)
    
                # figure out how to divide the night, 
                # given the new advice on how much time to spend
    
                n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
                    probs, times, hoursPerNight=hoursPerNight,
                    hoursAvailable=hoursOnTarget)
    
                best_slot = getHexObservations.now( 
                    n_slots, mapDirectory=outputDir, simNumber=trigger_id, 
                    maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot, 
                    exposure_list=exposure_length, filter_list=filter_list, 
                    trigger_type = trigger_type, skipJson =skipJson)
    else :
        econ_prob, econ_area, quality = 0,0, 1.0
        need_area = 11734.0 

    if not quick :
        if n_slots > 0 :
            # make observation plots
            print "We're going to do {} slots with best slot {}".format(n_slots, best_slot)
            n_plots = getHexObservations.makeObservingPlots(n_slots, trigger_id, best_slot, outputDir, outputDir, camera, allSky = False)
            string = "$(ls -v {}-observingPlot*)  {}_animate.gif".format(trigger_id, trigger_id)
            os.system("convert  -delay 40 -loop 0  " + string)
        else :
            n_plots = getHexObservations.nothingToObserveShowSomething(trigger_id, outputDir, outputDir)

    # Lets find out how well we did in covering Ligo probability
    sum_ligo_prob = \
        getHexObservations.how_well_did_we_do(skymap, trigger_id, outputDir, camera, resolution)
    return best_slot, n_slots, first_slot, econ_prob, econ_area, need_area, quality


#  sim,mjd,distance,snr, p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2015.txt", unpack=True);

#==== Big routine # 4:  collate information into a single file
#
# There is no such thing as vedi. It is veni, vidi, vici.
# This is concerned with the LIGO ancillary data
# which it is going to connect to....
#       the total probability for each day.
# and thus to a file. I don't know how the day maps are made....
#
# make the big files, the all-maxtimes files
def vedi2(sims, mjds, distances, do2015=True, doV=False) :
    import os.path
    import shutil
    savefile="all-maxtimes-{}.txt"
    snrFile = "/data/des30.a/data/annis/des-gw/ligo/{}_coinc.txt"
    checkFile = "/data/des30.a/data/annis/des-gw/ligo/check-{}.txt"
    if do2015:
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out/"
        if doV :
            data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2015-out-v/"
        savefile = savefile.format("2015")
        snrFile = snrFile.format("2015")
        snrSim, snrNet = np.genfromtxt(snrFile, unpack=True,skiprows=35, usecols=(0,3))
        checkFile = checkFile.format("2015")
        network = np.zeros(snrSim.size)
    else :
        data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out/"
        if doV :
            data_dir = "/data/des30.a/data/annis/des-gw/ligo/sims-2016-out-v/"
        savefile = savefile.format("2016")
        snrFile = snrFile.format("2016")
        snrSim, snrNet = np.genfromtxt(snrFile, unpack=True,skiprows=36, usecols=(0,3))
        network = np.genfromtxt(snrFile, unpack=True,skiprows=36, usecols=(2),dtype=("str"))
        checkFile = checkFile.format("2016")
    check_sim, check_ra, check_dec = np.genfromtxt(checkFile, unpack=True)
    csim,rapidarea = np.genfromtxt("../c2016.txt",unpack=True)
    totalProbs = dict()
    for day in range(0,10) :
        totalProbs[day]  = []
    snr = []
    hlv = []
    cra, cdec = [],[]
    carea =[]
    for sim, mjd, distance in zip(sims,mjds,distances) :
        print "sim, distance: ", sim, distance

        new_maxtimeFile = maxtimeFilename(sim, data_dir)
        maxtimes, maxprobs = np.genfromtxt(new_maxtimeFile, unpack=True)
        snr_ix = np.nonzero(snrSim == sim)
        if network[snr_ix] == "HL" : 
            hlv.append(0)
        elif network[snr_ix] == "HLV" : 
            hlv.append(1)
        elif network[snr_ix] == "HV" : 
            hlv.append(2)
        elif network[snr_ix] == "LV" : 
            hlv.append(3)
        else :
            hlv.append(4)
        snr.append( snrNet[snr_ix] )
        check_ix = np.nonzero(check_sim == sim)
        cra.append(check_ra[check_ix])
        cdec.append(check_dec[check_ix])
        check_ix = np.nonzero(csim == sim)
        carea.append(rapidarea[check_ix])

        for day in range(0,10) :
            simfile = data_dir+str(sim)+"-"+str(day)+"-map.hp"
            if not os.path.exists(simfile) : 
                totalProb = 0.0
            else :
                ligo = hp.read_map(simfile)
                simfile = data_dir+str(sim)+"-"+str(day)+"-probMap.hp"
                decam = hp.read_map(simfile)
                ligo = ligo/ligo.sum()
                totalProb = (decam*ligo).sum()
            totalProbs[day].append(totalProb)
    snr= np.array(snr)
    hlv= np.array(hlv).astype(int)
    cra = np.array(cra)
    cdec = np.array(cdec)
    carea = np.array(carea)
    for day in range(0,10) :
        totalProbs[day] = np.array(totalProbs[day])
    data = np.array([sims, mjds, distances, snr, hlv, cra, cdec, carea, totalProbs[0], \
        totalProbs[1], totalProbs[2], totalProbs[3], totalProbs[4], \
        totalProbs[5], totalProbs[6], totalProbs[7], totalProbs[8], \
        totalProbs[9]])
    np.savetxt(savefile, data.T, "%d %.5f %.0f %.1f %d %.6f %.5f %.3f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e")


# Get the saved maps for each day.
def readem (simNumber, day) :
    n=str(simNumber)+"-"+str(day)

    ra=hp.read_map(n+"-ra.hp"); 
    dec=hp.read_map(n+"-dec.hp"); 
    ha=hp.read_map(n+"-ha.hp"); 
    map=hp.read_map(n+"-map.hp");
    maglim=hp.read_map(n+"-maglim.hp");
    maglimg=hp.read_map(n+"-maglim-global.hp");
    prob=hp.read_map(n+"-prob.hp");
    probMap=hp.read_map(n+"-probMap.hp"); 
    hx=hp.read_map(n+"-hx.hp"); 
    hy=hp.read_map(n+"-hy.hp");
    return ra, dec, map, maglim, maglimg, prob, probMap, hx,hy

# Get the saved maps for each day and hour.
def readem2 (simNumber, day, hourno) :
    n=str(simNumber)+"-"+str(day)+"-"+str(hourno)

    ra=hp.read_map(n+"-ra.hp"); 
    dec=hp.read_map(n+"-dec.hp"); 
    ha=hp.read_map(n+"-ha.hp"); 
    map=hp.read_map(n+"-map.hp");
    maglim=hp.read_map(n+"-maglim.hp");
    prob=hp.read_map(n+"-prob.hp");
    probMap=hp.read_map(n+"-probMap.hp"); 
    hx=hp.read_map(n+"-hx.hp"); 
    hy=hp.read_map(n+"-hy.hp");
    return ra, dec, map, maglim, prob, probMap, hx,hy

# A clear cut save routine
def saven2 (maxtimes, maxprobs, simNumber, data_dir) :
    name = maxtimeFilename ( simNumber, data_dir)
    np.savetxt(name, np.array([maxtimes,maxprobs]).T, "%.5f %.5e")
    print "\t writing ", name

# A file finder routine
def maxtimeFilename ( simNumber, data_dir) :
    nameStem = data_dir + str(simNumber) 
    name = nameStem + "-maxtimes.txt"
    return name
    
#  for ten days,
#   find the maximum total probability for each day
#   return
def maximumProbabilityPerDay (totalProbs,times) :
    maxtimes = []
    maxprobs = []
    # for ten days
    for day in range(0,10) :
        # for times measured in days
        ix = (times >= day) & (times < day+(1.))
        max = totalProbs[ix].max()
        ix = np.nonzero((times >= day) & (times < day+(1.)) & (totalProbs == max))
        if totalProbs[ix].sum() == 0 :
            ix = ix[0][0]
        if totalProbs[ix].size > 1 :
            ix = ix[0][0]
        if type(ix) == np.int64 :
            t = times[ix]; tp = totalProbs[ix]
        else :
            t = times[ix][0]; tp = totalProbs[ix][0]
        maxtimes.append(t)
        maxprobs.append(tp)
    maxtimes =np.array(maxtimes)
    maxprobs = np.array(maxprobs)
    return maxtimes, maxprobs

    sim,mjd,distance,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2015.txt", unpack=True);
    ix = p0 > .01; p=p0[ix]; n=p0.size-p.size; ix=p0>.33; n2=p0[ix].size; plt.clf();a=plt.hist(p,bins=100,color="k"); plt.text(0.4,23,"day 0: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2))
    ix = p1 > .01; p=p1[ix]; n=p1.size-p.size;ix=p1>.33; n2=p1[ix].size; a=plt.hist(p,bins=100,color="r",histtype="step"); plt.text(0.4,21,"day 1: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="r")
    ix = p4 > .01; p=p4[ix]; n=p4.size-p.size;ix=p4>.33; n2=p4[ix].size; a=plt.hist(p,bins=100,color="g",histtype="step"); plt.text(0.4,19,"day 4: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="g")
    ix = p9 > .01; p=p9[ix]; n=p9.size-p.size;ix=p9>.33; n2=p9[ix].size; a=plt.hist(p,bins=100,color="b"); plt.text(0.4,17,"day 9: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="b")
    plt.ylim(0,28);plt.xlabel("probability of detection");plt.ylabel("N");plt.title("2015")
    plt.savefig("hist-2015.pdf")

    os.chdir("/data/des30.a/data/annis/des-gw/ligo/sims-2016-out")
    sim,mjd,distance,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = np.genfromtxt("all-maxtimes-2016.txt", unpack=True);
    ix = p0 > .01; p=p0[ix]; n=p0.size-p.size; ix=p0>.33; n2=p0[ix].size; plt.clf();a=plt.hist(p,bins=100,color="k"); plt.text(0.4,23,"day 0: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2))
    ix = p1 > .01; p=p1[ix]; n=p1.size-p.size;ix=p1>.33; n2=p1[ix].size; a=plt.hist(p,bins=100,color="r",histtype="step"); plt.text(0.4,21,"day 1: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="r")
    ix = p4 > .01; p=p4[ix]; n=p4.size-p.size;ix=p4>.33; n2=p4[ix].size; a=plt.hist(p,bins=100,color="g",histtype="step"); plt.text(0.4,19,"day 4: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="g")
    ix = p9 > .01; p=p9[ix]; n=p9.size-p.size;ix=p9>.33; n2=p9[ix].size; a=plt.hist(p,bins=100,color="b"); plt.text(0.4,17,"day 9: n with p < 0.01: {:3d}      n with p > 0.33: {:3d}".format(n,n2),color="b")
    plt.ylim(0,28);plt.xlabel("probability of detection");plt.ylabel("N");plt.title("2016")
    plt.savefig("hist-2016.pdf")
