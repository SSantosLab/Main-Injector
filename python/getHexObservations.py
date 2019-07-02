import matplotlib
matplotlib.use("Agg"); # matplotlib.use("TkAgg") important for off line image generation
import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import hp2np
import decam2hp
import jsonMaker
import obsSlots
import gw_map_configure
import pickle
import glob
#
#       the routine mapsAtTimeT.oneDayOfTotalProbability
#           breaks the night into slots of 32 minutes duration
#       the routine mapsAtTimeT.probabilityMapSaver
#           makes maps for each of these slots
#       the routine getHexObservations.observing
#           once told how many slots can be used
#           returns "hoursObserving" containing ra,dec,prob,time
#               of the hexes to observe each slot
#       the routine obsSlots.observingStats
#           places the ra,dec,id,prob,times per slot onto single  lists
#           tells sum(prob) on the list
#       the routine getHexObservations.observingPlot
#           needs to be told simNumber, "slot", data_dir, nslots
#           and then will plot the sky limitMag+ligo contour
#           and all hexes (in the night) on the map for that slot
#

#========================================================================
#
# main routines: these seven are called in Dillon's recycler.py
#   make_maps
#   make_divisions_of_time
#   make_jsons
#   economics
#   makeObservingPlots      
#   nothingToObserveShowSomething
#   readMaps   
#
#========================================================================

# ==== prep the observation by calculating necessary maps
#       distance is a problem- the question is whether we want a horizon
#       distance or the estimated distance from LIGO?
#       >>> distance = 75. ;# Mpc, as an estimated horizon distance
#
#    deltaTime = 0.0223  => 32 minute slots  for izzi 90sec exp, 4 hex/slot
#    deltaTime = 0.0446  => 62 minute slots  for izzi 90sec exp, 8 hex/slot
#    deltaTime = 0.0417  => 60 minute slots  for izz 90sec exp, 10 hex/slot
#    deltaTime = 0.0208  => 30 minute slots  for izz 90sec exp,  5 hex/slot
#       The first is flexible to observing conditions,
#       while the second is x2 faster as there are less maps to hexalate
#       The third shows a different tiling/filter complement
#
#   exposure_length is only used in the computation of the limiting mag map
#
def make_maps(gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) :
    import mapsAtTimeT
    import mags
    import modelRead
    import healpy as hp
    import os
    forcedistance = False

    skymap               = gw_map_trigger.skymap
    trigger_id           = gw_map_trigger.trigger_id
    trigger_type         = gw_map_trigger.trigger_type
    distance             = gw_map_trigger.distance
    days_since_burst     = gw_map_trigger.days_since_burst
    burst_mjd            = gw_map_trigger.burst_mjd

    camera               = gw_map_strategy.camera
    exposure_list        = gw_map_strategy.exposure_list
    filter_list          = gw_map_strategy.filter_list
    maxHexesPerSlot      = gw_map_strategy.maxHexesPerSlot
    overhead             = gw_map_strategy.overhead 

    resolution           = gw_map_control.resolution
    this_tiling          = gw_map_control.this_tiling
    reject_hexes         = gw_map_control.reject_hexes
    debug                = gw_map_control.debug
    data_dir             = gw_map_control.datadir
    snarf_mi_maps        = gw_map_control.snarf_mi_maps
    mi_map_dir           = gw_map_control.mi_map_dir

    print "burst_mjd = {:.2f} at a distance of {:.1f} Mpc".format(
        burst_mjd, distance)
    print "calculation starting at  {:.1f} days since burst\n".format(
         days_since_burst)
    start_mjd = burst_mjd + days_since_burst


    if snarf_mi_maps:
        if mi_map_dir == "/data/des41.a/data/desgw/O3FULL/Main-Injector/OUTPUT/O3REAL/":
            mi_map_dir = mi_map_dir + trigger_id 
        # get everything we need from a different directory
        print "copying from {}/ to {}/".format(mi_map_dir, data_dir)
        os.system("cp {}/*hp {}/".format(mi_map_dir, data_dir))
        os.system("cp {}/*probabilityTimeCache*txt {}/".format(mi_map_dir, data_dir))
        os.system("cp {}/*-hexVals.txt {}/".format(mi_map_dir, data_dir))
        # cp $other_dir/*hp .
        # cp $other_dir/probabilityTimeCache*txt .
        # cp $other_dir/*-hexVals.txt .
        return


    print "\t cleaning up"
    files = glob.glob(data_dir+"/*png"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*jpg"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*json"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*hp"); 
    for f in files: os.remove(f)
    files = glob.glob(data_dir+"/*txt"); 
    for f in files: os.remove(f)
    print "done cleaning up"

    exposure_list = np.array(exposure_list)
    filter_list = np.array(filter_list)
    ix = filter_list == "i"
    exposure_length = exposure_list[ix].sum()
    print "obs slots starting"

    answers = obsSlots.slotCalculations( burst_mjd, exposure_list, overhead, 
        hexesPerSlot=maxHexesPerSlot) 
    print "obs slots done"
    hoursPerNight = answers["hoursPerNight"] ;# in minutes
    slotDuration = answers["slotDuration"] ;# in minutes
    deltaTime = slotDuration/(60.*24.) ;# in days

        
    # ==== get the neutron star explosion models
    models = modelRead.getModels()

    # === prep the maps
    ra,dec,ligo=hp2np.hp2np(skymap, degrade=resolution, field=0)
    print "got map"
    ligo_dist, ligo_dist_sig, ligo_dist_norm  = \
        distance*np.ones(ra.size), np.zeros(ra.size), np.zeros(ra.size)
    if not forcedistance:
        try :
            junk,junk,ligo_dist =hp2np.hp2np(skymap, degrade=resolution, field=1)
            junk,junk,ligo_dist_sig =hp2np.hp2np(skymap, degrade=resolution, field=2)
            junk,junk,ligo_dist_norm =hp2np.hp2np(skymap, degrade=resolution, field=3)
        except:
            print "\t !!!!!!!! ------- no distance information in skymap ------ !!!!!!!!"

    obs = mags.observed(ra,dec,ligo, start_mjd+.3, verbose=False)
    obs.limitMag("i",exposure=exposure_length)
    print "finished setting up exposure calculation"

    # ==== calculate maps during a full night of observing
    probabilityTimesCache = os.path.join(data_dir,\
        "probabilityTimesCache_"+str(trigger_id)+".txt")
    probs,times,isdark = mapsAtTimeT.oneDayOfTotalProbability(
        obs, models, deltaTime, start_mjd, probabilityTimesCache,
        gw_map_trigger, gw_map_control)

    mapsAtTimeT.probabilityMapSaver (obs, models, times, probs, 
        gw_map_trigger, gw_map_strategy, gw_map_control)

    gw_map_results.probability_per_slot = probs
    gw_map_results.time_of_slot  = times
    gw_map_results.slotDuration  = slotDuration
    gw_map_results.hoursPerNight = hoursPerNight
    gw_map_results.isdark        = isdark
    return 

# ==== figure out what to observe
#
# Another fast routine
#   basic for this routine is how many hexes per slot
#
#   this_tiling and reject_hexes are inverse in thier effect:
#   only hexes within 36" of the ra,decs of the centers of this_tiling
#   are to be observed, while hexes wihtin within 36" of the ra,decs
#   of the centers of reject_hexes are skipped
#
#   if a doneFile is given, which should be
#   something with ra,dec in columns 1,2, like 
#       G184098-ra-dec-prob-mjd-slot.txt
#   then those ra,decs are interpreted as hex centers,
#   and hexes in the hexVals maps within 36" of the ra,decs
#   are removed- they are done, the question is what to observe next
#
def make_hexes( gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, 
    start_slot = -1, do_nslots=-1) :

    trigger_id = gw_map_trigger.trigger_id
    trigger_type  = gw_map_trigger.trigger_type
    maxHexesPerSlot  = gw_map_strategy.maxHexesPerSlot
    exposure_list  = gw_map_strategy.exposure_list
    filter_list  = gw_map_strategy.filter_list
    hoursAvailable  = gw_map_strategy.hoursAvailable
    data_dir = gw_map_control.datadir
    if start_slot != -1 or do_nslots != -1 :
        gw_map_control.start_slot = start_slot
        gw_map_control.do_nslots = do_nslots

    probs = gw_map_results.probability_per_slot 
    times = gw_map_results.time_of_slot 
    hoursPerNight = gw_map_results.hoursPerNight 
    if (hoursPerNight == False)  :
        reuse_results(data_dir, gw_map_trigger, gw_map_strategy, gw_map_results) 
        probs = gw_map_results.probability_per_slot 
        times = gw_map_results.time_of_slot 
        hoursPerNight = gw_map_results.hoursPerNight 

    n_slots, mapZero  = make_divisions_of_time (
        probs, times, hoursPerNight, hoursAvailable) 
    print "\t tonight is {:.2f} hours long with {} slots".format(hoursPerNight, n_slots)
        
    # if the number of slots is zero, nothing to observe or plot
    if n_slots == 0: 
        return 0
    # compute the observing schedule
    print "=============>>>>  observing"
    if start_slot != -1 or do_nslots != -1 :
        print "\t tonight we will use {} slots starting at {}".format(do_nslots, start_slot)
    hoursObserving=obsSlots.observing(
        trigger_id,n_slots, data_dir, mapZero=mapZero,
        maxHexesPerSlot = maxHexesPerSlot, do_nslots = do_nslots, start_slot=start_slot)

    # print stats to screen
    print "=============>>>>  observingStats"
    # save results to the record -- here is made the ra-dec-id- file
    writeObservingRecord(hoursObserving,   data_dir, gw_map_trigger, gw_map_control)
    ra,dec,id,prob,mjd,slotNumbers,islots = obsSlots.observingStats(hoursObserving, mapZero, do_nslots, start_slot)
    # Get slot number  of maximum probability
    maxProb_slot = obsSlots.maxProbabilitySlot(prob,slotNumbers)
    hoursObserving["maxSlot"] = maxProb_slot
    pickle.dump(hoursObserving, open("{}-hoursObserving.pickle".format(trigger_id),"w"))
    # if the length of ra is one and value zero, nothing to observe or plot
    if ra.size == 1 and ra[0] == 0: 
        return 0
    # shall we measure the total ligo probability covered?
    # Lets find out how well we did in covering Ligo probability
    try:
        sum_ligo_prob = how_well_did_we_do( gw_map_trigger, gw_map_strategy, gw_map_control)
    except:
        print '=============>>>> ZERO PROBABILITY!'
        sum_ligo_prob = 0.
        
    if do_nslots == -1 :
        cumulPlot(trigger_id, data_dir) 
        cumulPlot2(trigger_id, data_dir) 

    gw_map_results.n_slots = n_slots
    gw_map_results.first_slot = mapZero
    gw_map_results.best_slot = maxProb_slot
    gw_map_results.sum_ligo_prob = sum_ligo_prob

    return 

# Make the json files
def make_jsons(gw_map_trigger, gw_map_strategy, gw_map_control) :
    trigger_id    = gw_map_trigger.trigger_id
    trigger_type  = gw_map_trigger.trigger_type
    exposure_list = gw_map_strategy.exposure_list
    filter_list   = gw_map_strategy.filter_list
    propid        = gw_map_strategy.propid
    data_dir      = gw_map_control.datadir
    ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)
    print "\n=============>>>>  JSONs"
    turnObservingRecordIntoJSONs(
        ra,dec,id,prob,mjd,slotNum, trigger_id,
        exposure_list=exposure_list, filter_list=filter_list, 
            trigger_type=trigger_type, mapDirectory=data_dir, propid=propid) 


#
# ====== there are possibilities. Show them.
#
def makeGifs (gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results) :
    # make gif centered on hexes
    n_plots = makeObservingPlots(
        gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, allSky=False)

    if  gw_map_control.allSky == True :
    # make gif centered on ra=0,dec=0, all sky
        n_plots = makeObservingPlots(
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, allSky=True)

def makeObservingPlots( gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results, allSky=True) :
    trigger_id = gw_map_trigger.trigger_id
    camera = gw_map_strategy.camera

    data_dir = gw_map_control.datadir
    start_slot = gw_map_control.start_slot
    do_nslots = gw_map_control.do_nslots

    # we are doing a quick run, avoiding most calculations as they are already done and on disk
    if (gw_map_results.hoursPerNight == False)  :
        reuse_results(data_dir, gw_map_trigger, gw_map_strategy, gw_map_results, get_slots=True)
    n_slots = gw_map_results.n_slots
    first_slot = gw_map_results.first_slot
    best_slot = gw_map_results.best_slot
    probs =   gw_map_results.probability_per_slot 
    times =    gw_map_results.time_of_slot 
    isdark =    (gw_map_results.isdark).astype(int)

    if n_slots == 0:
        nothingToObserveShowSomething(trigger_id, data_dir)
        return

    print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    print "makeObservingPlots(",n_slots, trigger_id, best_slot,data_dir," )"
    print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
    print "We're going to do {} slots with best slot={}".format(np.nonzero(isdark)[0].size, best_slot)
    figure = plt.figure(1,figsize=(8.5*1.618,8.5))

    print "\t cleaning up old png files"
    files = glob.glob(data_dir+"/*png"); 
    for f in files: os.remove(f)

    # first, make the probability versus something plot
    ra,dec,id,prob,slotMjd,slotNumbers,dist = obsSlots.readObservingRecord(
        trigger_id, data_dir)

    # now make the hex observation plots
    counter = 1   ;# already made one
    #for i in np.unique(slotNumbers) :
    for i in range(isdark.size) :
        if not isdark[i]: continue
        
        observingPlot(figure,trigger_id,i,data_dir, n_slots, camera, allSky=allSky)
        name = str(trigger_id)+"-observingPlot-{}.png".format(i)
        plt.savefig(os.path.join(data_dir,name))
        counter += 1
        counter+= equalAreaPlot(figure,i,trigger_id,data_dir)

    label = ""
    if allSky == False: label="centered_"

    string = "$(ls -v {}-observingPlot*)  {}_{}animate.gif".format(data_dir+'/'+trigger_id, data_dir+'/'+trigger_id, label)
    print string
    os.system("convert  -delay 40 -loop 0  " + string)

    # return the number of plots made
    return counter
#
# ===== its a disaster, compute something
#
def nothingToObserveShowSomething(simNumber, data_dir) :
    import matplotlib.pyplot as plt
    print "nothingToObserveShowSomething",
    figure = plt.figure(1,figsize=(8.5*1.618,8.5))
    ra,dec,id,prob,mjd,slotNum,dist=obsSlots.readObservingRecord(simNumber, data_dir)
    ix = np.argmax(prob)
    try:
        slot = np.int(slotNum[ix])
        title = "slot {} hex maxProb {:.6f}, nothing to observe".format(slot, prob[ix])
    except:
        print "nothingToObserveShowSomething: fiasco."
        print "ra",ra
        print "dec",dec
        print "id",id
        print "prob",prob
        print "mjd",mjd
        print "slotNum", slotNum
        print "ix",ix
        print "ix index",np.nonzero(ix)
        print "nothingToObserveShowSomething: fiasco. attempting to go on"
        slot = 0
        #raise Exception("this failed once, but ran fine when I did it in the directory. NSF thing?")
        title = "slot {} hex maxProb {:.6f}, nothing to observe".format("nothing", 0.0)
    counter = equalAreaPlot(figure,slot,simNumber,data_dir,title) 
    return counter

def how_well_did_we_do(gw_map_trigger, gw_map_strategy, gw_map_control) :
    skymap = gw_map_trigger.skymap
    trigger_id = gw_map_trigger.trigger_id
    camera = gw_map_strategy.camera
    resolution = gw_map_control.resolution
    data_dir = gw_map_control.datadir

    ra = gw_map_trigger.ligo_ra
    dec = gw_map_trigger.ligo_dec
    ligo = gw_map_trigger.ligo
    name = os.path.join(data_dir, str(trigger_id) + "-ra-dec-id-prob-mjd-slot-dist.txt")
    raH, decH = np.genfromtxt(name, unpack=True, usecols=(0,1))
    treedata = decam2hp.buildtree(ra, dec, resolution, recompute=True) 
    tree = treedata[2] 
    sum = decam2hp.hexalateMap(ra,dec,ligo,tree, raH,decH,camera) 
    print "\nTotal Ligo probability covered by hexes observed: {:.3f}%".format(sum.sum()*100.)
    print "   (from decam2hp.hexalateMap of -ra-dec-id file)\n"
    return sum.sum()


# ========== do simple calculations on how to divide the night
#
#   one of the questions is how many hours to devote to observations
#       hoursAvailable,  another is the slot duration
#
# if the 1% cut isn't in place in mapsAtTimeT.oneDayOfTotalProbability
# then one expects circa 40-45 maps as there is about 2/hour
#   with the 1% cut, then one expects far fewer. Perhaps zero.
#
def make_divisions_of_time (
        probs, times, hoursPerNight= 10., hoursAvailable=6) :
    slotDuration = 30. # minutes
    if hoursAvailable > hoursPerNight:
        hoursAvailable = hoursPerNight
        
    # if the number of slots is zero, nothing to observe or plot
    if np.size(times) == 0 : return 0,0
    if probs.sum() < 1e-9 : return 0,0
    verbose = 0
    n_slots = obsSlots.findNSlots(hoursAvailable,slotDuration=slotDuration)
    n_maps = times.size
    if verbose: print n_slots, n_maps
    if n_maps == n_slots : 
        mapZero = 0
    elif n_maps < n_slots : 
        mapZero = 0
        n_slots = n_maps
    elif n_maps > n_slots :
        mapZero = obsSlots.findStartMap ( probs, times, n_slots )
    else :
        raise Exception ("no possible way to get here")
    print "=============>>>>  make_divisions_of_time:"
    print "\t n_maps = {}, n_slots = {}, mapZero = {}, prob_max = {:.6}".format(
        n_maps, n_slots, mapZero, probs.max())
    return n_slots, mapZero

# make_hexes computes slot information, so get_slots = false
# make_observingPlots needs slot information, so get_slots = true
def reuse_results(data_dir, gw_map_trigger,gw_map_strategy, gw_map_results, get_slots=False) :
        trigger_id = gw_map_trigger.trigger_id
        probabilityTimesCache = os.path.join(data_dir,\
        "probabilityTimesCache_"+str(trigger_id)+".txt")
        print "=============>>>> Reuse results via reuse_results"
        print "\t using probabilities, times, and maps from", probabilityTimesCache
        if os.stat(probabilityTimesCache).st_size == 0 :
            probs, times = np.array([0,]),np.array([0,])
            print "\t got nothing- we're skipping this one"
        else :
            data = np.genfromtxt(probabilityTimesCache, unpack=True)
            probs, times,isdark = data[0],data[1],data[2]
        gw_map_results.probability_per_slot = probs
        gw_map_results.time_of_slot = times
        gw_map_results.isdark = isdark

        exposure_list   = gw_map_strategy.exposure_list
        burst_mjd       = gw_map_trigger.burst_mjd 
        maxHexesPerSlot = gw_map_strategy.maxHexesPerSlot
        overhead        = gw_map_strategy.overhead
        answers = obsSlots.slotCalculations(burst_mjd, exposure_list, overhead,
            maxHexesPerSlot)
        hoursPerNight = answers["hoursPerNight"] ;# in minutes
        slotDuration = answers["slotDuration"] ;# in minutesk
        gw_map_results.slotDuration = slotDuration
        gw_map_results.hoursPerNight = hoursPerNight
        # but be aware that make_divisions of time has a hardwired slotDuration

        if get_slots:
            hoursObserving = pickle.load(open("{}-hoursObserving.pickle".format(trigger_id),"r"))
            gw_map_results.n_slots = hoursObserving["nslots"]
            gw_map_results.first_slot = hoursObserving["mapZero"]
            gw_map_results.best_slot = hoursObserving["maxSlot"]

# ===== The economics analysis
#
#   area_left is th enumber of hexes we have left to observe this season
#   days_left is the number of days left in the season
#   rate is the effective rate of triggers
#       p_gw is that for which the table cumul_table_pgw50.txt was  made.
#
def economics (simNumber, best_slot, mapDirectory, 
        area_left=200., days_left=60., rate=1/30.,  p_gw = 0.10) :
    import healpy as hp
    import cumul
    import des_optimization
    import os
    gw_data_dir = os.environ["DESGW_DATA_DIR"]
    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(mapDirectory, simNumber, best_slot)

    area_bar_p,area_bar = np.genfromtxt(
        gw_data_dir+"/area_bar_table.txt",unpack=True)
    avge_cumu_area,avge_cumu = np.genfromtxt(
        gw_data_dir+"/cumul_table_pgw10.txt",unpack=True)

    obsProb = ligo*prob
    nsides = hp.get_nside(obsProb)
    # max area viewable by Blanco at one time is 11734. sq-degrees
    max_area=11734.
    area, cum_prob  = cumul.area(ra,dec,obsProb, p_gw, nsides, max_area=max_area)
    area_to_cover_p_gw = area
    #print avge_cumu_area
    #print area
    ix = np.searchsorted(avge_cumu_area, area)
    if ix >= avge_cumu_area.size :
        fraction_of_sims_better_than_this_trigger = 1.0
    else :
        fraction_of_sims_better_than_this_trigger = avge_cumu[ix]

    prob, N_max = des_optimization.evaluate_average_event(
        area_left, days_left, rate, avge_cumu, avge_cumu_area, area_bar, area_bar_p)

    if fraction_of_sims_better_than_this_trigger < 1./N_max :
        area, cum_prob = cumul.area(ra,dec,obsProb, prob, nsides, max_area=max_area)
        if area>area_left:
            print "\t maxing out area: \t {:.3f} -> ".format( cum_prob),
            cum_prob = cumul.probability_covered(ra,dec,obsProb, area_left, nsides, max_area=max_area)
            print "{:.3f}".format(cum_prob)
            area=area_left
    else :
        print "\t ignore event"
        area = 0
        prob = 0

    probability_covered = cum_prob
    quality = fraction_of_sims_better_than_this_trigger 
    return probability_covered, area, area_to_cover_p_gw, quality

#
# no, no, no, we actually can see something: lets see the best plots
#
#   raMap, decMap, ligoMap, maglimMap, probMap, haMap, xMap,yMap, hxMap,hyMap = readMaps(
#   ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = readMaps(
def readMaps(mapDir, simNumber, slot) :
    import healpy as hp
    # get the maps for a reasonable slot
    name = os.path.join(mapDir, str(simNumber) + "-"+str(slot))
    print "\t reading ",name+"-ra.hp  & etc"
    raMap     =hp.read_map(name+"-ra.hp", verbose=False);
    decMap    =hp.read_map(name+"-dec.hp", verbose=False);
    haMap     =hp.read_map(name+"-ha.hp", verbose=False);
    xMap      =hp.read_map(name+"-x.hp", verbose=False);
    yMap      =hp.read_map(name+"-y.hp", verbose=False);
    hxMap     =hp.read_map(name+"-hx.hp", verbose=False);
    hyMap     =hp.read_map(name+"-hy.hp", verbose=False);
    ligoMap   =hp.read_map(name+"-map.hp", verbose=False);
    maglimMap =hp.read_map(name+"-maglim.hp", verbose=False);
    probMap   =hp.read_map(name+"-probMap.hp", verbose=False);
    haMap=haMap/(2*np.pi/360.)
    raMap=raMap/(2*np.pi/360.)
    decMap=decMap/(2*np.pi/360.)
    return raMap, decMap, ligoMap, maglimMap, probMap, \
        haMap, xMap, yMap, hxMap, hyMap

#========================================================================
# 
# support routines
# 
#========================================================================
# for economics analysis
def time_cost_per_hex (nvisits, overhead, exposure_length) :
    tot_exptime = (np.array(overhead)+np.array(exposure_length)).sum
    time_cost_per_hex = nvisits * tot_exptime #sec
    return time_cost_per_hex
    
# for economics analysis
def area_left (area_per_hex, time_budget, time_cost_per_hex) :
    area_per_hex * (time_budget * 3600)/(time_cost_per_hex)
    return area_per_hex
#time_cost_per_hex = nvisits * nexposures * (overhead + exposure_length) #sec
#area_left =  area_per_hex * (time_budget * 3600)/(time_cost_per_hex)

# place holder for the code brought from desisurvey...
def hoursPerNight (mjd) :
    import mags
    night,sunset,sunrise = mags.findNightDuration(mjd)
    night = night*24.
    return night


#===================================================================
#
    
def turnObservingRecordIntoJSONs(
        ra,dec,id,prob,mjd,slotNumbers, simNumber, 
        exposure_list, filter_list, trigger_type, mapDirectory, propid) :
    seqtot =  ra.size
    seqzero = 0

    # write slot json files
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        slotMJD = mjd[ix][0]  ;# just get first mjd in this slot
        tmpname, name = jsonUTCName(slot, slotMJD, simNumber, mapDirectory)
        jsonMaker.writeJson(ra[ix],dec[ix],id[ix],
            simNumber, seqzero, seqtot, exposureList= exposure_list, propid=propid,
            filterList= filter_list, trigger_type=trigger_type, jsonFilename=tmpname)

        desJson(tmpname, name, mapDirectory) 
        seqzero += ra[ix].size
        
# verbose can be 0, 1=info, 2=debug
def desJson(tmpname, name, data_dir, verbose = 1) :
    import os
    import logging
    from collections import defaultdict
    import gwwide
    gw_data_dir          = os.environ["DESGW_DATA_DIR"]
    des_json  = gw_data_dir + "all_wide_sispi_queue.json"
    log_levels = defaultdict(lambda:logging.WARNING)
    log_levels[1] = logging.INFO
    log_levels[2] = logging.DEBUG
    logging.basicConfig(filename=data_dir+"json-conversion.log",
        format='%(asctime)s %(message)s', level=verbose)

    logging.info("Begin")
    gwwide.file_gwwide(tmpname, des_json, name)

def jsonUTCName (slot, mjd, simNumber, mapDirectory) :
    time = utcFromMjd(mjd)
    tmpname, name = jsonName(slot, time, simNumber, mapDirectory)
    return tmpname, name

def utcFromMjd (mjd) :
    year,month,day,hour,minute = utc_time_from_mjd(mjd)
    time = "UTC-{}-{:02d}-{:02d}-{:02d}:{:02d}:00".format(year,month,day,hour,minute)
    return time

def utc_time_from_mjd (mjd) :
    from pyslalib import slalib
    date = slalib.sla_djcl(mjd)
    year = np.int(date[0])
    month= np.int(date[1])
    day = np.int(date[2])
    hour = np.int(date[3]*24.)
    minute = np.int( (date[3]*24.-hour)*60.  )
    return year, month, day, hour, minute

def jsonName (slot, utcString, simNumber, mapDirectory) :
    slot = "-{}-".format(np.int(slot))
    tmpname = os.path.join(mapDirectory, str(simNumber) + slot + utcString + "-tmp.json")
    name = os.path.join(mapDirectory, str(simNumber) + slot + utcString + ".json")
    return tmpname, name

def jsonFromRaDecFile(radecfile, nslots, slotZero, 
        hexesPerSlot, simNumber, mjdList, trigger_type, data_dir) :
    ra,dec = np.genfromtxt(radecfile, unpack=True,usecols=(0,1),comments="#",propid="propid")

    seqtot =  ra.size
    seqzero = 0

    # instead, just reorder the ra,dec before feeding to this routine
    #ix = np.argsort(ra)

    counter = 0
    slot = slotZero
    slotRa = np.array([])
    slotDec = np.array([])
    for i in range(0,ra.size) :
        slotRa = np.append(slotRa, ra[i])
        slotDec = np.append(slotDec, dec[i])
        counter += 1
        if counter == hexesPerSlot :
            tmpname, name = jsonName(slot, mjdList[slot-slotZero], 
                simNumber,data_dir)
            jsonMaker.writeJson(slotRa,slotDec, 
                simNumber, seqzero+(hexesPerSlot*(slot-slotZero)), 
                seqtot, trigger_type, jsonFilename=tmpname, propid=propid)
            desJson(tmpname, name, data_dir) 
            counter = 0
            slot += 1
            slotRa = np.array([])
            slotDec = np.array([])
    if counter > 0 :
        tmpname, name = jsonName(slot, mjdList[slot-slotZero], 
            simNumber,data_dir)
        jsonMaker.writeJson(slotRa,slotDec, 
            simNumber, seqzero+(hexesPerSlot*(slot-slotZero)), 
            seqtot, trigger_type, jsonFilename=tmpname, propid=propid)
        desJson(tmpname, name, data_dir) 
        
# ==================================
# plotting 
# ==================================

def cumulPlot(trigger_id, data_dir) :
    from scipy.interpolate import interp1d
    ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)
    u_slotNum = np.unique(slotNum)
    u_dur = np.unique(mjd)
    duration = (u_dur[1]-u_dur[0])*24*60.

    new_x, new_y, new_mjd, new_hour = np.array([]), np.array([]), np.array([]), np.array([])
    for u in u_slotNum :
        ix = slotNum == u
        y = prob[ix].sum()
        new_x = np.append(new_x, u)
        new_y = np.append(new_y, y)
        new_mjd = np.append(new_mjd, mjd[ix].min())
        year,month,day,hour,minute = utc_time_from_mjd(mjd[ix].min())
        new_hour = np.append(new_hour, hour+minute/60.)

    plt.close()
    fig, ax1 = plt.subplots() ; 
    ax2 = ax1.twinx()
    ax3 = ax1.twiny()

    ax1.scatter(new_x,new_y,c="k", s=15, zorder=10)
    ax1.plot(new_x,new_y,c="b", zorder=9)
    ax1.set_xlabel("slot number")
    ax1.set_ylabel("LIGO probability in slot")
    ax1.set_ylim(0,new_y.max()*1.1)
    ax1.grid(alpha=0.2)

    plotProb, plotSlot,plotN = np.array([]), np.array([]), np.array([])
    for i in np.unique(slotNum) :
        ix = slotNum == i
        if prob[ix].sum() > 0 :
            plotN = np.append(plotN, prob[ix].size)
            plotSlot = np.append(plotSlot,i)
            plotProb = np.append(plotProb,100.*prob[ix].sum())
    ax1.text(0.99,0.02,"total probability: {:5.1f}%".format(prob.sum()*100.),
        transform = ax1.transAxes,   horizontalalignment='right',
        verticalalignment='center')
    avghex = str( np.round(plotN.mean(),1) )
    ax1.text(0.99,0.07,"slot duration: {:5.1f}$^m$".format(duration),
        transform = ax1.transAxes,   horizontalalignment='right',
        verticalalignment='center')
    ax1.text(0.99,0.12,"hexes per slot: {}".format(avghex),
        transform = ax1.transAxes,   horizontalalignment='right',
        verticalalignment='center')

    ix = slotNum == u_slotNum[0]
    year,month,day,hour,minute = utc_time_from_mjd(mjd[ix].min())
    new_cy = np.cumsum(new_y)
    ax2.scatter(new_x,new_cy,color="k",s=15,zorder=10)
    ax2.plot(new_x,new_cy,color="g",zorder=9)
    ax2.tick_params(axis='y', labelcolor="g")
    ax2.set_ylabel("cumulative LIGO probability",color="g")
    ax2.set_ylim(0,new_cy.max())
    fig.tight_layout()
    qtiles = new_cy/new_cy[-1]
    interp = interp1d(qtiles, new_x, fill_value="extrapolate")
    for q in [0.25, 0.5, 0.75] :
        ax2.plot([interp(q),interp(q)], [0, new_cy.max()],alpha=0.3,c="r", ls="dotted")
        ax2.text(interp(q), new_cy.max()*0.95, "{:2d}%".format(int(q*100)), 
            verticalalignment='bottom', alpha=0.3, color="r")

    mjd_h = (new_mjd - new_mjd[0])*24.
    interp = interp1d(mjd_h,new_x)
    x_ticks = interp( np.arange(0, mjd_h.max(), 1) )

    mjd_24 = (new_mjd - new_mjd[0])*24. + new_hour[0]
    ix = mjd_24 > 24
    mjd_24[ix] = mjd_24[ix]-24.0
    interp = interp1d(mjd_h, mjd_24)
    hours_labels = interp( np.arange(0, mjd_h.max(), 1) )
    labels = []
    for i in range(hours_labels.size) :
        labels.append( "{:.1f}".format(hours_labels[i]) )
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticks(x_ticks)
    ax3.set_xticklabels(labels)
    ix = slotNum == u_slotNum[0]
    year,month,day,hour,minute = utc_time_from_mjd(mjd[ix].min())
    ix = slotNum == u_slotNum[-1]
    eyear,emonth,eday,ehour,eminute = utc_time_from_mjd(mjd[ix].max())
    ax3.set_xlabel("UT hour, starting {}/{}/{}, ending {}/{}".format(
        year, month, day, emonth, eday))
    print "\t writing {}-slot-probabilities.png".format(trigger_id)
    plt.savefig("{}-slot-probabilities.png".format(trigger_id))



def cumulPlot2(trigger_id, data_dir) :
    from scipy.interpolate import interp1d
    ra,dec,id,prob,mjd,slotNum,dist = obsSlots.readObservingRecord(trigger_id,data_dir)
    ix = np.argsort(prob)[::-1]
    ra,dec,id,prob,mjd,slotNum,dist = \
        ra[ix],dec[ix],id[ix],prob[ix],mjd[ix],slotNum[ix],dist[ix]  
    
    plt.close()
    fig, ax1 = plt.subplots() ; 
    ax2 = ax1.twinx()

    x = range(0,prob.size)
    ax1.scatter(x,prob,c="b",s=15)
    ax1.plot(x,prob,c="k",zorder=100)
    ax1.set_xlabel("hex count,  max to min")
    ax1.set_ylabel("LIGO probability in slot")
    ax1.set_ylim(0,prob.max()*1.1)
    ax1.grid(alpha=0.2)

    new_cy = np.cumsum(prob)
    ax2.scatter(x,new_cy,color="g",s=15, zorder=90)
    ax2.plot(x,new_cy,c="k",zorder=101)
    ax2.tick_params(axis='y', labelcolor="g")
    ax2.set_ylabel("cumulative LIGO probability",color="g")
    ax2.set_ylim(0,new_cy.max())
    fig.tight_layout()

    plt.title("Hexes ordered by probability")
    qtiles = new_cy/new_cy[-1]
    interp = interp1d(qtiles, x, fill_value="extrapolate")
    for q in [0.25, 0.5, 0.75] :
        ax2.plot([interp(q),interp(q)], [0, new_cy.max()],alpha=0.3,c="r", ls="dotted")
        ax2.text(interp(q), new_cy.max()*0.95, "{:2d}%".format(int(q*100)), 
            verticalalignment='bottom', alpha=0.3, color="r")
    print "\t writing {}-hex-probabilities.png".format(trigger_id)
    plt.savefig("{}-hex-probabilities.png".format(trigger_id))

def equalAreaPlot(figure,slot,simNumber,data_dir, title="") :
    import matplotlib.pyplot as plt
    from equalArea import mcplot
    from equalArea import mcbryde
    import insideDesFootprint

    ra, dec, ligo, maglim, prob, ha, x,y, hx,hy = \
        readMaps(data_dir, simNumber, slot)
    # x,y are the mcbryde projection of ra, dec
    # hx,hy are the mcbryde projection of ha, dec
    ra, dec = x, y

    # des footprint
    # plots the DES footprint on the maps
    desra, desdec = insideDesFootprint.getFootprintRaDec()
    desx, desy = mcbryde.mcbryde(desra, desdec)


    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-ligo-eq.png")
    print "making ",name,
    plt.clf();mcplot.plot(ra,dec,ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-maglim-eq.png")
    print name,
    plt.clf();mcplot.plot(ra,dec,maglim,vmin=17);
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-prob-eq.png")
    print name,
    plt.clf();mcplot.plot(ra,dec,prob)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)

    name = os.path.join(data_dir,str(simNumber)+"-"+str(slot)+"-probXligo-eq.png")
    print name
    plt.clf();mcplot.plot(ra,dec,prob*ligo)
    plt.plot(desx,desy,color="w")
    plt.xlabel("RA");plt.ylabel("Dec")
    plt.title(title)
    plt.savefig(name)
    # return the number of plots made
    return 4 

# modify mcbryde to have alpha=center of plot
#   "slot" is roughly hour during the night at which to make plot
def observingPlot(figure, simNumber, slot, data_dir, nslots, camera, allSky=False) :
    import plotMapAndHex

    # get the planned observations
    ra,dec,id,prob,mjd,slotNumbers,dist = obsSlots.readObservingRecord(simNumber, data_dir)
    ix = slotNumbers == slot
    if np.any(ix):
        the_mjd = mjd[ix][0]
        time = utcFromMjd(the_mjd)
    else :
        time = ""
    
    title = " Slot {}    {} ".format(slot, time)
    title = title + "      {}".format(simNumber)


# this is useful to debug the plots
    #print "making plotMapAndHex.mapAndHex(figure, ", simNumber, ",", slot, ",", data_dir, ",", nslots, ",ra,dec,", camera, title,"allSky=",allSky,") "

    d=plotMapAndHex.mapAndHex(figure, simNumber, slot, data_dir, nslots, ra, dec, \
        camera, title, slots=slotNumbers, allSky=allSky)
    return d

def writeObservingRecord(slotsObserving, data_dir, gw_map_trigger, gw_map_control) :
    trigger_id = gw_map_trigger.trigger_id
    just_sort_by_ra = gw_map_control.just_sort_by_ra

    name = os.path.join(data_dir, str(trigger_id) + "-ra-dec-id-prob-mjd-slot-dist.txt")
    ra,dec,id,prob,mjd,slotNum,islot = obsSlots.slotsObservingToNpArrays(slotsObserving)
    if just_sort_by_ra :
        # rearrange inside the slots if desired
        ix = np.argsort(ra)
        ra,dec,id,prob,islot = \
            ra[ix],dec[ix],id[ix],prob[ix],islot [ix]

    map_distance = gw_map_trigger.ligo_dist
    nside = hp.npix2nside(map_distance.size)
    #ang2pix(nside,theta,phi,nest=False)
    pix_num = hp.ang2pix(nside,ra,dec, lonlat=True)
    dist = map_distance[pix_num]
    #fixing the distance
    #dist = 60.
    fd = open(name,'w')
    fd.write("# ra, dec, id, prob, mjd, slotNum, dist\n")
    unique_slots = np.unique(slotNum)
    for slot in unique_slots:
        ix = slot==slotNum
        iy = np.argsort(ra[ix])
        # sort by ra inside the slot
        if (np.any(ra[ix] > 90) and np.any(ra[ix] < -90) ) :
            dummy_ra = ra[ix]
            iy = dummy_ra < 0
            dummy_ra[iy] = dummy_ra[iy] + 360
            iy = np.argsort(dummy_ra)
        else :
            iy = np.argsort(ra[ix])
        for r,d,i,p,m,s,di in zip(
                ra[ix][iy], dec[ix][iy], id[ix][iy],
                prob[ix][iy], mjd[ix][iy], slotNum[ix][iy], dist[ix][iy]):
                #ra[ix], dec[ix], id[ix],
                #prob[ix], mjd[ix], slotNum[ix], dist[ix]):
            fd.write("{:.6f} {:.5f} {:s} {:.7f} {:.4f} {:.1f} {:.2f}\n".format(r,d,i,p,m,s,di))
    fd.close()
    return ra,dec,id,prob,mjd,slotNum,dist

