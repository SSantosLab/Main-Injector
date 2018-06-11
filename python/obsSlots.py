import numpy as np
import os

# ========== do simple calculations on how to divide the night
#
#   one of the questions is how many hours to devote to observations
#       hoursAvailable,  another is the slot duration
#
# if the 1% cut isn't in place in mapsAtTimeT.oneDayOfTotalProbability
# then one expects circa 40-45 maps as there is about 2/hour
#   with the 1% cut, then one expects far fewer. Perhaps zero.

#===============================================================================
# These calculations used to be spread over hither and yon.
# Bring them together.
#
#    deltaTime = 0.0223  => 32 minute slots  for izzi 90sec exp, 4 hex/slot
#    deltaTime = 0.0446  => 62 minute slots  for izzi 90sec exp, 8 hex/slot
#    deltaTime = 0.0417  => 60 minute slots  for izz 90sec exp, 10 hex/slot
#    deltaTime = 0.0208  => 30 minute slots  for izz 90sec exp,  5 hex/slot
#       The first is flexible to observing conditions,
#       while the second is x2 faster as there are less maps to hexalate
#       The third shows a different tiling/filter complement
#
# So-
#
# Let us redfine this: a slot is the length of time it takes
# to do 6 hexes to completion. That is usually somewhere between 30 minutes
# and one hour, so close to the original defintion, and by force is an
# even number of hexes. Ok. Use n=6 for the forcing definition
#
# This is based on the current NS observing strategy, izz at 90s each.
# If instead we move to a BH strategy that is 1xi at 90s, then
# one could still use 6 hexes/slot but that means the number of slots rises by x3.
# Instead, what if we used 18 hexes/slot for the BH case.
#
# Ok, then, code:
#
def slotCalculations(mjd, exposure_lengths, overhead, hexesPerSlot = 6) :
    tot_exptime = (np.array(overhead)+np.array(exposure_lengths)).sum()
    slot_time = tot_exptime*hexesPerSlot
    slot_duration = slot_time/60. ;# in minutes
    from getHexObservations import hoursPerNight
    hoursAvailable = hoursPerNight(mjd)
    answers = dict()
    answers["slotDuration"] = slot_duration
    answers["hoursPerNight"] = hoursAvailable
    return answers

# find the number of slots per night
def findNSlots(hoursAvailable, slotDuration=32.) :
    verbose = 0
    if verbose:
        print hoursAvailable
        print hoursAvailable*60./slotDuration, round(hoursAvailable*60./slotDuration)
        print int(round(hoursAvailable*60./slotDuration))
    nslots = int(round(hoursAvailable*60./slotDuration))   ;# 32 minutes/slot
    return nslots

# ok, I designed observing for the case
#   where the number of slots in a night
#   was equal to the number of maps made.
# This is unrealistic for two reasons:
#   the night could be longer than the time desgw wishes to allocate
#   the object could be up for less than the time desgw could allocate
# so:
    # possibilities
    # n_maps = n_slots     design case
    # n_maps < n_slots     set n_slots = n_maps
    # n_maps > n_slots     pick highest contiguous set of n_slot maps
# the first two are easy. 
# the third runs into the naming convention of the maps
def findStartMap ( probs, times, n_slots ) :
    n_maps = times.size
    mapNums = np.arange(0,n_maps)
    # sum the probability in each allowed nslot range
    n_mapsToTProb = np.array([])
    n_mapStart = np.array([])
    for map in range(0,n_maps-n_slots) :
        ix = (mapNums >= map) & (mapNums < map+n_slots)
        n_mapsToTProb = np.append(n_mapsToTProb, probs[ix].sum() )
        n_mapStart = np.append(n_mapStart, map)
    minToMax = np.argsort(n_mapsToTProb)
    bestStart = n_mapStart[minToMax[-1]]
    bestStart = int(bestStart)
    return bestStart


# Load all n -hexVals files,
#   pick the highest probability one, 
#   put it into one of the n slots, unless that slots is maxed out
#   remove that hex from all the hexVals lists
#   do it again, untill all n time slots are full.
#       maxHexesPerSlot=4 comes from 32 minute duration slots
#       and 8 minutes/hex (izzi 2 min/image)
#
#   if we do zzi at 2 mins/image then 4 min/hex + 2 min/hex2 = 6 mins
#   call it 60 minute slots  and 10 hexes/slot
def observing(sim, nslots, data_dir, 
        maxHexesPerSlot = 4, mapZero = 0, verbose=0) :
    # prep the observing lists
    observingSlots = np.arange(0,nslots)
    slotsObserving = dict()
    slotsObserving["nslots"] = nslots
    for i in observingSlots :
        slotsObserving[i] = 0
        slotsObserving[i,"ra"]   = np.array([])
        slotsObserving[i,"dec"]  = np.array([])
        slotsObserving[i,"id"]  = np.array([])
        slotsObserving[i,"prob"] = np.array([])
        slotsObserving[i,"mjd"] = np.array([])
        slotsObserving[i,"slotNum"] = np.array([]) ;#non-zero prob slots
        slotsObserving[i,"islot"] = np.array([]) ;# range(0,nslots)

    # read in the hexelated probability data
    hexData = dict()
    for i in observingSlots :
        map_i = i + mapZero
        raHexen, decHexen, idHexen, hexVal, rank, mjd, slotNum = \
           loadHexalatedProbabilities( sim, map_i, data_dir)
        islot = i*np.ones(raHexen.size)
        print "\t", map_i, "map size= {};".format(raHexen.size), 

        impossible = 1e-5
        impossible = 1e-7
        ix = np.nonzero(hexVal < impossible)
        raHexen, decHexen, idHexen, hexVal, mjd, slotNum, islot  = \
            np.delete(raHexen, ix), \
            np.delete(decHexen, ix), \
            np.delete(idHexen, ix), \
            np.delete(hexVal, ix), \
            np.delete(mjd, ix) , \
            np.delete(slotNum, ix), \
            np.delete(islot, ix)
        print " n hexes >{} probability=".format(str(impossible)),
        print "{:4d};".format(raHexen.size),
        print "  sum prob= {:7.4f} %".format( 100*hexVal.sum())
        hexData[i] = raHexen, decHexen, idHexen, hexVal, mjd, slotNum, islot
        #print np.sort(hexVal), hexVal.sum(), 100.*hexVal.sum(),"%"

    # start the search for all max probabilities
    # we'll assume the list is less than 40,000 long, the n-sq-degrees/sky
    for n in range(0,40000) :
        # search for a single max probabilities
        maxRa, maxDec, maxId, maxProb, maxMjd, maxSlotNum, maxIslot  = \
            findMaxProbOfAllHexes(hexData, observingSlots, n, verbose) 
        maxData = maxRa,maxDec,maxId, maxProb,maxMjd,maxSlotNum, maxIslot

        # we've found the maximum probability on the lists, 
        # so add it to the obs lists # unless not possible. 
        # If the latter, delete it from that slot
        slot = maxIslot
        if slotsObserving[slot] < maxHexesPerSlot : 
            # it is possible to make the observation, 
            # put it onto the observing lists
            slotsObserving = addObsToSlot (slotsObserving, maxData, slot)
            if verbose >= 1: print n, "slot of max:",slot
        else :
            # but if this slot of observing is full, it is not possible 
            # to make the observation,
            # so move on AFTER deleting it from the list
            hexData = deleteHexFromSlot (hexData, slot, maxProb) 
        #if verbose >= 2: 
        #   if n > 7: raise Exception("jack")
    
        # perform the necessary bookkeeping, 
        # eliminating this hex from future observing
        hexData = deleteHexFromAllSlots (
            hexData, observingSlots, maxRa, maxDec, verbose, n) 

        # do some summary statistics
        sumHexes = 0
        sumObs = 0
        for i in range(0,nslots) :
            sumHexes += hexData[i][0].size
            sumObs += slotsObserving[i]

        if verbose >= 2: 
            print "sumHexes =", sumHexes, 
            print "   slots left=", len(observingSlots),
            print "   slots=",observingSlots,
            print "   n_obs=",
            for i in observingSlots:
                print slotsObserving[i],
            print "   sum prob= ",
            for i in observingSlots:
                print " {:8.6f}".format( slotsObserving[i,"prob"].sum()) ,
            print ""

        # eliminate full observing slots
        observingSlots = eliminateFullObservingSlots(\
            hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) 

        # check and see if we are done
        # two conditions: observing is full, or candidates empty
        if (len(observingSlots)==0) | (sumHexes == 0) :
            print "\t======================================== "
            if verbose >= 1: 
                print "n slots =", len(observingSlots)," == 0?"
                print "sumHexes = ", sumHexes, "==? 0"
            print "\tnumber of hexes observed = ", sumObs
            print "\t======================================== "
            return slotsObserving 

        # otherwise, go back and do it again
    
    # we've done everything on the lists, we can observe it all,
    # return this great success that will never be reached.
    return slotsObserving
    
#
# examine the statistics of the observing lists
#
def observingStats( slotsObserving ) :
    nslots = slotsObserving["nslots"]
    for i in range(0,nslots) :
        print "\t",i, 
        #print "slotnum={} ".format( slotsObserving[i,"slotNum"]),
        print "n obs= {}".format( slotsObserving[i,"ra"].size), 
        print "  sum prob= {:7.4f} %".format( 100*slotsObserving[i,"prob"].sum())
    ra,dec,id,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 

    print "\tobservingStats:  ",
    print "observable prob_tot = {:.1f}%".format(100.*prob.sum())
    return ra,dec,id,prob,mjd,slotNum,islot

def observingRecord(slotsObserving, simNumber, data_dir) :
    name = os.path.join(data_dir, str(simNumber) + "-ra-dec-id-prob-mjd-slot.txt")
    ra,dec,id,prob,mjd,slotNum,islot = slotsObservingToNpArrays(slotsObserving) 
    data = np.array([ra, dec, id, prob, mjd, slotNum]).T
    f = open(name,'w')
    unique_slots = np.unique(slotNum)
    for slot in unique_slots:
        ix = slot==slotNum
        iy = np.argsort(ra[ix])
        for r,d,i,p,m,s in zip(
                ra[ix][iy], dec[ix][iy], id[ix][iy], 
                prob[ix][iy], mjd[ix][iy], slotNum[ix][iy]):
            f.write("{:.6f} {:.5f} {:s} {:.7f} {:.4f} {:.1f}\n".format(r,d,i,p,m,s))
    f.close()
    #np.savetxt(name, data, "%.6f %.5f %s %.6f %.4f %d")
    return ra,dec,id,prob,mjd,slotNum

#     ra,dec,id,prob,mjd,slotNum,islot = readObservingRecord(simNumber, data_dir)
def readObservingRecord(simNumber, data_dir) :
    import os
    name = os.path.join(data_dir, str(simNumber) + "-ra-dec-id-prob-mjd-slot.txt")
    if not os.path.exists(name) or os.stat(name).st_size == 0 :
        ra,dec,id,prob,mjd,slotNum = \
            np.array(0),np.array(0),np.array("0"), \
            np.array(0),np.array(0),np.array(0)
    else :
        ra,dec,prob,mjd,slotNum = np.genfromtxt(name,unpack=True,comments="#",usecols=(0,1,3,4,5))
        id = np.genfromtxt(name,unpack=True,comments="#", usecols=(2),dtype="str")
    return ra,dec,id,prob,mjd,slotNum

def slotsObservingToNpArrays(slotsObserving) :
    nslots = slotsObserving["nslots"]
    ra = np.array([])
    dec = np.array([])
    id = np.array([])
    prob = np.array([])
    mjd = np.array([])
    slotNum = np.array([])
    islot = np.array([])
    for i in range(0,nslots) :
        ra = np.append(ra, slotsObserving[i,"ra"])
        dec = np.append(dec, slotsObserving[i,"dec"])
        id = np.append(id, slotsObserving[i,"id"])
        prob = np.append(prob, slotsObserving[i,"prob"])
        mjd = np.append(mjd, slotsObserving[i,"mjd"])
        slotNum = np.append(slotNum, slotsObserving[i,"slotNum"])
        islot = np.append(islot, slotsObserving[i,"islot"])
    return ra,dec,id,prob,mjd,slotNum,islot
    

#
# search for the single highest probability hex over all of the possible hexes
# in the hexData slots 
#
def findMaxProbOfAllHexes(hexData, observingSlots, n="", verbose = 0) :
    maxProb = -1
    for i in observingSlots :
        data = hexData[i]
        hexRa     = data[0]
        hexDec    = data[1]
        hexId     = data[2]
        hexVal    = data[3]
        hexMjd    = data[4]
        hexMyslot = data[5]
        if hexVal.size == 0: continue
        if verbose >= 2: 
            if i == 2: print n,"====",i, "hexSize =",hexRa.size
        # now check for max prob
        newProb = hexVal.max()
        if verbose >= 4: print n,i, maxProb, ">?", newProb, "     n=",hexVal.size
        if newProb > maxProb :
            if verbose >= 1: print n,"==== new max", i, "       ",newProb , ">", maxProb
            ix = hexVal == newProb
            maxRa     = hexRa[ix]
            maxDec    = hexDec[ix]
            maxId     = hexId[ix]
            maxVal    = hexVal[ix]
            maxMjd    = hexMjd[ix]
            maxProb   = newProb
            maxSlot   = hexMyslot[ix]
            islot = i
    #print "observingSlots", observingSlots
    if maxProb == -1 : 
        maxRa, maxDec, maxId, maxVal, maxMjd, maxSlot, islot = \
            0,0,0,0,0,0,0
        #raise Exception("no max probability found")
    return maxRa, maxDec, maxId, maxVal, maxMjd, maxSlot, islot

# we've found a hex,slot that can be observed so add it the the observing lists
def addObsToSlot (slotsObserving, maxData, slot) :
    maxRa  = maxData[0]
    maxDec = maxData[1]
    maxId  = maxData[2]
    maxVal = maxData[3]
    maxMjd = maxData[4]
    maxSlotNum = maxData[5]
    maxIslot = maxData[6]
    slotsObserving[slot,"ra"]   =  np.append(slotsObserving[slot,"ra"], maxRa)
    slotsObserving[slot,"dec"]  =  np.append(slotsObserving[slot,"dec"], maxDec)
    slotsObserving[slot,"id"]   =  np.append(slotsObserving[slot,"id"], maxId)
    slotsObserving[slot,"prob"] =  np.append(slotsObserving[slot,"prob"], maxVal)
    slotsObserving[slot,"mjd"]  =  np.append(slotsObserving[slot,"mjd"], maxMjd)
    slotsObserving[slot,"slotNum"]   =  np.append(slotsObserving[slot,"slotNum"], maxSlotNum)
    slotsObserving[slot,"islot"] =  np.append(slotsObserving[slot,"islot"], maxIslot)
    slotsObserving[slot] += 1
    return slotsObserving
# there can be no more observing in this slot, so this hex,slotj
# is impossible, delete it from the list.
def deleteHexFromSlot (hexData, slot, maxProb) :
    hexRa, hexDec, hexId, hexVal, hexMjd, hexSlotNum, hexIslot = hexData[slot] 
    ix = np.nonzero(hexVal == maxProb) 
    hexRa  = np.delete(hexRa, ix)
    hexDec = np.delete(hexDec, ix)
    hexId = np.delete(hexId, ix)
    hexVal = np.delete(hexVal, ix)
    hexMjd = np.delete(hexMjd, ix)
    hexSlotNum = np.delete(hexSlotNum, ix)
    hexIslot   = np.delete(hexIslot, ix)
    hexData[slot] = hexRa, hexDec, hexId, hexVal, hexMjd, hexSlotNum, hexIslot
    return hexData
# the hex,slot has made it onto an observing list, so remove the hex
# from all hex lists ( rm hex,*)
def deleteHexFromAllSlots (hexData, observingSlots, maxRa, maxDec, verbose=0, n="") :
    for i in observingSlots:
        data = hexData[i]
        hexRa  = data[0]
        hexDec = data[1]
        hexId  = data[2]
        hexVal = data[3]
        hexMjd = data[4]
        hexSlotNum = data[5]
        hexIslot   = data[6]
        ix = np.nonzero((hexRa == maxRa) & (hexDec == maxDec))
        if verbose >=4 : print ix, hexRa, maxRa
        if verbose >=2 : 
            ixs = np.shape(ix)[1]
            print n,"bookkeeping",i,"  nHex=",hexRa.size, "   ix.size", ixs, 
        hexRa  = np.delete(hexRa, ix)
        hexDec = np.delete(hexDec, ix)
        hexId  = np.delete(hexId, ix)
        hexVal = np.delete(hexVal, ix)
        hexMjd = np.delete(hexMjd, ix)
        hexSlotNum = np.delete(hexSlotNum, ix)
        hexIslot   = np.delete(hexIslot, ix)
        hexData[i] = hexRa, hexDec, hexId, hexVal, hexMjd, hexSlotNum, hexIslot

        if verbose >= 2: 
            print "\t after delete nHex=",hexRa.size,
            if ixs > 0 : print "index = ",ix[0][0]
            else: print ""
    return hexData

# it is useful to remove full observing slots from further processing,
# though they survive to be returned in the final lists
def eliminateFullObservingSlots(
        hexData, slotsObserving, observingSlots, maxHexesPerSlot, verbose) :
    full = []
    for i in range(0,len(observingSlots)):
        # either we've observed as much as we can, or there is nothing to see
        slot = observingSlots[i]
        if (slotsObserving[slot] >= maxHexesPerSlot) | (hexData[slot][0].size ==0): 
            full.append(i)
    if verbose >= 2: 
        print "hiding full observing slot ", 
        for f in full: print observingSlots[f],
        print ""
    observingSlots = np.delete(observingSlots, full)
    return observingSlots


def maxProbabilitySlot(prob,slotNumbers) :
    # find slot with the maximum probability
    maxProb = -1; maxProb_slot = -1
    for slot in np.unique(slotNumbers) :
        ix = slotNumbers == slot
        sumProb = prob[ix].sum()
        if sumProb > maxProb :
            maxProb = sumProb
            maxProb_slot = slot
    maxProb_slot = np.int(maxProb_slot)
    return maxProb_slot

#===================================================================
#
# Read in all of the hexalated probability files
#   hex inclusion and removal depends on hex id
#
def loadHexalatedProbabilities(sim, slot, data_dir) :
    import hexalate
    nameStem = os.path.join(data_dir, str(sim) + "-{}".format(str(slot)))
    name = nameStem + "-hexVals.txt"
    if os.path.isfile(name) :
        raHexen, decHexen, hexVal, rank, mjd = np.genfromtxt(name, unpack=True, 
            delimiter=",", usecols=(0,1,3,4,5))
        idHexen = np.genfromtxt(name, unpack=True, delimiter=",", usecols=(2), dtype="str")
        slots = np.ones(raHexen.size)*slot
    else :
        raHexen, decHexen, hexVal, rank, mjd, idHexen, slots = \
            7*[np.zeros(0)]

    return raHexen, decHexen, idHexen, hexVal, rank, mjd, slots

