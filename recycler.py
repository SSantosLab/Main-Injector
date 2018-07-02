import os
import sys, getopt, traceback
import numpy as np
#import triggerpages2 as tp
import triggerpagesfinal as tp
import getHexObservations
import subprocess
import datetime
import yaml
import obsSlots
import time
#import jobmanager
import pytz
from threading import Thread
from copy import copy
#sys.path.append("/data/des41.a/data/desgw/")


class event:
    def __init__(self, skymap_filename, outfolder, trigger_id, mjd, config):
        
        #if xml exists then just continue
        #else read in skymap_url.txt and get xml 
        #like so wget --auth-no-challenge https://gracedb.ligo.org/apibasic/events/G268556/files/G268556-4-Initial.xml -O G268556-4-Initial.xml 

        #read in xml and set params

        oldoutfolder = copy(outfolder)
        self.oldoutfolder = oldoutfolder
        skymap_newoutfolder = skymap_filename.split('/')[-1].split('.')[0]
        self.skymap_newoutfolder = skymap_newoutfolder
        #print skymap_filename.strip()+'.processing'
        #raw_input()
        os.system('touch '+skymap_filename.strip()+'.processing')
        os.system('mkdir '+outfolder+'/'+skymap_newoutfolder)
        outfolder = outfolder+'/'+skymap_newoutfolder+'/'
        os.system('cp '+skymap_filename.strip()+' '+outfolder)
        os.system('cp '+oldoutfolder+'/'+trigger_id +'_params.npz '+outfolder)
        print 'cp '+skymap_filename.strip()+' '+outfolder
        print outfolder
        #print skymap_filename
        #print trigger_id
        #raw_input()
        #self.skymap = os.path.join(outfolder, 'lalinference.fits.gz')
        #if not os.path.exists(self.skymap):
        #    self.skymap = os.path.join(outfolder, 'bayestar.fits.gz')
        self.skymap = skymap_filename
        #print 'skymappp'*10
        print self.skymap
        #self.skymap = os.path.join(outfolder, config['default_map_name'])
        self.outfolder = outfolder
        self.trigger_id = trigger_id
        self.mjd = mjd
        self.config = config

        season_start_date = datetime.datetime.strptime(config["start_of_season_date"], "%m/%d/%Y")
        now = datetime.datetime.now()
        
        self.now = now
        if config["force_recycler_mjd"]:
            self.recycler_mjd = config["recycler_mjd"]
        else:
            self.recycler_mjd = self.getmjd(now)
            #self.recycler_mjd = 57981.25
            #print 'hereeee'*100
            print self.recycler_mjd
            print self.mjd
            print now
#raw_input()
            #self.recycler_mjd = config["start_of_season"] + (now - season_start_date).days


        # Setup website directories
        self.mapspath = os.path.join(outfolder, "maps/")
        if not os.path.exists(self.mapspath):
            os.makedirs(self.mapspath)
        self.imagespath = "./DES_GW_Website/Triggers/" + trigger_id + "/"+outfolder.split('/')[-1]
        if not os.path.exists(self.imagespath):
            os.makedirs(self.imagespath)
        if not os.path.exists(self.imagespath+'/images'):
            os.makedirs(self.imagespath+'/images')

        #self.website_jsonpath = "./DES_GW_Website/Triggers/" + trigger_id + "/"
        #self.website_imagespath = "./DES_GW_Website/Triggers/" + trigger_id + "/images/"
        self.website_imagespath = self.imagespath+'/'+self.skymap_newoutfolder
        self.website_jsonpath = self.imagespath
        if not os.path.exists(self.website_imagespath):
            os.makedirs(self.website_imagespath)

        self.event_paramfile = os.path.join(outfolder, trigger_id + '_params.npz')
        print self.event_paramfile
        #raw_input()
        self.weHaveParamFile = True
        try:
            self.event_params = np.load(self.event_paramfile)
        except:
            self.event_params = {}
            self.weHaveParamFile = False
        os.system('cp recycler.yaml ' + os.path.join(outfolder, 'strategy.yaml'))
        print '***** Copied recycler.yaml to ' + os.path.join(outfolder,
                                                              'strategy.yaml') + ' for future reference *****'
        '''
        krbdir = '/usr/krb5/bin'
        ticket_cache = '/var/keytab/desgw.keytab'
        pid = os.getpid()
        krb_cache = '/tmp/krb5cc_desgw_%s' % pid
        os.environ['KRB5CCNAME']='FILE:%s' % krb_cache
        principal = 'desgw/des/des41.fnal.gov@FNAL.GOV'
        kinit_cmd = '%s/kinit -A -c %s -k -t %s %s' % (krbdir,krb_cache,ticket_cache,principal)
        os.system(kinit_cmd)
        '''
        os.system('kinit -k -t /var/keytab/desgw.keytab desgw/des/des41.fnal.gov@FNAL.GOV')


# Let's guess that mapMaker is the counterpart to recyc.mainInjector from
# desgw-maps. 

    def mapMaker(self, trigger_id, skymap, config):
        import os
        import yaml
        import getHexObservations


        overhead = config["overhead"]
        #nvisits = config["nvisits"]
        area_per_hex = config["area_per_hex"]
        start_of_season = config["start_of_season"]
        end_of_season = config["end_of_season"]
        events_observed = config["events_observed"]
        skipAll = config["skipAll"]
        mjd = self.mjd
        #print mjd, self.event_params['MJD']
        # if self.event_params['MJD'] == 'NAN':
        #     self.event_params['MJD'] = str(mjd)
        #raw_input()
        outputDir = self.outfolder
        mapDir = self.mapspath
        recycler_mjd = self.recycler_mjd

        start_days_since_burst = self.recycler_mjd - self.mjd
        #start_days_since_burst = 1.


        #recycler_mjd = 57773

        if self.skymap is None:
            self.skymap = os.path.join(outputDir,'lalinference.fits.gz')

        # If distance is not set in config use xml distance
        # if config["force_distance"]:
        #     distance = config["distance"]
        # else:
        #     if self.weHaveParamFile:
        #         distance = self.event_params["MaxDistance"]
        #     else:
        #         print 'THERE IS NO PARAMFILE, HARDCODING THE DISTANCE TO THE CONFIG DIST.'
        #         distance = config["distance"]

        
        eventtype = self.event_params['boc']


        try:
            probhasns = self.event_params['probhasns']
        except:
            probhasns = 0. #for old maps...

        if config['forceProbHasNS']: probhasns = config['probHasNS']

        self.probhasns = probhasns
        gethexobstype = None
        
        #print 'eventtype',eventtype
        if eventtype == 'Burst':
            gethexobstype = 'BH'
            self.distance = 1.
        elif eventtype == 'CBC':
            #print 'probhasns'*100
            print 'PROB HAS NS',probhasns
            if probhasns > config['probHasNS_threshold']:
                gethexobstype = 'NS'
                self.distance = -999
            else:
                gethexobstype = 'BH'
                self.distance = 1.
        else: #we dont know what we're looking at... do default obs for lightcurve
            print 'WE DONT KNOW WHAT WERE LOOKING AT!'*5
            gethexobstype = 'BH'
            self.distance = 1.



        if gethexobstype == 'BH':
            filter_list = config["exposure_filter_BH"]
            maxHexesPerSlot = config["maxHexesPerSlot_BH"]
            exposure_length = config["exposure_length_BH"]
            hoursAvailable = config["time_budget_for_BH"]
            tiling_list = config["exposure_tiling_BH"]
        else:
            filter_list = config["exposure_filter_NS"]
            maxHexesPerSlot = config["maxHexesPerSlot_NS"]
            exposure_length = config["exposure_length_NS"]
            hoursAvailable = config["time_budget_for_NS"]
            tiling_list = config["exposure_tiling_NS"]

        exposure_length = np.array(exposure_length)
        self.exposure_length = exposure_length
        self.time_budget = hoursAvailable

        if config["force_distance"]:
            self.distance = config["distance"]

        #self.distance = distance

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)


        self.gethexobstype = gethexobstype
        #print 'triggertype'*100
        print 'TRIGGER TYPE:',self.gethexobstype
        # make the maps
        #try:
        #where = 'getHexObservations'
        #line = '103'
            #try:
            #    probs, times, slotDuration, hoursPerNight = getHexObservations.prepare(
            #        skymap, mjd, trigger_id, outputDir, mapDir, distance=distance,
            #        exposure_list=exposure_length, filter_list=filter_list,
            #        overhead=overhead, maxHexesPerSlot=maxHexesPerSlot, skipAll=skipAll)
            #except ValueError:

        # print 'skymap',self.skymap
        #
        #
        # print 'skymap',self.skymap
        # print 'distance',self.distance
        # print 'gethexobstype',gethexobstype
        # print 'start_days_since_burst',start_days_since_burst
        # print 'exposure_length',exposure_length
        # print 'filter_list',filter_list
        # print 'resolution',config['resolution']
        # print 'halfnight',config['ishalfnight']
        # print 'firsthalf', config['isfirsthalf']
        # print 'overhead',overhead
        # print 'maxHexesPerSlot',maxHexesPerSlot
        # print 'skipAll',skipAll

        #raw_input()
        probs, times, slotDuration, hoursPerNight = getHexObservations.prepare(
                    self.skymap, trigger_id, outputDir, mapDir, distance=self.distance,
                    trigger_type=gethexobstype,start_days_since_burst=start_days_since_burst,
                    exposure_list=exposure_length, filter_list=filter_list,resolution=config['resolution'],
                    halfNight=config['ishalfnight'], firstHalf=config['isfirsthalf'],
                    #isCustomDark=config['isCustomDark'],customDarkIndices=config['customDarkSlots'],
                    overhead=overhead, maxHexesPerSlot=maxHexesPerSlot, skipAll=skipAll)
            # figure out how to divide the night
            # where = 'getHexObservations.contemplateTheDivisionsOfTime()'
            # line = '102'

        # print 'probs',probs
        # print 'times', times
        # print 'slotDuration', slotDuration
        # print 'hoursPerNight', hoursPerNight
        #raw_input()
        print probs,times,slotDuration, hoursPerNight
        #raw_input('probs times')
        n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
                probs, times, hoursPerNight=hoursPerNight,
                hoursAvailable=hoursAvailable)
        print n_slots, first_slot
        #raw_input('contemplate')
            # compute the best observations
            # where = 'getHexObservations.now()'
            # line = '109'
        best_slot = getHexObservations.now(
                n_slots, mapDirectory=mapDir, simNumber=trigger_id,
                maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot,
                exposure_list=exposure_length, filter_list=filter_list,
                #tiling_list = tiling_list,
                trigger_type=gethexobstype, skipJson=config['skipjson'])
        # except:
        #     try:
        #         print 'skymap', self.skymap
        #         self.skymap = os.path.join(outputDir,'bayestar.fits.gz')
        #
        #         probs, times, slotDuration, hoursPerNight = getHexObservations.prepare(
        #             self.skymap, mjd, trigger_id, outputDir, mapDir, distance=distance,
        #             exposure_list=exposure_length, filter_list=filter_list,
        #             overhead=overhead, maxHexesPerSlot=maxHexesPerSlot, skipAll=skipAll)
        #         # figure out how to divide the night
        #         # where = 'getHexObservations.contemplateTheDivisionsOfTime()'
        #         # line = '102'
        #         n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
        #             probs, times, hoursPerNight=hoursPerNight,
        #             hoursAvailable=hoursAvailable)
        #
        #         # compute the best observations
        #         # where = 'getHexObservations.now()'
        #         # line = '109'
        #         best_slot = getHexObservations.now(
        #             n_slots, mapDirectory=mapDir, simNumber=trigger_id,
        #             maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot,
        #             exposure_list=exposure_length, filter_list=filter_list,
        #             skipJson=True)
        #     except:
        #         try:
        #             print 'skymap', self.skymap
        #             self.skymap = os.path.join(outputDir, 'lalinference.fits.gz')
        #
        #             probs, times, slotDuration, hoursPerNight = getHexObservations.prepare(
        #                 self.skymap, mjd, trigger_id, outputDir, mapDir, distance=distance,
        #                 exposure_list=exposure_length, filter_list=filter_list,
        #                 overhead=overhead, maxHexesPerSlot=maxHexesPerSlot, skipAll=skipAll)
        #             # figure out how to divide the night
        #             # where = 'getHexObservations.contemplateTheDivisionsOfTime()'
        #             # line = '102'
        #             n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
        #                 probs, times, hoursPerNight=hoursPerNight,
        #                 hoursAvailable=hoursAvailable)
        #
        #             # compute the best observations
        #             # where = 'getHexObservations.now()'
        #             # line = '109'
        #             best_slot = getHexObservations.now(
        #                 n_slots, mapDirectory=mapDir, simNumber=trigger_id,
        #                 maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot,
        #                 exposure_list=exposure_length, filter_list=filter_list,
        #                 skipJson=True)
        #         except:
        #             e = sys.exc_info()
        #             trace = traceback.format_exc(sys.exc_info())
        #             print trace
        #             self.send_processing_error(e, where, line, trace)
        #             sys.exit()

        skipecon = True


        if not skipecon:
            if n_slots > 0:
                print "================ N_SLOTS > 0 =================== "
                #   area_left is th enumber of hexes we have left to observe this season
                #   T_left is the number of days left in the season
                #   rate is the effective rate of triggers
                #
                # in seconds
                time_cost_per_hex = nvisits * np.sum(overhead + exposure_length)
                area_left = area_per_hex * \
                            (hoursAvailable * 3600) / (time_cost_per_hex)
                time_left = end_of_season - start_of_season
                rate = len(events_observed) / (recycler_mjd - start_of_season)

                # do Hsun-yu Chen's
                try:
                    where = 'getHexObservations.economics()'
                    line = '136'
                    econ_prob, econ_area, need_area, quality = \
                        getHexObservations.economics(trigger_id,
                                                     best_slot, mapDirectory=mapDir,
                                                     area_left=area_left, days_left=time_left,
                                                     rate=rate)

                    hoursOnTarget = (econ_area / area_per_hex) * (time_cost_per_hex / 3600.)

                    # figure out how to divide the night,
                    # given the new advice on how much time to spend
                    where = 'getHexObservations.contemplateTheDivisionsOfTime()'
                    line = '148'
                    n_slots, first_slot = \
                        getHexObservations.contemplateTheDivisionsOfTime(
                            probs, times, hoursPerNight=hoursPerNight,
                            hoursAvailable=hoursOnTarget)

                    where = 'getHexObservations.now()'
                    line = '156'
                    best_slot = getHexObservations.now(
                        n_slots, mapDirectory=mapDir, simNumber=trigger_id,
                        maxHexesPerSlot=maxHexesPerSlot, mapZero=first_slot,
                        exposure_list=exposure_length, filter_list=filter_list,
                        skipJson=False)
                except:
                    e = sys.exc_info()
                    trace = traceback.format_exc(sys.exc_info())
                    print trace
                    self.send_processing_error(e, where, line, trace)
                    sys.exit()
        else:
            econ_prob = 0
            econ_area = 0
            #best_slot = 0
            #print 'setting best slot to zero'*10
            #raw_input()
            need_area = 11734.0
            quality = 1.0

        # make observation plots
        # try:
        #     where = 'getHexObservations.makeObservingPlots()'
        #     line = '176'
        #     print '888' * 20
        #     print n_slots, trigger_id, best_slot, outputDir, mapDir
        #     print '888' * 20
        #     if not config['skipPlots']:
        #         n_plots = getHexObservations.makeObservingPlots(
        #             n_slots, trigger_id, best_slot, outputDir, mapDir, allSky=True )
        #     #string = "$(ls -v {}-observingPlot*)"
        # except:
        #     e = sys.exc_info()
        #     trace = traceback.format_exc(sys.exc_info())
        #     print trace
        #     self.send_processing_error(e, where, line, trace)
        #     sys.exit()


        try:
            self.sumligoprob = getHexObservations.how_well_did_we_do(
                self.skymap, trigger_id, mapDir)
        except:
            e = sys.exc_info()
            exc_type, exc_obj, exc_tb = e[0],e[1],e[2]
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            where = fname
            line = exc_tb.tb_lineno
            trace = traceback.format_exc(sys.exc_info())
            print trace
            self.send_processing_error(e, where, line, trace)
            sys.exit()

        self.best_slot = best_slot
        self.n_slots = n_slots
        self.first_slot = first_slot
        self.econ_prob = econ_prob
        self.econ_area = econ_area
        self.need_area = need_area
        self.quality = quality


        np.savez(os.path.join(self.outfolder, 'mapmaker_results.npz')
                 , best_slot=best_slot
                 , n_slots=n_slots
                 , first_slot=first_slot
                 , econ_prob=econ_prob
                 , econ_area=econ_area
                 , need_area=need_area
                 , quality=quality
                 )

        ra, dec, id, self.prob, mjd, slotNum = \
            obsSlots.readObservingRecord(self.trigger_id, mapDir)

        integrated_prob = np.sum(self.prob)
        print '-'*20+'>','LIGO PROB: %.3f \tLIGO X DES PROB: %.3f' % (self.sumligoprob,integrated_prob)
        #raw_input('checking comparison of probs!!!!'*10)


        self.outputDir = outputDir
        self.mapDir = mapDir

        self.weHaveParamFile = True

        if self.weHaveParamFile:
            np.savez(self.event_paramfile,
                     MJD=self.mjd,
                     ETA=self.event_params['ETA'],
                     FAR=self.event_params['FAR'],
                     ChirpMass=self.event_params['ChirpMass'],
                     MaxDistance=self.event_params['MaxDistance'],
                     DESXLIGO_prob=integrated_prob,
                     LIGO_prob=self.sumligoprob,
                     M1=self.event_params['M1'],
                     M2=self.event_params['M2'],
                     nHexes=self.prob.size,
                     time_processed=self.now.strftime("%H:%M %B %d, %Y "),
                     boc=self.event_params['boc'],
                     CentralFreq=self.event_params['CentralFreq'],
                     best_slot=self.best_slot,
                     n_slots=self.n_slots,
                     first_slot=self.first_slot,
                     econ_prob=self.econ_prob,
                     econ_area=self.econ_area,
                     need_area=self.need_area,
                     quality=self.quality,
                     codeDistance=self.distance,
                     exposure_times=exposure_length,
                     exposure_filter=filter_list,
                     hours=self.time_budget,
                     nvisits=-999,#config['nvisits'],
                     mapname='NAN',
                     filename=self.skymap,
                     gethexobstype=self.gethexobstype,
                     probhasns=self.probhasns
                     )

        else:
            np.savez(self.event_paramfile,
                     MJD='NAN',
                     ETA='NAN',
                     FAR='NAN',
                     ChirpMass='NAN',
                     MaxDistance='NAN',
                     DESXLIGO_prob=integrated_prob,
                     LIGO_prob=self.sumligoprob,
                     M1='NAN',
                     M2='NAN',
                     nHexes=self.prob.size,
                     time_processed=self.now.strftime("%H:%M %B %d, %Y "),
                     boc='NAN',
                     CentralFreq='NAN',
                     best_slot=self.best_slot,
                     n_slots=self.n_slots,
                     first_slot=self.first_slot,
                     econ_prob=self.econ_prob,
                     econ_area=self.econ_area,
                     need_area=self.need_area,
                     quality=self.quality,
                     codeDistance=self.distance,
                     exposure_times=exposure_length,
                     exposure_filter=filter_list,
                     hours=self.time_budget,
                     nvisits=-999,#config['nvisits'],
                     mapname='NAN',
                     filename=self.skymap,
                     gethexobstype=self.gethexobstype,
                     probhasns = 'NAN'
                     )

    def getContours(self, config):
        import matplotlib.pyplot as plt

        #if exposure_length is None:
        #    exposure_length = config["exposure_length"]
        exposure_length= self.exposure_length
        image_dir = self.website_imagespath
        map_dir = self.mapspath
        
        if self.n_slots<1:
            counter = getHexObservations.nothingToObserveShowSomething(self.trigger_id, self.outfolder, self.mapspath)
            # iname = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
            # oname = self.trigger_id + "-observingPlot.gif"
            # os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
            iname = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
            oname = self.trigger_id + "-probabilityPlot.png"
            os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
        #if self.n_slots > 0:
        if True:
            # print 'Converting Observing Plots to .gif'
            # os.system('convert $(for ((a=0; a<50; a++)); do printf -- "-delay 50 '+os.path.join(map_dir,self.trigger_id)+'-observingPlot-%s.png " $a; done;) '+os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif')
            # #os.system('convert -delay 70 -loop 0 '+os.path.join(map_dir,self.trigger_id)+'-observingPlot-*.png '+
            # #          os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif')
            # os.system('cp '+os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif '+ image_dir)
            iname = self.trigger_id + "-" + str(self.best_slot) + "-maglim-eq.png"
            oname = self.trigger_id + "_limitingMagMap.png"
            os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
            iname = self.trigger_id + "-" + str(self.best_slot) + "-prob-eq.png"
            oname = self.trigger_id + "_sourceProbMap.png"
            os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
            iname = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
            oname = self.trigger_id + "_LIGO.png"
            os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
            iname = self.trigger_id + "-" + str(self.best_slot) + "-probXligo-eq.png"
            oname = self.trigger_id + "_sourceProbxLIGO.png"
            os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
            # DESGW observation map
            inname = self.trigger_id + "-observingPlot-{}.png".format(
                self.best_slot)
            outname = self.trigger_id + "-observingPlot.png"
            os.system('cp ' + os.path.join(map_dir, inname) + ' ' + os.path.join(image_dir, outname))
            # probability plot
            name = self.trigger_id + "-probabilityPlot.png"
            os.system('cp ' + os.path.join(self.outfolder, name) + ' ' + image_dir)
            #raw_input('getting contours stopped')
        # oname = 'observingPlots.gif'
        # giffile = os.path.join(self.outfolder,iname)+' '+ os.path.join(image_dir,oname)
        # oname = self.trigger_id+'-observingPlot-*.png'
        # pngs = os.path.join(self.outfolder,iname)+' '+ os.path.join(image_dir,oname)
        # os.system('convert -delay 10 -loop 0 '+pngs+' '+giffile)

        # else:
        #     # there is nothing to observe, make default plots
        #     try:
        #         where = 'getHexObservations.nothingToObserveShowSomething()'
        #         line = '240'
        #
        #         ra = [-999]
        #         dec = [-999]
        #         ligo = -999
        #         maglim = [-999]
        #         probMap = [-999]
        #     except:
        #         e = sys.exc_info()
        #         trace = traceback.format_exc(sys.exc_info())
        #         print trace
        #         self.send_processing_error(e, where, line, trace)
        #         #sys.exit()
        #         ra = [-999]
        #         dec = [-999]
        #         ligo = -999
        #         maglim = [-999]
        #         probMap = [-999]
        #
        #
        #     print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
        #     print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
        #     print "faking it with getHexObservations.nothingToObserveShowSomething("
        #     print self.skymap, self.mjd, exposure_length
        #     print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
        #     print "================ >>>>>>>>>>>>>>>>>>>>> =================== "
        #
        #     figure = plt.figure(1,figsize=(8.5*1.618,8.5))
        #     plt.figure(figsize=(8.5*1.618,8.5))
        #     # computing limiting mag
        #     # plot as ra,dec map
        #     plt.clf()
        #     print 'ra', ra
        #     print 'dec', dec
        #     print 'maglim', maglim
        #
        #     plt.hist(maglim, bins=np.arange(-11, -7, .2))
        #     plt.xlabel('maglim')
        #     plt.ylabel('counts')
        #     name = self.trigger_id + "_limitingMagHist.png"
        #     plt.savefig(os.path.join(self.outfolder, name))
        #
        #     # plt.hexbin( ra, dec, maglim, vmin=15)
        #     plt.hexbin(ra, dec, maglim)
        #     plt.colorbar()
        #     plt.xlabel('RA')
        #     plt.ylabel('DEC')
        #     name = self.trigger_id + "_limitingMagMap.png"
        #     plt.savefig(os.path.join(self.outfolder, name))
        #     os.system('cp ' + os.path.join(self.outfolder, name) + ' ' + image_dir)
        #
        #     # Calculate source probability map
        #     plt.clf()
        #     plt.hexbin(ra, dec, probMap, )
        #     plt.colorbar()
        #     plt.xlabel('RA')
        #     plt.ylabel('DEC')
        #     name = self.trigger_id + "_sourceProbMap.png"
        #     plt.savefig(os.path.join(self.outfolder, name))
        #     os.system('cp ' + os.path.join(self.outfolder, name) + ' ' + image_dir)
        #
        #     # DES Source Prob Map x Ligo Sky Map
        #     plt.clf()
        #     plt.hexbin(ra, dec, probMap * ligo)
        #     plt.colorbar()
        #     plt.xlabel('RA')
        #     plt.ylabel('DEC')
        #     name = self.trigger_id + "_sourceProbxLIGO.png"
        #     plt.savefig(os.path.join(self.outfolder, name))
        #     os.system('cp ' + os.path.join(self.outfolder, name) + ' ' + image_dir)
        #
        #     plt.clf()
        #     plt.hexbin(ra, dec, ligo)
        #     plt.colorbar()
        #     plt.xlabel('RA')
        #     plt.ylabel('DEC')
        #     name = self.trigger_id + "_LIGO.png"
        #     plt.savefig(os.path.join(self.outfolder, name))
        #     os.system('cp ' + os.path.join(self.outfolder, name) + ' ' + image_dir)

        return

    def makeJSON(self, config):

        mapmakerresults = np.load(os.path.join(self.outfolder, 'mapmaker_results.npz'))

        self.best_slot = mapmakerresults['best_slot']
        self.n_slots = mapmakerresults['n_slots']
        self.first_slot = mapmakerresults['first_slot']
        self.econ_prob = mapmakerresults['econ_prob']
        self.econ_area = mapmakerresults['econ_area']
        self.need_area = mapmakerresults['need_area']
        self.quality = mapmakerresults['quality']

        # DESGW json file (to be files once that is done)
        json_dir = self.website_jsonpath
        map_dir = self.mapspath
        jsonname = self.trigger_id + "_"+ self.skymap_newoutfolder +"_JSON.zip"
        jsonFile = os.path.join(map_dir, jsonname)
        jsonfilelistld = os.listdir(map_dir)
        jsonfilelist = []
        for f in jsonfilelistld:
            if '-tmp' in f:
                os.remove(os.path.join(map_dir, f))
            elif '.json' in f:
                jsonfilelist.append(f)

        if self.n_slots > 0:
            # get statistics
            ra, dec, id, self.prob, mjd, slotNum = \
                obsSlots.readObservingRecord(self.trigger_id, map_dir)

            # adding integrated probability to paramfile
            integrated_prob = np.sum(self.prob)
            nHexes = str(self.prob.size)
        else:
            integrated_prob = 0
            nHexes = str(0)

        from time import gmtime, strftime
        timeprocessed = strftime("%H:%M:%S GMT \t %b %d, %Y", gmtime())

        #exptimes = ', '.join(map(str, config['exposure_length']))
        #expf = ', '.join(map(str, config['exposure_filter']))

        try:
            boc = self.event_params['boc']
        except:
            boc = 'NA'

        # Copy json file to web server for public download
        if not os.path.exists(jsonFile):
            if integrated_prob == 0:
                print "zero probability, thus no jsonFile at ", jsonFile
            else:
                # try:
                os.chmod(self.mapspath, 0o777)
                for js in os.listdir(self.mapspath):
                    os.chmod(os.path.join(self.mapspath,js), 0o777)

                os.system('zip -j ' + jsonFile + ' ' + self.mapspath + '/*0.json')
                # except:
                #    print "no jsonFiles at ", jsonFile
        else:
            os.remove(jsonFile)
            os.system('zip -j ' + jsonFile + ' ' + self.mapspath + '/*0.json')
            os.system('cp ' + jsonFile + ' ' + self.website_jsonpath)
        return jsonfilelist

    def send_nonurgent_Email(self,sendtexts=False):
        import smtplib
        from email.mime.text import MIMEText

        text = 'DESGW Webpage Created. See \nhttp://des-ops.fnal.gov:8080/desgw/Triggers/' + self.trigger_id + '/' + self.trigger_id + '_' +self.skymap_newoutfolder+ '_trigger.html'

        # Create a text/plain message
        msg = MIMEText(text)

        # me == the sender's email address
        # you == the recipient's email address
        me = 'automated-desGW@fnal.gov'
        if self.config['sendEmailsToEveryone']:
            yous = ['djbrout@gmail.com', 'marcelle@fnal.gov', 'annis@fnal.gov']
        else:
            yous = ['djbrout@gmail.com']

        if sendtexts:
            t = ['7737578495@msg.fi.google.com', '3017883369@mms.att.net', '6173357963@mms.att.net',
                 '2153008763@mms.att.net', '6307654596@tmomail.net']
            yous.extend(t)


        msg['Subject'] = 'DESGW Webpage Created for ' + self.trigger_id + ' Map: '+self.skymap_newoutfolder
        msg['From'] = me
        for you in yous:
            msg['To'] = you

            s = smtplib.SMTP('localhost')
            s.sendmail(me, [you], msg.as_string())
            s.quit()
        print 'Email sent...'
        return

    def send_processing_error(self, error, where, line, trace):
        import smtplib
        from email.mime.text import MIMEText

        message = 'Processing Failed for Trigger ' + str(self.trigger_id) + '\n\nFunction: ' + str(
            where) + '\n\nLine ' + str(line) + ' of recycler.py\n\nError: ' + str(error) + '\n\n'
        message += '-' * 60
        message += '\n'
        message += trace
        message += '\n'
        message += '-' * 60
        message += '\n'

        # Create a text/plain message                                                                                                                                                                      
        msg = MIMEText(message)

        me = 'automated-desGW@fnal.gov'
        if self.config['sendEmailsToEveryone']:
            yous = ['djbrout@gmail.com', 'marcelle@fnal.gov', 'annis@fnal.gov']
        else:
            yous = ['djbrout@gmail.com']
        msg['Subject'] = 'Trigger ' + self.trigger_id + ' '+self.skymap_newoutfolder+ ' Processing FAILED!'
        msg['From'] = me
        for you in yous:
            msg['To'] = you
            s = smtplib.SMTP('localhost')
            s.sendmail(me, [you], msg.as_string())
            s.quit()
        print 'Email sent...'
        return

    def updateTriggerIndex(self, real_or_sim=None):
        if real_or_sim == 'real':
            fff = './DES_GW_Website/real-trigger_list.txt'
        if real_or_sim == 'sim':
            fff = './DES_GW_Website/test-trigger_list.txt'
        l = open(fff, 'r')
        lines = l.readlines()
        l.close()
        a = open(fff, 'a')
        triggers = []
        for line in lines:
            triggers.append(line.split(' ')[0])
        if not self.trigger_id in np.unique(triggers):
            a.write(self.trigger_id + ' ' + self.outfolder + '\n')
        a.close()
        tp.make_index_page('./DES_GW_Website', real_or_sim=real_or_sim)
        return

    def make_cumulative_probs(self):
        print ['python', './sims_study/cumulative_plots.py', '-d',
               '/data/des41.a/data/desgw/maininjector/sims_study/data', '-p', self.outfolder, '-e', self.trigger_id,
               '-f', os.path.join(self.outfolder, 'maps', self.trigger_id + '-ra-dec-id-prob-mjd-slot.txt')]
        subprocess.call(['python', './sims_study/cumulative_plots.py', '-d',
                         '/data/des41.a/data/desgw/maininjector/sims_study/data', '-p', self.outfolder, '-e',
                         self.trigger_id, '-f',
                         os.path.join(self.outfolder, 'maps', self.trigger_id + '-ra-dec-id-prob-mjd-slot.txt')])
        os.system('scp ' + os.path.join(self.outfolder,
                                        self.trigger_id + '-and-sim-cumprobs.png') + ' ./DES_GW_Website/Triggers/' + self.trigger_id + '/images/')

    def updateWebpage(self,real_or_sim):
        os.system('scp -r DES_GW_Website/Triggers/'+self.trigger_id+' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/')
        os.system('scp DES_GW_Website/* codemanager@desweb.fnal.gov:/des_web/www/html/desgw/')
        tp.makeNewPage(os.path.join(self.oldoutfolder, self.trigger_id +'_'+ self.skymap_newoutfolder+ '_trigger.html'), self.trigger_id,self.event_paramfile,self.skymap_newoutfolder,real_or_sim=real_or_sim)
        os.system('scp -r ' + os.path.join(self.oldoutfolder,
                                           self.trigger_id +'_' + self.skymap_newoutfolder + '_trigger.html') + ' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + self.trigger_id + '/')
        os.system('scp -r ' + os.path.join(self.oldoutfolder,
                                           self.trigger_id +'_'+ self.skymap_newoutfolder + '_trigger.html') + ' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + self.trigger_id + '/'+self.trigger_id +'_trigger.html')

        os.system('scp -r ' + self.outfolder + ' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + self.trigger_id + '/')
        os.system('cp ' + self.outfolder+'/'+self.skymap_newoutfolder + 'recycler.log ' + self.website_jsonpath)
        #os.system('scp ' + self.oldoutfolder + '/*/*.html ' + ' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + self.trigger_id + '/')
        return

    def getmjd(self,datet):
        mjd_epoch = datetime.datetime(1858, 11, 17)
        print 'FIX ME UTC OR NOT?!?'
        mjdd = datet-mjd_epoch
        mjd = 5./24. + mjdd.total_seconds() / 3600. / 24.
        return mjd

    def mjd_to_datetime(self,mjd):
        mjd_epoch = datetime.datetime(1858, 11, 17, tzinfo=pytz.utc)
        d = mjd_epoch + datetime.timedelta(mjd)
        return d

    def makeObservingPlots(self):
        try:
            where = 'getHexObservations.makeObservingPlots()'
            line = '776'
            if not self.config['skipPlots']:
                n_plots = getHexObservations.makeObservingPlots(
                    self.n_slots, self.trigger_id, self.best_slot, self.outputDir, self.mapDir, allSky=True )

                image_dir = self.website_imagespath
                map_dir = self.mapspath

                if self.n_slots < 1:
                    counter = getHexObservations.nothingToObserveShowSomething(self.trigger_id, self.outfolder,
                                                                               self.mapspath)
                    iname = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
                    oname = self.trigger_id + "-observingPlot.gif"
                    os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
                    iname = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
                    oname = self.trigger_id + "-probabilityPlot.png"
                    os.system('cp ' + os.path.join(self.outfolder, iname) + ' ' + os.path.join(image_dir, oname))
                # if self.n_slots > 0:
                if True:
                    print 'Converting Observing Plots to .gif'
                    os.system('convert $(for ((a=0; a<250; a++)); do printf -- "-delay 50 ' + os.path.join(map_dir,
                                        self.trigger_id) + '-observingPlot-%s.png " $a; done;) ' + os.path.join(
                        map_dir, self.trigger_id) + '-observingPlot.gif')
                    # os.system('convert -delay 70 -loop 0 '+os.path.join(map_dir,self.trigger_id)+'-observingPlot-*.png '+
                    #          os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif')
                    os.system('cp ' + os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif ' + image_dir)
            #string = "$(ls -v {}-observingPlot*)"
        except:
            e = sys.exc_info()
            trace = traceback.format_exc(sys.exc_info())
            print trace
            self.send_processing_error(e, where, line, trace)
            sys.exit()

if __name__ == "__main__":

    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "tp:tid:mjd:exp:sky",
            longopts=["triggerpath=", "triggerid=", "mjd=", "exposure_length=", "skymapfilename="])

    except getopt.GetoptError as err:
        print str(err)
        print "Error : incorrect option or missing argument."
        print __doc__
        sys.exit(1)

    # Read in config
    with open("recycler.yaml", "r") as f:
        config = yaml.safe_load(f);
    # Set defaults to config
    trigger_path = config["trigger_path"]

    real_or_sim = config["real_or_sim"]

    if config["skymap_filename"] == 'Default':
        skymap_filename = None
    else:
        skymap_filename = config["skymap_filename"]

    trigger_ids = [config["trigger_id"]]

    force_mjd = config["force_mjd"]

    #exposure_length = config["exposure_length"]

    # Override defaults with command line arguments
    # THESE NOT GUARANTEED TO WORK EVER SINCE WE SWITCHED TO YAML

    dontwrap = False
    for o,a in opt:
        print 'Option'
        print o
        print a
        print '-----'
        if o in ["-tp","--triggerpath"]:
            trigger_path = str(a)
        elif o in ["-tid","--triggerid"]:
            trigger_ids = [str(a)]
            dontwrap = True
        elif o in ["-mjd","--mjd"]:
            mjd = float(a)
        #elif o in ["-exp","--exposure_length"]:
        #    exposure_length = float(a)
        elif o in ["-hours","--hours_available"]:
            hours_available = float(a)
        elif o in ["-sky","--skymapfilename"]:
            skymap_filename = str(a)
        else:
            print "Warning: option", o, "with argument", a, "is not recognized"

    # Clear bad triggers, only used for wrapping all triggers...
    badtriggers = open('badtriggers.txt', 'w')
    badtriggers.close()



    ####### BIG MONEY NO WHAMMIES ###############################################
    if config["wrap_all_triggers"]:
        if not dontwrap:
            trigger_ids = os.listdir(trigger_path)
            trigger_ids = trigger_ids[2:]
    for trigger_id in trigger_ids:
        if force_mjd:
            mjd = config["mjd"]
        else:
            try:
                mjd = open(os.path.join(trigger_path, trigger_id, trigger_id + '_eventMJD.txt'), 'r').read()
            except:
                mjd = '99999'
        if skymap_filename is None:
            try:
            #if True:
                #mapname = open(os.path.join(trigger_path,
                #                            trigger_id,
                #                            config['default_map_name']), 'r').read()
                #skymap_filename = os.path.join(trigger_path,
                #                               trigger_id, config['default_map_name'])
                #print os.path.join(trigger_path, trigger_id,'default_skymap.txt')
                #print os.path.join(trigger_path, trigger_id,'default_skymap.txt').read()
                skymap_filename = os.path.join(trigger_path, trigger_id,
                                               open(os.path.join(trigger_path, trigger_id,
                                                                 'default_skymap.txt'),'r').read())
            except:
               badtriggers = open('badtriggers.txt', 'a')
               badtriggers.write(trigger_id + '\n')
               print 'Could not find skymap url file'

        if 'bayestar' in skymap_filename:
            print 'bayestar' * 500
            #print 'waiting 20 seconds for lalinference map otherwise compute using bayestar...'
            #time.sleep(20)
            #if os.path.exists(os.path.join(trigger_path,trigger_id) + '/wehavelal'):
            #    print 'bayestar skipped because we have lalinference map'
            #    sys.exit()
            #else:
            #    print 'running bayestar, never recieved lalinference map'
            #    pass

        try:
            try:
                mjd = float(mjd)
            except:
                badtriggers = open('badtriggers.txt', 'a')
                badtriggers.write(trigger_id + '\n')
                print 'WARNING: Could not convert mjd to float. Trigger: ' + trigger_id + ' flagged as bad.'
# here is where the object is made, and parts of it are filed in
            e = event(skymap_filename,
                      os.path.join(trigger_path,
                                   trigger_id),
                      trigger_id, mjd, config)

# e has variables and code assocaiated with it. The mapMaker is called "e" or "self"
            e.mapMaker(trigger_id, skymap_filename, config)
            e.getContours(config)
            jsonfilelist = e.makeJSON(config)
            e.make_cumulative_probs()
            e.updateTriggerIndex(real_or_sim=real_or_sim)
            e.updateWebpage(real_or_sim)
            e.send_nonurgent_Email()
            e.makeObservingPlots()
            e.getContours(config)
            e.updateWebpage(real_or_sim)

            # ISREALTRIGGER = True
            # eventmngr = Thread(target=jobmanager.eventmanager, args=(trigger_id, jsonfilelist,os.path.join(trigger_path,trigger_id),
            #                                                 os.path.join(trigger_path, trigger_id, 'maps'),ISREALTRIGGER,trigger_path))
            # eventmngr.start()

            #e.send_nonurgent_Email()
            #eventmngr.join()

        except KeyError:
            print "Unexpected error:", sys.exc_info()
            badtriggers = open('badtriggers.txt', 'a')
            badtriggers.write(trigger_id + '\n')
    #############################################################################

    print 'Done'
