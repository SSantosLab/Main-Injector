import os
import glob
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
    def __init__(self, skymap_filename, master_dir, trigger_id, mjd, config):

        '''
        event:

        master_dir:  the directory containing all trigger subdirectories
        trigger_id: the name of the trigger event that comes from LIGO
        mjd:        the modified julian date of the event (in FITS header, among other places)
        config:     the config filename with the parameters that controls the routine
        '''
        
        # set up the working directories
        self.modify_filesystem(skymap_filename, master_dir, trigger_id, mjd) 
        work_area = self.work_area 

        # read config file
        if config["force_recycler_mjd"]:
            self.recycler_mjd = config["recycler_mjd"]
        else:
            self.recycler_mjd = self.getmjd(datetime.datetime.now())

        self.event_paramfile = os.path.join(work_area, trigger_id + '_params.npz')
        print(self.event_paramfile)
        #raw_input()
        self.weHaveParamFile = True
        try:
            self.event_params = np.load(self.event_paramfile)
        except:
            self.event_params = {}
            self.weHaveParamFile = False
        yaml_dir = os.path.join(work_area, 'strategy.yaml')
        os.system('cp recycler.yaml ' + yaml_dir)
        print('***** Copied recycler.yaml to ' + yaml_dir + ' for future reference *****')
        os.system('kinit -k -t /var/keytab/desgw.keytab desgw/des/des41.fnal.gov@FNAL.GOV')


    def modify_filesystem(self, skymap_filename, master_dir, trigger_id, mjd) :

        # master_dir is the directory holding all the trigger directories
        # trigger_dir is the directory holding everything to do with a single trigger
        # work_area is master_dir + trigger_dir
            # work_area is derived
            # trigger_dir is derived from the path of the skymap_file

        skymap_filename = skymap_filename.strip()
        trigger_dir = skymap_filename.split('/')[-1].split('.')[0]
        self.trigger_dir = trigger_dir
        os.system('touch '+skymap_filename+'.processing')
        os.system('mkdir '+master_dir+'/'+ trigger_dir)

        self.master_dir = master_dir
        work_area = master_dir+'/'+ trigger_dir +'/'

        os.system('cp '+skymap_filename.strip()+' '+work_area)
        os.system('cp '+self.master_dir+'/'+trigger_id +'_params.npz '+work_area)
        print('cp '+skymap_filename.strip()+' '+work_area)
        print(work_area)
        self.skymap = skymap_filename
        print(self.skymap)

        # Setup website directories
        website = "./DES_GW_Website/"
        if not os.path.exists(website) :
            os.mkdir(website)

        self.mapspath = os.path.join(work_area, "maps/")
        if not os.path.exists(self.mapspath):
            os.makedirs(self.mapspath)
        self.imagespath = website + "Triggers/" + trigger_id + "/"+work_area.split('/')[-1]
        if not os.path.exists(self.imagespath):
            os.makedirs(self.imagespath)
        if not os.path.exists(self.imagespath+'/images'):
            os.makedirs(self.imagespath+'/images')

        self.website_imagespath = self.imagespath+'/'+self.trigger_dir
        self.website_jsonpath = self.imagespath
        if not os.path.exists(self.website_imagespath):
            os.makedirs(self.website_imagespath)

        self.master_dir = master_dir
        self.work_area = work_area
        self.trigger_id = trigger_id
        self.mjd = mjd
        self.config = config
        self.website=website

# Let's guess that mapMaker is the counterpart to recyc.mainInjector from
# desgw-maps. 

    def mapMaker(self, trigger_id, skymap, config):
        import os
        import yaml
        import getHexObservations

        # debug
        debug = config["debug"]

        # camera
        camera   = config["camera"]

       #resolution
        resolution = config["resolution"]
        
        overhead = config["overhead"]
        #nvisits = config["nvisits"]
        area_per_hex = config["area_per_hex"]
        start_of_season = config["start_of_season"]
        end_of_season = config["end_of_season"]
        events_observed = config["events_observed"]
        skipAll = config["skipAll"]
        mjd = self.mjd
        outputDir = self.work_area
        mapDir = self.mapspath
        recycler_mjd = self.recycler_mjd

        start_days_since_burst = self.recycler_mjd - self.mjd

        #Set parameters if the camera is HSC
        if camera == "hsc":
            overhead = 20.
            area_per_hex = 1.5

        if self.skymap is None:
            self.skymap = os.path.join(outputDir,'lalinference.fits.gz')

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
            print('PROB HAS NS',probhasns)
            if probhasns > config['probHasNS_threshold']:
                gethexobstype = 'NS'
                self.distance = -999
            else:
                gethexobstype = 'BH'
                self.distance = 1.
        else: #we dont know what we're looking at... do default obs for lightcurve
            print('WE DONT KNOW WHAT WERE LOOKING AT!'*5)
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
        self.camera = camera

        if config["force_distance"]:
            self.distance = config["distance"]

        #self.distance = distance

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)


        self.gethexobstype = gethexobstype
        #print 'triggertype'*100
        print('TRIGGER TYPE:',self.gethexobstype)

        #raw_input()
        probs, times, slotDuration, hoursPerNight = getHexObservations.prepare(
                    self.skymap, trigger_id, outputDir, mapDir, distance=self.distance,
                    trigger_type=gethexobstype,start_days_since_burst=start_days_since_burst,
                    exposure_list=exposure_length, filter_list=filter_list,resolution=resolution, debug=debug, camera=camera,
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
        print(probs,times,slotDuration, hoursPerNight)
        #raw_input('probs times')
        n_slots, first_slot = getHexObservations.contemplateTheDivisionsOfTime(
                probs, times, hoursPerNight=hoursPerNight,
                hoursAvailable=hoursAvailable)
        print(n_slots, first_slot)
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
        skipecon = True


        if not skipecon:
            if n_slots > 0:
                print("================ N_SLOTS > 0 =================== ")
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
                    # do Hsun-yu Chen's 
                    print("======================================>>>>>>>>>>>>>>>>>>")
                    print(" economics ")
                    print("getHexObservations.economics (", trigger_id, ",",\
                        best_slot, ", mapDirectory= \"",outputDir, "\" ,",\
                        "area_left=",area_left, ", days_left=",time_left, ",rate=",rate,") ")
                    print("======================================>>>>>>>>>>>>>>>>>>")
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
                except:
                    e = sys.exc_info()
                    exc_type, exc_obj, exc_tb = e[0],e[1],e[2]
                    where = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    line = exc_tb.tb_lineno
                    trace = traceback.format_exc(e)
                    print(trace)
                    self.send_processing_error(e, where, line, trace)
                    sys.exit()
        else:
            econ_prob = 0
            econ_area = 0
            need_area = 11734.0
            quality = 1.0

        try:
            self.sumligoprob = getHexObservations.how_well_did_we_do(
                self.skymap, trigger_id, mapDir, camera, resolution)
        except:
            e = sys.exc_info()
            exc_type, exc_obj, exc_tb = e[0],e[1],e[2]
            where = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            line = exc_tb.tb_lineno
            trace = traceback.format_exc(e)
            print(trace)
            self.send_processing_error(e, where, line, trace)
            sys.exit()

        self.best_slot = best_slot
        self.n_slots = n_slots
        self.first_slot = first_slot
        self.econ_prob = econ_prob
        self.econ_area = econ_area
        self.need_area = need_area
        self.quality = quality


        np.savez(os.path.join(self.work_area, 'mapmaker_results.npz')
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
        print('-'*20+'>','LIGO PROB: %.3f \tLIGO X DES PROB: %.3f' % (self.sumligoprob,integrated_prob))
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
                     time_processed=self.recycler_mjd,
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
                     time_processed=self.recycler_mjd,
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

        bestslot_name = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
        cp_string = os.path.join(self.work_area, bestslot_name) + ' ' + image_dir +"/"
        trigger_id =  self.trigger_id 
        trigger_best_slot =  trigger_id + "-" + str(self.best_slot) 
        
        if self.n_slots<1:
            counter = getHexObservations.nothingToObserveShowSomething(trigger_id, self.work_area, self.mapspath)
            bestslot_name = trigger_best_slot + "-ligo-eq.png"
            oname = trigger_id + "-probabilityPlot.png"
            os.system('cp ' + cp_string + oname)
        if True:
            bestslot_name = trigger_best_slot + "-maglim-eq.png"
            oname = trigger_id + "_limitingMagMap.png"
            os.system('cp ' + cp_string + oname)
            bestslot_name = trigger_best_slot + "-prob-eq.png"
            oname = trigger_id + "_sourceProbMap.png"
            os.system('cp ' + cp_string + oname)
            bestslot_name = trigger_best_slot + "-ligo-eq.png"
            oname = trigger_id + "_LIGO.png"
            os.system('cp ' + cp_string + oname)
            bestslot_name = trigger_best_slot + "-probXligo-eq.png"
            oname = trigger_id + "_sourceProbxLIGO.png"
            os.system('cp ' + cp_string + oname)
            # DESGW observation map
            os.system('cp ' + cp_string + oname)
            # probability plot
            name = trigger_id + "-probabilityPlot.png"
            os.system('cp ' + os.path.join(self.work_area, name) + ' ' + image_dir)
            #raw_input('getting contours stopped')

        return

    def makeJSON(self, config):

        mapmakerresults = np.load(os.path.join(self.work_area, 'mapmaker_results.npz'))

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
        jsonname = self.trigger_id + "_"+ self.trigger_dir +"_JSON.zip"
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
                print("zero probability, thus no jsonFile at ", jsonFile)
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

        text = 'DESGW Webpage Created. See \nhttp://des-ops.fnal.gov:8080/desgw/Triggers/' + self.trigger_id + '/' + self.trigger_id + '_' +self.trigger_dir + '_trigger.html'

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


        msg['Subject'] = 'DESGW Webpage Created for ' + self.trigger_id + ' Map: '+self.trigger_dir
        msg['From'] = me
        for you in yous:
            msg['To'] = you

            s = smtplib.SMTP('localhost')
            s.sendmail(me, [you], msg.as_string())
            s.quit()
        print('Email sent...')
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
        msg['Subject'] = 'Trigger ' + self.trigger_id + ' '+self.trigger_dir+ ' Processing FAILED!'
        msg['From'] = me
        for you in yous:
            msg['To'] = you
            s = smtplib.SMTP('localhost')
            s.sendmail(me, [you], msg.as_string())
            s.quit()
        print('Email sent...')
        return

    def updateTriggerIndex(self, real_or_sim=None):
        website = self.website
        if real_or_sim == 'real':
            fff = website + 'real-trigger_list.txt'
        if real_or_sim == 'sim':
            fff = website + 'test-trigger_list.txt'

        if not os.path.exists(fff) :
            lines = []
        else :
            l = open(fff, 'r')
            lines = l.readlines()
            l.close()

        a = open(fff, 'a')
        if lines == [] :
            a.write(self.trigger_id + ' ' + self.work_area + '\n')
        else  :
            triggers = []
            for line in lines:
                triggers.append(line.split(' ')[0])
            if not self.trigger_id in np.unique(triggers):
                a.write(self.trigger_id + ' ' + self.work_area + '\n')
        a.close()
        tp.make_index_page(website, real_or_sim=real_or_sim)
        return

    def make_cumulative_probs(self):
        GW_website_dir = os.path.join(self.website, '/Triggers/')
        sim_study_dir = '/data/des41.a/data/desgw/maininjector/sims_study/data'
        radecfile = os.path.join(self.work_area, 'maps', self.trigger_id + '-ra-dec-id-prob-mjd-slot.txt')
        cumprobs_file = os.path.join(self.work_area, self.trigger_id + '-and-sim-cumprobs.png') 
        print(['python', './python/cumulative_plots.py', '-d',
               sim_study_dir, '-p', self.work_area, '-e', self.trigger_id,
               '-f', radecfile])
        subprocess.call(['python', './python/cumulative_plots.py', '-d',
               sim_study_dir, '-p', self.work_area, '-e', self.trigger_id, 
                '-f',radecfile])
        os.system('scp '+ cumprobs_file + ' ' + GW_website_dir + self.trigger_id + '/images/')

        
    def updateWebpage(self,real_or_sim):
        trigger_id = self.trigger_id
        trigger_dir = self.trigger_dir
        GW_website_dir = self.website
        desweb = "codemanager@desweb.fnal.gov:/des_web/www/html/desgw/"
        GW_website_dir_t = GW_website_dir + "Triggers/"
        desweb_t = desweb + "Triggers/"
        desweb_t2 = desweb + "Triggers/" + trigger_id 
        trigger_html = os.path.join(self.master_dir, trigger_id +'_' + 
            trigger_dir + '_trigger.html') 

        os.system('scp -r ' + GW_website_dir_t + self.trigger_id+ ' ' + desweb_t)
        os.system('scp ' + GW_website_dir + '/* ' + desweb)
        tp.makeNewPage(self.master_dir,
            os.path.join(self.master_dir, trigger_id +'_'+ trigger_dir+ '_trigger.html'), 
            trigger_id,self.event_paramfile,trigger_dir, real_or_sim=real_or_sim)
        os.system('scp -r ' + trigger_html + ' ' + desweb_t2 + "/")
        os.system('scp -r ' + trigger_html + ' ' + desweb_t2 +'_trigger.html')

        #os.system('cp ' + self.work_area +  trigger_dir + '/'+ 'recycler.log ' + self.website_jsonpath)
        return

    def getmjd(self,datet):
        mjd_epoch = datetime.datetime(1858, 11, 17)
        print('FIX ME UTC OR NOT?!?')
        mjdd = datet-mjd_epoch
        mjd = 5./24. + mjdd.total_seconds() / 3600. / 24.
        return mjd

    def mjd_to_datetime(self,mjd):
        mjd_epoch = datetime.datetime(1858, 11, 17, tzinfo=pytz.utc)
        d = mjd_epoch + datetime.timedelta(mjd)
        return d

    def makeObservingPlots(self):
        try:
            if not self.config['skipPlots']:
                n_plots = getHexObservations.makeObservingPlots(
                    self.n_slots, self.trigger_id, self.best_slot, self.outputDir, self.mapDir, self.camera, allSky=True )

                image_dir = self.website_imagespath
                map_dir = self.mapspath

                bestslot_name = self.trigger_id + "-" + str(self.best_slot) + "-ligo-eq.png"
                if self.n_slots < 1:
                    counter = getHexObservations.nothingToObserveShowSomething(
                        self.trigger_id, self.work_area, self.mapspath)
                    oname = self.trigger_id + "-observingPlot.gif"
                    os.system('cp ' + os.path.join(self.work_area, bestslot_name) + ' ' + os.path.join(image_dir, oname))
                    oname = self.trigger_id + "-probabilityPlot.png"
                    os.system('cp ' + os.path.join(self.work_area, bestslot_name) + ' ' + os.path.join(image_dir, oname))
                # if self.n_slots > 0:
                if True:
                    print('Converting Observing Plots to .gif')
                    files=np.array(glob.glob(os.path.join(map_dir, self.trigger_id)+'-observingPlot-*.png'))
                    split=[i.split('-', 2)[2] for i in files]
                    number=[i.split('.', 1)[0] for i in split]
                    f=np.array(number).astype(np.int)
                    maximum=str(np.max(f))
                    minimum=str(np.min(f))
                    os.system('convert $(for ((a='+minimum+'; a<='+maximum+'; a++)); do printf -- "-delay 50 ' + os.path.join(map_dir,
                                        self.trigger_id) + '-observingPlot-%s.png " $a; done;) ' + os.path.join(
                        map_dir, self.trigger_id) + '-observingPlot.gif')
                    # os.system('convert -delay 70 -loop 0 '+os.path.join(map_dir,self.trigger_id)+'-observingPlot-*.png '+
                    #          os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif')
                    os.system('cp ' + os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif ' + image_dir)
            #string = "$(ls -v {}-observingPlot*)"


        except:
            e = sys.exc_info()
            exc_type, exc_obj, exc_tb = e[0],e[1],e[2]
            where = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            line = exc_tb.tb_lineno
            trace = traceback.format_exc(e)
            print(trace)
            self.send_processing_error(e, where, line, trace)
            sys.exit()

if __name__ == "__main__":

    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "tp:tid:mjd:exp:sky",
            longopts=["triggerpath=", "triggerid=", "mjd=", "exposure_length=", "skymapfilename="])

    except getopt.GetoptError as err:
        print(str(err))
        print("Error : incorrect option or missing argument.")
        print(__doc__)
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
        print('Option')
        print(o)
        print(a)
        print('-----')
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
            print("Warning: option", o, "with argument", a, "is not recognized")

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
               print('Could not find skymap url file')

        if 'bayestar' in skymap_filename:
            print('bayestar' * 50)

        try:
            try:
                mjd = float(mjd)
            except:
                badtriggers = open('badtriggers.txt', 'a')
                badtriggers.write(trigger_id + '\n')
                print('WARNING: Could not convert mjd to float. Trigger: ' + trigger_id + ' flagged as bad.')
# here is where the object is made, and parts of it are filed in
            master_dir = os.path.join(trigger_path, trigger_id)
            e = event(skymap_filename, master_dir, trigger_id, mjd, config)

# e has variables and code assocaiated with it. The mapMaker is called "e" or "self"

            e.mapMaker(trigger_id, skymap_filename, config)
            e.getContours(config)
            jsonfilelist = e.makeJSON(config)
            e.make_cumulative_probs()
            e.updateTriggerIndex(real_or_sim=real_or_sim) # generates the homepage 
            e.updateWebpage(real_or_sim) #make a blank page with the basic info that is available
            e.makeObservingPlots()
            e.getContours(config)
            e.send_nonurgent_Email()
            e.updateWebpage(real_or_sim)

        except KeyError:
            print("Unexpected error:", sys.exc_info())
            badtriggers = open('badtriggers.txt', 'a')
            badtriggers.write(trigger_id + '\n')
    #############################################################################

    print('Done')
