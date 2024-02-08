# To-do:
# Strategy code integration
# How to pick models
# Need to find code that calculates teff based on moon
# Main injector does this for every pixel throughout the night (images in gif)
# Need to steal the bit that does the quick teff calc for the zenith for all filters during the time we are going to observer
# Make simple version of map making – get_sky_info
# Takes the ra and dec of event
# Calculates time of zenith (highest point in sky)
# Calculates where the moon is at time of zenith
# Gets teff
# Best  time to observe?
# MI workflow
# Run get_sky_info
# Andre will look into this
# Jim says look at simplicity.py
# Use teff and time of obs as input to run strategy code
# Returns exposure time, exposure areas/prob level, filters


import datetime
import os
import traceback
import sys
import yaml
import glob

import numpy as np
import matplotlib.pyplot as plt

import obsSlots
import observations
# import send_texts_and_emails
import mags
import trigger_pages as tp
import gw_map_configure

#
# calc produces information about a single nights observing.
#   the information about what to observe is in the *ra-dec-*txt file that strategy produces
#   the "ra_dec_file_dir" tells the code where that lives
#   tigger_id tells it what the name of the *ra-dec*.txt file is
# The expTime is used only in the limiting magnitude calculation, not in what to observe when
#


class Event:
    def __init__(self, skymap_filename: str, master_dir: str, trigger_id: int,
                 mjd: float, config: str, official, hasrem):
        """
        Event base blass for Main-Injector.

        Attributes:
        -----------
        master_dir: str
            the directory containing all trigger subdirectories
        trigger_id: str
            the name of the trigger event that comes from LIGO
        mjd: float
            the modified julian date of the event (in FITS header, among other places)
        config: str
            the config filename with the parameters that controls the routine
        official:

        hasrem: bool

        """

        
        # set up the working directories
        self.modify_filesystem(
            skymap_filename, master_dir, trigger_id, mjd, hasrem)
        work_area = self.work_area
        self.trigger_id = trigger_id
        # read config file
        if config["force_recycler_mjd"]:
            self.recycler_mjd = config["recycler_mjd"]
        else:
            self.recycler_mjd = self.getmjd(datetime.datetime.now())

        self.event_paramfile = os.path.join(
            work_area, trigger_id + '_params.npz')
        #self.event_paramfile = os.path.join(work_area,'..', trigger_id + '_params.npz')

        # print(self.event_paramfile)
        # asdf
        # raw_input()
        self.weHaveParamFile = True
        try:
            self.event_params = np.load(self.event_paramfile)
        except:
            self.event_params = {
                'ETA': 'NAN',
                'FAR': 'NAN',
                'ChirpMass': 'NAN',
                'MaxDistance': 'NAN',
                'M1': 'NAN',
                'M2': 'NAN',
                'boc': 'NAN',
                'CentralFreq': 'NAN'
            }

            self.weHaveParamFile = False

        # asdf
        yaml_dir = os.path.join(work_area, 'strategy.yaml')
        os.system('cp recycler.yaml ' + yaml_dir)
        print('***** Copied recycler.yaml to ' +
              yaml_dir + ' for future reference *****')
        #os.system('kinit -k -t /var/keytab/desgw.keytab desgw/des/des41.fnal.gov@FNAL.GOV')
        self.official = official

    def modify_filesystem(self, skymap_filename, master_dir, trigger_id, mjd, hasrem):
        """
        Parses the skymap_filename and master_dir to the correct format.

        Parameters:
        -----------
        skymap_filename: str
            filename path og an given skymap
        master_dir: str

        """
        #skymap_filename, master_dir, trigger_id, mjd, config, official, hasrem

        # master_dir is the directory holding all the trigger directories
        # trigger_dir is the directory holding everything to do with a single trigger
        # work_area is master_dir + trigger_dir
        # work_area is derived
        # trigger_dir is derived from the path of the skymap_file

        skymap_filename = skymap_filename.strip()
        trigger_dir = skymap_filename.split('/')[-1].split('.')[0]
        self.trigger_dir = trigger_dir
        os.system('touch '+skymap_filename+'.processing')
        if not os.path.exists(master_dir+'/' + trigger_dir):
            os.system('mkdir '+master_dir+'/' + trigger_dir)

        self.master_dir = master_dir
        if hasrem:
            work_area = master_dir+'/' + trigger_dir + '/hasrem/'
        else:
            work_area = master_dir+'/' + trigger_dir + '/norem/'

        os.system('cp '+skymap_filename.strip()+' '+work_area)
        os.system('cp '+self.master_dir+'/' +
                  trigger_id + '_params.npz '+work_area)

        self.skymap = skymap_filename

        # Setup website directories
        website = os.environ["DESGW_WEB"]
        if not os.path.exists(website):
            os.mkdir(website)

        self.mapspath = os.path.join(work_area, "maps/")
        if not os.path.exists(self.mapspath):
            os.makedirs(self.mapspath)

        self.imagespath = website + "Triggers/" + \
            trigger_id + "/"+work_area.split('/')[-1]
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
        self.website = website

# Let's guess that mapMaker is the counterpart to recyc.mainInjector from
# desgw-maps.

    def read_config_file(self, trigger_id, skymap, config, hasrem,
                         snarf_mi_maps=False, start_slot=-1,
                         do_nslots=-1,  mi_map_dir="./"):

        # debug
        debug = config["debug"]
        camera = config["camera"]
        resolution = float(config["resolution"])
        overhead = config["overhead"]
        allSky = config['allSky']
        area_per_hex = config["area_per_hex"]
        start_of_season = config["start_of_season"]
        end_of_season = config["end_of_season"]
        events_observed = config["events_observed"]
        skipAll = config["skipAll"]
        mjd = self.mjd
        outputDir = self.work_area
        mapDir = self.mapspath
        recycler_mjd = self.recycler_mjd
        kasen_fraction = config['kasen_fraction']
        debug = config["debug"]
        camera = config["camera"]
        resolution = config["resolution"]
        do_make_maps = config["do_make_maps"]
        do_make_hexes = config["do_make_hexes"]
        do_make_jsons = config["do_make_jsons"]
        do_make_gifs = config["do_make_gifs"]
        days_since_burst = config["days_since_burst"]

        self.distance = 1.

    def set_strategy(self, config: str) -> None:
        """
        Strategy setup for observation plans

        Parameters:
        -----------

        config: configuration filepath
        """

        try:
            days_since_burst = config["days_since_burst"]
        except:
            pass
    # strategy
        exposure_length_rem = config["exposure_length_Rem"]
        filter_list_rem = config["exposure_filter_Rem"]
        maxHexesPerSlot_rem = config["maxHexesPerSlot_Rem"]
        exposure_length_bh = config["exposure_length_BH"]
        filter_list_bh = config["exposure_filter_BH"]
        maxHexesPerSlot_bh = config["maxHexesPerSlot_BH"]

        # ag added
        exposure_tiling_rem = config["exposure_tiling_Rem"]
        exposure_tiling_bh = config["exposure_tiling_BH"]
        max_number_of_hexes_to_do = config["max_number_of_hexes_to_do"]
        hoursAvailable = 20.
        self.time_budget = hoursAvailable

        #hasremnant = self.event_params['hasremnant']

        if self.hasrem:
            trigger_type = 'hasrem'
        else:
            trigger_type = 'norem'

    # configure strategy for the event type
        if trigger_type == "hasrem":
            exposure_length = exposure_length_rem
            filter_list = filter_list_rem
            maxHexesPerSlot = maxHexesPerSlot_rem
            tiling_list = exposure_tiling_rem
            propid = config['propid_Rem']
        elif trigger_type == "norem":
            exposure_length = exposure_length_bh
            filter_list = filter_list_bh
            maxHexesPerSlot = maxHexesPerSlot_bh
            tiling_list = exposure_tiling_bh
            propid = config['propid_BH']

        else:
            raise Exception(
                "trigger_type={}  ! Can only compute BH or Rem".format(trigger_type))
        exposure_length = np.array(exposure_length)

        gif_resolution = config['gif_resolution']

        gw_map_control = gw_map_configure.control(resolution,
                                                  outputDir,
                                                  debug,
                                                  allSky=allSky,
                                                  snarf_mi_maps=snarf_mi_maps,
                                                  mi_map_dir=mi_map_dir,
                                                  gif_resolution=gif_resolution)

        gw_map_trigger = gw_map_configure.Trigger(skymap,
                                                  trigger_id,
                                                  trigger_type,
                                                  resolution,
                                                  days_since_burst=days_since_burst)
        # ag test jul 30
        use_teff = 1.0
        gw_map_strategy = gw_map_configure.Strategy(camera,
                                                    exposure_length,
                                                    filter_list,
                                                    tiling_list,
                                                    maxHexesPerSlot,
                                                    hoursAvailable,
                                                    propid,
                                                    max_number_of_hexes_to_do,
                                                    kasen_fraction,
                                                    use_teff)
        # strat need max number of hexes and tiling list

        gw_map_results = gw_map_configure.results()

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        if do_make_maps:
            # make the computationally expensive maps of everything
            observations.make_maps(
                gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)

        if do_make_hexes:
            # compute the best observations
            observations.make_hexes(
                gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results,
                start_slot=start_slot, do_nslots=do_nslots)
            # if start_slot = -1, do_nslots = -1, then do whole night, as if called like:
            #    observations.make_hexes(
            #        gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)
            # the web page maker version of main injector should default to -1, -1

        if do_make_jsons:
            # make the jsons
            observations.make_jsons(
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)

        if do_make_gifs:
            observations.makeGifs(
            gw_map_trigger, gw_map_strategy, gw_map_control, gw_map_results)

        allSky = config['allSky']

        try:
            eventtype = self.event_params['boc']
        except:
            eventtype = None

        try:
            probhasns = self.event_params['probhasns']
        except:
            probhasns = 0.  # for old maps...

        ra, dec, id, self.prob, mjd, slotNum, dist = \
            obsSlots.readObservingRecord(self.trigger_id, mapDir)

        self.slotNum = slotNum

        integrated_prob = np.sum(self.prob)
        try:
            print('-'*20+'>', 'LIGO PROB: %.3f \tLIGO X DES PROB: %.3f' %
                  (gw_map_results.sum_ligo_prob, integrated_prob))
        except:
            pass

        self.best_slot = gw_map_results.best_slot
        self.n_slots = gw_map_results.n_slots
        self.first_slot = gw_map_results.first_slot
        self.exposure_length = exposure_length
#        if do_make_maps:
        if 1 == 1:
            np.savez(self.event_paramfile,
                     MJD=self.mjd,
                     ETA=self.event_params['ETA'],
                     FAR=self.event_params['FAR'],
                     ChirpMass=self.event_params['ChirpMass'],
                     MaxDistance=self.event_params['MaxDistance'],
                     DESXLIGO_prob=integrated_prob,
                     LIGO_prob=0.0,  # AG changed gw_map_results.sum_ligo_prob -> 0.0 because sum_ligo_prob doesn't make sense here??? sept 6 2022
                     M1=self.event_params['M1'],
                     M2=self.event_params['M2'],
                     nHexes=self.prob.size,
                     time_processed=self.recycler_mjd,
                     boc=self.event_params['boc'],
                     CentralFreq=self.event_params['CentralFreq'],
                     best_slot=gw_map_results.best_slot,
                     n_slots=gw_map_results.n_slots,
                     first_slot=gw_map_results.first_slot,
                     econ_prob=0,  # self.econ_prob,
                     econ_area=0,  # self.econ_area,
                     need_area=0,  # self.need_area,
                     quality=0,  # self.quality,
                     codeDistance=self.distance,
                     exposure_times=exposure_length,
                     exposure_filter=filter_list,
                     hours=self.time_budget,
                     nvisits=-999,  # config['nvisits'],
                     mapname='NAN',
                     filename=self.skymap,
                     gethexobstype=trigger_type,
                     # probhasns=self.probhasns
                     probhasns=probhasns
                     )

        map_dir = mapDir
        jsonname = self.trigger_id + "_" + self.trigger_dir + "_JSON.zip"
        jsonFile = os.path.join(map_dir, jsonname)
        jsonfilelistld = os.listdir(map_dir)
        jsonfilelist = []
        for f in jsonfilelistld:
            if '-tmp' in f:
                os.remove(os.path.join(map_dir, f))
            elif '.json' in f:
                jsonfilelist.append(f)

        os.system(f"zip -j {jsonFile} {self.mapspath}/*0.json")
        os.system(f"cp {jsonFile} {self.website_jsonpath}")

        os.system(f"cp + {os.path.join(map_dir, self.trigger_id)}\
            _centered_animate.gif {self.website_imagespath}")
        
        os.system(f"cp + {os.path.join(map_dir, self.trigger_id)}\
            _animate.gif {self.website_imagespath}")
        
        os.system(f"cp + {os.path.join(map_dir, self.trigger_id)}\
            *.png {self.website_imagespath}")
        

    def pp(self):
        """
        Probability plot for Event
        """
        plt.clf()
        plt.plot(self.slotNum, self.prob, label='Total Prob %.3f' %
                 np.sum(self.prob))
        plt.scatter(self.slotNum, self.prob)
        plt.xlabel('Slot Number')
        plt.ylabel('Probability Per Slot')
        plt.title('decam*ligo')
        plt.legend()
        name = self.trigger_id + "-probabilityPlot.png"
        plt.savefig(os.path.join(self.mapspath, name))
        plt.clf()

    def getContours(self, config):
        import matplotlib.pyplot as plt

        # if exposure_length is None:
        #    exposure_length = config["exposure_length"]
        exposure_length = self.exposure_length
        image_dir = self.website_imagespath
        map_dir = self.mapspath

        bestslot_name = self.trigger_id + "-" + \
            str(self.best_slot) + "-ligo-eq.png"
        cp_string = os.path.join(
            self.work_area, bestslot_name) + ' ' + image_dir + "/"
        trigger_id = self.trigger_id
        trigger_best_slot = trigger_id + "-" + str(self.best_slot)

        if True:
            bestslot_name = trigger_best_slot + "-maglim-eq.png"
            cp_string = os.path.join(
                map_dir, bestslot_name) + ' ' + image_dir + "/"
            oname = trigger_id + "_limitingMagMap.png"
            os.system('cp ' + cp_string + oname)
            bestslot_name = trigger_best_slot + "-prob-eq.png"
            cp_string = os.path.join(
                map_dir, bestslot_name) + ' ' + image_dir + "/"
            oname = trigger_id + "_sourceProbMap.png"
            os.system('cp ' + cp_string + oname)
            bestslot_name = trigger_best_slot + "-ligo-eq.png"
            cp_string = os.path.join(
                map_dir, bestslot_name) + ' ' + image_dir + "/"
            oname = trigger_id + "_LIGO.png"
            os.system('cp ' + cp_string + oname)
            bestslot_name = trigger_best_slot + "-probXligo-eq.png"
            cp_string = os.path.join(
                map_dir, bestslot_name) + ' ' + image_dir + "/"
            oname = trigger_id + "_sourceProbxLIGO.png"
            os.system('cp ' + cp_string + oname)
            # DESGW observation map
            os.system('cp ' + cp_string + oname)
            # probability plot
            self.makeProbabilityPlot()
            name = trigger_id + "-probabilityPlot.png"
            os.system('cp ' + os.path.join(map_dir, name) + ' ' + image_dir)
            #raw_input('getting contours stopped')

        return

    def makeJSON(self, config):

        mapmakerresults = np.load(os.path.join(
            self.work_area, 'mapmaker_results.npz'))

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
        jsonname = self.trigger_id + "_" + self.trigger_dir + "_JSON.zip"
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
            ra, dec, id, self.prob, mjd, slotNum, dist = \
                obsSlots.readObservingRecord(self.trigger_id, map_dir)
            self.slotNum = slotNum
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
                    os.chmod(os.path.join(self.mapspath, js), 0o777)

                os.system('zip -j ' + jsonFile + ' ' +
                          self.mapspath + '/*0.json')
                # except:
                #    print "no jsonFiles at ", jsonFile
        else:
            os.remove(jsonFile)
            os.system('zip -j ' + jsonFile + ' ' + self.mapspath + '/*0.json')
            os.system('cp ' + jsonFile + ' ' + self.website_jsonpath)
        return jsonfilelist

    def send_nonurgent_Email(self, sendtexts=False):
        text = f"""
        DESGW Webpage Created for REAL event.
        See http://des-ops.fnal.gov:8080/desgw/Triggers/\
            {self.trigger_id}/{self.trigger_id}_{self.trigger_dir}_trigger.html
        DO NOT REPLY TO THIS THREAD, NOT ALL USERS WILL SEE YOUR RESPONSE.
        """
        subject = f"""
        DESGW Webpage Created for REAL event {self.trigger_id}\
        Map: {self.trigger_dir} NOREPLY
        """

        send_texts_and_emails.send(subject, text, official=self.official)
        print('Email sent...')
        return

    def send_processing_error(self, error, where, line, trace):
        import smtplib
        from email.mime.text import MIMEText

        message = 'Processing Failed for REAL Trigger ' + str(self.trigger_id) + '\n\nFunction: ' + str(
            where) + '\n\nLine ' + str(line) + ' of recycler.py\n\nError: ' + str(error) + '\n\n'
        message += '-' * 60
        message += '\n'
        message += trace
        message += '\n'
        message += '-' * 60
        message += '\n'

        subject = 'REAL Trigger ' + self.trigger_id + \
            ' '+self.trigger_dir + ' Processing FAILED!'
        send_texts_and_emails.send(subject, message, official=self.official)
        print('Email sent...')
        return

    def updateTriggerIndex(self, real_or_sim=None):
        website = self.website
        if real_or_sim == 'real':
            fff = website + 'real-trigger_list.txt'
        if real_or_sim == 'sim':
            fff = website + 'test-trigger_list.txt'

        if not os.path.exists(fff):
            lines = []
        else:
            l = open(fff, 'r')
            lines = l.readlines()
            l.close()

        a = open(fff, 'a')
        if lines == []:
            a.write(self.trigger_id + ' ' + self.work_area + '\n')
        else:
            triggers = []
            for line in lines:
                triggers.append(line.split(' ')[0])
            if not self.trigger_id in np.unique(triggers):
                a.write(self.trigger_id + ' ' + self.work_area + '\n')
        a.close()
        tp.make_index_page(website, real_or_sim=real_or_sim)
        return

    def updateWebpage(self, real_or_sim):
        trigger_id = self.trigger_id
        trigger_dir = self.trigger_dir
        GW_website_dir = self.website
        desweb = "codemanager@desweb.fnal.gov:/des_web/www/html/desgw/"
        GW_website_dir_t = GW_website_dir + "Triggers/"
        desweb_t = desweb + "Triggers/"
        desweb_t2 = desweb + "Triggers/" + trigger_id
        trigger_html = os.path.join(self.master_dir, trigger_id + '_' +
                                    trigger_dir + '_trigger.html')

        print('scp -r ' + GW_website_dir_t + self.trigger_id + ' ' + desweb_t)
        print('scp ' + GW_website_dir + '/* ' + desweb)
        # asdf
        os.system('cp '+trigger_html+' '+GW_website_dir)

        #os.system('scp -r ' + GW_website_dir_t + self.trigger_id + ' ' + desweb_t)
        #os.system('scp ' + GW_website_dir + '/* ' + desweb)
        # master_dir,outfilename,trigger_id,event_paramfile,mapfolder,processing_param_file=None,real_or_sim='real',secondtimearound=False
        print(self.master_dir)
        print(os.path.join(self.master_dir, trigger_id +
              '_' + trigger_dir + '_trigger.html'))
        #os.system('cp '+trigger_html+' '+GW_website_dir)
        print(trigger_id, self.event_paramfile, trigger_dir)

        tp.makeNewPage(self.master_dir, os.path.join(self.master_dir, trigger_id + '_' + trigger_dir +
                       '_trigger.html'), trigger_id, self.event_paramfile, trigger_dir, real_or_sim=real_or_sim)
        print('here1')
        print('scp -r ' + trigger_html + ' ' + desweb_t2 + "/")
        #os.system('scp -r ' + trigger_html + ' ' + desweb_t2 + "/")
        print('here2')
        #os.system('scp -r ' + trigger_html + ' ' + desweb_t2 + '_trigger.html')
        print('here3')
        #os.system('cp ' + self.master_dir + '/' + trigger_dir + '/' + trigger_dir + '_recycler.log ' + self.website_jsonpath)
        return

    def make_obs_plots(self):
        try:
            if not self.config['skipPlots']:
                # n_plots = observations.makeObservingPlots(
                #    self.n_slots, self.trigger_id, self.best_slot, self.outputDir, self.mapDir, self.camera, allSky=True )

                image_dir = self.website_imagespath
                map_dir = self.mapspath

                bestslot_name = self.trigger_id + "-" + \
                    str(self.best_slot) + "-ligo-eq.png"
                if self.n_slots < 1:
                    # counter = observations.nothingToObserveShowSomething(
                    #    self.trigger_id, self.work_area, self.mapspath)
                    oname = self.trigger_id + "-observingPlot.gif"
                    os.system('cp ' + os.path.join(self.work_area,
                              bestslot_name) + ' ' + os.path.join(image_dir, oname))
                    oname = self.trigger_id + "-probabilityPlot.png"
                    os.system('cp ' + os.path.join(self.work_area,
                              bestslot_name) + ' ' + os.path.join(image_dir, oname))
                # if self.n_slots > 0:
                if True:
                    print('Converting Observing Plots to .gif')
                    files = np.array(glob.glob(os.path.join(
                        map_dir, self.trigger_id)+'-observingPlot-*.png'))
                    split = [i.split('-', 2)[2] for i in files]
                    number = [i.split('.', 1)[0] for i in split]
                    f = np.array(number).astype(np.int)
                    maximum = str(np.max(f))
                    minimum = str(np.min(f))
                    os.system('convert $(for ((a='+minimum+'; a<='+maximum+'; a++)); do printf -- "-delay 50 ' + os.path.join(map_dir,
                                                                                                                              self.trigger_id) + '-observingPlot-%s.png " $a; done;) ' + os.path.join(
                        map_dir, self.trigger_id) + '-observingPlot.gif')
                    # os.system('convert -delay 70 -loop 0 '+os.path.join(map_dir,self.trigger_id)+'-observingPlot-*.png '+
                    #          os.path.join(map_dir, self.trigger_id) + '-observingPlot.gif')
                    os.system('cp ' + os.path.join(map_dir,
                              self.trigger_id) + '-observingPlot.gif ' + image_dir)
            #string = "$(ls -v {}-observingPlot*)"

        except:
            e = sys.exc_info()
            exc_type, exc_obj, exc_tb = e[0], e[1], e[2]
            where = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            line = exc_tb.tb_lineno
            trace = traceback.format_exc(e)
            print(trace)
            self.send_processing_error(e, where, line, trace)
            sys.exit()


def calc(trigger_id, ra_dec_file_dir, mjd, expTime, filter, best=False, camera="decam"):
    # do dark siren stuff
    #ra,dec,prob,filter = np.genfromtxt(trigger_id,unpack=True)
    #filter = np.genfromtxt(trigger_id,unpack=True,dtype="str",usecols=3)
    #filter = "g"

    # get the hexes
    ra, dec, id, prob, obs_mjd, slotNum, dist = obsSlots.readObservingRecord(
        trigger_id, ra_dec_file_dir)
    # find the highest prob hex
    best_ix = np.argmax(prob)

    # find night statistics
    night, sunset, sunrise = mags.findNightDuration(mjd, camera)
    night = np.float(night)*24.
    sunset = np.float(sunset)
    sunrise = np.float(sunrise)

    # will work every hour from sunset to sunrise
    mjd_list = np.arange(sunset, sunrise+1./24., 1./24.)

    # calculate the limiting magnitudes
    limit_mag = []
    best_limit_mag = []
    moon_sep = []
    moon_phase = []
    obs = mags.observed(ra, dec, prob, sunset, do_maps=False, verbose=False)

    for mjd in mjd_list:
        obs.resetTime(mjd)
        obs.limitMag(filter, exposure=expTime)
        limit_mag.append(obs.maglim)
        best_limit_mag.append(obs.maglim[best_ix])
        moon_sep.append(obs.moonSep[0]*360./2/np.pi)
        moon_phase.append(obs.moonPhase)
    limit_mag = np.vstack(limit_mag)
    ix = limit_mag < 0
    limit_mag[ix] = 0
    moon_sep = np.array(moon_sep)
    moon_phase = np.array(moon_phase)
    # moon phase: where 0 = full, 90 equals half, and  180 = new
    # convert to %full
    moon_phase = ((180.-moon_phase)/180.)*100.
    bins = np.arange(21, 25.5, 0.5)

    # now print the answers
    print("{:15s}  {:10s}".format("MJD".rjust(15), "best".rjust(10))),
    print(" {:10s}".format("moon sep".rjust(10))),
    print(" {:10s}".format("moon phase".rjust(10))),
    print("      hexes w/ limiting mag in "),
    print("")
    print("{:15s}  {:10s}".format("days".rjust(15), "hex".rjust(10))),
    print(" {:10s}".format("degrees".rjust(10))),
    print("  {:10s}".format("% full ".rjust(10))),
    for i in range(bins.size-1):
        print("{:4.1f}-".format(bins[i])),
    print("{:4.1f}".format(bins[i])),
    print("")
    for i in range(0, mjd_list.size):
        print("{:15.4f} ".format(mjd_list[i]))
        print("{:10.2f} ".format(best_limit_mag[i]))
        print("{:10.2f} ".format(moon_sep[i]))
        print("{:10.0f}  ".format(moon_phase[i])),
        counts = np.histogram(limit_mag[i], bins)
        # print "\n counts ", np.shape(counts), counts[0], counts[1]
        for j in range(bins.size-1):
            if counts[0][j] != 0:
                print("{:3d}  ".format(counts[0][j])),
            else:
                print("     "),
        print("")
    return obs
