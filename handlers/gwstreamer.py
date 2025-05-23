from yaml.loader import SafeLoader
from base64 import b64decode
from datetime import datetime
from threading import Timer
from gcn_kafka import Consumer
from io import BytesIO
from astropy.table import Table
from argparse import ArgumentParser
from astropy.time import Time
from astropy.io import fits
from pathlib import Path
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.bayestar import rasterize
from .slack import SlackBot
from .emails import EmailBot
import glob
import time
import subprocess
import os
os.environ["API_BASE_URL"] = "https://desgw-api-physics-soares-santos-flaskapi.apps.gnosis.lsa.umich.edu/api/v0.1/"
import numpy as np
import astropy_healpix as ah
import pprint
import json
import yaml
from .short_latency_plots import make_plots_initial, make_alert_skymap,make_agn_plot
import sys
sys.path.insert(0, '/data/des70.a/data/desgw/O4/Main-Injector-O4b/desgw_db_writer')
import desgw_db_writer.api as DESGWApi
import pickle

def elapsedTimeString(start):
    elapsed = int(time.time() - start)
    return "{}h {:02d}m {:02d}s".format(elapsed//(60*60), elapsed//60%60, elapsed%60)


def chirp_mass(distmean, diststd, area90, area50):
    test_event = [[distmean, diststd, area90, area50]] 

    with open('/data/des70.a/data/desgw/O4/Main-Injector-O4b/models/model.dat', 'rb') as f:
        clf = pickle.load(f)

    pdf = clf.predict_proba(test_event)[0]
    if clf.missing_bins.shape[0]:
        for c in clf.missing_bins:
            pdf = np.insert(pdf, c-1, 0.)
    expv = np.sum(pdf*clf.midpoints)
    std = np.sqrt(np.sum(pdf*clf.midpoints*clf.midpoints) - expv**2)

    return expv, std
class GWStreamer():

    def __init__(self, mode):

        mi_dir = os.path.dirname(os.path.abspath(__file__))
        root = os.path.dirname(mi_dir)

        self.mode = mode
        self.OUTPUT_PATH = None
        self.OUTPUT_TRIGGER = None
        self.SKYMAP_OUTPUT = None
        self._ROOT = root
        self.email_bot = EmailBot(mode=mode)
        self.slack_bot = SlackBot(mode=mode)
        self.slack_bot2 = SlackBot(mode="dummy")
        self.weather_protocol_CTIO = "https://noirlab.edu/science/index.php/observing-noirlab/observing-ctio/cerro-tololo/bad-weather-protocol-at-ctio"
        self.weather_forecast = "https://www.wunderground.com/forecast/cl/la-serena/ICOQUIMB2"
        self.weather_CTIO_currently = "https://noirlab.edu/science/observing-noirlab/weather-webcams/cerro-tololo/environmental-conditions"
        self.website_base_url = "https://des-ops.fnal.gov:8082/desgw-new/#/trigger/"
        self.FAR_threshold = 100 # in yr^-1
        self.api = DESGWApi.DESGWApi()

    def _get_max_prob(self, record: dict) -> tuple:
        """Returns Event Source and max probability."""
        Classification = {
            record['event']['classification']['BNS']: 'BNS',
            record['event']['classification']['NSBH']: 'NSBH',
            record['event']['classification']['BBH']: 'BBH',
            record['event']['classification']['Terrestrial']: 'Terrestrial'
        }
    
        classification_probs = [p[0] for p in Classification.items()]
        EVENT_PROB = classification_probs[np.argmax(classification_probs)]
        source = Classification[EVENT_PROB]
    
        return source, EVENT_PROB

    def _format_message(self,
                        trigger_id: str,
                        record: dict,
                        skymap_plot_link: str = None,
                        moon_plot_link: str = None,
                        retraction: bool = False) -> tuple:
        """
        Format Subject and Text for messages.
        Returns a tuple in the format (subject, text).
        
        Parameters:
        -----------

        trigger_id: int
            Trigger identifier to a env reported by LVK.

        record: dict
            GCN alert content.

        retraction: bool (Default: False)
            A Retraction content to know if an alert is a real
            astrophysical event or not.

        Returns:
        --------
            subject, text : tuple(str, str)

        """

        trigger_id = record['superevent_id']

        if retraction:
            subject = f"Retraction for {str(trigger_id)}"
            text = ''
            self.email_bot.send_email(subject, text)
            return subject, text


        FAR = record['event']['far']
        FAR = round(1./float(FAR)/60./60./24./365., 2)
        FAR = f'1 per {FAR} Years'

        source, prob_source = self._get_max_prob(record)

        text = f"*Superevent ID*: {record['superevent_id']}\n"+\
            f"*MJD*: {record['MJD']}\n"+\
            f'*Classification*: {source} ({prob_source})\n' +\
            f"BNS: {record['event']['classification']['BNS']:.3f}\n"+\
            f"NSBH: {record['event']['classification']['NSBH']:.3f}\n"+\
            f"BBH: {record['event']['classification']['BBH']:.3f}\n"+\
            f"*FAR*: {FAR}\n"+\
            f"*Alert Type*: {record['alert_type']}\n\n"+\
            f"Event Time: {record['event']['time']} \n"+\
            f"Alert Time: {record['time_created']}\n"+\
            f"Detectors: {record['event']['instruments']}\n"+\
            f"Terrestrial: {record['event']['classification']['Terrestrial']}\n"+\
            f"Has Remmant: {record['event']['properties']['HasRemnant']}\n"+\
            f"*DISTMEAN*: {record['event']['distmean']:.2f} Mpc\n"+\
            f"*DISTSIGMA*: {record['event']['distsigma']:.2} Mpc\n\n"+\
            f"Has Mass Gap: {record['event']['properties']['HasMassGap']}\n\n"+\
            f"GraceDB Link: {record['urls']['gracedb']}\n"+\
            f"Skymap Link: https://gracedb.ligo.org/apiweb/superevents/{trigger_id}/files/bayestar.png\n"+\
            f"SLIP Moon plot: https://des-ops.fnal.gov:8082/desgw-new/{trigger_id}/initial_plots/Moon.png\n"+\
            f"SLIP Initial Skymap: https://des-ops.fnal.gov:8082/desgw-new/{trigger_id}/initial_plots/initial_skymap.png"
        
        subject = ""
        return subject, text

    def flatten_skymap(self,
                       input_skymap: str, 
                       output_skymap: str, 
                       nside: int = 1024) -> None:
        """Flattens skymap
        
        Parameters:
        -----------
        input : str
            input skymap fullpath name.
        output: str
            output skymap fullpath name.

        nside : int
            HEALPix NSIDE.

        Returns:
        --------
            flatten skymap.
        """
        
        os.system("ligo-skymap-flatten --nside {} {} {}".format(nside, input_skymap, output_skymap))
        
        return
    
    def process_external_coinc(self,
                               trigger_id: str,
                               alert_type: str,
                               external_coinc: dict) -> None:
        """Process external_coinc from gcn alert"""

        subject = f'{alert_type} {trigger_id} have external coincidence!'

        COINC_FAR = external_coinc['time_coincidence_far']
        COINC_FAR = round(1./float(COINC_FAR)/60./60./24./365., 2)
        COINC_FAR = f'{COINC_FAR} Years'

        COINC_TIME_SKY_FAR = external_coinc['time_sky_position_coincidence_far']
        COINC_TIME_SKY_FAR = round(1./float(COINC_TIME_SKY_FAR)/60./60./24./365., 2)
        COINC_TIME_SKY_FAR = f'1 per {COINC_TIME_SKY_FAR} Years'

        text = f"GCN Notice: {external_coinc['gcn_notice_id']}\n" +\
            f"IVORN: {external_coinc['ivorn']}\n" +\
            f"Observatory: {external_coinc['observatory']}\n" +\
            f"Search: {external_coinc['search']}\n" +\
            f"Time difference between detecctions: {external_coinc['time_difference']}\n" +\
            f"Time coincidence FAR: {COINC_FAR} Years\n" +\
            f"FAR using Timing and Sky position: {COINC_TIME_SKY_FAR}\n"
            
        combined_skymap_str = external_coinc.get('combined_skymap', '')
        if combined_skymap_str != '':
            skymap_bytes = b64decode(combined_skymap_str)
            skymap = Table.read(BytesIO(skymap_bytes))
            combined_skymap_output = os.path.join(self.OUTPUT_TRIGGER,
                                                  'combined_skymap.multiorder.fits')
            combined_skymap_output_flatten = os.path.join(self.OUTPUT_TRIGGER,
                                                          'combined_skymap.fits.gz')
            
            skymap.write(combined_skymap_output, overwrite=True)        
            self.flatten_skymap(combined_skymap_output,
                                combined_skymap_output_flatten)
                    
    def handle(self, gcn_alert: dict) -> None:
        """
        Parsers gcn kafka notice
        """
        
        t0 = time.time()
        record = gcn_alert
        alert_type = record['alert_type']
        trigger_id = record['superevent_id']
        print('Handling trigger {}, alert type {}'.format(trigger_id, alert_type), flush=True)

        if self.mode == 'observation':
            if record['superevent_id'][0] != 'S':
                print('MI running in Observing Mode; unofficial trigger discarded', flush=True)
                return
            is_mock = False
            self.OUTPUT_PATH = os.path.join(self._ROOT,
                                            "OUTPUT/O4REAL")
        else:
            is_mock=True
                
        if (self.mode == 'test') or (self.mode == 'mock'):
            if record['superevent_id'][0] != 'M':
                print('MI running in a testing mode; non-testing trigger discarded', flush=True)
                return
            
            self.OUTPUT_PATH = os.path.join(self._ROOT,
                                "OUTPUT/TESTING")

        if record['alert_type'] == 'RETRACTION':
            subject, text = self._format_message(trigger_id=trigger_id,
                                                record=record,
                                                retraction=True)
            
            self.slack_bot.post_message(subject, text)
            return None

        if record['event']['group'] != 'CBC':
            print('Non-CBC event discarded', flush=True)
            self.slack_bot2.post_message(subject="", text="Non-CBC event {}, discarded".format(trigger_id))
            return

        if record['event']['pipeline'] == 'CWB':
            print('Coherent waveburst search event discarded', flush=True)
            self.slack_bot2.post_message(subject="", text="Coherent waveburst search event {}, discarded".format(trigger_id))
            return
        
        if alert_type == 'PRELIMINARY':
            self.OUTPUT_TRIGGER = os.path.join(self.OUTPUT_PATH,
                                               trigger_id,
                                               f'{alert_type}_0')
            self.alert_num =f'{alert_type}_0' 
            
            if os.path.isdir(self.OUTPUT_TRIGGER):
                self.OUTPUT_TRIGGER = self.OUTPUT_TRIGGER.replace('_0','_1')
                self.alert_num =f'{alert_type}_1'  

        else:
            self.OUTPUT_TRIGGER = os.path.join(self.OUTPUT_PATH,
                                               trigger_id,
                                               alert_type)
            self.alert_num =alert_type
              
        if not os.path.exists(self.OUTPUT_TRIGGER):
            os.makedirs(self.OUTPUT_TRIGGER)
        
        with open(f'{self.OUTPUT_TRIGGER}/{trigger_id}.json', 'w') as jsonfile:
            json.dump(record, jsonfile)
            os.chmod(f'{self.OUTPUT_TRIGGER}/{trigger_id}.json', 0o0777)

    # Adding trigger call 
        
        print('Posting to website...',flush=True)

        trigger_data = {
                        "trigger_label":trigger_id,
                        "mock":is_mock,
                        "season":"1599",
                        "mjd":float(Time(record['event']['time']).mjd),
                        "event_datetime":str(record['time_created']), # record['event']['time']
                        "detectors":str(record['event']['instruments'])[1:-1],
                        "lvc_event_url":record['urls']['gracedb']}

        
        print("Trigger data to be posted to website",flush=True)
        print("",flush=True)
        for key, val in zip(trigger_data.keys(),trigger_data.values()):
            print("Key:",key,flush=True)
            print("Value:",val,flush=True) 
            print("Value datatype:",type(val),flush=True)  
            print("",flush=True)
        
        if self.alert_num=="PRELIMINARY_0":
            self.api.add_trigger(trigger_data)
            self.firstAlert=True
        else:
            self.api.update_trigger(trigger_data)
            self.firstAlert=False

        weather_text = "*Current weather at CTIO*: {} \n*CTIO Weather protocol*: {} \n*CTIO Weather forecast*: {}".format(self.weather_CTIO_currently,self.weather_protocol_CTIO,self.weather_forecast)

        print('Handling Trigger...', flush=True)
        skymap_str = record.get('event', {}).pop('skymap')
        if skymap_str:
            skymap_bytes = b64decode(skymap_str)
            skymap = Table.read(BytesIO(skymap_bytes))
        
        DISTANCE = skymap.meta["DISTMEAN"]
        DISTANCE_SIGMA = skymap.meta["DISTSTD"]
        record['event']['distmean'] = DISTANCE
        record['event']['distsigma'] = DISTANCE_SIGMA

        if record['external_coinc'] != None:
            external_coinc = record['external_coinc']
            self.process_external_coinc(trigger_id=trigger_id,
                                        alert_type=alert_type,
                                        external_coinc=external_coinc)
            record['external_coinc'].pop('combined_skymap', None)

        pprint.pprint(record)
        OUTPUT_SKYMAP = os.path.join(self.OUTPUT_TRIGGER,
                                     f'bayestar.multiorder.fits')
        
        OUTPUT_FLATTEN = OUTPUT_SKYMAP.replace('.multiorder.fits','.fits.gz')

        if not os.path.isfile(OUTPUT_SKYMAP):
            skymap.write(OUTPUT_SKYMAP, overwrite=True)
            os.chmod(OUTPUT_SKYMAP, 0o0777)
            print('Skymap saved at '+OUTPUT_FLATTEN)
        
        self.flatten_skymap(OUTPUT_SKYMAP, OUTPUT_FLATTEN)
        os.chmod(OUTPUT_FLATTEN, 0o0777)

        nt = Time(record['event']['time'])
        trigger_mjd = round(nt.mjd, 4)

        record['MJD'] = str(trigger_mjd)
        record['url'] = record['urls']['gracedb']
        
        
        event_paramfile = os.path.join(self.OUTPUT_TRIGGER,
                                       f"{trigger_id}_params.npz")
        np.savez(event_paramfile, record)
        os.chmod(event_paramfile, 0o0777)
        print('Parameters saved at '+event_paramfile, flush=True)
        FAR = record['event']['far']
        FAR = round(1./float(FAR)/60./60./24./365., 2)
        source, EVENT_PROB = self._get_max_prob(record)

        if source == 'BBH':
            if FAR < self.FAR_threshold or not all(x in record['event']['instruments'] for x in ['H1', 'L1', 'V1']):
                print("BBH event does not pass FAR and detector cut, discarding")

                self.slack_bot2.post_message(subject="", text="BBH event <{}|{}> does not pass FAR and detector cut, discarding \n\n FAR: {} [yr^-1]\n\n FAR threshold: {} [yr^-1]\n\n Detectors: {}".format(record['url'],trigger_id,FAR,self.FAR_threshold,record['event']['instruments']))

                return
        
        if source == 'Terrestrial':
            print("Terrestrial event, discarding")
            self.slack_bot2.post_message(subject="", text="Terrestrial event {}, discarding".format(trigger_id))
            return
            
        info_text = "New GCN received, starting strategy on event: *{}* \n\n*Website page for this event*: {}\n\n :milky_way: *Weather report* :milky_way:\n\n{}\n\n".format(trigger_id,self.website_base_url+trigger_id,weather_text)
 
        area50, area90, maxprob_ra, maxprob_dec, maxprob_dist, maxprob_distsigma, levels, nside, prob = make_alert_skymap(f'{self.OUTPUT_TRIGGER}/bayestar.fits.gz') # the halpix map path

        mass_chirp,chirp_mass_std = chirp_mass(DISTANCE,DISTANCE_SIGMA,area90,area50)
        
        print('Plotting...', flush=True)
        plots_path = Path(os.path.join(self.OUTPUT_TRIGGER, "initial_plots"))
        
        plots_path.mkdir(parents=True, exist_ok=True)
        skymap_plot, moon_plot = make_plots_initial(OUTPUT_FLATTEN, plots_path.as_posix(),str(mass_chirp),DISTANCE,DISTANCE_SIGMA,self.slack_bot)
        server_dir = os.path.join("/des_web","www","html","desgw-new",f"{trigger_id}")
        os.system("ssh -k codemanager@desweb.fnal.gov 'mkdir -p {}'".format(server_dir))
        desweb = f"codemanager@desweb.fnal.gov:{server_dir}"
        os.system(f"rsync -a {plots_path} {desweb}")
        
        subject, text = self._format_message(trigger_id=trigger_id,
                                        record=record,
                                        retraction=False)
        

        
        self.slack_bot.post_message(subject=subject, text=info_text + text)
        self.slack_bot.post_image(skymap_plot,"Skymap - {}".format(trigger_id),"Skymap for {}".format(trigger_id))
        self.slack_bot.post_image(moon_plot,"MoonPlot - {}".format(trigger_id),"MoonPlot for {}".format(trigger_id))
        
        
        self.email_bot = EmailBot(mode=self.mode)
        self.email_bot.send_email(subject=subject,text=text)
        
        OUTPUT_IMAGE = OUTPUT_FLATTEN.replace('bayestar.fits.gz', 'bayestar.png')
        ## Here is where the BBH-AGN observability code should be implmented
        if source=="BBH":
            vis_plot = make_agn_plot(plots_path,trigger_id,mass_chirp,DISTANCE,DISTANCE_SIGMA)
            while not os.path.isfile(vis_plot):
                time.sleep(1)
            self.slack_bot.post_image(vis_plot,"AGN Visibility plot - {}".format(trigger_id),"AGN Visibility plot for {}".format(trigger_id))
        trigger_data = {
                "trigger_label":trigger_id,
                "date":str(record['time_created']),
                "type":str(alert_type),
                # "ligo_prob":float(EVENT_PROB),
                "far":float(FAR),
                "distance":float(DISTANCE),
                # "n_hexes":,
                # "econ_prob":,
                # "econ_area":,
                # "need_area":,
                # "quality":,
                # "exp_time":,
                # "filter":,
                # "hours":,
                # "n_visits":,
                # "n_slots":,
                # "b_slot":,
                "prob_region_50":float(area50),
                "prob_region_90":float(area90),
                # "prob_coverage":,
                # "snr":,
                "chirp_mass":float(mass_chirp),
                # "component_mass_1":,
                # "component_mass_2":,
                # "season":,
                # "prob_vs_slot_plot":,
                # "centered_gif_plot":,
                # "ligo_prob_contour_plot":,
                # "des_prob_vs_ligo_prob_plot":,
                # "des_limit_mag_map":,
                # "des_limit_mag_map_src":,
                # "highest_prob_json":,
                # "low_tt_json":,
                # "log_link":,
                # "strategy_table":,
                "initial_skymap":skymap_plot,
                # "final_skymap":,
                # "airmass":,
                # "cumulative_hex_prob":,
                # "galaxies_plot_initial":,
                # "galaxies_plot_final":,
                # "galaxy_percentage_file":,
                "moon":moon_plot
                }

        
        print("Trigger data to be posted to website",flush=True)
        print("",flush=True)
        for key, val in zip(trigger_data.keys(),trigger_data.values()):
            print("Key:",key,flush=True)
            print("Value:",val,flush=True) 
            print("",flush=True)

        if self.firstAlert:
            self.api.add_trigger_by_day(trigger_data)
        else:
            self.api.update_trigger_by_day(trigger_data)
        
        creationTime = str(record['time_created'])
        print('Passing event to Recycler. Handler took '+elapsedTimeString(t0), flush=True)
        root_dir = os.environ["ROOT_DIR"]
        recycler = 'python ' +\
                    f'{root_dir}/recycler.py ' +\
                    f'--trigger-id {trigger_id} ' +\
                    f'--skymap {self.OUTPUT_TRIGGER}/bayestar.fits.gz ' +\
                    f'--event {source} ' +\
                    f'--official ' +\
                    f'--ltt ' +\
                    f'--alertNum {self.alert_num} ' +\
                    f'--creationTime {creationTime}'
                
### TESTING CHANGES HERE ####### BIG MONEY NO WHAMMIES ### PLEASE CHANGE #####
#        skymap_location = '/data/des80.a/data/eliseke/Main-Injector/new_test_skymaps/o4_hlv_bns/348.fits.gz'
 #       recycler = 'python ' +\
  #                 f'{root_dir}/recycler.py ' +\
   #                f'--trigger-id {trigger_id} ' +\
    #               f'--skymap {skymap_location} ' +\
     #              f'--event {source} ' +\
      #             f'--ltt'


        subprocess.Popen(recycler,shell=True)
