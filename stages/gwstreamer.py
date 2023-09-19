import os
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
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.bayestar import rasterize
from .slack import SlackBot
from .communication import EmailBot
import glob
import subprocess
import os
import numpy as np
import numpy as np
import astropy_healpix as ah
import pprint
import json
import yaml

class GWStreamer():

    def __init__(self,
                 mode):

        mi_dir = os.path.dirname(os.path.abspath(__file__))
        root = os.path.dirname(mi_dir)

        self.mode = mode
        self._ROOT = root
        self.OUTPUT_PATH = self._set_output_path()
        self.OUTPUT_TRIGGER = self._set_output_trigger()
        self.SKYMAP_OUTPUT = f'{self.OUTPUT_TRIGGER}/bayestar.fits.gz'
        self.email_bot = EmailBot(mode=mode)
        self.slack_bot = SlackBot(mode=mode)

    def _set_output_path(self):

        if self.mode == 'observation':
            self.OUTPUT_PATH = os.path.join(self._ROOT,
                                            "OUTPUT/O4REAL")
        
        if (self.mode == 'test') or (self.mode == 'mock'):
            self.OUTPUT_PATH = os.path.join(self._ROOT,
                                "OUTPUT/TESTING")

        return self.OUTPUT_PATH
    
    def _set_output_trigger(self):

        alert_type = self._gcn_alert['alert_type']
        trigger_id = self._gcn_alert['superevent_id']

        if alert_type == 'PRELIMINARY':
            self.OUTPUT_TRIGGER = os.path.join(self.OUTPUT_PATH,
                                                trigger_id,
                                                f'{alert_type}_0')
            
        if os.path.isdir(self.OUTPUT_TRIGGER):
            self.OUTPUT_TRIGGER = self.OUTPUT_TRIGGER.replace('_0','_1')

        else:
            self.OUTPUT_TRIGGER = os.path.join(self.OUTPUT_PATH,
                                                trigger_id,
                                                alert_type)

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
            subject = f"Rectraction for {str(trigger_id)}"
            text = ''
            self.email_bot.send_email(subject, text)
            return subject, text


        FAR = record['event']['far']
        FAR = round(1./float(FAR)/60./60./24./365., 2)
        FAR = f'1 per {FAR} Years'

        text = f"Alert Type: {record['alert_type']}\n"+\
            f"Superevent ID: {record['superevent_id']}\n"+\
            f"Event Time: {record['event']['time']} \n"+\
            f"Alert Time: {record['time_created']}\n"+\
            f"MJD: {record['MJD']}\n"+\
            f"FAR: {FAR}\n"+\
            f"Detectors: {record['event']['instruments']}\n"+\
            f"BNS: {record['event']['classification']['BNS']:.3f}\n"+\
            f"NSBH: {record['event']['classification']['NSBH']:.3f}\n"+\
            f"BBH: {record['event']['classification']['BBH']:.3f}\n"+\
            f"Terrestrial: {record['event']['classification']['Terrestrial']}\n"+\
            f"Has Remmant: {record['event']['properties']['HasRemnant']}\n"+\
            f"DISTMEAN: {record['event']['distmean']:.2f} Mpc\n"+\
            f"DISTSIGMA: {record['event']['distsigma']:.2} Mpc\n"+\
            f"Has Mass Gap: {record['event']['properties']['HasMassGap']}\n"+\
            f"GraceDB Link: {record['urls']['gracedb']}\n"+\
            f"Skymap Link: https://gracedb.ligo.org/apiweb/superevents/{trigger_id}/files/bayestar.png"
        
        source, prob_source = self._get_max_prob(record)

        subject = f'Trigger {trigger_id} MJD: {record["MJD"]} '+\
            f'Classification: {source} ({prob_source}) FAR: {FAR} '+\
            f'Alert Type {record["alert_type"]}.'

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

        hdu = fits.open(input_skymap)
        order = ah.nside_to_level(nside)
        table = read_sky_map(hdu, moc=True)
        table = rasterize(table, order=order)
        write_sky_map(output_skymap, table, nest=True)
        
        return
    
    def make_trigger_data(self, record: dict) -> dict:
        """parses the record, returning a dict with
        the appropriate format for the website."""

        columns = ["trigger_label",
                   "type",
                   "ligo_prob",
                   "far",
                   "distance",
                   "mjd",
                   "event_datetime",
                   "mock",
                   "image_url"]
        
        trigger_id = record['superevent_id']
        alert_type, prob = self._get_max_prob(record)
        far = record['event']['far']
        distance = record['event']['distmean']
        mjd = float(Time(record['event']['time']).mjd)
        event_datetime = str(record['event']['time'])
        mock = True
        #image_url = (f'https://des-ops.fnal.gov:8082/desgw/Triggers/'
                     #f'{trigger_id}/bayestar/{trigger_id}_animate.gif')
        image_url = 'https://des-ops.fnal.gov:8082/desgw/Triggers/S200128d/bayestar/S200128d_animate.gif'
        values = [trigger_id,alert_type,prob,far,distance,mjd,
                  event_datetime,mock,image_url]
        trigger_data = {key: value for key, value in zip(columns,values)}

        return trigger_data

    def make_trigger_by_day(self, record: dict) -> dict:
        pass
    
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
            
        combined_skymap_str = external_coinc.get('combined_skymap', {})
        skymap_bytes = b64decode(combined_skymap_str)
        skymap = Table.read(BytesIO(skymap_bytes))
        combined_skymap_output = os.path.join(self.OUTPUT_TRIGGER,
                                              'combined_skymap.multiorder.fits')
        combined_skymap_output_flatten = os.path.join(self.OUTPUT_TRIGGER,
                                                      'combined_skymap.fits.gz')
        
        skymap.write(combined_skymap_output, overwrite=True)        
        self.flatten_skymap(combined_skymap_output,
                            combined_skymap_output_flatten)
        self.flatten_skymap(combined_skymap_output,
                            combined_skymap_output_flatten)
                    
    def handle(self, gcn_alert: dict) -> None:
        """
        Parsers gcn kafka notice
        """

        record = gcn_alert
        alert_type = record['alert_type']
        trigger_id = record['superevent_id']

        if self.mode == 'observation':
            if record['superevent_id'][0] != 'S':
                return
                
        if (self.mode == 'test') or (self.mode == 'mock'):
            if record['superevent_id'][0] != 'M':
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
            return
        
            
        if not os.path.exists(self.OUTPUT_TRIGGER):
            os.makedirs(self.OUTPUT_TRIGGER)


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
            record['external_coinc'].pop('combined_skymap')

        with open(f'{self.OUTPUT_TRIGGER}/{trigger_id}.json', 'w') as jsonfile:
            json.dump(record, jsonfile)

        pprint.pprint(record)
        OUTPUT_SKYMAP = os.path.join(self.OUTPUT_TRIGGER,
                                     f'bayestar.multiorder.fits')
        
        OUTPUT_FLATTEN = OUTPUT_SKYMAP.replace('.multiorder.fits','.fits.gz')

        if not os.path.isfile(OUTPUT_SKYMAP):
            skymap.write(OUTPUT_SKYMAP, overwrite=True)
        
        self.flatten_skymap(OUTPUT_SKYMAP, OUTPUT_FLATTEN)

        nt = Time(record['event']['time'])
        trigger_mjd = round(nt.mjd, 4)

        record['MJD'] = str(trigger_mjd)
        record['url'] = record['urls']['gracedb']
        
        
        event_paramfile = os.path.join(self.OUTPUT_TRIGGER,
                                       f"{trigger_id}_params.npz")
        np.savez(event_paramfile, record)
        FAR = record['event']['far']
        FAR = round(1./float(FAR)/60./60./24./365., 2)
        source, EVENT_PROB = self._get_max_prob(record)

        if FAR < 1000.0:
            if source == 'BBH':
                return
        
        if source == 'Terrestrial':
            return
        
        subject, text = self._format_message(trigger_id=trigger_id,
                                             record=record,
                                             retraction=False)

        #self.email_bot = EmailBot(mode=self.mode)
        #self.email_bot.send_email(subject=subject,text=text)
        #self.slack_bot = SlackBot(mode=self.mode)
        #self.slack_bot.post_message(subject=subject,text=text)
        
        recycler = os.environ["ROOT_DIR"]
        recycler = (f'python {os.path.join(recycler, "recycler.py")} '
                    f'--trigger-id {trigger_id} '
                    f'--skymap {self.OUTPUT_TRIGGER}/bayestar.fits.gz '
                    f'--event {source} '
                    '--official')
        
        subprocess.Popen(recycler,shell=True)
