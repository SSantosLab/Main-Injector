"""
Python script to listen from LVK alerts trough gcn
"""
from base64 import b64decode
from datetime import datetime
from threading import Timer
from typing import Optional, Union
from gcn_kafka import Consumer
from io import BytesIO
from astropy.table import Table
from pathlib import Path
from argparse import ArgumentParser
from subprocess import run, STDOUT, PIPE
from astropy.time import Time
from astropy.io import fits
from ligo.skymap.io.fits import read_sky_map, write_sky_map
from ligo.skymap.bayestar import rasterize

import logging as log
import subprocess
import os
import sys
import numpy as np
import healpy as hp
import smtplib
import numpy as np
import astropy_healpix as ah
import pprint
import send_texts_and_emails
import json

FORMAT = '%(asctime)s %(message)s'
log.basicConfig(format=FORMAT)

def send_first_trigger_email(trigger_id: int,
                             event_params: dict,
                             retraction: bool = False,
                             mode: str = 'test') -> None:
    """
    Send an email / SLACK message with trigger alert from LVK.

    Parameters:
    -----------

    trigger_id: int
        Trigger identifier to a env reported by LVK
    far: float
        Distance to a event from Earth.
    mapname: str (Default: 'NA')
        Map name of the event.
    retraction: bool (Default: False)
        A Retraction content to know if an alert is a real
        astrophysical event or not.
    event_params: float (Optional, Default: None)

    Returns:
    --------

    None
    """

    
    if mode == 'observation':
        plus = 'REAL'
        official=True

    if mode == 'test':
        plus = 'FAKE'
        official=False
    
    if retraction:
        subject = f"Rectraction for {str(trigger_id)}"
        text = ''
        send_texts_and_emails.postToSLACK(
            subject, text, official=official, atchannel=False)
        send_texts_and_emails.send(subject, text, official=official)
        print("Retraction notice sent, returning to listening")
        return


    plus = ''
    far = event_params['FAR']
    classfication_scores = [
        ('BBH', event_params['BBH']),
        ('BNS', event_params['BNS']),
        ('NSBH', event_params['NSBH']),
        ('terrestrial', event_params['terrestrial'])
    ]

    max_prob = np.argmax(np.array([p[1] for p in classfication_scores]))
    EVENT_KIND = classfication_scores[max_prob][0]
    EVENT_PROB = classfication_scores[max_prob][1]

    text = f"""\
        Trigger {trigger_id}
        HasRemnant: {event_params['hasremnant']}
        Alert Type: {event_params['alerttype']}
        FAR: {far}
        URL: {event_params['url']}
        Classification: {EVENT_KIND}: {EVENT_PROB}
        MJD: {event_params['MJD']}
        Group: {event_params['boc']}
        DISTMEAN: {event_params['DISTMEAN']:.2f} Mpc
        DISTSIGMA: {event_params['DISTSIGMA']:.2} Mpc
        """
    
    # if external_coinc:
    #     text += f"""\
    #         This alert have external coincidence!
    #         """
    #     msg = ''
    #     for key in external_coinc.keys():
    #         msg = f'{key}: {external_coinc[key]}\n'

    #     text += msg


    subject = f'{plus} Trigger {trigger_id} FAR: {far}. '

    send_texts_and_emails.postToSLACK(
        subject,
        text,
        official=official,
        atchannel=True
    )

    log.info('Trigger email sent...')
    send_texts_and_emails.send(subject=subject,text=text,official=official)

def sendFailedEmail(trigger_id: int, message: str = 'FAILED') -> None:
    """
    Send an email / SLACK message with trigger alert from LVK.
    By default, send with 'FAILED' message.

    Parameters:
    -----------

    trigger_id: int
        Trigger identifier to a env reported by LVK
    message: stor (Default: 'FAILED')

    Returns:
    --------

    None
    """

    plus = ''
    if mode == 'observation':
        plus = 'REAL'
    if mode == 'test':
        plus = 'FAKE'

    text = message

    subject = f"{plus} Trigger {trigger_id} FAILED!"

    send_texts_and_emails.send(subject, text, atchannel=True)


def flatten_skymap(input: str, 
                   output: str, 
                   nside: int = 512, 
                   overwrite: bool = True) -> None:
    """Flattens skymap"""

    hdu = fits.open(input)
    order = ah.nside_to_level(nside)
    table = read_sky_map(hdu, moc=True)
    table = rasterize(table, order=order)
    write_sky_map(output, table, nest=True)

def process_kafka_gcn(payload: dict, mode: str = 'test') -> None:
    """
    Parsers gcn kafka notice
    """

    log.info('GOT GCN LIGO EVENT')
    payload = json.loads(payload)
    trigger_id = payload['superevent_id']

    if mode == 'test':
        if payload['superevent_id'][0] != 'M':
            return
        
        OUTPUT_PATH = "OUTPUT/TESTING"

    elif mode == 'observation':
        if payload['superevent_id'][0] != 'S':
            return
    
        OUTPUT_PATH = "OUTPUT/04REAL"

    if payload['alert_type'] == 'RETRACTION':
        event_params = None
        log.info(payload['superevent_id'], 'was retracted')
        send_first_trigger_email(trigger_id=trigger_id,
                                 event_params=event_params,
                                 retraction=True,
                                 mode=mode)
        return    

    if payload['event']['group'] != 'CBC':
        return


    
    OUTPUT_TRIGGER = os.path.join(OUTPUT_PATH, trigger_id)
    OUTPUT_MOC = os.path.join(OUTPUT_TRIGGER, 'moc')
    OUTPUT_COMBINED = os.path.join(OUTPUT_TRIGGER, 'combined_skymap')

    if not os.path.exists(OUTPUT_TRIGGER):
        os.makedirs(OUTPUT_TRIGGER)

    if not os.path.exists(OUTPUT_MOC):
        os.makedirs(OUTPUT_MOC)

    if not os.path.exists(OUTPUT_COMBINED):
        os.makedirs(OUTPUT_COMBINED)

    skymap_str = payload.get('event', {}).pop('skymap')

    # if payload.get('external_coinc') is not None:

    #     try:
    #         combined_skymap = payload.get('external_coinc', {}).pop('combined_skymap')
    #         # external_coinc = payload.get('external_coinc', {})

    #         skymap_bytes = b64decode(combined_skymap)
    #         skymap = Table.read(BytesIO(skymap_bytes))
    #         OUTPUT_COMBINED_SKYMAP = os.path.join(OUTPUT_COMBINED,
    #                                  'bayestar_combined_moc.fits.gz',)
        
    #         if not os.path.isfile(OUTPUT_COMBINED_SKYMAP):
    #             skymap.write(OUTPUT_COMBINED_SKYMAP, overwrite=True)
    #             flatten_skymap(OUTPUT_COMBINED_SKYMAP,
    #                        f'{OUTPUT_COMBINED}/bayestar_combined.fits.gz')
    #     except:
    #         external_coinc = {}
        

    if skymap_str:
        skymap_bytes = b64decode(skymap_str)
        skymap = Table.read(BytesIO(skymap_bytes))
        OUTPUT_SKYMAP = os.path.join(OUTPUT_MOC,
                                     'bayestar_moc.fits.gz',)
        
        if not os.path.isfile(OUTPUT_SKYMAP):
            skymap.write(OUTPUT_SKYMAP, overwrite=True)
            flatten_skymap(OUTPUT_SKYMAP, f'{OUTPUT_TRIGGER}/bayestar.fits.gz')
        
    DISTANCE = skymap.meta["DISTMEAN"]
    DISTANCE_SIGMA = skymap.meta["DISTSTD"]

    
    pprint.pprint(payload)

    try:
        alerttype = payload['alert_type']
    except:
        alerttype = 'N/A'
    try:
        hasremnant = payload['event']['properties']['HasRemnant']
    except:
        hasremnant = -9
    try:
        terrestrial = payload['event']['classification']['Terrestrial']
    except:
        terrestrial = -9
    try:
        massgap = payload['event']['properties']['HasMassGap']
    except:
        massgap = -9
    try:
        BNS = payload['event']['classification']['BNS']
        BBH = payload['event']['classification']['BBH']
        NSBH =payload['event']['classification']['NSBH']
    except:
        BNS = -9
        BBH = -9
        NSBH = -9


    log.info(f'Trigger outpath: {OUTPUT_TRIGGER}')

    event_paramfile = os.path.join(OUTPUT_TRIGGER, f"{trigger_id}_params.npz")
    event_params = {}

    
    nt = Time(payload['event']['time'])
    trigger_mjd = round(nt.mjd, 4)

    event_params['MJD'] = str(trigger_mjd)
    event_params['alerttype'] = alerttype
    event_params['terrestrial'] = terrestrial
    event_params['hasremnant'] = hasremnant
    event_params['massgap'] = massgap
    event_params['BBH'] = BBH
    event_params['BNS'] = BNS
    event_params['NSBH'] = NSBH
    event_params['url'] = payload['urls']['gracedb']
    try:
        FAR = payload['event']['far']
        event_params['FAR'] = f'{round(1./float(FAR)/60./60./24./365., 2)} Years'
    except:
        event_params['FAR'] = '-999.'

    try:
        event_params['probhasns'] = payload['event']['classfication']['BNS']
    except:
        event_params['probhasns'] = '0.'
    try:
        event_params['CentralFreq'] = payload['event']['central_frequency']
    except:
        event_params['CentralFreq'] = '-999.'
    try:
        event_params['boc'] = payload['event']['group']
    except:
        event_params['boc'] = 'Not Available'
    try:
        event_params['DISTMEAN'] = DISTANCE
    except:
        event_params['DISTMEAN'] = '60.0 Mpc'
    try:
        event_params['DISTSIGMA'] = DISTANCE_SIGMA
    except:
        event_params['DISTSIGMA'] = '0'

    event_params['ETA'] = '-999.'  # WHAT ETA Means?
    event_params['ChirpMass'] = '-999.'
    event_params['time_processed'] = '-999'
    event_params['integrated_prob'] = '-9.999'
    event_params['M1'] = '-999'
    event_params['M2'] = '-999'
    event_params['nHexes'] = '-999'

    send_first_trigger_email(trigger_id=trigger_id,
                             event_params=event_params,
                             retraction=False,
                             mode=mode)

    

    
    log.info(f"Trigger ID: {trigger_id}")
    log.info('saving event paramfile:', event_paramfile)

    np.savez(event_paramfile,
             MJD=event_params['MJD'],
             ETA=event_params['ETA'],
             FAR=event_params['FAR'],
             ChirpMass=event_params['ChirpMass'],
             DSSTMEAN=event_params['DISTMEAN'],
             DISTSIGMA=event_params['DISTSIGMA'],
             integrated_prob=event_params['integrated_prob'],
             M1=event_params['M1'],
             M2=event_params['M2'],
             nHexes=event_params['nHexes'],
             CentralFreq=event_params['CentralFreq'],
             time_processed=event_params['time_processed'],
             boc=event_params['boc'],
             mapname='bayestar.fits.gz',
             hasremnant=event_params['hasremnant'],
             massgap=event_params['massgap'],
             terrestrial=event_params['terrestrial'],
             alerttype=event_params['alerttype'],
             BBH=event_params['BBH'],
             BNS=event_params['BNS'],
             NSBH=event_params['NSBH']
             )


    args_rem = ['python',
                'recycler.py', 
                f'--skymapfilename={OUTPUT_TRIGGER}/bayestar.fits.gz',
                f'--triggerpath={OUTPUT_PATH}',
                f'--triggerid={trigger_id}',
                f'--mjd='+str(trigger_mjd),
                '--official',
                '--hasrem']
    
    args_norem = ['python',
                'recycler.py', 
                f'--skymapfilename={OUTPUT_TRIGGER}/bayestar.fits.gz',
                f'--triggerpath={OUTPUT_PATH}',
                f'--triggerid={trigger_id}',
                f'--mjd='+str(trigger_mjd),
                '--official',
                '--norem']

    try:
        os.mkdir(os.path.join(OUTPUT_TRIGGER, 'hasrem', 'bayestar'))
    except:
        pass
    try:
        os.mkdir(os.path.join(OUTPUT_TRIGGER, 'norem', 'bayestar'))
    except:
        pass

    hasrem_log = open(f'{OUTPUT_TRIGGER}/recycler_rem.log', 'w')
    norem_log = open(f'{OUTPUT_TRIGGER}/recycler_norem.log', 'w')

    hasrem = subprocess.Popen(args_rem,
                              stdout=hasrem_log,
                              stderr=hasrem_log,
                              text=True)
    
    norem = subprocess.Popen(args_norem,
                             stdout=norem_log,
                             stderr=norem_log,
                             text=True)

    hasrem_log.close()
    norem_log.close()
    # Need to send an email here saying analysis code was fired
    log.info('fired off job!')
    log.info(f'See log here: {OUTPUT_TRIGGER}')

def imAliveEmail():
    pass

    # import smtplib
    # from email.mime.text import MIMEText

    # text = 'MainInjector Is Alive'
    # msg = MIMEText(text)

    # me = 'imAlive-desGW@fnal.gov'
    # if config.sendEveryoneEmails:
    #     you = config.allemails
    # else:
    #     you = ['djbrout@gmail.com']

    # for y in you:
    #     msg['Subject'] = 'MainInjector is Alive'
    #     msg['From'] = me
    #     msg['To'] = y

    #     s = smtplib.SMTP('localhost')
    #     s.sendmail(me, y, msg.as_string())
    #     s.quit()
    # print('Im alive email sent...')
    # Timer(43200, imAliveEmail).start()

def kinit():
    os.system(
        'kinit -k -t /var/keytab/desgw.keytab desgw/des/des41.fnal.gov@FNAL.GOV'
    )
    Timer(43000, kinit).start()


if __name__ == "__main__":


    parser = ArgumentParser()
    parser.add_argument('--mode',
                        '-m',
                        default='test')
    parser.add_argument('--offline',
                        action='store_true',
                        help='Activates offline testing.')
    
    args = parser.parse_args()

    mode = args.mode
    offline = args.offline

    if offline:
        with open('MS181101ab-preliminary.json') as f:
            record = f.read()

        process_kafka_gcn(record)

    if mode == 'test' and not offline:
        consumer = Consumer(client_id='44cs6etajmhmbvc2k8opnlj0e7',
                            client_secret='1e7rbf0ore0lfn8446c6f0cmblq51n4povhuq1t3nrbcqiunq4sq')
        consumer.subscribe(['igwn.gwalert'])

        while True:
            for message in consumer.consume(timeout=1):
                value = message.value()
                process_kafka_gcn(value)

    try:
        if mode == 'observation':
            consumer = Consumer(client_id='44cs6etajmhmbvc2k8opnlj0e7',
                                client_secret='1e7rbf0ore0lfn8446c6f0cmblq51n4povhuq1t3nrbcqiunq4sq')
            consumer.subscribe(['igwn.gwalert'])

            while True:
                for message in consumer.consume(timeout=1):
                    value = message.value()
                    process_kafka_gcn(value,mode=mode)
    except:
        send_texts_and_emails.send_emergencial_email()

