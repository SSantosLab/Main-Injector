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

import logging
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

def send_first_trigger_email(trigger_id: int,
                             far: float,
                             mapname: str = 'NA',
                             retraction: Union[str, int] = 0,
                             event_params: Optional[float] = None,
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
    retraction: int (Default: 0)
        A Retraction content to know if an alert is a real
        astrophysical event or not.
    event_params: float (Optional, Default: None)

    Returns:
    --------

    None
    """

    plus = ''
    if mode == 'observation':
        plus = 'REAL'
    if mode == 'test':
        plus = 'FAKE'

    if retraction == False:
        print(f"AG TEST: print retraction: {str(retraction)}")
    
    # if retraction == 'Retraction':
    #     subject = f"Rectraction for {str(trigger_id)}"
    #     text = ''
    #     send_texts_and_emails.postToSLACK(
    #         subject, text, official=official, atchannel=False)
    #     print("Retraction notice sent, exiting")
    #     sys.exit()

    if event_params is None:
        text = f"""\
            Trigger: {trigger_id}
            Alert Type: {retraction}
            FAR: {str(far)}
            Map: {mapname}
            URL: https://gracedb.ligo.org/superevents/{trigger_id}
            view analysis has begun, please hold tight for a DESGW webpage
            which will be located here shortly:
            http://des-ops.fnal.gov:8080/desgw/Triggers/{trigger_id}/{trigger_id}_trigger.html
            DO NOT REPLY TO THIS THREAD, NOT ALL USERS WILL SEE YOUR RESPONSE.
            """
    else:
        text = f"""\
            Trigger {trigger_id}
            HasRemnant: {str(event_params['hasremnant'])}
            Alert Type: '{str(retraction)}
            FAR: {str(far)}
            Map: {mapname}
            URL: https://gracedb.ligo.org/superevents/{trigger_id}
            view Analysis has begun, please hold tight for a DESGW webpage
            which will be located here shortly:
            http://des-ops.fnal.gov:8080/desgw/Triggers/trigger_id/{trigger_id}_trigger.html
            DO NOT REPLY TO THIS THREAD, NOT ALL USERS WILL SEE YOUR RESPONSE.
            """

    subject = plus+' Trigger '+trigger_id + \
        ' FAR: '+str(far)+' Map: '+mapname+' NOREPLY'

    send_texts_and_emails.postToSLACK(
        subject, text, official=official, atchannel=True)

    print_('Trigger email sent...')

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


def flatten_skymap(input: str, output: str, nside: int = 512) -> None:
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

    print('GOT GCN LIGO EVENT')
    payload = json.loads(payload)
    
    if mode == 'test':
        if payload['superevent_id'][0] != 'M':
            return
        
        OUTPUT_PATH = "OUTPUT/TESTING"

    elif mode == 'observation':
        if payload['superevent_id'][0] != 'S':
            return
    
        OUTPUT_PATH = "OUTPUT/04REAL"

    if payload['alert_type'] == 'RETRACTION':
        print(payload['superevent_id'], 'was retracted')
        return

    if payload['event']['group'] != 'CBC':
        return


    trigger_id = payload['superevent_id']
    os.makedirs(os.path.join(OUTPUT_PATH, trigger_id))

    skymap_str = payload.get('event', {}).pop('skymap')

    if skymap_str:
        skymap_bytes = b64decode(skymap_str)
        skymap = Table.read(BytesIO(skymap_bytes))
        OUTPUT_SKYMAP = os.path.join(OUTPUT_PATH,
                                     trigger_id,
                                     'bayestar_moc.fits')
        
        skymap.write(OUTPUT_SKYMAP)
        flatten_skymap(OUTPUT_SKYMAP, f'{OUTPUT_PATH}/{trigger_id}/bayestar.fits.gz')
        
    DISTANCE = skymap.meta["DISTMEAN"]
    DISTANCE_ERR = skymap.meta["DISTSTD"]

    
    pprint.pprint(payload)

    try:
        alerttype = payload['alert_type']
    except:
        alerttype = 'N/A'
    try:
        hasremnant = payload['properties']['HasRemnant']
    except:
        hasremnant = -9
    try:
        terrestrial = payload['classification']['Terrestrial']
    except:
        terrestrial = -9
    try:
        massgap = payload['properties']['HasMassGap']
    except:
        massgap = -9
    try:
        BNS = payload['classfication']['BNS']
        BBH = payload['classfication']['BBH']
        NSBH =payload['classfication']['NSBH']
    except:
        BNS = -9
        BBH = -9
        NSBH = -9


    print('Trigger outpath')
    outfolder = os.path.join(os.path.abspath(OUTPUT_PATH), trigger_id)
    print(outfolder)
    
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    os.chmod(outfolder, 777)

    event_paramfile = os.path.join(outfolder, f"{trigger_id}_params.npz")
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

    try:
        event_params['FAR'] = str(
            round(1./float(payload['event']['FAR'])/60./60./24./365., 2))+' Years'
    except:
        event_params['FAR'] = '-999.'

    try:
        event_params['probhasns'] = payload['classfication']['BNS']
    except:
        event_params['probhasns'] = '0.'

    event_params['time_processed'] = '-999'
    event_params['ChirpMass'] = '-999.'
    event_params['CentralFreq'] = '-999.'
    event_params['integrated_prob'] = '-9.999'
    event_params['M1'] = '-999'
    event_params['M2'] = '-999'
    event_params['nHexes'] = '-999'

    # try:
    #     send_first_trigger_email(trigger_id,
    #                              event_params['FAR'],
    #                              mapname=trigger_id,
    #                              retraction=alerttype,
    #                              event_params=event_params)
    # except:
    #     send_first_trigger_email(trigger_id,
    #                              event_params['FAR'],
    #                              retraction=alerttype)

    print(f"Trigger ID: {trigger_id}")

    # Read sky map
    print('saving event paramfile', event_paramfile)

    # Keep this
    np.savez(event_paramfile,
             MJD=event_params['MJD'],
             #ETA=event_params['ETA'],
             FAR=event_params['FAR'],
             ChirpMass=event_params['ChirpMass'],
             MaxDistance=event_params['MaxDistance'],
             integrated_prob=event_params['integrated_prob'],
             M1=event_params['M1'],
             M2=event_params['M2'],
             nHexes=event_params['nHexes'],
             CentralFreq=event_params['CentralFreq'],
             time_processed=event_params['time_processed'],
             #boc=event_params['boc'],
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
                f'--skymapfilename={OUTPUT_PATH}/{trigger_id}/bayestar.fitz.gz',
                f'--triggerpath={OUTPUT_PATH}',
                f'--triggerid={trigger_id}',
                f'--mjd='+str(trigger_mjd),
                '--official',
                '--hasrem']
    
    args_norem = ['python',
                'recycler.py', 
                f'--skymapfilename={OUTPUT_PATH}/{trigger_id}/bayestar.fitz.gz',
                f'--triggerpath={OUTPUT_PATH}',
                f'--triggerid={trigger_id}',
                f'--mjd='+str(trigger_mjd),
                '--official',
                '--norem']

    try:
        os.mkdir(os.path.join(OUTPUT_PATH, trigger_id, 'hasrem/', 'bayestar'))
    except:
        print('path exists')
    try:
        os.mkdir(os.path.join(OUTPUT_PATH, trigger_id, 'norem', 'bayestar'))
    except:
        print('path exists')

    hasrem = subprocess.Popen(args_rem, stdout=PIPE, stderr=STDOUT, text=True)
    norem = subprocess.Popen(args_norem, stdout=PIPE, stderr=STDOUT, text=True)

    hasrem_output = os.path.join(OUTPUT_PATH,trigger_id,'hasrem')
    norem_output = os.path.join(OUTPUT_PATH,trigger_id,'norem')
    
    with open(f'{hasrem_output}/hasrem.log') as f:
        f.write(hasrem.stdout)
    
    with open(f'{norem_output}/norem.log') as f:
        f.write(norem.stdout)

    # Need to send an email here saying analysis code was fired


def imAliveEmail():

    import smtplib
    from email.mime.text import MIMEText

    text = 'MainInjector Is Alive'
    msg = MIMEText(text)

    me = 'imAlive-desGW@fnal.gov'
    if config.sendEveryoneEmails:
        you = config.allemails
    else:
        you = ['djbrout@gmail.com']

    for y in you:
        msg['Subject'] = 'MainInjector is Alive'
        msg['From'] = me
        msg['To'] = y

        s = smtplib.SMTP('localhost')
        s.sendmail(me, y, msg.as_string())
        s.quit()
    print('Im alive email sent...')
    Timer(43200, imAliveEmail).start()

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
    
    args = parser.parse_args()

    mode = args.mode

    if mode == 'test':
        consumer = Consumer(client_id='44cs6etajmhmbvc2k8opnlj0e7',
                            client_secret='1e7rbf0ore0lfn8446c6f0cmblq51n4povhuq1t3nrbcqiunq4sq')
        consumer.subscribe(['igwn.gwalert'])

        while True:
            for message in consumer.consume(timeout=3):
                value = message.value()
                process_kafka_gcn(value)
