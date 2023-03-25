"""
Python script to listen from LVK alerts trough gcn
"""

from array import array
from datetime import datetime
from threading import Timer
import getopt
import logging
import subprocess
import os
import sys

from typing import Optional, Union, Tuple

import numpy as np
import healpy as hp
import smtplib
#from email.mime.text import MIMEText

#import send_texts_and_emails
# import checkevent_config as config


from gcn_kafka import Consumer
from io import BytesIO
import json
from astropy.table import table
import numpy as np
import pprint

## if name == main part

consumer = Consumer(client_id='44cs6etajmhmbvc2k8opnlj0e7',
                    client_secret='1e7rbf0ore0lfn8446c6f0cmblq51n4povhuq1t3nrbcqiunq4sq')

consumer.subscribe(['igwn.gwalert'])

while True:
    for message in consumer.consume(timeout=3):
        value = message.value()
        get_skymap(value)

sys.path.append(os.environ("GCN_DIR"))

global official

official = False
atchannel = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
logger.addHandler(logging.FileHandler('checkevent.log', 'a'))
print_ = logger.info


def send_first_trigger_email(trigger_id: int,
                             far: float,
                             mapname: str = 'NA',
                             retraction: Union[str, int] = 0,
                             event_params: Optional[float] = None) -> None:
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
    if config.mode.lower() == 'observation':
        plus = 'REAL'
    if config.mode.lower() == 'mdc':
        plus = 'FAKE'
    # if retraction == False:
    print(f"AG TEST: print retraction: {str(retraction)}")
    if retraction == 'Retraction':
        subject = f"Rectraction for {str(trigger_id)}"
        text = ''
        send_texts_and_emails.postToSLACK(
            subject, text, official=official, atchannel=False)
        print("Retraction notice sent, exiting")
        sys.exit()

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
    if config.mode.lower() == 'observation':
        plus = 'REAL'
    if config.mode.lower() == 'MDC':
        plus = 'FAKE'

    text = message

    #    you.extend(t)

    # for y in you:
    subject = plus+' Trigger '+trigger_id+' FAILED!'

    # global official
    send_texts_and_emails.send(
        subject, text, official=official, atchannel=True)

    #    msg['From'] = me
    #    msg['To'] = y

    #    s = smtplib.SMTP('localhost')
    #    s.sendmail(me, y, msg.as_string())
    #    s.quit()


    
def get_skymap(skymap_url: str, skymap_path: str):
    """
    Look up URL of sky map in VOEvent XML document,
    download sky map, and parse FITS file.

    Parameters:
    -----------
    skymap_url: str
        URL of skymap in VOEvent XML format.
    skymap_path: str
        payh of skymap destination local

    Returns:
    --------
    Tuple[array, astropy.fits.header, str]
      Returns a Tuple containing an array, header and a skymap's name.
    """
    mapname = skymap_url.split('/')[-1]
    skymap_filename = os.path.join(skymap_path, mapname)
    # try:
    if True:
        print(f"wget {skymap_url} -O {skymap_filename}")
        os.system(f"wget {skymap_url} -O {skymap_filename}")
        # os.system('wget  --auth-no-challenge '+skymap_url+' -O '+skymap_filename)
    # except:
    #    print 'excepted 1'
    #    sendFailedEmail(
    # skymap_url.split('/')[-1].split('.')[0],
    # message=f"""FAILED TO OBTAIN SKYMAP WITH WGET
    #            wget --auth-no-challenge {skymap_url} -O {skymap_filename}
    #         """
    #    raise ValueError
    try:
        skymap, header = hp.read_map(skymap_filename, h=True, verbose=False)
    except:
        print_(f"failed to read skymap {skymap_filename}")
        sys.exit()
    header = dict(header)

    a = open(os.path.join(skymap_path, 'default_skymap.txt'), 'w')
    a.write(mapname)
    a.close()
    # Done!
    return skymap, header, mapname

# Function to call every time a GCN is received.
# Run only for notices of type LVC_INITIAL or LVC_UPDATE.

def flatten_skymap(skymap: str or Path) -> None:
    """
    Flattens moc skymap
    """

    from subprocess import run
    flatten_name = os.path.basename(skymap)
    if '.fits.gz' in flatten_name:
        flatten_name = flatten_name.rstrip('.fits.gz')
    
    if '.fits' in flatten_name:
        flatten_name = flatten_name.rstrip('.fits')
    
    run('ligo-skymap-flatten s+{skymap} {skymap}_flatten.fits.gz', shell=True)


def process_kafka_gcn(payload: dict) -> None:
    """
    Parsers gcn kafka notice
    """

    print('GOT GCN LIGO EVENT')
    pprint.pprint(payload)


    # for real events
    # if payload['superevent_id'][0] != 'S':
    #     return

    # for mock alerts
    if payload['superevent_id'][0] != 'M':
        return

    if payload['event']['group'] != 'CBC':
        return

    SKYMAP_NAME = payload['superevent_id'][0]
    # Parse sky map
    skymap_str = payload.get('event', {}).pop('skymap')

    if skymap_str:
        skymap_bytes = b64decode(skymap_str)
        skymap = Table.read(BytesIO(skymap_bytes))
        skymap.write(f"{SKYMAP_NAME}.fits")

        # level, ipix = ah.uniq_to_level_ipix(
        #     skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
        # )
        # ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level),
        #                                order='nested')

        DISTANCE = skymap.meta["DISTMEAN"]
        DISTANCE_ERR = skymap.meta["DISTSTD"]

    print('Got LIGO VOEvent!!!!!!!')

    print('payload:')

    payload_without_skymap_string = {
        key: payload[key] for key in paylaod.keys() if key != "skymap"
    }

    pprint(payload_without_skymap_string)

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

    if payload['alert_type'] == 'RETRACTION':
        send_first_trigger_email(trigger_id, 'NA', retraction=alerttype)
        return


    print('Trigger outpath')
    outfolder = os.path.join(config.trigger_outpath, trigger_id)
    print(outfolder)

    

    skymap_filename = os.path.join(outfolder, skymap_url.split('/')[-1])
    # skymap_filename = os.path.join(config.trigger_outpath,trigger_id,skymap_name)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    os.chmod(outfolder, 0o777)
    # Numpy file to save all parameters for future use in webpage.
    event_paramfile = os.path.join(outfolder, f"{trigger_id}_params.npz")
    # Dictionary containing all pertinant parameters
    event_params = {}

    try:

        trigger_mjd_rounded = float(params['Trigger_TJD'])+40000.
        trigger_sod = float(params['Trigger_SOD'])
        seconds_in_day = 86400.
        trigger_mjd = round(trigger_mjd_rounded +
                            trigger_sod/seconds_in_day, 4)
    except:
        # IF NO MJD WE CAN SPECIFY THE MJD OF RECEIVED TRIGGER
        trigger_mjd = 56555.
        from astropy.time import Time
        nt = Time.now()
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
        event_params['ETA'] = str(params['Eta'])
    except:
        event_params['ETA'] = '-999.'
    try:
        event_params['FAR'] = str(
            round(1./float(params['FAR'])/60./60./24./365., 2))+' Years'
    except:
        event_params['FAR'] = '-999.'
    try:
        event_params['ChirpMass'] = str(
            (root.find('./What/Param[@name='"'ChirpMass'"']')
                .attrib['value'])
        )+' '+str(
            root.find('./What/Param[@name="ChirpMass"]').attrib['unit']
        )
    except:
        event_params['ChirpMass'] = '-999.'
    try:
        event_params['MaxDistance'] = str(params['Distance'])
        if float(event_params['MaxDistance'].split()[0]) < 0.:
            event_params['MaxDistance'] = '60.0 Mpc'
    except:
        event_params['MaxDistance'] = '60. Mpc'
    try:
        event_params['boc'] = str(params['Group'])
    except:
        event_params['boc'] = 'Not Available'
    try:
        event_params['CentralFreq'] = str(
            (root.find("./What/Param[@name='CentralFreq']")
                .attrib['value'])
        )+' '+str(
            root.find("./What/Param[@name='CentralFreq']").attrib['unit']
        )
    except:
        event_params['CentralFreq'] = '-999.'
    try:
        event_params['probhasns'] = str(params['BNS'])
    except:
        event_params['probhasns'] = '0.'

    event_params['time_processed'] = '-999'

    # set to defualt
    event_params['integrated_prob'] = '-9.999'
    event_params['M1'] = '-999'
    event_params['M2'] = '-999'
    event_params['nHexes'] = '-999'

    # if config.mode.lower() == 'observation':

    try:
        send_first_trigger_email(trigger_id, event_params['FAR'],
                                 mapname=skymap_url.split(
                                     '/')[-1].split('.')[0],
                                 retraction=alerttype, event_params=event_params)
    except:
        send_first_trigger_email(trigger_id, event_params['FAR'],
                                 retraction=alerttype)

    print(f"Trigger ID HEERERERERE {trigger_id}")
    # save payload to file

    # if not dontwritepayload:
    open(os.path.join(outfolder, trigger_id+'_payload.xml'),
         'w').write(str(payload))
    # save url to file
    open(os.path.join(outfolder, trigger_id+'_skymapURL.txt'),
         'w').write(skymap_url)
    # save mjd to file
    open(os.path.join(outfolder, trigger_id+'_eventMJD.txt'),
         'w').write(str(trigger_mjd))

    # Read sky map
    skymap, header, mapname = get_skymap(skymap_url, outfolder)

    print('saving event paramfile', event_paramfile)

    # Keep this
    np.savez(event_paramfile,
             MJD=event_params['MJD'],
             ETA=event_params['ETA'],
             FAR=event_params['FAR'],
             ChirpMass=event_params['ChirpMass'],
             MaxDistance=event_params['MaxDistance'],
             integrated_prob=event_params['integrated_prob'],
             M1=event_params['M1'],
             M2=event_params['M2'],
             nHexes=event_params['nHexes'],
             CentralFreq=event_params['CentralFreq'],
             time_processed=event_params['time_processed'],
             boc=event_params['boc'],
             mapname=mapname,
             hasremnant=event_params['hasremnant'],
             massgap=event_params['massgap'],
             terrestrial=event_params['terrestrial'],
             alerttype=event_params['alerttype'],
             BBH=event_params['BBH'],
             BNS=event_params['BNS'],
             NSBH=event_params['NSBH']
             )

    # Fire off analysis code
    # if skymap_url.split('/')[-1] == 'bayestar.fits.gz':
    if official:
        args_rem = ['python', 'recycler.py', '--skymapfilename='+skymap_filename, '--triggerpath=' +
                    config.trigger_outpath, '--triggerid='+trigger_id, '--mjd='+str(trigger_mjd), '--official', '--hasrem']
        args_norem = ['python', 'recycler.py', '--skymapfilename='+skymap_filename, '--triggerpath=' +
                      config.trigger_outpath, '--triggerid='+trigger_id, '--mjd='+str(trigger_mjd), '--official', '--norem']

        print_('ARGSSSSSSSSSSSSSSSSSSSSS')
        print_(args)
        ppp = ''
        for argrem, argno in zip(args_rem, args_norem):
            try:
                ppp += argrem+' '+argno+' '
            except:
                pass
    else:
        args = ['python', 'recycler.py', '--skymapfilename='+skymap_filename, '--triggerpath=' +
                config.trigger_outpath, '--triggerid='+trigger_id, '--mjd='+str(trigger_mjd)]
        print_('ARGSSSSSSSSSSSSSSSSSSSSS')

        print_(args)
        ppp = ''
        for arg in args:
            try:
                ppp += arg+' '
            except:
                pass
    print_(ppp)
    # os.mkdir(os.path.join(config.trigger_outpath,trigger_id))
    if official:
        try:
            os.mkdir(os.path.join(config.trigger_outpath, trigger_id +
                     'hasrem/', skymap_filename.split('/')[-1].split('.')[0]))
        except:
            print('path exists')
        try:
            os.mkdir(os.path.join(config.trigger_outpath, trigger_id +
                     'norem/', skymap_filename.split('/')[-1].split('.')[0]))
        except:
            print('path exists')

    f_rem = open(os.path.join(config.trigger_outpath, trigger_id, skymap_filename.split(
        '/')[-1].split('.')[0], skymap_filename.split('/')[-1].split('.')[0]+'_recycler_rem.log'), "w")
    f_norem = open(os.path.join(config.trigger_outpath, trigger_id, skymap_filename.split(
        '/')[-1].split('.')[0], skymap_filename.split('/')[-1].split('.')[0]+'_recycler_norem.log'), "w")

    subprocess.Popen(args_rem, stdout=f_rem, stderr=f_rem)
    subprocess.Popen(args_norem, stdout=f_norem, stderr=f_norem)
    f_rem.close()
    f_norem.close()

    # Need to send an email here saying analysis code was fired

    print('Finished downloading, fired off job')
    # ,skymap_filename.split('/')[-1].split('.')[0],skymap_filename.split('/')[-1].split('.')[0]+'_recycler.log'))
    print('See log here: ' + os.path.join(config.trigger_outpath, trigger_id))


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


def imAliveHTML():
    current_time = datetime.utcnow().strftime("%a, %d %b %Y %H:%M:%S")

    f = open('imalive.html', 'w')
    f.write(
        """
        <!doctype html>
        <html>
        <head>
            <title>Is Main Injector Alive</title>
        </head>
        <body>
            <p>Most recent Main Injector timestamp was: </p>
            <p><break><strong>{current_time} GMT</strong></p>
            <break>If this is >10 minutes older than now:
            <strong><p id="datetime"></p></strong>
            <break>
            <p>Then something could be wrong and please notify Dillon Brout:
            dbrout@physics.upenn.edu and +1 215-300-8763</p>
        <script>
            var dt = new Date(); document.getElementById("datetime").innerHTML = dt.toUTCString();
        </script>
        </body>
        </html>
        """
    )
    f.close()
    Timer(500, imAliveHTML).start()
    #os.system('scp imalive.html codemanager@desweb.fnal.gov:/des_web/www/html/desgw/')

    return


def kinit():
    os.system(
        'kinit -k -t /var/keytab/desgw.keytab desgw/des/des41.fnal.gov@FNAL.GOV'
    )
    Timer(43000, kinit).start()


if __name__ == "__main__":

    
    try:
        args = sys.argv[1:]
        opt, arg = getopt.getopt(
            args, "pl:lm",
            longopts=["payload=", "ligomap=", "official", ])

    except getopt.GetoptError as err:
        print_(str(err))
        print_("Error : incorrect option or missing argument.")
        print_(__doc__)
        sys.exit(1)

    # global official

    official = False
    # try:
    # runnow=False
    for o, a in opt:
        print_(o)
        if o in ['--payload']:
            if '--official' in np.array(opt):
                official = True
            payloadpath = a
            payload = open(payloadpath).readlines()
            import xml.etree.ElementTree as ET
            tree = ET.parse(payloadpath)
            root = tree.getroot()
            process_gcn(payload, root)
            runnow = True
            sys.exit()
        if o in ['--official']:
            official = True
            print_('Official!!!'*10)
            if config.mode.lower() != 'observation':
                print_('This is official, but mode is not observation!!!'*5)
        else:
            print_('Not Official!'*10)

    if config.mode.lower() == 'test':
        import xml.etree.ElementTree as ET
        tree = ET.parse('MS181101ab-1-Preliminary.xml')
        root = tree.getroot()
        payload = open('MS181101ab-1-Preliminary.xml', 'rb').read()
        process_gcn(payload, root)
        sys.exit()
    if config.mode.lower() == 'mdc':
        pass
    elif config.mode.lower() == 'observation':
        pass
    else:
        KeyError(
            """
            checkevent.py Mode must be set to either test or observation.
            Exiting...
            """
        )

    # except:
    # if not runnow:
# Start timer - use threading to say I'm Alive
    print('Started Threading')
    # imAliveEmail()
    imAliveHTML()
    kinit()
# Listen for GCNs until the program is interrupted
# (killed or interrupted with control-C).

    print('Listening...')
    gcn.listen(port=8096, handler=process_gcn)

# IF YOU END UP HERE THEN SEND AN EMAIL AND REBOOT
