import os
import urlparse
import matplotlib.pyplot as plt
import numpy as np; import os
import hp2np
import tempfile
import shutil
import sys
import glob
import gcn
import gcn.handlers
import gcn.notice_types
import requests
import healpy as hp
import subprocess
#from astropy.time import Time
import time
import checkevent_config as config

def sendFirstTriggerEmail(trigger_id,far,mapname='NA'):
    import smtplib
    from email.mime.text import MIMEText

    text = 'Trigger '+trigger_id+'\nFAR: '+str(far)+'\nMap: '+mapname+'\nURL: https://gracedb.ligo.org/events/view/'+trigger_id+' \nAnalysis has begun, please hold tight for a DESGW webpage...'
    msg = MIMEText(text)

    me = 'automated-desGW@fnal.gov'
    if config.sendEveryoneEmails:
        you = config.allemails
    else:
        you = ['djbrout@gmail.com']

    if config.sendtexts:
        t = config.alltexts
        you.extend(t)
    else:
        you.extend(['2153008763@mms.att.net'])

    for y in you:
        try:
            msg['Subject'] =  'Trigger '+trigger_id+' FAR: '+str(far)+' Map: '+mapname
            msg['From'] = me
            msg['To'] = y
            
            s = smtplib.SMTP('localhost')
            s.sendmail(me, y, msg.as_string())
            s.quit()
        except:
            print 'could not send alert! '*100

    #if config.sendtexts:
    #    os.system('curl http://textbelt.com/text -d number=2153008763 -d "message=New Trigger FAR:'+str(far)+'"')
    print 'Trigger email sent...'

def sendFailedEmail(trigger_id,message='FAILED'):
    import smtplib
    from email.mime.text import MIMEText

    text = message
    msg = MIMEText(text)

    me = 'automated-desGW@fnal.gov'
    if config.sendEveryoneEmails:
        you = config.allemails
    else:
        you = ['djbrout@gmail.com']

    if config.sendtexts:
        t = config.alltexts
        you.extend(t)

    for y in you:
        msg['Subject'] =  'Trigger '+trigger_id+' FAILED!'
        msg['From'] = me
        msg['To'] = y

        s = smtplib.SMTP('localhost')
        s.sendmail(me, y, msg.as_string())
        s.quit()

def get_skymap(skymap_url,skymap_path):
    """
    Look up URL of sky map in VOEvent XML document,
    download sky map, and parse FITS file.
    """
    mapname = skymap_url.split('/')[-1]
    skymap_filename = os.path.join(skymap_path,mapname)
    #try:
    if True:
        print 'wget'+' --auth-no-challenge '+skymap_url+' -O '+skymap_filename 
        os.system('wget  --auth-no-challenge '+skymap_url+' -O '+skymap_filename)
    #except:
    #    print 'excepted 1'
    #    sendFailedEmail(skymap_url.split('/')[-1].split('.')[0],message='FAILED TO OBTAIN SKYMAP WITH WGET\n'+'wget'+' --auth-no-challenge '+skymap_url+' -O '+skymap_filename)
    #    raise ValueError
    try:
        skymap, header = hp.read_map(skymap_filename, h=True, verbose=False)
    except:
        print 'failed to read skymap '+skymap_filename
        sys.exit()
    header = dict(header)
    
    a = open(os.path.join(skymap_path,'default_skymap.txt'),'w')
    a.write(mapname)
    a.close()
    # Done!
    return skymap, header, mapname

# Function to call every time a GCN is received.
# Run only for notices of type LVC_INITIAL or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)
def process_gcn(payload, root):
    # Print the alert
    #import checkevent_config as config

    # Respond only to 'test' events.
    # VERY IMPORTANT! Replce with the following line of code
    # to respond to only real 'observation' events. DO SO IN CONFIG FILE
    print payload
    print 'GOT GCN LIGO EVENT'

    if root.attrib['role'] != config.mode.lower():
        print 'This event was not of type '+str(config.mode.upper())
        return #This can be changed in the config file
    

    trigger_id = str(root.find("./What/Param[@name='GraceID']").attrib['value'])

    sendFirstTriggerEmail(trigger_id,'NA')


    # Respond only to 'CBC' events. Change 'CBC' to "Burst' to respond to only
    # unmodeled burst events.
    event_type = root.find("./What/Param[@name='Group']").attrib['value']
    #print 'event_type'*50
    #print event_type
    #return
    #if event_type.lower() != 'CBC':
    #    if event_type.lower() != 'Burst':
    #        print 'not cbc or burst'
    #        return
    print('Got LIGO VOEvent!!!!!!!')
    print(payload)
    #print(root.attrib.items())


    # Read out integer notice type (note: not doing anythin with this right now)
    try:
        notice_type = int(root.find("./What/Param[@name='Packet_Type']").attrib['value'])
    except:
        notice_type = 'Not Available'
    #Creating unique folder to save trigger data
    #if event_type == 'CBC':
    skymap_url = root.find(
        "./What/Param[@name='SKYMAP_URL_FITS_BASIC']").attrib['value']

    #    #skymap_name = trigger_id+'_bayestar.fits'
    #if event_type == 'Burst':
    #    skymap_url = 'https://gracedb.ligo.org/events/'+str(trigger_id)+'/files/skyprobcc_cWB_complete.fits'
    #    #skymap_name = trigger_id+'_skyprobcc_cWB_complete.fits'

    print('Trigger outpath')
    outfolder = os.path.join(config.trigger_outpath,trigger_id)
    print(outfolder)

    if 'lalinference' in skymap_url:
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        os.system('touch '+outfolder+'/wehavelal')

    skymap_filename = os.path.join(outfolder,skymap_url.split('/')[-1])
    #skymap_filename = os.path.join(config.trigger_outpath,trigger_id,skymap_name)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    os.chmod(outfolder, 0o777)
    # Numpy file to save all parameters for future use in webpage.
    event_paramfile = os.path.join(outfolder,trigger_id+'_params.npz')
    # Dictionary containing all pertinant parameters
    event_params = {}

    try:
        trigger_mjd_rounded = float(root.find("./What/Param[@name='Trigger_TJD']").attrib['value'])+40000.
        trigger_sod = float(root.find("./What/Param[@name='Trigger_SOD']").attrib['value'])
        seconds_in_day = 86400.
        trigger_mjd = round(trigger_mjd_rounded + trigger_sod/seconds_in_day,4)    
    except:
        #IF NO MJD WE CAN SPECIFY THE MJD OF RECEIVED TRIGGER
        trigger_mjd = 56555.

    event_params['MJD'] = str(trigger_mjd)
    try:
        event_params['ETA'] = str(root.find("./What/Param[@name='Eta']").attrib['value'])
    except:
        event_params['ETA'] = '-999.'
    try:
        event_params['FAR'] = str(round(1./float(root.find("./What/Param[@name='FAR']").attrib['value'])/60./60./24.,1))+' Days'
    except:
        event_params['FAR'] = '-999.'
    try:
        event_params['ChirpMass'] = str(root.find("./What/Param[@name='ChirpMass']").attrib['value'])+' '+str(root.find("./What/Param[@name='ChirpMass']").attrib['unit'])
    except:
        event_params['ChirpMass'] = '-999.'
    try:
        event_params['MaxDistance'] = str(root.find("./What/Param[@name='MaxDistance']").attrib['value'])+' '+str(root.find("./What/Param[@name='MaxDistance']").attrib['unit'])
        if event_params['MaxDistance'] < 0.:
            event_params['MaxDistance'] = '60.0 Mpc'
    except:
        event_params['MaxDistance'] = '60. Mpc'
    try:
        event_params['boc'] = str(root.find("./What/Param[@name='Group']").attrib['value'])
    except:
        event_params['boc'] = 'Not Available'
    try:
        event_params['CentralFreq'] = str(root.find("./What/Param[@name='CentralFreq']").attrib['value'])+' '+str(root.find("./What/Param[@name='CentralFreq']").attrib['unit'])
    except:
        event_params['CentralFreq'] = '-999.'
    try:
        event_params['probhasns'] = str(root.find("./What/Param[@name='ProbHasNS']").attrib['value']) + ' ' + str(
            root.find("./What/Param[@name='ProbHasNS']").attrib['unit'])
    except:
        event_params['probhasns'] = '0.'


    event_params['time_processed'] = '-999'

    #set to defualt
    event_params['integrated_prob'] = '-9.999'
    event_params['M1'] = '-999'
    event_params['M2'] = '-999'
    event_params['nHexes'] = '-999'


    #if config.mode.lower() == 'observation':
    try:
        sendFirstTriggerEmail(trigger_id,event_params['FAR'],mapname=skymap_url.split('/')[-1].split('.')[0])
    except:
        sendFirstTriggerEmail(trigger_id,event_params['FAR'])
    print 'Trigger ID HEERERERERE '+trigger_id
    #save payload to file
    open(os.path.join(outfolder,trigger_id+'_payload.xml'), 'w').write(payload)
    #save url to file
    open(os.path.join(outfolder,trigger_id+'_skymapURL.txt'), 'w').write(skymap_url)
    #save mjd to file
    open(os.path.join(outfolder,trigger_id+'_eventMJD.txt'), 'w').write(str(trigger_mjd))
    
    # Read sky map
    skymap, header, mapname = get_skymap(skymap_url,outfolder)


    np.savez(event_paramfile,
             MJD=event_params['MJD'],
             ETA=event_params['ETA'],
             FAR=event_params['FAR'],
             ChirpMass=event_params['ChirpMass'],
             MaxDistance=event_params['MaxDistance'],
             integrated_prob=event_params['integrated_prob'],
             M1 = event_params['M1'],
             M2 = event_params['M2'],
             nHexes = event_params['nHexes'],
             CentralFreq = event_params['CentralFreq'],
             time_processed = event_params['time_processed'],
             boc = event_params['boc'],
             mapname= mapname
             )

    
    #Fire off analysis code    
    #if skymap_url.split('/')[-1] == 'bayestar.fits.gz':
    args = ['python', 'recycler.py','--skymapfilename='+skymap_filename, '--triggerpath='+config.trigger_outpath, '--triggerid='+trigger_id, '--mjd='+str(trigger_mjd)]    
    print 'ARGSSSSSSSSSSSSSSSSSSSSS'
    print args
    #os.mkdir(os.path.join(config.trigger_outpath,trigger_id))
    os.mkdir(os.path.join(config.trigger_outpath,trigger_id,skymap_filename.split('/')[-1].split('.')[0]))
    f = open(os.path.join(config.trigger_outpath,trigger_id,skymap_filename.split('/')[-1].split('.')[0],skymap_filename.split('/')[-1].split('.')[0]+'recycler.log'), "w")
    subprocess.Popen(args,stdout=f,stderr=f)
    f.close()
    #Need to send an email here saying analysis code was fired
    
    print 'Finished downloading, fired off job, listening again...'

from threading import Timer
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
    print 'Im alive email sent...'
    Timer(43200,imAliveEmail).start()

    return


if config.mode.lower() == 'test':
    pass
elif config.mode.lower() == 'observation':
    pass
else:
    KeyError('checkevent_config.py Mode must be set to either test or observation.\nExiting...')

import logging
# Set up logger
logging.basicConfig(level=logging.INFO)

#Start timer - use threading to say I'm Alive
print 'Started Threading'
imAliveEmail()

# Listen for GCNs until the program is interrupted
# (killed or interrupted with control-C).
print 'Listening...'
gcn.listen(host='68.169.57.253', port=8096, handler=process_gcn)
#gcn.listen(host='68.169.57.253', port=8096, handler=process_gcn,im_alive_filename='/data/des41.a/data/desgw/maininjector/imalivetest.txt')

#IF YOU END UP HERE THEN SEND AN EMAIL AND REBOOT
