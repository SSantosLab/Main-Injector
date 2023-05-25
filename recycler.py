
import sys
import getopt
import os
import yaml
from subprocess import run, Popen
import OneRing
import send_texts_and_emails as send
import pandas as pd
import datetime
import logging as log
import multiprocessing

FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
log.basicConfig(format=FORMAT)

try:
    args = sys.argv[1:]
    opt, arg = getopt.getopt(
        args, "tp:tid:mjd:exp:sky",
        longopts=["triggerpath=",
                  "triggerid=",
                  "mjd=",
                  "event=",
                  "exposure_length=",
                  "official",
                  "skymapfilename=",
                  "hasrem",
                  "norem"])

except getopt.GetoptError as err:
    print(str(err))
    print("Error : incorrect option or missing argument.")
    print(__doc__)
    sys.exit(1)

# Read in config
with open(
    os.path.join(os.environ["ROOT_DIR"], "recycler.yaml"), "r"
) as f:
    config = yaml.safe_load(f)
# Set defaults to config
trigger_path = config["trigger_path"]

real_or_sim = config["real_or_sim"]

official = False

if config["skymap_filename"] == 'Default':
    skymap_filename = None
else:
    skymap_filename = config["skymap_filename"]

trigger_ids = [config["trigger_id"]]

force_mjd = config["force_mjd"]
resolution = config["resolution"]

#exposure_length = config["exposure_length"]

# Override defaults with command line arguments
# THESE NOT GUARANTEED TO WORK EVER SINCE WE SWITCHED TO YAML
hasrem = False
#    norem = False

dontwrap = False
for o, a in opt:
    print('Option')
    print(o)
    print(a)
    print('-----')
    if o in ["-tp", "--triggerpath"]:
        trigger_path = str(a)
    elif o in ["-tid", "--triggerid"]:
        trigger_ids = [str(a)]
        dontwrap = True
    elif o in ["-mjd", "--mjd"]:
        mjd = float(a)
    # elif o in ["-exp","--exposure_length"]:
    #    exposure_length = float(a)
    elif o in ["-hours", "--hours_available"]:
        hours_available = float(a)
    elif o in ["-sky", "--skymapfilename"]:
        skymap_filename = str(a)
    elif o in ['--hasrem']:
        hasrem = True  # str(a)
        print("HASREM ", hasrem)
    elif o in ['--official']:
        official = True
    elif o in ['--hasrem']:
        hasrem = str(a)
    elif o in ['--norem']:
       hasrem = False
    elif o in ['--event']:
        event = str(a)

    else:
        print("Warning: option", o, "with argument", a, "is not recognized")

# Clear bad triggers, only used for wrapping all triggers...
badtriggers = open('badtriggers.txt', 'w')
badtriggers.close()

trigger_path = trigger_path.rstrip("/")
trigger_id = trigger_ids[0]
log.info(f'Running strategy for {trigger_ids[0]}')

def run_strategy_and_onering(skymap_filename,
                             mjd,
                             trigger_path,
                             trigger_id,
                             sky_condition: str = 'moony',
                             event: str = 'BNS'):

    if event == 'BNS':
        kn_type = 'blue'
    elif event == 'NSBH':
        kn_type = 'red'
    else:
        kn_type = None

    if kn_type is None:
        MSG = 'The current event is either BBH or Terrestrial.' +\
              'We don\'t have a strategy for those kind of events.' +\
              'Hence, not running strategy.'
        log.info(MSG)
        return
    
    mjd = str(mjd).replace('.','')
    current_time = mjd #datetime.datetime.now().strftime('%Y%m%d%H%M')
    cmd = 'python ' +\
        'python/knlc/kn_strategy_sims_MI.py '+\
        f'--input {trigger_path}'+'/'+f'{trigger_id} '+\
        f'--output {trigger_path}'+'/'+f'{trigger_id} '+\
        f'--teff-type {sky_condition} ' +\
        f'--kn-type {kn_type} ' +\
        f'--time {current_time}'
    
    path = os.path.join(f'{trigger_path}/{trigger_ids[0]}')
    output_log = os.path.join(path, f'{sky_condition}_{kn_type}_strategy.log')
    strategy_log = open(output_log, 'w')

    run(cmd,
        shell=True,
        stdout=strategy_log,
        stderr=strategy_log,
        text=True)
    
    log.info(f'Strategy for {trigger_id} done!')
    log.info(f'See log report in {output_log}')
    strategy_file = f'bayestar_{sky_condition}_{kn_type}_{current_time}' +\
                    '_allconfig.csv'
    
    strategy_log.close()
    strategy = os.path.join(path, strategy_file)
    
    df = pd.read_csv(strategy, header=1)
    df.sort_values(by='Deprob1', ascending=False, inplace=True)
    optimal_strategy = df.iloc[0]
    outer, inner, filt, exposure_outer, exposure_inner = optimal_strategy[1:6]
    json_output = os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                               trigger_path,
			                   trigger_id,
                               f"des-gw_{current_time}_{sky_condition}.json"
                               )
    
    OneRing.run_or(
        skymap_filename,
        outer,
        inner,
        filt[0],
        exposure_inner,
        exposure_outer,
        mjd,
        resolution=resolution,
        jsonFilename=json_output
    )
    subject = f'Strategy for event {event}'
    text = f"""\
        Outer region coverage: {outer}
        Inner region coverage: {inner}
        Filter Combination: {filt}
        Exposure for Outer Region: {exposure_outer}
        Exposure for Inner Region: {exposure_inner}
        Json file path: {json_output}
        """
    send.postToSLACK(subject=subject, text=text, official=True, atchannel=True)

moony_strategy = multiprocessing.Process(target=run_strategy_and_onering,
                                        args=(skymap_filename,
                                            mjd,
                                            trigger_path,
                                            trigger_id,
                                            'moony',
                                            event,))

notmoony_strategy = multiprocessing.Process(target=run_strategy_and_onering,
                                            args=(skymap_filename,
                                                mjd,
                                                trigger_path,
                                                trigger_id,
                                                'notmoony',
                                                event,))

moony_strategy.start()
notmoony_strategy.start()
####### BIG MONEY NO WHAMMIES ###############################################
# if config["wrap_all_triggers"]:
#     if not dontwrap:
#         trigger_ids = os.listdir(trigger_path)
#         trigger_ids = trigger_ids[2:]
# for trigger_id in trigger_ids:
#     if force_mjd:
#         mjd = config["mjd"]
#     else:
#         try:
#             mjd = open(os.path.join(trigger_path, trigger_id,
#                         trigger_id + '_eventMJD.txt'), 'r').read()
#         except:
#             mjd = '99999'
#     if skymap_filename is None:
#         try:
#             # if True:
#             # mapname = open(os.path.join(trigger_path,
#             #                            trigger_id,
#             #                            config['default_map_name']), 'r').read()
#             # skymap_filename = os.path.join(trigger_path,
#             #                               trigger_id, config['default_map_name'])
#             # print os.path.join(trigger_path, trigger_id,'default_skymap.txt')
#             # print os.path.join(trigger_path, trigger_id,'default_skymap.txt').read()
#             skymap_filename = os.path.join(trigger_path, trigger_id,
#                                             open(os.path.join(trigger_path, trigger_id,
#                                                                 'default_skymap.txt'), 'r').read())
#         except:
#             badtriggers = open('badtriggers.txt', 'a')
#             badtriggers.write(trigger_id + '\n')
#             print('Could not find skymap url file')

#     if 'bayestar' in skymap_filename:
#         print('bayestar' * 50)

# #        try:
#     if 1 == 1:
#         try:
#             mjd = float(mjd)
#         except:
#             badtriggers = open('badtriggers.txt', 'a')
#             badtriggers.write(trigger_id + '\n')
#             print('WARNING: Could not convert mjd to float. Trigger: ' +
#                     trigger_id + ' flagged as bad.')
# # here is where the object is made, and parts of it are filed in
#         master_dir = os.path.join(trigger_path, trigger_id)
        
#         e = event.Event(skymap_filename,
#                         master_dir,
#                         trigger_id,
#                         mjd,
#                         config,
#                         official,
#                         hasrem)

# # e has variables and code assocaiated with it. The mapMaker is called "e" or "self"

#         e.mapMaker(trigger_id, skymap_filename, config, hasrem) # work end to end
#         # e.getContours(config)  # work
#         #e.makeObservingPlots()  # not working
#         # jsonfilelist = e.makeJSON(config)  # note working
#         # e.make_cumulative_probs()
#         # os.system('cp '+e.event_paramfile+' '+master_dir)
#         # generates the homepage
#         # e.updateTriggerIndex(real_or_sim=real_or_sim)
#         # make a blank page with the basic info that is available
#         # e.updateWebpage(real_or_sim)
        
#         e.send_nonurgent_Email()
        

#        except KeyError:
#            print("Unexpected error:", sys.exc_info())
#            badtriggers = open('badtriggers.txt', 'a')
#            badtriggers.write(trigger_id + '\n')
#############################################################################
