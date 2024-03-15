import os
import numpy as np
import OneRing
import pandas as pd
import multiprocessing
from subprocess import run
from argparse import ArgumentParser
from astropy.io import fits
from handlers.slack import SlackBot
import time

parser = ArgumentParser()
parser.add_argument('--trigger-id',
                    type=str,
                    help='Superevent ID. example: S230615az')
parser.add_argument('--skymap',
                    type=str,
                    help='Path to superevent\'s skymap.')
parser.add_argument('--event',
                    type=str,
                    help='Source of GW Sginal. Can be BNS, NSBH or BBH.')
parser.add_argument('--max-hex-count',
                    type=float,
                    default=None,
                    nargs='?',
                    help='Limit the number of hexes in json file. Default is None')
parser.add_argument('--max-hex-time',
                    type=float,
                    default=None,
                    nargs='?',
                    help='Limit the number of hexes in json file based on time in seconds. Default is None.')

parser.add_argument('--official',
                    action='store_true',
                    default=False,
                    help='If official, starts recycler for an official event.')

#NEW ADDITION -- LEAST TELESCOPE TIME STRATEGY?#
parser.add_argument('--ltt',
                    action='store_true',
                    default=False,
                    help='If true, uses least telescope time strategy.')

args = parser.parse_args()
timer_start = time.perf_counter()
print('Settings for Recycler:')

arguments = vars(args)
for key,value in arguments.items():
    print('---------------')
    print(f'{key}: {value}')
    print()
    print('---------------')

trigger_id = args.trigger_id
skymap = args.skymap
event = args.event
max_hex_time = args.max_hex_time
max_hex_count = args.max_hex_count
official = args.official
least_telescope = args.ltt

with fits.open(skymap) as f:
    header = f[1].header
    if header['ORDERING'] == 'NUNIQ':
        skymap_flatten_path = skymap.replace(skymap.split('/')[-1], skymap.split('/')[-1].partition('.')[0] + "_flatten.fits.gz")
        os.system('ligo-skymap-flatten --nside 1024 {} {}'.format(skymap, skymap_flatten_path))
        with fits.open(skymap_flatten_path) as f_flat:
            header = f_flat[1].header
        skymap = skymap_flatten_path 

mjd = header['MJD-OBS']

def run_strategy_and_onering(skymap_filename,
                             trigger_id,
                             mjd,
                             sky_condition: str = 'moony',
                             event: str = 'BNS'):
    if event != 'BBH':
        if event == 'BNS':
            kn_type = 'blue'

        if event == 'NSBH':
            kn_type = 'red'
            
        #ELISE TOOK OUT THIS LINE. DECIMAL IS NEEDED FOR DEPENDENT CODE.
        outname = str(mjd).replace('.','')
        print(f'mjd in runstrategy code: {mjd}')
        current_time = mjd

        input_dir = os.path.dirname(skymap)
        output_dir = input_dir
        cmd = 'python ' +\
            'python/knlc/kn_strategy_sims_MI.py '+\
            f'--input-skymap {skymap} '+\
            f'--output {output_dir} '+\
            f'--teff-type {sky_condition} ' +\
            f'--kn-type {kn_type} ' +\
            f'--time {outname}'
            
        output_log = os.path.join(output_dir,
                                f'{sky_condition}_{kn_type}_strategy.log')
        strategy_log = open(output_log, 'w')

        print('strategy started!')
        run(cmd,
            shell=True,
            stdout=strategy_log,
            stderr=strategy_log,
            text=True)
            
        strategy_file = f'bayestar_{sky_condition}_{kn_type}_{outname}' +\
                        '_allconfig.csv'
            
        strategy_log.close()
        strategy = os.path.join(output_dir, strategy_file)
        df = pd.read_csv(strategy, header=1)

        if not least_telescope:
            df.sort_values(by='Detprob1', ascending=False, inplace=True)
            optimal_strategy = df.iloc[0]

#     print(df)
        if least_telescope:
## NEW LINES: THESE ARE FOR LEAST TELESCOPE TIME. WRITTEN BY JOHNNY NOT FULLY UNDER ##
            #read in the strategy csv file
            #determine time delays. note: ask Johnny what specifically this does. he's said Observation01 = time delays after merger on first pass, and Observation02 = time delays after merger on second pass.
            strategy_time_delays=np.add(-1*np.array(df["Observation01"].values),df["Observation02"].values) 
            strategy_time_delays=np.add(-1*np.array(df["Observation01"].values),df["Observation02"].values) 

            #limit the strategy time delay. this somehow does least telescope time
            df=df[np.logical_or(strategy_time_delays > 0.6,strategy_time_delays < 0.4) ]

            total_telescope_time=np.add(df["Telescope_time01"].values,df["Telescope_time02"].values)
            total_telescope_time=np.add(df["Telescope_time01"].values,df["Telescope_time02"].values)
            prob_all=df["Detection Probability"].values#top["Detection Probability"].values
            prob_top_test=max(prob_all)-0.01#-1.5
            all_df_top=df[df["Detection Probability"]>prob_top_test].copy().reset_index(drop=True)
            prob_top=max(df["Detection Probability"].values)
            
            ltt_config=[0.05,0.10,0.15]
            df_ltt=df[df["Detection Probability"].values > (prob_top-(ltt_config[1]*prob_top))]
    
            df.sort_values(by='Detprob1', ascending=False, inplace=True)
            optimal_strategy = df_ltt.iloc[0]
#     print(optimal_strategy)

    outer = optimal_strategy['Region Coverage']
    inner = optimal_strategy['Region Coverage_deep']
    filt = optimal_strategy['Filter_comb']
    exposure_outer1 = optimal_strategy['Exposure01']
    exposure_inner1 = optimal_strategy['Exposure01_deep']
    exposure_outer2 = optimal_strategy['Exposure02']
    exposure_inner2 = optimal_strategy['Exposure02_deep']
    json_output = os.path.join(output_dir, f"des-gw_{outname}_{sky_condition}.json")

    df.assign(json_output=json_output)
    # else:
    #     # Default Strategy for BBHs
    #     outer, inner, filt = 0.9, 0.5, 'i'
    #     exposure_inner, exposure_outer = 90, 90

    exposure_inner = [exposure_inner1, exposure_inner2]
    exposure_outer = [exposure_outer1, exposure_outer2]
        
    print(f'OneRing inputs: skymap: {skymap}, outer: {outer}, inner: {inner}, filt: {filt}, exp_out: {exposure_outer}, exposure_inner: {exposure_inner}, mjd:{mjd}')
    #run updated onering!
    OneRing.run_or(
        skymap,
        outer,
        inner,
        filt,
        exposure_inner,
        exposure_outer,
        mjd,
        resolution=64,
        jsonFilename=json_output
    )

    subject = f'Strategy for event {event}'
    text = f'Outer region coverage: {outer}\n' +\
        f'Inner region coverage: {inner}\n' +\
        f'Filter: {filt}\n' +\
        f'Exposure for Outer Region: {exposure_outer}\n' +\
        f'Exposure for Inner Region: {exposure_inner}\n' +\
        f'Json file path: {json_output}\n'
    
    if official:
        mode = 'observation'
    else:
        mode = 'test'

    slack_bot = SlackBot(mode=mode)

    slack_bot.post_message(subject=subject,
                           text=text)

moony_strategy = multiprocessing.Process(target=run_strategy_and_onering,
                                        args=(skymap,
                                            trigger_id,
                                            mjd,
                                            'moony',
                                            event,))

notmoony_strategy = multiprocessing.Process(target=run_strategy_and_onering,
                                            args=(skymap,
                                                trigger_id,
                                                mjd,
                                                'notmoony',
                                                event,))

moony_strategy.start()
notmoony_strategy.start()

timer_end = time.perf_counter()
print(f"Finished recycler in {timer_end - timer_start:0.4f} seconds")

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
