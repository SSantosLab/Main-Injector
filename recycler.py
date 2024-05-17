import os
os.environ["API_BASE_URL"] = "https://desgw-api-physics-soares-santos-flaskapi.apps.gnosis.lsa.umich.edu/api/v0.1/"
import numpy as np
import OneRing
import pandas as pd
import multiprocessing
from subprocess import run
from argparse import ArgumentParser
from astropy.io import fits
from handlers.slack import SlackBot
import time
import sys
from datetime import datetime
sys.path.insert(0, '/data/des70.a/data/desgw/O4/Main-Injector-O4b/desgw_db_writer')
import desgw_db_writer.api as DESGWApi

def elapsedTimeString(start):
    elapsed = int(time.time() - start)
    return "{}h {:02d}m {:02d}s".format(elapsed//(60*60), elapsed//60%60, elapsed%60)

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

parser.add_argument('--alertNum',
                    type=str,
                    default="PRELIMINARY_0",
                    help='Stores the preliminary iterator')

#NEW ADDITION -- LEAST TELESCOPE TIME STRATEGY?#
parser.add_argument('--ltt',
                    action='store_true',
                    default=False,
                    help='If true, uses least telescope time strategy.')

parser.add_argument('--creationTime',
                    type=str,
                    default=str(datetime.now()),
                    help='The date the GCN was created.')

args = parser.parse_args()
t0 = time.time()
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
alertNum = args.alertNum
least_telescope = args.ltt
creationTime = args.creationTime
desgw =  DESGWApi.DESGWApi() # Update the base URL of the base API

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

        print('strategy started!', flush=True)
        run(cmd,
            shell=True,
            stdout=strategy_log,
            stderr=strategy_log,
            text=True)
            
        strategy_file = f'bayestar_{sky_condition}_{kn_type}_{outname}' +\
                        '_allconfig.csv'
        os.chmod(strategy_file, 0o0777)
            
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
    
            df_ltt.sort_values(by='Detprob1', ascending=False, inplace=True)
            optimal_strategy = df_ltt.iloc[0]
#     print(optimal_strategy)

    outer = optimal_strategy['Region Coverage']
    inner = optimal_strategy['Region Coverage_deep']
    filt = optimal_strategy['Filter_comb']
    exposure_outer1 = optimal_strategy['Exposure01']
    exposure_inner1 = optimal_strategy['Exposure01_deep']
    exposure_outer2 = optimal_strategy['Exposure02']
    exposure_inner2 = optimal_strategy['Exposure02_deep']
    detP_1 = optimal_strategy['Detprob1']
    detP_2 = optimal_strategy['Detprob2']
    json_output = os.path.join(output_dir, f"des-gw_{outname}_{sky_condition}.json")

    df.assign(json_output=json_output)
    # else:
    #     # Default Strategy for BBHs
    #     outer, inner, filt = 0.9, 0.5, 'i'
    #     exposure_inner, exposure_outer = 90, 90

    exposure_inner = [exposure_inner1, exposure_inner2]
    exposure_outer = [exposure_outer1, exposure_outer2]
    detP = [detP_1,detP_2]

    print(f'OneRing inputs: skymap: {skymap}, outer: {outer}, inner: {inner}, filt: {filt}, exp_out: {exposure_outer}, exposure_inner: {exposure_inner}, mjd:{mjd}', flush=True)
    #run updated onering!
    local_prob, disco_prob = OneRing.run_or(skymap,
                                            outer,
                                            inner,
                                            filt,
                                            exposure_inner,
                                            exposure_outer,
                                            mjd,
                                            detP,
                                            creationTime,
                                            resolution=64,
                                            jsonFilename=json_output)
    

    local_prob = round(local_prob * 100,1)
    disco_prob = round(disco_prob,1)


    trigger_data = {"trigger_label":trigger_id,
                    "date": creationTime, # Add the date from the gcn
                    'ligo_prob':local_prob,
                    "exp_time":[exposure_inner,exposure_outer].__str__(),
                    "filter":filt,
                    "prob_coverage":disco_prob,
                    "low_tt_json":json_output,
                    "log_link":output_log,
                    "strategy_table":strategy_file
                    }
    
    print("Trigger data to be posted to website",flush=True)
    print("",flush=True)
    for key, val in zip(trigger_data.keys(),trigger_data.values()):
        print("Key:",key,flush=True)
        print("Value:",val,flush=True) 
        print("",flush=True)

    if alertNum=="PRELIMINARY_0":
        desgw.add_trigger_by_day(trigger_data)
    else:
        desgw.update_trigger_by_day(trigger_data)
                                 
    subject = ""
    text = f'*Strategy for event: {trigger_id}* \n\n' +\
        f'*Assumptions* \n' +\
        f'Sky Conditions: {sky_condition}\n' +\
        f'Event type: {event}\n' +\
        f'Light curve model: {kn_type}\n' +\
        f'Localization cumulative probability: {local_prob}%\n' +\
        f'Discovery cumulative probability: {disco_prob}%\n\n' +\
        f'*Optimal Strategy Parameters*\n\n' +\
        f'Outer region coverage: {outer}\n' +\
        f'Inner region coverage: {inner}\n' +\
        f'Filter: {filt} \n' +\
        f'Exposure for Outer Region (seconds): {exposure_outer}\n' +\
        f'Exposure for Inner Region (seconds): {exposure_inner}\n' +\
        f'Json file path: `{json_output}`'
    
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

print('Strategy begun. Recycler exiting after '+elapsedTimeString(t0), flush=True)

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


    # trigger_data = { # for add_trigger_by_day
    #                 "trigger_label":trigger_id,
    #                 "date": str(datetime.now()),
                    # "type":,
                    # "ligo_prob":,
                    # "far":,
                    # "distance":,
                    # "n_hexes":, # total number of hexes - this is in OneRing
                    # "econ_prob":, # not yet incorporated
                    # "econ_area":, # not yet incorporated
                    # "need_area":,  # not yet incorporated
                    # "quality":, # not yet incorporated
                    # "exp_time":[exposure_inner,exposure_outer].__str__(),
                    # "filter":filt,
                    # "hours":, # this is in OneRing
                    # "n_visits":, # number of visits to a hex - this is in OneRing, but probably not important
                    # "n_slots":, # defunct
                    # "b_slot":, # defunct
                    # "prob_region_50":,
                    # "prob_region_90":,
                    # "prob_coverage":disco_prob,
                    # "snr":, # recycler? "An outdated concept" says Isaac
                    # "chirp_mass":, # in gwstreamer
                    # "component_mass_1":, # huh recycler
                    # "component_mass_2":, # huh
                    # "season":, # DONT CARE
                    # "prob_vs_slot_plot":, # defunct
                    # "centered_gif_plot":, # Post OneRing, Isaac probably has code for it - testPlotSkymaps jupyter notebook
                    # "ligo_prob_contour_plot":, # This is one of the SLIPS
                    # "des_prob_vs_ligo_prob_plot":, # OneRing allegedly??
                    # "des_limit_mag_map":, # defunct? maybe calculated in OneRing? No, awesomeness functions uses mags.py to calculate some stufffffff
                    # "des_limit_mag_map_src":, # defunct??
                    # "highest_prob_json":, # we don't use this so defunct
                    # "low_tt_json":json_output,
                    # "log_link":, # OneRing? No, recycler
                    # "strategy_table":, # recycler
                    # "initial_skymap":, # recycler
                    # "final_skymap":, # recycler
                    # "airmass":, # not yet incorporated??? lol
                    # "cumulative_hex_prob":, # This was made somewhere but idk where it is anymore 
                    # "galaxies_plot_initial":, # not yet incorporated
                    # "galaxies_plot_final":, # not yet incorporated
                    # "galaxy_percentage_file":, # not yet incorporated
                    # "moon": # moonplot
                    # }
