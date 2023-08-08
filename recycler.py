import os
import OneRing
import pandas as pd
import multiprocessing
import yaml
from datetime import datetime
from subprocess import run
from argparse import ArgumentParser
from astropy.io import fits
from handlers.slack import SlackBot
import handlers.observations as gwplot
from handlers.gwstreamer import GWStreamer
from loguru import logger
from api import DESGWApi
def parser():
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
    return parser

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
            
        mjd = str(mjd).replace('.','')
        current_time = mjd

        input_dir = os.path.dirname(skymap)
        output_dir = input_dir
        cmd = 'python ' +\
            'python/knlc/kn_strategy_sims_MI.py '+\
            f'--input-skymap {skymap} '+\
            f'--output {output_dir} '+\
            f'--teff-type {sky_condition} ' +\
            f'--kn-type {kn_type} ' +\
            f'--time {current_time}'
            
        output_log = os.path.join(output_dir,
                                f'{sky_condition}_{kn_type}_strategy.log')
        strategy_log = open(output_log, 'w')

        print('strategy started!')
        run(cmd,
            shell=True,
            stdout=strategy_log,
            stderr=strategy_log,
            text=True)
            
        strategy_file = f'bayestar_{sky_condition}_{kn_type}_{current_time}' +\
                        '_allconfig.csv'
            
        strategy_log.close()
        strategy = os.path.join(output_dir, strategy_file)
            
        df = pd.read_csv(strategy, header=1)
        df.sort_values(by='Detprob1', ascending=False, inplace=True)
        optimal_strategy = df.iloc[0]
        outer, inner, filt, exposure_outer, exposure_inner = optimal_strategy[1:6]
        json_output = os.path.join(output_dir, 
                                f"des-gw_{current_time}_{sky_condition}.json")
        
        df.assign(json_output=json_output)
        filt = filt[0]

    else:
        # Default Strategy for BBHs
        outer, inner, filt = 0.9, 0.5, 'i'
        exposure_inner, exposure_outer = 90, 90
    
    OneRing.run_or(
        skymap,
        outer,
        inner,
        filt,
        exposure_inner,
        exposure_outer,
        mjd,
        resolution=64,
        jsonFilename=json_output,
        max_hex_count=max_hex_count,
        max_hex_time=max_hex_time
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

if __name__ == '__main__':

    GW = GWStreamer(mode='observation')
    DESGW = DESGWApi()
    p = parser()
    args = p.parse_args()
    trigger_path = GW.OUTPUT_TRIGGER
    trigger_id = args.trigger_id
    output_alert = os.path.join(trigger_path, trigger_id)
    logger.add(f"{output_alert}.log", level="DEBUG", backtrace=True, diagnose=True)
    logger.add(lambda message: print(message, end=''), level="DEBUG", backtrace=True, diagnose=True)

    logger.info('Settings for Recycler:')
    arguments = vars(args)
    for key,value in arguments.items():
        logger.info(f'{key}: {value}')

    logger.info('Creating recycler.yaml config file')
    with open(os.path.join(trigger_path,'recyler.yaml'), 'w') as f:
        f.write(yaml.dump(arguments))
    
    skymap = args.skymap
    event = args.event
    max_hex_time = args.max_hex_time
    max_hex_count = args.max_hex_count
    official = args.official
    with fits.open(skymap) as f:
        header = f[1].header

    mjd = header['MJD-OBS']

    logger.info('Making Initial plots')
    gwplot.make_plots_initial(url=output_alert,
                              name=trigger_id)
    
    alert = f'{output_alert}/{trigger_id}.json'
    alert_data = GW.make_trigger_data(alert)
    DESGW.add_trigger(alert)

    # moony_strategy = multiprocessing.Process(target=run_strategy_and_onering,
    #                                         args=(skymap,
    #                                             trigger_id,
    #                                             mjd,
    #                                             'moony',
    #                                             event,))

    # notmoony_strategy = multiprocessing.Process(target=run_strategy_and_onering,
    #                                             args=(skymap,
    #                                                 trigger_id,
    #                                                 mjd,
    #                                                 'notmoony',
    #                                                 event,))
    
    # logger.info('Firing off moony and notmoony strategy')
    # moony_strategy.start()
    # notmoony_strategy.start()
