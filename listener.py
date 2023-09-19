"""
Python script to listen from LVK alerts trough gcn.
"""

import glob
import traceback
import json
from datetime import datetime
import pytz
from argparse import ArgumentParser
import yaml
from yaml.loader import SafeLoader
from gcn_kafka import Consumer
from stages.communication import Email
from subprocess import Popen
from loguru import logger
import os

def parser() -> ArgumentParser:

    parser = ArgumentParser()
    parser.add_argument('--mode',
                        default='observation',
                        type=str,
                        help='Chose `mode` for listen to gcn alerts. ' +\
                        'Default mode is `observation`. ' +\
                        'For streamed mock alerts, use `test` mode. ' +\
                        'For real streamed alerts, use `observation` mode.',
                        choices=['test', 'observation'])

    return parser

def create_directories(*args: str) -> None:
    """Setup default directories to store raw alerts as json and logs file."""

    for dirpath in args:
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)

    return

def insert_mode_into_config_file(*configs_file: str,
                               mode: str) -> None:
    """Load mode option into yaml config file.
    
    config_file (str):
        path to config yaml file. Can be more than one.
    mode (bool):
        Listener mode.
    """
    
    for config_file in configs_file:
        with open(config_file, 'r', encoding='utf-8') as f:
            conf = yaml.load(f, Loader=SafeLoader)
        conf['mode'] = mode

        with open(config_file, 'w', encoding='utf-8') as f:
            conf = yaml.dump(conf, f)

    return

def process_raw_alert(alert: dict,
                      alert_path: str,
                      mode: str = 'observation') -> None:
    """preprocess initial raw alert, saving content as json file.
    Return path to saved json file.
    
    Parameters:
    -----------
    alert (dict):
        alert content as dictionary.

    alert_path (str):
        output directory to save alert_path.

    mode (str):
        'test' to preprocess initial mock raw alert.
        'observation' to preprocess initial real raw alert from LVK.
    
    Returns:
        None.
    """

    event_id = alert['superevent_id']
    if mode == 'test':
        if event_id[0] != 'M':
            return
    
    if mode == 'observation':
        if event_id[0] != 'S':
            return
        
    alert_name = os.path.join(alert_path, f'alert{event_id}.json')
    with open(alert_name, 'w') as f:
        json.dump(alert, f)

    logger.info(f"alert content from {event_id} saved in {alert_name}.")
    return alert_name

def listen() -> None:
    """Listens to LVK alerts trough Kafka streaming service."""

    ROOT_DIR = os.environ['ROOT_DIR']
    options = parser().parse_args()
    config_path = os.path.join(ROOT_DIR, 'configs')
    log_path = os.path.join(ROOT_DIR, 'logs')
    raw_alerts_path = os.path.join(ROOT_DIR,'raw_alerts')

    create_directories(config_path, log_path, raw_alerts_path)

    communication_config = os.path.join(config_path, 'communications.yaml')
    gcn_config = os.path.join(config_path, 'gcn_credentials.yaml')
    listener_log = os.path.join(log_path,'listener.log')

    insert_mode_into_config_file(communication_config,
                                 gcn_config,
                                 mode=options.mode)
    
    logger.add(listener_log, level="INFO", rotation="00:00")

    with open(gcn_config, 'r', encoding='utf-8') as f:
        gcn = yaml.load(f, Loader=SafeLoader)

        consumer = Consumer(client_id=gcn['client_id'],
                            client_secret=gcn['client_secret'])
        consumer.subscribe(['igwn.gwalert'])

        try:
            while True:
                for message in consumer.consume(timeout=1):
                    gcn_alert = json.loads(message.value())
                    date = (datetime
                            .now(pytz.utc)
                            .strftime('%H:%M:%S (UT) - %y/%m/%d')
                    )
                    logger.info(f"New alert from {message.topic()} at {date}")

                    alert_name = process_raw_alert(gcn_alert,
                                                    alert_path=raw_alerts_path,
                                                    mode=options.mode,
                                                    )

                # cmd = (f'python recycler.py --alert {alert_name}')

        except Exception as e:
            print(e)
            email_bot = Email.from_config_file(config_path=communication_config)
            email_bot.send_email(subject='Listener went down, see traceback',
                                text=traceback.format_exc(),
                                emergency=True)

if __name__ == "__main__":
    listen()