"""
Python script to listen from LVK alerts trough gcn.
"""

import glob
import traceback
import json
from argparse import ArgumentParser
import yaml
from yaml.loader import SafeLoader
from gcn_kafka import Consumer
from handlers.gwstreamer import GWStreamer
from handlers.emails import EmailBot
from subprocess import Popen
    
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--mode',
                        default='test',
                        help='Chose `mode` for listen to gcn alerts. ' +\
                        'Default mode is `test`. ' +\
                        'For offline testing, use `test` mode. ' +\
                        'For streamed mock alerts, use `mock` mode. ' +\
                        'For real streamed alerts, use `observation` mode.',
                        choices=['test', 'observation'])
    parser.add_argument('--test-file',
                        default=None,
                        help='Make a test run with a test json file')

    args = parser.parse_args()
    mode = args.mode
    test_file = args.mode
    gw_streamer = GWStreamer(mode=mode)
    email_bot = EmailBot(mode=mode)

    with open('configs/gcn_credentials.yaml', 'r', encoding='utf-8') as f:
        gcn = yaml.load(f, Loader=SafeLoader)

    if test_file != None:
        test_alert = json.load(f)
        gw_streamer.handle(test_alert)

    else:    
        if mode == 'test':
            fake_alert_list = glob.glob('test_hexes/MS181101ab*preliminary*.json')
            for fake_alert in fake_alert_list:
                with open(fake_alert, 'r', encoding='utf-8') as f:
                    gcn_fake_alert = json.load(f)

                gw_streamer.handle(gcn_fake_alert)

        if mode == 'observation':
            try:
                consumer = Consumer(client_id=gcn['client_id'],
                                    client_secret=gcn['client_secret'])
                consumer.subscribe(['igwn.gwalert'])

                while True:
                    for message in consumer.consume(timeout=10):
                        gcn_alert = json.loads(message.value())
                        gw_streamer.handle(gcn_alert)

            except Exception as e:
                print(e)
                email_bot.send_email(subject='Listener went down, see traceback',
                                    text=traceback.format_exc(),
                                    emergency=True)
