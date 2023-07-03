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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--mode',
                        default='test',
                        help='Chose `mode` for listen to gcn alerts. ' +\
                        'Default mode is `test`. ' +\
                        'For offline testing, use `test` mode. ' +\
                        'For streamed mock alerts, use `mock` mode. ' +\
                        'For real streamed alerts, use `observation` mode.',
                        choices=['test', 'mock', 'observation'])

    args = parser.parse_args()
    mode = args.mode

    gw_streamer = GWStreamer(mode=mode)
    email_bot = EmailBot(mode=mode)

    with open('configs/gcn_credentials.yaml', 'r', encoding='utf-8') as f:
        gcn = yaml.load(f, Loader=SafeLoader)

    if mode == 'test':
        fake_alert_list = glob.glob('test_hexes/MS181101ab*preliminary*.json')
        for fake_alert in fake_alert_list:
            with open(fake_alert, 'r', encoding='utf-8') as f:
                gcn_fake_alert = json.load(f)

            gw_streamer.handle(gcn_fake_alert)

    else:
        try:
            consumer = Consumer(client_id=gcn['client_id'],
                                client_secret=gcn['client_secret'])
            consumer.subscribe(['igwn.gwalert'])

            while True:
                for message in consumer.consume(timeout=1):
                    gcn_alert = message.value()
                    gw_streamer.handle(gcn_alert)

        except Exception as e:
            print(e)
            # EmailBot.send_email(subject='Listener went down, see traceback',
            #                     text=traceback.format_exc(),
            #                     emergency=True)
