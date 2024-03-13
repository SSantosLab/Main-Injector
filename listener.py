"""
Python script to listen from LVK alerts trough gcn.
"""

import glob
import time
import traceback
import json
from argparse import ArgumentParser
import yaml
from yaml.loader import SafeLoader
from gcn_kafka import Consumer
from handlers.gwstreamer import GWStreamer
from handlers.emails import EmailBot
from test_hexes.mock_bayestar_event import makeBayestarMock

def elapsedTimeString(start):
    elapsed = int(time.time() - start)
    return "{}h {:02d}m {:02d}s".format(elapsed//(60*60), elapsed//60%60, elapsed%60)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--mode',
                        default='test',
                        help='Chose `mode` for listen to gcn alerts. ' +\
                        'Default mode is `test`. ' +\
                        'For offline testing, use `test` mode. ' +\
                        'For streamed mock alerts, use `mock` mode. ' +\
                        'For simulated BAYESTAR mock alerts, use `mock-bayestar` mode.' +\
                        'For real streamed alerts, use `observation` mode.',
                        choices=['test', 'mock', 'mock-bayestar', 'observation'])

    args = parser.parse_args()
    mode = args.mode
    start_time = time.time()
    print("Listener Activating")
    print('Running Listener in {} mode...'.format(mode))
    
    gw_streamer = GWStreamer(mode=mode)
    email_bot = EmailBot(mode=mode)

    with open('configs/gcn_credentials.yaml', 'r', encoding='utf-8') as f:
        gcn = yaml.load(f, Loader=SafeLoader)

    if mode == 'test':
        print('Reading test event...')
        fake_alert_list = glob.glob('test_hexes/MS181101ab*preliminary*.json')
        print('Passing event to Handler - Listener took '+elapsedTimeString(start_time))
        for fake_alert in fake_alert_list:
            with open(fake_alert, 'r', encoding='utf-8') as f:
                gcn_fake_alert = json.load(f)

            gw_streamer.handle(gcn_fake_alert)
    	
    elif mode == 'mock-bayestar':
        print('Simulating BAYESTAR event...')
    	fake_alert = makeBayestarMock()
    	with open(fake_alert, 'r', encoding='utf-8') as f:
    	    gcn_fake_alert = json.load(f)
    	
    	print('Passing event to Handler - Listener took '+elapsedTimeString(start_time))
    	gw_streamer = GWStreamer(mode='mock')
    	gw_streamer.handle(gcn_fake_alert)
    	
    else:
        try:
            consumer = Consumer(client_id=gcn['client_id'],
                                client_secret=gcn['client_secret'])
            consumer.subscribe(['igwn.gwalert'])

            while True:
                for message in consumer.consume(timeout=10):
                    print('Trigger Received...')
                    gcn_alert = json.loads(message.value())
                    print('Passing event to Handler.')
                    gw_streamer.handle(gcn_alert)

        except Exception as e:
            print(e)
            email_bot.send_email(subject='Listener went down, see traceback',
                                text=traceback.format_exc(),
                                emergency=True)
    print('Main Injector exiting after '+elapsedTimeString(start_time))
