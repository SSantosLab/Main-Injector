"""
Python CLI Tool to make a config file.
"""
import argparse
import yaml
import os

def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--trigger-path', default='./', nargs='?', help='trigger-path.')
    parser.add_argument('--trigger-id', help='TriggerID of LVK alert.')
    parser.add_argument('--type', default='sim', nargs='?', help='type of running: real or sim.')
    parser.add_argument('--debug', default='False', nargs='?', help='Set to debug mode.')
    parser.add_argument('--camera', default='decam', nargs='?', help='camera used for observation strategy.')
    parser.add_argument('--resolution', default=64, nargs='?', help='Camera resolution')
    parser.add_argument('--gif-resolution', default=0.1, nargs='?')
    parser.add_argument('--make-maps', default=True ,nargs='?',)
    parser.add_argument('--make-hex',default=True, nargs='?',)
    parser.add_argument('--make-json',default=True, nargs='?',)
    parser.add_argument('--make-gif',default=True, nargs='?',)
    parser.add_argument('--area-per-hex',default=3.0, nargs='?', help="Area per hex in sq deg.")
    parser.add_argument('--overhead', default=30, nargs='?',)
    parser.add_argument('--all-sky', default=True, nargs='?',)
    parser.add_argument('--season-start', default='', nargs='?',)
    parser.add_argument('--season-end',default=True, nargs='?',)
    parser.add_argument('--distance',default=True, nargs='?',)
    parser.add_argument('--force-distance', default=True, nargs='?',)
    parser.add_argument('--map-name', default='bayestar.fits',  nargs='?',)
    parser.add_argument('--output', default='./', help='Output path.')
    parser.add_argument('--skip-all',  nargs='?',)
    parser.add_argument('--skip-plots',  nargs='?',)
    parser.add_argument('--ns-threshold', default=.5,  nargs='?',)
    parser.add_argument('--observed-event', nargs='?',)
    parser.add_argument('--force-mjd', nargs='?',)
    parser.add_argument('--mjd', nargs='?',)
    parser.add_argument('--rem-strategy', nargs='?',)

    return parser

def main():

    p = parser()
    args = p.parse_args()
    
    
    destination = os.path.join(args.output, 'recycler.yaml')
    with open(destination, 'w') as outfile:
        yaml.dump(vars(args), outfile, default_flow_style=False)

main()