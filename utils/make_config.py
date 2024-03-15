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
    parser.add_argument('--do-make-maps', default=True ,nargs='?',)
    parser.add_argument('--do-make-hexes',default=True, nargs='?',)
    parser.add_argument('--do-make-jsons',default=False, nargs='?',)
    parser.add_argument('--do-make-gifs',default=True, nargs='?',)
    parser.add_argument('--area-per-hex',default=3.0, nargs='?', help="Area per hex in sq deg.")
    parser.add_argument('--overhead', default=30, nargs='?',)
    parser.add_argument('--allSky', default=True, nargs='?',)
    parser.add_argument('--season-start', default='', nargs='?',)
    parser.add_argument('--season-end',default=True, nargs='?',)
    parser.add_argument('--distance',default=True, nargs='?',)
    parser.add_argument('--force-distance', default=True, nargs='?',)
    parser.add_argument('--map-name', default='bayestar.fits',  nargs='?',)
    parser.add_argument('--output', default='./', help='Output path.')
    parser.add_argument('--skipAll',  nargs='?',)
    parser.add_argument('--skipPlots',  nargs='?',)
    parser.add_argument('--ns-threshold', default=.5,  nargs='?',)
    parser.add_argument('--observed-event', nargs='?',)
    parser.add_argument('--force-mjd', nargs='?',)
    parser.add_argument('--mjd', nargs='?',)
    parser.add_argument('--rem-strategy', nargs='?',)
    parser.add_argument('--real-or-sim', nargs='?', default='sim')
    parser.add_argument('--skymap-filename', nargs='?', default='Default')
    parser.add_argument('--wrap-all-triggers', nargs='?', default=False)
    parser.add_argument('--force-recycler-mjd', nargs='?',default=False)
    parser.add_argument('--start-of-season', nargs='?',default=57444)
    parser.add_argument('--end-of-season', nargs='?',default=57996)
    parser.add_argument('--events-observed', nargs='?',default=['GW150914' , 'GW151226'])
    parser.add_argument('--kasen-fraction', nargs='?',default=50)
    parser.add_argument('--exposure-length-Rem', nargs='?',default=[ 60., 90. ])
    parser.add_argument('--exposure-filter-Rem', nargs='?',default=[ 'i', 'z'])
    parser.add_argument('--exposure-tiling-Rem', nargs='?',default=[ 0, 5])
    parser.add_argument('--maxHexesPerSlot-Rem', nargs='?',default=6)
    parser.add_argument('--exposure-length-BH', nargs='?',default=[ 90.])
    parser.add_argument('--exposure-filter-BH', nargs='?',default=['r'])
    parser.add_argument('--exposure-tiling-BH', nargs='?',default=[ 10])
    parser.add_argument('--maxHexesPerSlot-BH', nargs='?',default=18)
    parser.add_argument('--propid-BH', nargs='?',default='propideBH')
    parser.add_argument('--propid-Rem', nargs='?',default='propideRem')
    parser.add_argument('--max-number-of-hexes-to-do', nargs='?',default=10000)
    return parser

def main():

    p = parser()
    args = p.parse_args()
    
    
    destination = os.path.join(args.output, 'recycler.yaml')
    with open(destination, 'w', encoding='utf8') as outfile:
        yaml.dump(vars(args), outfile, default_flow_style=False)

main()