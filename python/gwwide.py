#!/usr/bin/env python
"""Convert a native GW followup json obs queue to a DES wide compatible one
"""

import sys
import json
import logging

from argparse import ArgumentParser
from collections import defaultdict
from functools import partial
from math import sin, cos, radians, degrees, acos
from copy import deepcopy

# constants

max_match_angle = 3.0/60.0
camera_radius = 1.3

# exception classes
# interface functions

def gwwide(wide_queue, gw_queue):
    """Replace entries in a GW quere with nearby entries in the DES wide survey

    :Parameters:
        - `wide_queue`: the DES wide survey queue (list of dicts)
        - `gw_queue`: the naieve GW queue (list of dicts)

    @returns: a queue (list of dicts) with entries replaced
    """
    fixed_queue = [fix_obs(gw_obs, wide_queue) for gw_obs in gw_queue]
    return fixed_queue

# classes
# internal functions & classes

def angle(obs1, obs2):
    """Angular distance between two queue entries (deg)

    :Parameters:
        - `obs1`: queue entry 1
        - `obs2`: queue entry 2

    @returns: angular distance in degress
    """

    ra1, decl1 = obs1['RA'], obs1['dec']
    ra2, decl2 = obs2['RA'], obs2['dec'] 
    if abs(ra1-ra2)<(1.0/(60*60*10000)) and abs(decl1-decl2)<(1.0/(60*60*10000)):
        return 0.0

    sep = degrees( acos( sin(radians(decl1))*sin(radians(decl2))
                         + cos(radians(decl1))*cos(radians(decl2))*cos(radians(ra1-ra2)) ) )
    return sep

def wide_covered(gw_obs, wide_queue):
    """Return whether a queue entry is covered by wide survey hexes

    Assume it is covered if the pointing is covered by hexes
    in all 10 DES wide tilings.

    :Parameters:
        - `gw_obs`: the queue entry to check
        - `wide_queue`: the list of all DES wide queue entrues

    @returns: true iff an obs is covered by wide survey hexes
    """
    # If we are within a camera radius of a pointing in every tiling, we are in the footprint
    for tiling_id in range(1,11):
        wide_in_tiling = [w for w in wide_queue 
                          if w['tiling_id']==tiling_id 
                          and w['filter']==gw_obs['filter']]
        nearest_in_tiling = min(wide_in_tiling, 
                                key=partial(angle, gw_obs))
        angle_to_nearest = angle(gw_obs, nearest_in_tiling)
        if angle_to_nearest > camera_radius:
            return False
    return True
    
def fix_obs(gw_obs, all_filter_wide_queue):
    """Adjust an exposure queue entry to match a wide field entry, if close

    :Parameters:
        - `gw_obs`: the GW obs queue entry
        - `all_filter_wide`: a list of all DES wide survey queue entrues

    @returns: a queue entry, adjusted if necessary
    """
    wide_queue = [w for w in all_filter_wide_queue 
                  if w['filter']==gw_obs['filter']]

    logging.debug("GW exposure:   %3.4f, %3.4f", gw_obs['RA'], gw_obs['dec'])

    nearest_wide_obs = min(wide_queue, key=partial(angle, gw_obs))

    nearest_angle = angle(gw_obs, nearest_wide_obs)
    logging.debug("wide exposure: %3.4f, %3.4f (tiling %d) is %3.4f degrees away",
                  nearest_wide_obs['RA'], nearest_wide_obs['dec'], 
                  nearest_wide_obs['tiling_id'],
                  nearest_angle)

    if not 'exptime' in gw_obs:
        gw_obs['exptime'] = gw_obs['expTime']
    if nearest_angle > max_match_angle or gw_obs['exptime']!=90:
        logging.warning("Exposure at %(RA)3.4f %(dec)3.4f not matched", gw_obs)
        for kw in ['count','seqtot','seqnum']:
            gw_obs[kw] = int(gw_obs[kw])

        return gw_obs

    fixed_obs = deepcopy(nearest_wide_obs)
    for kw in ['count','seqid','seqnum','seqtot','note','comment','propid']:
        fixed_obs[kw] = gw_obs[kw]

    for kw in ['count','seqtot','seqnum']:
        fixed_obs[kw] = int(fixed_obs[kw])

    return fixed_obs

def file_gwwide(gw_fname, wide_fname, fixed_fname):
    """Adjust exposure queue from a file to match a wide field entries, if close

    :Parameters:
        - `gw_fname`: the GW obs queue json fiel
        - `wide_fname`: the json file with all DES wide field queue entries
        - `fixed_fnamd`: the json obs queue file into which to write entries

    """

    with open(wide_fname, 'r') as f:
        wide_queue = json.load(f)
    logging.info("Loaded %d wide survey entries from %s", len(wide_queue), wide_fname)

    print 'gw_fname',gw_fname
    with open(gw_fname, 'r') as f:
        gw_queue = json.load(f)
    logging.info("Loaded %d GW followup entries from %s", len(gw_queue), gw_fname)

    fixed_queue = gwwide(wide_queue, gw_queue)
    logging.info("Finished correcting queue entries")

    with open(fixed_fname, 'w') as f:
        json.dump(fixed_queue, f, indent=4)
    logging.info("Saved %d GW followup entries in %s", len(fixed_queue), fixed_fname)

    return 0

def main():

    parser = ArgumentParser(
        description='Fix a GW followup script to match wide-survey exposures')    
    parser.add_argument('gw_script', help='Filename of GW obs script')
    parser.add_argument('wide_script', help='Filename of wide obs script')
    parser.add_argument('fixed_script', help='Filename we write fixed script to')
    parser.add_argument('-v', '--verbose', action="count", 
                        help="be verbose")

    args = parser.parse_args()

    log_levels = defaultdict(lambda:logging.WARNING)
    log_levels[1] = logging.INFO
    log_levels[2] = logging.DEBUG
    logging.basicConfig(format='%(asctime)s %(message)s',
                        level=log_levels[args.verbose])

    logging.info("Begin")
    
    status = file_gwwide(args.gw_script, args.wide_script, args.fixed_script)
    return status


if __name__ == '__main__':
    status = main()
    sys.exit(status)    

