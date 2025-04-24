import sys
sys.path.insert(1, '/data/des70.a/data/desgw/O4/Main-Injector-O4b/python/')

import OneRing

mjd = 60762.23659722

event_id = 'S250328ae'
event_type = 'UPDATE'

lp, dp = OneRing.run_or(f'OUTPUT/O4REAL/{event_id}/{event_type}/bayestar.fits.gz', 0.90, 0.50, 'ii', [90, 90], [90, 90], mjd, [1, 1], 1, 1, resolution=64, jsonFilename=f'OUTPUT/O4REAL/{event_id}/{event_type}/{event_id}_manual.json')
