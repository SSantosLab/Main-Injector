import os, io, string, random, datetime, base64, json
from astropy.time import Time

def makeBayestarMock():
    t = Time(str(datetime.datetime.utcnow()), format='iso', scale='utc')
    trigger_id = 'MS'+t.isot[2:4]+t.isot[5:7]+t.isot[8:10]+''.join(random.choices(string.ascii_lowercase, k=3))
    outdir = "test_hexes/mock_simulations/"+trigger_id+"/"
    os.makedirs(outdir, exist_ok=True)
    sourcefile = "test_hexes/host_galaxies.txt"
    
    os.system('test_hexes/bayestar_injection.sh {} {} {} {}'.format(t.gps-1, t.gps-0.5, sourcefile, outdir))
    
    if os.path.isfile(outdir+'0.fits'):
        with open(outdir+"0.fits", "rb") as skymap_binary:
            skymap_bytes = base64.b64encode(skymap_binary.read())
    else:
        print('No .fits file: BAYESTAR failure or event not detected')
    
    alert = {
        "alert_type": "PRELIMINARY",
        "time_created": t.isot,
        "superevent_id": trigger_id,
        "urls": {
            "gracedb": "https://example.org/superevents/{}/view/".format(trigger_id)
        },
        "event": {
            "time": t.isot,
            "far": 9.11069936486e-14,
            "significant": True,
            "instruments": [
                "H1",
                "L1",
                "V1"
            ],
            "group": "CBC",
            "pipeline": "gstlal",
            "search": "MDC",
            "properties": {
                "HasNS": 0.95,
                "HasRemnant": 0.91,
                "HasMassGap": 0.01
            },
            "classification": {
                "BNS": 0.95,
                "NSBH": 0.01,
                "BBH": 0.03,
                "Terrestrial": 0.0
            },
            "duration": None,
            "central_frequency": None,
            "skymap": skymap_bytes.decode("utf-8")
        },
        "external_coinc": None
    }
    
    json_path = outdir+trigger_id+"-preliminary.json"
    with open(json_path, "w") as outfile: 
        json.dump(alert, outfile)
    return json_path