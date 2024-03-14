import numpy as np

# Support:
#
# How to use JSONs
#   import json
#   from pprint import pprint
#   json_data=open('json_data')
#   data = json.load(json_data)
#   pprint(data)
#   json_data.close()
#
# What is format of JSON for CTIO?
#   http://dns.ctio.noao.edu/noao/content/Preparing-Observing-Scripts
#   {"count": 3,
#       "expType": "object",
#       "object": "Sgr",
#       "filter": "g",
#       "RA": 281.0,
#       "dec": -29.866,
#       "expTime": 60
#       },

#
# This code wants a list ra,decs and will write a JSON file
# to cause the Blanco to observe them.
#

# Nov13, 2022 simplified by Jim Annis, will break use in the main_injector,
# but this is what is really needed

def writeJson(
    ra,
    dec,
    exposureTimeList,
    flt="r",
    trigger_id="LIGO/Virgo",
    trigger_type="bright",
    propid='propid',
    skymap="bayestar.fits",
    jsonFilename="des-gw.json"
):

    fd = open(jsonFilename, "w")
    fd.write("[\n")

    nHexes = ra.size
    for i in range(0, nHexes):
        exp = exposureTimeList[i]
        tra = ra[i]
        tdec = dec[i]
        filt = flt[i]
        # comment = "DESGW: LIGO {} event {}: {} of {}, hex {} ".format(
        #    trigger_type, seqid, seqnum, seqtot, id[i], 9)
        comment = "{} strategy {} on {}: image {}, filter {},".format(
            trigger_id, trigger_type, skymap, i+1, filt)
        object = trigger_id

        fd.write("{")
        fd.write(" \"expType\" : \"object\",\n")
        fd.write("  \"object\" : \"{}\",\n".format(object))
        fd.write("  \"expTime\" : {:d},\n".format(int(exp)))
        fd.write("  \"wait\" : \"False\",\n")
        fd.write("  \"filter\" : \"{}\",\n".format(filt))
        fd.write("  \"program\" : \"des gw\",\n")
        fd.write("  \"RA\" : {:.6f},\n".format(tra))
        fd.write("  \"dec\" : {:.5f},\n".format(tdec))
        fd.write("  \"propid\" : \"{}\",\n".format(propid))
        fd.write("  \"comment\" : \"{}\",\n".format(comment))
        fd.write("  \"HEX\" : -999 \n")
        # note lack of comma for end
        fd.write("}")
        if (i == nHexes-1):
            pass
        else:
            fd.write(",")
        fd.write("\n")
    fd.write("]\n")
    fd.close()

