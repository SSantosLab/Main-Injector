import os

trigger_id = 123456
retraction = 'False'
mapname = 'bayestar'
far = '140mpc'

text = f"""\
    Trigger {trigger_id}\n
    Alert Type: {retraction}
    FAR: {str(far)}
    Map: {mapname}
    URL: https://gracedb.ligo.org/superevents/{trigger_id}
    view analysis has begun, please hold tight for a DESGW webpage.
    URL for DESGW webpage: http://des-ops.fnal.gov:8080/desgw/Triggers/\
{trigger_id}/{trigger_id}_trigger.html
    DO NOT REPLY TO THIS THREAD, NOT ALL USERS WILL SEE YOUR RESPONSE.
"""

print(text)