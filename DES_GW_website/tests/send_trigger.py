import json
import os
import sys
sys.path.insert(0, '/data/des70.a/data/desgw/O4/Main-Injector-O4b/desgw_db_writer')
import desgw_db_writer.api as DESGWApi
os.environ["API_BASE_URL"] = "https://desgw-api-physics-soares-santos-flaskapi.apps.gnosis.lsa.umich.edu/api/v0.1/"

js = {"trigger_label": "S240413p", 
"mock":False, 
"season": "1599", 
"mjd": 60413.0978125, 
"event_datetime": "2024-04-13 02:20:33", 
"detectors": "H1, L1, V1", 
"lvc_event_url": "https://gracedb.ligo.org/superevents/S240413p/view/"}
# Now js is in the format that we need it to be

api = DESGWApi.DESGWApi()

print(api.add_trigger(js))

# print(api.__dir__(),"\n",js)

