import json
import os
import sys
sys.path.insert(0, '/data/des70.a/data/desgw/O4/Main-Injector-O4b/desgw_db_writer')
import desgw_db_writer.api as DESGWApi
os.environ["API_BASE_URL"] = "https://desgw-api-physics-soares-santos-flaskapi.apps.gnosis.lsa.umich.edu/api/v0.1/"


print(DESGWApi.test_connection())