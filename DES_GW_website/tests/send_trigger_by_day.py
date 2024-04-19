import json
import os
import sys
sys.path.insert(0, '/data/des70.a/data/desgw/O4/Main-Injector-O4b/desgw_db_writer')
import desgw_db_writer.api as DESGWApi
os.environ["API_BASE_URL"] = "https://desgw-api-physics-soares-santos-flaskapi.apps.gnosis.lsa.umich.edu/api/v0.1/"

js = {
"trigger_label":"S240413p",
"date":"2024-04-18 12:00:00",
"type":"BBH",
"ligo_prob":float(0),
"far":float(0),
"distance":float(0),
"n_hexes":int(0),
"econ_prob":float(0),
"econ_area":float(0),
"need_area":float(0),
"quality":float(0),
"exp_time":"",
"filter":"",
"hours":float(0),
"n_visits":int(0),
"n_slots":int(0),
"b_slot":int(0),
"prob_region_50":float(0),
"prob_region_90":float(0),
"prob_coverage":float(0),
"snr":float(0),
"chirp_mass":float(0),
"component_mass_1":float(0),
"component_mass_2":float(0),
"season":"",
"prob_vs_slot_plot":"",
"centered_gif_plot":"",
"ligo_prob_contour_plot":"",
"des_prob_vs_ligo_prob_plot":"",
"des_limit_mag_map":"",
"des_limit_mag_map_src":"",
"highest_prob_json":"",
"low_tt_json":"",
"log_link":"",
"strategy_table":"",
"initial_skymap":"",
"final_skymap":"",
"airmass":"",
"cumulative_hex_prob":"",
"galaxies_plot_initial":"",
"galaxies_plot_final":"",
"galaxy_percentage_file":"",
"moon":""}
# Now js is in the format that we need it to be

api = DESGWApi.DESGWApi()

newjs = json.dumps(js)

print(api.add_trigger_by_day(js))

# print(api.__dir__(),"\n",js)

