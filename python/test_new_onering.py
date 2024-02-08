import numpy as np
import hp2np 
import hexalate
import decam2hp
import jsonMaker
from os import getenv
import pandas as pd
import sys
import  matplotlib.pyplot as plt
import hex_functions
from hex_functions import get_hexinfo 
import OneRing


skymap = "/data/des70.a/data/the_main_injectors/get_realistics_exp_for_main_injector/skymaps/134.fits.gz"
mjd = 55961
strategy = f'134lowdist152.0_lowarea168.0_sortTT.csv'
df = pd.read_csv(strategy)
df.sort_values(by='Deprob1', ascending=False, inplace=True)
optimal_strategy = df.iloc[0]
outer = optimal_strategy['Region Coverage']
inner = optimal_strategy['Region Coverage_deep']
filt = optimal_strategy['Filter_comb'][0]
exposure_outer = optimal_strategy['Exposure01']
exposure_inner = optimal_strategy['Exposure01_deep']

OneRing.run_or(skymap, outer, inner, filt, exposure_outer, exposure_inner, 55961, jsonFilename="test_afactor.json")
