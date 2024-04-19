import handlers.short_latency_plots as slp
import numpy as np
target_coords = [165.72,11.64]
arr = ["2024-04-29"]
for day in arr:
    slp.moon_airmass("S240413p", day, target_coords)
