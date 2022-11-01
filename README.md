# Main Injector Docs

Main Injector is a python package designed by the SoaresSantosLab team in order to retrieve alerts from LVC GCN and work with those alerts to search for Electromagnetic Counterparts of Gravitational Wave Events

## Folder Structure:

```
src
└── main-injector
    ├── __init__.py
    ├── configs
    │   ├── TESTINGrecycler.yaml
    │   ├── configs.yaml
    │   ├── resimulator.yaml
    │   └── strategy.yaml
    ├── events
    │   └── checkevent.py
    ├── gwhelper
    │   ├── dark_siren_hexes.py
    │   ├── gw_map_configure.py
    │   └── gwwide.py
    ├── hex
    │   ├── desHexen.py
    │   ├── hexalate.py
    │   └── observations.py
    ├── kasen
    │   └── kasen_modelspace.py
    ├── knlc
    │   ├── README.md
    │   ├── __init__.py
    │   ├── data
    │   │   └── grouped_photometry.csv
    │   ├── dotunnel.sh
    │   ├── example.sh
    │   ├── kn_brightness_estimate.py
    │   ├── kn_exp_plot.png
    │   ├── kn_mag_plot.png
    │   ├── kn_report.txt
    │   ├── launchjupyter.sh
    │   └── setup.sh
    ├── main.py
    ├── plots
    │   ├── cumulative_plots.py
    │   ├── mcplot.py
    │   └── plotMapAndHex.py
    ├── sky
    │   ├── atmosphere.py
    │   ├── dust.py
    │   ├── info.py
    │   ├── seeingModel.py
    │   ├── simplicity.py
    │   ├── sky_model.py
    │   └── telescope.py
    ├── skymaps
    │   ├── inside_footprint.py
    │   ├── mags.py
    │   ├── mapsAtTimeT.py
    │   ├── mcbryde.py
    │   ├── modelRead.py
    │   └── sourceProb.py
    ├── test.py
    ├── trigger
    │   └── trigger_pages.py
    └── utils
        ├── all_maps.py
        ├── checkevent_config.py
        ├── cumul.py
        ├── decam2hp.py
        ├── des_optimization.py
        ├── distance.py
        ├── hp2np.py
        ├── jsonMaker.py
        ├── make_recycler_config.py
        ├── obsSlots.py
        ├── rotate.py
        └── send_texts_and_emails.p
```

### 1. Setup Main-Injector

To setup Main-Injector, we need to
    
1. Manually install pyslalib. Inside pyslalib folder located in main_injector/pyslalib, do 
    ```bash
    python setup.py install
    ```

2. source SOURCEME file

3. Run test
    ```bash
    python main_injector/recycler.py --skymapfilename=OUTPUT/TESTING/S190814bv/bayestar.fits.gz --triggerpath=OUTPUT/TESTING/ --triggerid=S190814bv
    ```