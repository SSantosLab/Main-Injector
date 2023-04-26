# Main Injector Docs

Main Injector is a python package designed by the SoaresSantosLab team in order to retrieve alerts from LVC GCN and work with those alerts to search for Electromagnetic Counterparts of Gravitational Wave Events

## Setup Main-Injector

To setup Main-Injector, we need to:
    
1. Edit SOURCEME file

    the SOURCEME file sets custom env variables that are necessary in order to recycler.py work fully. One can add specific modifications for your hostname. (If you don't know your hostname, type `echo $HOSTNAME` in a bash shell) 

2. source SOURCEME file.
    - ROOT_DIR: current Main-Injector location.
    - DATA_DIR: path to data directory containing files used by Main-Injector.
    - WEB: path to DES_GW_website.
    - PYTHONPATH: Location to custom Main-Injector scripts.
    - KNLC_ROOT: Root location of knlc folder. Don't need to inlude knlc directory itself.
3. Run test:
    ```bash
    python recycler.py --skymapfilename=OUTPUT/TESTING/MS190327o/bayestar.fits.gz --triggerpath=OUTPUT/TESTING/ --triggerid=MS190327o
    ```

