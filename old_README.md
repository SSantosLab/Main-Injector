# Main-Injector

As of June 2018, run as user gw in the gw account. 

We'd like to be able to run it in test mode as user XXX (say alenon, annis, brout, marcelle). This entails having the directory structure in the gw account being listed in the maininjector yaml configuration file.

```unix
% git clone git@github.com:SSantosLab/Main-Injector.git
% cd Main-Injector
```

You'll need to run the setup script each time you login
```
% source SOURCEME
```

Then you can run the recycler (this is automatically in test mode)
```
% python recycler.py
```

Example:
Before running tests:

Check inputs to recycler.yaml, specifically:

real_or_sim: 'sim' in test mode

kasen_fraction: default is 50. This means that if the apparent magnitude will not recover 50% of the Kasen models, the code will say there is no probability of detecting this event

days_since_burst: 0.5 is the default

exposure_length_Rem: defualt [60., 90.]

exposure_filter_Rem: default [ 'i', 'z']

exposure_length_BH : [ 90., ] #sec / exposure_filter_BH : [ 'r', ]

tilings - exposure_tiling_Rem and exposure_tiling_BH defualt [ 0, 1 ]

```
python recycler.py --skymapfilename=OUTPUT/TESTING/S190814bv/bayestar.fits.gz --triggerpath=OUTPUT/TESTING/ --triggerid=S190814bv
```
--skymapfilename: path to skymap. Code will read if its bayestar or LAL from filename

--triggerpath: path to output dir

--triggerid: name of event

optional: --hasrem: this will run recycler in "bright"/"has remnant". Omission of this argument will run the code in "dark"/"no remnant" mode


## Getting extra package

get pyslalib, python over industrial grade spherical astronomy code slalib
https://github.com/scottransom/pyslalib
```
wget  https://github.com/scottransom/pyslalib/archive/master.zip
unzip master.zip
cd pyslalib-master
make
python setup.py install --home=$WORK_DIR/python-home-stash/
python test/test_slalib.py
PYTHONPATH=$PYTHONPATH:$WORK_DIR/python-home-stash/ ;export PYTHONPATH
```