
import sys
import getopt
import os
import yaml

from main_injector.events import event
try:
    args = sys.argv[1:]
    opt, arg = getopt.getopt(
        args, "tp:tid:mjd:exp:sky",
        longopts=["triggerpath=", "triggerid=", "mjd=", "exposure_length=", "official", "skymapfilename=", "hasrem"])

except getopt.GetoptError as err:
    print(str(err))
    print("Error : incorrect option or missing argument.")
    print(__doc__)
    sys.exit(1)

# Read in config
with open(
    os.path.join(os.environ["DESGW_CONFIG"], "recycler.yaml"), "r"
) as f:
    config = yaml.safe_load(f)
# Set defaults to config
trigger_path = config["trigger_path"]

real_or_sim = config["real_or_sim"]

official = False

if config["skymap_filename"] == 'Default':
    skymap_filename = None
else:
    skymap_filename = config["skymap_filename"]

trigger_ids = [config["trigger_id"]]

force_mjd = config["force_mjd"]

#exposure_length = config["exposure_length"]

# Override defaults with command line arguments
# THESE NOT GUARANTEED TO WORK EVER SINCE WE SWITCHED TO YAML
hasrem = False
#    norem = False

dontwrap = False
for o, a in opt:
    print('Option')
    print(o)
    print(a)
    print('-----')
    if o in ["-tp", "--triggerpath"]:
        trigger_path = str(a)
    elif o in ["-tid", "--triggerid"]:
        trigger_ids = [str(a)]
        dontwrap = True
    elif o in ["-mjd", "--mjd"]:
        mjd = float(a)
    # elif o in ["-exp","--exposure_length"]:
    #    exposure_length = float(a)
    elif o in ["-hours", "--hours_available"]:
        hours_available = float(a)
    elif o in ["-sky", "--skymapfilename"]:
        skymap_filename = str(a)
    elif o in ['--hasrem']:
        hasrem = True  # str(a)
        print("HASREM ", hasrem)
    elif o in ['--official']:
        official = True
#        elif o in ['--hasrem']:
#            hasrem = str(a)
    # elif o in ['--norem']:
    #    hasrem = False

    else:
        print("Warning: option", o, "with argument", a, "is not recognized")

# Clear bad triggers, only used for wrapping all triggers...
badtriggers = open('badtriggers.txt', 'w')
badtriggers.close()

####### BIG MONEY NO WHAMMIES ###############################################
if config["wrap_all_triggers"]:
    if not dontwrap:
        trigger_ids = os.listdir(trigger_path)
        trigger_ids = trigger_ids[2:]
for trigger_id in trigger_ids:
    if force_mjd:
        mjd = config["mjd"]
    else:
        try:
            mjd = open(os.path.join(trigger_path, trigger_id,
                        trigger_id + '_eventMJD.txt'), 'r').read()
        except:
            mjd = '99999'
    if skymap_filename is None:
        try:
            # if True:
            # mapname = open(os.path.join(trigger_path,
            #                            trigger_id,
            #                            config['default_map_name']), 'r').read()
            # skymap_filename = os.path.join(trigger_path,
            #                               trigger_id, config['default_map_name'])
            # print os.path.join(trigger_path, trigger_id,'default_skymap.txt')
            # print os.path.join(trigger_path, trigger_id,'default_skymap.txt').read()
            skymap_filename = os.path.join(trigger_path, trigger_id,
                                            open(os.path.join(trigger_path, trigger_id,
                                                                'default_skymap.txt'), 'r').read())
        except:
            badtriggers = open('badtriggers.txt', 'a')
            badtriggers.write(trigger_id + '\n')
            print('Could not find skymap url file')

    if 'bayestar' in skymap_filename:
        print('bayestar' * 50)

#        try:
    if 1 == 1:
        try:
            mjd = float(mjd)
        except:
            badtriggers = open('badtriggers.txt', 'a')
            badtriggers.write(trigger_id + '\n')
            print('WARNING: Could not convert mjd to float. Trigger: ' +
                    trigger_id + ' flagged as bad.')
# here is where the object is made, and parts of it are filed in
        master_dir = os.path.join(trigger_path, trigger_id)
        
        e = event.Event(skymap_filename,
                        master_dir,
                        trigger_id,
                        mjd,
                        config,
                        official,
                        hasrem)

# e has variables and code assocaiated with it. The mapMaker is called "e" or "self"

        e.mapMaker(trigger_id, skymap_filename, config, hasrem)
        e.getContours(config)
        #jsonfilelist = e.makeJSON(config)
        # e.make_cumulative_probs()
        os.system('cp '+e.event_paramfile+' '+master_dir)
        # generates the homepage
        e.updateTriggerIndex(real_or_sim=real_or_sim)
        # make a blank page with the basic info that is available
        e.updateWebpage(real_or_sim)
        # e.makeObservingPlots()
        # e.getContours(config)
        # e.send_nonurgent_Email()
        # e.updateWebpage(real_or_sim)

#        except KeyError:
#            print("Unexpected error:", sys.exc_info())
#            badtriggers = open('badtriggers.txt', 'a')
#            badtriggers.write(trigger_id + '\n')
#############################################################################

print('Done')
