import os
import sys
import multiprocessing
import pandas as pd
from argparse import ArgumentParser
from datetime import datetime
from loguru import logger
from astropy.table import Table
from astropy.time import Time
from stages.gw import GW
from stages.communication import Email, Slack
from stages.plots import make_plots_initial, make_plots_post_strat
from stages.api import DESGWApi
from stages.rank_galaxies import find_galaxy_list
from pathlib import Path
#from python.OneRing import run_or
from subprocess import run

def run_strategy(
        skymap_location: str or os.PathLike,
        time: float,
        output_dir: str or os.PathLike,
        sky_condition: str = "moony",
        event: str = "BNS",
        queue: multiprocessing.Queue = None
):
    """
    Run strategy code and return output csv file location.
    
    Parameters:
    -----------
    
    skymap_location: str or os.Pathlike
        input skymap location (flatten, constant nside).
    time: float
        Modified Julian Date (MJD) of event.
    output_dir: str or os.Pathlike
        output folder location for strategy
    sky_condition: str
        Estimate for sky conditions during observing night.
        Options are 'moony' and 'notmoony'.
    event: str
        Classification for event. can be either 'BNS' or 'NSBH'.
    """

    if event == 'BNS':
        kn_type = 'blue'

    if event == 'NSBH':
        kn_type = 'red'

    if event == 'BBH':
        logger.warning(f"Trying to run strategy for {event} kind. Skipping!")
        sys.exit()

    mjd = str(time).replace('.','')
    strategy_filename = f"bayestar_{sky_condition}_{kn_type}_{mjd}_allconfig.csv"
    strategy_output = os.path.join(output_dir, strategy_filename)
    cmd = f"python " +\
        "python/knlc/kn_strategy_sims_MI.py " +\
        f"--input-skymap {skymap_location} " +\
        f"--teff-type {sky_condition} " +\
        f"--kn-type {kn_type} " +\
        f"--time {mjd} " +\
        f"--output-strategy {strategy_output}"
        
    output_log = os.path.join(
        output_dir,
        f'{sky_condition}_{kn_type}_strategy.log'
    )

    strategy_log = open(output_log, 'w')
    logger.info(f'Using {kn_type} kn for {event} event with teff {sky_condition}')

    run(cmd,
        shell=True,
        stdout=strategy_log,
        stderr=strategy_log,
        text=True)
        
    strategy_log.close()

    return strategy_output


def parser() -> ArgumentParser:

    p = ArgumentParser()
    p.add_argument("-f", "--file",
                   help="gw alert content json file location.",
                   default=None)
    
    p.add_argument("-s", "--stages",
                   help="Stages to run for main injector. Default is all.",
                   default="all",
                   type=str,
                   nargs="+",
                   choices=["all",
                            "initial-plots",
                            "handle",
                            "rank-galaxies",
                            "add-trigger",
                            "add-trigger-by-day-initial",
                            "strategy",
                            "post-plots",
                            "onering"])
    
    p.add_argument("--skip-message",
                   help="Use this option to not sent notification to people.",
                   action="store_true",
                   default=False)
    
    p.add_argument("--max-hex-time",
                   type=float,
                   default=None,
                   nargs="?",
                   help="Limit the number of hexes in json file based on time in seconds. Default is None.")
    
    p.add_argument("--max-hex-count",
                   type=float,
                   default=None,
                   nargs="?",
                   help="Limit the number of hexes in json file in number of hexes. Default is None.")
    return p

# Entrypoint
if __name__ == "__main__":

    p = parser()
    options = p.parse_args()
    gw = GW.from_alert(options.file)
    alert_type = gw.alert_type
    gw_id = gw.superevent_id
    ROOT_DIR = Path(os.environ.get("ROOT_DIR"))
    DATA_DIR = ROOT_DIR/"data"
    LOG_DIR = ROOT_DIR/"logs"
    listener_log = LOG_DIR/"recyler.log"
    logger.add(listener_log, level="INFO", rotation="00:00")

    '''if gw_id[0] == "S": #previouslt was 'alert_type[0]'
        EVENT_DIR = ROOT_DIR/"OUTPUT"/"O4REAL_test"/gw_id/f'{alert_type}_0' #This will throw an error since there is PRELIMINARY_0 and PRELIMINARY_1
        #EVENT_DIR = ROOT_DIR/"OUTPUT"/"O4REAL"/gw_id/"PRELIMINARY_0" #This will throw an error since there is PRELIMINARY_0 and PRELIMINARY_1
    else:
        EVENT_DIR = ROOT_DIR/"OUTPUT"/"TESTING_test"/gw_id/alert_type'''

    if alert_type == 'PRELIMINARY':
        alert_type_codemanager = f'{alert_type}_0'
        if gw_id[0] == "S": #previouslt was 'alert_type[0]'
            EVENT_DIR = ROOT_DIR/"OUTPUT"/"O4REAL"/gw_id/f'{alert_type}_0'
            is_mock = False
        else:
            EVENT_DIR = ROOT_DIR/"OUTPUT"/"TESTING"/gw_id/f'{alert_type}_0'
            is_mock = True

        if os.path.isdir(EVENT_DIR):
            alert_type_codemanager = alert_type_codemanager.replace('_0', '_1')
            EVENT_DIR = EVENT_DIR.replace('_0','_1')
    else:
        alert_type_codemanager = alert_type
        if gw_id[0] == "S": #previouslt was 'alert_type[0]'
            EVENT_DIR = ROOT_DIR/"OUTPUT"/"O4REAL"/gw_id/alert_type
            is_mock = False
           
        else:
            EVENT_DIR = ROOT_DIR/"OUTPUT"/"TESTING"/gw_id/alert_type
            is_mock = True

    
    
    #conf = ROOT_DIR/"configs"/"communications.yaml" 
    
    #email = Email.from_config_file(conf)
    #slack = Slack.from_config_file(conf)
    
    moc_skymap = EVENT_DIR/"bayestar.multiorder.fits"
    flatten_skymap = Path(moc_skymap.as_posix().replace(".multiorder.fits", ".fits.gz"))
    plots_path = EVENT_DIR/"initial_plots" 

    far = round(1/float(gw.event.far)/60/60/24/365, 2)
    source, event_prob = gw.find_source_and_prob()

    if ("all" in options.stages) or ("initial-plots" in options.stages):
        logger.info(f"In stage initial-plots for {alert_type} {gw_id} alert.")
        #plots_path = EVENT_DIR/"initial_plots"
        #plots_path.mkdir(parents=True, exist_ok=True)
        plots_path = Path(os.path.join(EVENT_DIR, "initial_data"))
        plots_path.mkdir(parents=True, exist_ok=True)

        skymap_plot, moon_plot, area50, area90 = make_plots_initial(flatten_skymap.as_posix(), plots_path.as_posix())

        logger.info(f"rsync data from alert {gw.superevent_id} to website")
        desgw = DESGWApi(os.environ.get("API_BASE_URL"))

        # Make a directory within the server hosting the webpage corresponding to the gw_id
        server_dir = os.path.join("/des_web","www","html","desgw-new",f"{gw_id}", f"{alert_type_codemanager}")
        os.system("ssh -k codemanager@desweb.fnal.gov 'mkdir -p {}'".format(server_dir))
        # Define a variable "desweb" to be the directory to store our plots/files in 
        desweb = f"codemanager@desweb.fnal.gov:{server_dir}"
        # rsync the relevant plots to the website server 
        os.system(f"rsync -a {plots_path} {desweb}")

        logger.info(f"Initial plots for {alert_type} {gw_id} stored at {plots_path}.")

    if ("all" in options.stages) or ("handle" in options.stages):
        logger.info(f"Handling Event {gw_id} with alert type {alert_type}")

        if not moc_skymap.parent.exists():
            moc_skymap.parent.mkdir(parents=True, exist_ok=True)

        gw.retrieve_skymaps(moc_path=moc_skymap,
                            flatten_path=flatten_skymap,
                            nside=1024)
        
        logger.info(f"skymaps stored at {EVENT_DIR}.")
        
        skymap = Table.read(moc_skymap)
        distmean = skymap.meta["DISTMEAN"]
        distsigma = skymap.meta["DISTSTD"]
        mjd = round(Time(gw.event.time).mjd, 4)
        gracedb = "https://gracedb.ligo.org"
        link = os.path.join(gracedb,
                            "apiweb",
                            "superevents",
                            gw_id,
                            "files",
                            "bayestar.png")
        
        subject = f"Trigger {gw_id} MJD: {mjd} "+\
            f"Classification: {source} ({event_prob}) FAR: {far} "+\
            f"Alert Type {gw.alert_type}."
        
        message = (f"Alert Type: {gw.alert_type}\n"
            f"Superevent ID: {gw.superevent_id}\n"
            f"Event Time: {gw.event.time} \n"
            f"Alert Time: {gw.time_created}\n"
            f"MJD: {mjd}\n"
            f"FAR: {far}\n"
            f"Detectors: {gw.event.instruments}\n"
            f"BNS: {gw.event.classification.BNS:.3f}\n"
            f"NSBH: {gw.event.classification.NSBH:.3f}\n"
            f"BBH: {gw.event.classification.BBH:.3f}\n"
            f"Terrestrial: {gw.event.classification.Terrestrial}\n"
            f"Has Remmant: {gw.event.properties.HasRemnant}\n"
            f"DISTMEAN: {distmean:.2f} Mpc\n"
            f"DISTSIGMA: {distsigma:.2} Mpc\n"
            f"Has Mass Gap: {gw.event.properties.HasMassGap}\n"
            f"GraceDB Link: {gw.urls.gracedb}\n"
            f"Skymap Link: {link}\n"
            f"SLIP Moon plot: https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/initial_plots/Moon.png\n"
            f"SLIP Initial Skymap: https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/initial_plots/initial_skymap.png")
        
        if source == "BBH" and (far < 1000):
            logger.info("BBH with low FAR, skipping event!")
            sys.exit()

        if source == "Terrestrial":
            logger.info("Terrestrial source, skipping event!")
            sys.exit()

        if gw.event.group != "CBC":
            logger.info("gw alert without classification, skipping event!")
            sys.exit()

        if not options.skip_message:

            if (source == 'BNS') or (source == 'NSBH'):
                subject = '@channel ' + subject
            
            email.send_message(subject=subject, text=message)
            slack.send_message(subject=subject,text=message)
            logger.info(f"Messages about {alert_type} {gw_id} event sent!")
    

    if ("all" in options.stages) or ("rank-galaxies" in options.stages):
        logger.info(f"Ranking galaxies for {alert_type} {gw_id} alert.")
        catalog_path = EVENT_DIR/"initial_data"
        galaxy_list = find_galaxy_list(
            map_path=flatten_skymap.as_posix(),
            event_name=gw.superevent_id,
            galax = ROOT_DIR/"data"/"franken_gals_array.npy",
            out_path=catalog_path/"ranked_galaxies_list.csv"
        )

        os.system(f"rsync -a {galaxy_list} {desweb}")
        logger.info(f"Ranked galaxies list stored at {catalog_path}")

    if ("all" in options.stages) or ("add-trigger" in options.stages):

        logger.info(f"Pushing data from alert {gw.superevent_id} to website")
        desgw = DESGWApi(os.environ.get("API_BASE_URL"))

        season = -9  # Change to official season later

        trigger_data = {
            "trigger_label": gw_id,
            "mjd": float(mjd),
            "event_datetime": Time(gw.event.time).strftime("%Y-%m-%d %H:%M:%S"),
            "mock": is_mock, # True for mock event, False for real event.
            "detectors":  ' '.join(str(detector) for detector in gw.event.instruments), #list of the detectors
            "lvc_event_url": gw.urls.gracedb, 
            "season": season
            }

        desgw.add_trigger(trigger_data = trigger_data)

    if ("all" in options.stages) or ("add-trigger-by-day-initial" in options.stages):

        logger.info(f"rsync data from alert {gw.superevent_id} to website")
        desgw = DESGWApi(os.environ.get("API_BASE_URL"))

        
        trigger_data = {
            "type": source,
            "ligo_prob": event_prob,
            "far": far,
            "distance": distmean,
            "galaxy_percentage_file": f"https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/ranked_galaxies_list.csv", #output from galaxy ranking file. (csv filepath)
            "initial_skymap": f"https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/initial_skymap.png", # output initial skymap plot filepath
            "moon": f"https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/Moon.png",
            "season": "-9",
            "prob_region_50": area50,
            "prob_region_90": area90,
            "prob_coverage": 'Need to figure out what this is',
            "snr": 'Need to determine how to find this from FAR',
            "chirp_mass": 'Awaiting Isaacs code',
            "component_mass": 'Updated afer LVC releases these values'
        }
        logger.info(f"Pushing initial data from alert {gw.superevent_id} to website")
        desgw.add_trigger_by_day(trigger_data = trigger_data)

    if ("all" in options.stages) or ("strategy" in options.stages):

        strategies = multiprocessing.Queue()
        processes = []

        print(os.path.dirname(flatten_skymap))
        for teff_kind in ['moony', 'notmoony']:
            process = multiprocessing.Process(
                target=run_strategy,
                args=(
                    flatten_skymap,
                    mjd,
                    os.path.dirname(flatten_skymap),
                    teff_kind,
                    source,
                    strategies
                )
            )
            processes.append(process)
            process.start()

        for process in processes:
            process.join()

        result_strategy = []
        while not strategies.empty():
            result = strategies.get()
            result_strategy.append(result)
        
        print(result_strategy)

    if ("all" in options.stages) or ("post-plots" in options.stages):
        logger.info(f"In stage post-plots for {alert_type} {gw_id} alert.")
        airmass_hexes, skymap_obs_hexes, cum_hex_prob = make_plots_post_strat(flatten_skymap.as_posix(), plots_path.as_posix(), 'CHOOSE_JSON')
        logger.info(f"Post strategy plots for {alert_type} {gw_id} stored at {plots_path}.")

        # Make a directory within the server hosting the webpage corresponding to the gw_id
        server_dir = os.path.join("/des_web","www","html","desgw-new",f"{gw_id}")
        # Define a variable "desweb" to be the directory to store our plots/files in 
        desweb = "codemanager@desweb.fnal.gov:{}".format(server_dir)
        # scp the relevant plots to the website server 
        os.system(f"rsync -a {plots_path} {desweb}")


    if ("all" in options.stages) or ("add-trigger-by-day-final" in options.stages):

        logger.info(f"scp data from alert {gw.superevent_id} to website")

        #season = -9  # Change to official season later
        trigger_data = {
            "date": datetime.now()
            "n_hexes":
            "econ_prob":
            "econ_area":
            "need_area":
            "quality":
            "exp_time":
            "filter":
            "hours":
            "n_visits":
            "n_slots":
            "b_slot":
            "prob_vs_slot_prob":
            "centered_gif_plot":
            "ligo_prob_contour_plot":
            "des_prob_vs_ligo_prob_plot":
            "des_limit_mag_map":
            "des_limit_mag_map_src":
            "json_link":
            "log_link":
            "strategy_table":
            "final_skymap": f"https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/skymap_obs_hexes.png", 
            "airmass": f"https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/airmass_hexes.png", 
            "cumulative_hex_prob": f"https://des-ops.fnal.gov:8082/desgw-new/{gw_id}/{alert_type_codemanager}/cum_hex_prob.png"
            }

        desgw.add_trigger(trigger_data = trigger_data)





    # if ("all" in options.stages) or ("onering" in options.stages):
        
    #     for strategy in result_strategy:

    #         strategy = pd.read_csv()
    #         run_or(
    #             skymap=flatten_skymap,
    #             probArea_inner=0.5,
    #             probArea_outer=0.9,
    #             flt='i',
    #             expTime_inner=90,
    #             expTime_outer=90,
    #             mjd=60209.0859,
    #             hexFile='data/all-sky-hexCenters-decam.txt'
    #         )
    #     print('done')
    #     # run_onering()