import os
import sys
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
                            "handle",
                            "initial-plots",
                            "rank-galaxies",
                            "add-trigger",
                            "strategy",
                            "post-plots"
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

    if alert_type[0] == "S":
        EVENT_DIR = ROOT_DIR/"OUTPUT"/"O4REAL"/gw_id/alert_type
    else:
        EVENT_DIR = ROOT_DIR/"OUTPUT"/"TESTING"/gw_id/alert_type
    
    conf = ROOT_DIR/"configs"/"communications.yaml" 
    
    email = Email.from_config_file(conf)
    slack = Slack.from_config_file(conf)
    
    moc_skymap = EVENT_DIR/"bayestar.multiorder.fits"
    flatten_skymap = Path(moc_skymap.as_posix().replace(".multiorder.fits", ".fits.gz"))

    far = round(1/float(gw.event.far)/60/60/24/365, 2)
    source, event_prob = gw.find_source_and_prob()

    if "all" or "handle" in options.stages:
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
        gracedb = Path("https://gracedb.ligo.org")
        link = gracedb/"apiweb"/"superevents"/gw_id/"files"/"bayestar.png"
        
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
            f"Skymap Link: {link}")
        
        # Thresholds
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
            email.send_message(subject=subject, text=message)
            slack.send_message(subject=subject,text=message)
            logger.info(f"Messages about {alert_type} {gw_id} event ent!")
    
    if ("all" in options.stages) or ("initial-plots" in options.stages):
        logger.info(f"In stage initial-plots for {alert_type} {gw_id}")
        plots_path = flatten_skymap.parent

        make_plots_initial(url=flatten_skymap.as_posix(),
                           name=(plots_path/gw_id).as_posix())

        logger.info(f"Initial plots for {alert_type} {gw_id} stored at {plots_path}.")

    if ("all" in options.stages) or ("rank-galaxies" in options.stages):
        logger.info(f"Ranking galaxies for {alert_type} {gw_id} alert.")
        find_galaxy_list(
            map_path=flatten_skymap,
            event_name=gw.superevent_id,
            galax = ROOT_DIR/"data"/"franken_gals_array.npy",
            out_path=EVENT_DIR/"ranked_galaxies_list.csv"
        )

        logger.info(f"Ranked galaxies list stored at {os.path.dirname(flatten_skymap)}")

    if ("all" in options.stages) or ("add-trigger" in options.stages):

        logger.info("In stage add trigger")
        desgw = DESGWApi(os.environ.get("API_BASE_URL"))

        season = 1000  # Change to official season later
        server_dir = os.path.join("/des_web","www","html","desgw-new",f"dp{season}",f"{gw_id}")
        # Make a directory within the server hosting the webpage corresponding to the season and trigger_id
        os.system("ssh -k codemanager@desweb.fnal.gov 'mkdir -p {}'".format(server_dir))
        # Define a variable "desweb" to be the directory to store our plots/files in 
        desweb = "codemanager@desweb.fnal.gov:{}".format(server_dir)
        # Assuming we store all the relevant plots in our working directory, in a directory called "MI_plots", we just "scp" the entire directory into the location on the server 
        os.system(f"scp -r {plots_path}*.png {desweb}")
        trigger_data = {
            "trigger_label": gw_id,
            "type": source,
            "ligo_prob": event_prob,
            "far": far,
            "distance": distmean,
            "mjd": float(mjd),
            "event_datetime": Time(gw.event.time).strftime("%Y-%m-%d %H:%M:%S"),
            "mock": True, # True for mock event, False for real event.
            # "galaxy_percentage_file": f"https://des-ops.fnal.gov:8082/desgw-new/dp{season}/{gw_id}/{galaxy_name}", #output from galaxy ranking file. (csv filepath)
            "initial_skymap": f"https://des-ops.fnal.gov:8082/desgw-new/dp{season}/{gw_id}/S230518h_initial_skymap.png", # output initial skymap plot filepath
            "moon": f"https://des-ops.fnal.gov:8082/desgw-new/dp{season}/{gw_id}/{gw_id}_Moon.png",
            "season": "1000"
        }

        desgw.add_trigger(trigger_data = trigger_data)

    if ("all" in options.stages) or ("strategy" in options.stages):
        pass

    if ("all" in options.stages) or ("onering" in options.stages):
        pass
        # run_onering()