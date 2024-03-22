# Main Injector Docs

Main Injector is a python package designed by the SoaresSantosLab team in order to retrieve alerts from LVC GCN and work with those alerts to search for Electromagnetic Counterparts of Gravitational Wave Events. A full standard operating procedure can be found [here](https://seanmacb.notion.site/Main-Injector-Standard-OperAting-Procedure-MISOAP-d9795b62a2644bae8daf4f2c990b40cd). It includes information about procedures, data structures, FAQ, and more. It is a living document, so it will likely have the most up to date information about `Main-Injector`. For an older procedure, continue reading.


## Cloning `Main-Injector`

When cloning `Main-Injector`, there are three files that must be created for it to work properly.
1. `DESGW_O4_PEOPLE.txt` - this file is a .txt with two tab spacing, comma separated, containinig the name, email, and phone numbers of people to be alerted through email about GCN alerts. The header for the file is "#NAME,          EMAIL,          PHONE"  
2. `configs/slack_token.txt` - this file is a .csv with no spaces, with two values: "webhook,channel_name", where webhook is the Webhook URL set up within the Slack API for the bot to post in whatever channel you would like it to, and channel_name is the name of the channel in the form "#Channel-name" 
3. `configs/slack_image_creds.txt` - this file is a .csv with no spaces and two values: "OAuth_token,channel_id", where OAuth_token is the Bot User OAuth Token for the slackbot, and the channel_id is the ID of the slack channel you would like the bot to post images to. channel_id takes the form "C12A34BC5DEF".

## Running `Main-Injector`

Find the most up to date procedure [here](https://seanmacb.notion.site/Main-Injector-Standard-OperAting-Procedure-MISOAP-d9795b62a2644bae8daf4f2c990b40cd).
    
1. Edit SOURCEME file

    the SOURCEME file sets custom env variables that are necessary in order to recycler.py work fully. One can add specific modifications for your hostname. (If you don't know your hostname, type `echo $HOSTNAME` in a bash shell) 

2. source SOURCEME file.
    - ROOT_DIR: current Main-Injector root directory.
    - DATA_DIR: path to data directory containing files used by Main-Injector.
    - WEB: path to DES_GW_website.
    - PYTHONPATH: Location to custom Main-Injector scripts.
    - KNLC_ROOT: Root location of knlc folder. Don't need to inlude knlc directory itself.

The current workflow of Main-Injector is:

1. listener.py "listen" for LVK alerts from GW Events.

2. Once an alert is retrieved, it's sent to gwstreamer.py, responsable for processing it's content and send important information to slack, emails and cellphones. gwstreamer.py also and initiates recycler.py trough a child process. At this point, recyler.py and it's chield process is independent from gwstreamer.py and listener.py.

3. in recyler.py, the setup for strategy and onering is made, running two instances for 'moony' and 'notmoony' sky conditions.

3. Run tests:

    1. From listener.py
        ```bash
        python listener.py --mode test
        ```
    
    2. From recyler.py
    ```bash
    python recycler.py --triger-id MS181101ab --skymap OUTPUT/TESTING/MS181101ab/PRELIMINARY_0/bayestar.fits.gz --event BNS --official
    ```

    4. Run Strategy:
    ```bash
    python python/knlc/kn_strategy.py --input-skymap OUPUT/TESTING/MS181101ab/PRELIMINARY_0/bayestar.fits.gz --ouput OUPUT/TESTING/MS181101ab/PRELIMINARY_0 -teff-type moony --kn-type blue --time 584239324
    ```
    

