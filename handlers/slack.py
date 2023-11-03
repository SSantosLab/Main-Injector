import os
import requests
import base64
import yaml
from yaml.loader import SafeLoader

class SlackBot():
    """Slack Bot Class to handle interationcs with Slack trough IFTTT"""

    def __init__(self, mode:bool = 'test') -> None:

        ROOT = os.path.abspath(__file__)
        ROOT = os.path.dirname(ROOT)
        ROOT = os.path.dirname(ROOT)
        self.CONFIG = os.path.join(ROOT, 'configs','communications.yaml')
        with open(self.CONFIG) as f:
            slack_config = yaml.load(f, Loader=SafeLoader)

        self.mode = mode
        self._link = slack_config['ifttt_link']

    def post_image(self, impath: str) -> None:
        """
        Post_image method. Post a image to slack trough POST request.

        Parameters:
        ----------
        impath: str
            image path.

        """

        post = {}
        post["value1"] = ''
        post["value2"] = open(impath,'rb').read()

        requests.post(self._link, data=post)


    def post_message(self, subject: str, text: str) -> None:
        """Post text to slack"""

        if self.mode == 'test':
            subject = 'OFFLINE TESTING NO PANICKING IS REQUIRED ' + subject

        if self.mode == 'mock':
            subject = 'GCN STREAM TEST ' + subject
            
        if self.mode == 'observation':
            subject = '@channel ' + subject
        
        post = {}
        post["value1"] = subject
        post["value2"] = text
        requests.post(self._link, data=post)
