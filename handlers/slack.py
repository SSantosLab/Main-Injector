import os
import requests
import base64
import yaml
from yaml.loader import SafeLoader
import numpy as np
import json
import logging
from slack_sdk import WebClient

class SlackBot():
    """
    Slack Bot Class to handle interationcs with Slack

    contact @seanmacb if you have questions about the slackbot
    """

    def __init__(self, mode:bool = 'test') -> None:

        ROOT = os.path.abspath(__file__)
        ROOT = os.path.dirname(ROOT)
        ROOT = os.path.dirname(ROOT)
        self.CONFIG = os.path.join(ROOT, 'configs','slack_token.txt') # need to make a new file called 'slack_token.txt' with the oauth token, and include it in the .gitignore
        with open(self.CONFIG) as f:
            self.token,self.channel = np.loadtxt(self.CONFIG,dtype=str,delimiter=",",comments="\0") # Open and log the token and channel
        self.mode = mode
        
        self.image_config = os.path.join(ROOT,"configs","slack_image_creds.txt")

        with open(self.image_config) as f:
            self.image_token,self.image_channel =  np.loadtxt(self.image_config,dtype=str,delimiter=",")
        
        
            
    def post_image(self, impath: str,title:str,comment:str) -> None:
        """
        Post_image method. Post a image to slack trough POST request.

        Parameters:
        ----------
        impath: str
            image path.

        """

        client = WebClient(self.image_token)

        # new_file = client.files_upload_v2(title="TestFile",filename=impath,content="Placeholder Text")

        # file_url = new_file.get("file").get("permalink")
        return client.files_upload_v2(channel=self.image_channel,file=impath,title=title,initial_comment=comment)

        # maybe use the below?
        # upload_and_then_share_file = client.files_upload_v2(channel="C123456789",title="Test text data",filename="test.txt",content="Hi there! This is a text file!",initial_comment="Here is the file:")

        # post = {}
        # post["value1"] = ''
        # post["value2"] = open(impath,'rb').read()

        # requests.post(self._link, data=post)

        # my_file = {'file' : (impath, open(impath, 'rb'), 'png')}

        # payload={"filename":impath,"token":self.token,"channels":self.channel}

        # dicto = {'token':self.token,'filename': impath, 'channels': self.channel}
            
        # return requests.post("https://slack.com/api/files.upload", params=payload, files=my_file)


    def post_message(self, subject: str, text: str) -> None:
        """Post text to slack"""

        if self.mode == 'test':
            subject = '*OFFLINE TESTING NO PANICKING IS REQUIRED* \n\n' + subject

        elif self.mode == 'mock':
            subject = '*GCN STREAM TEST* \n\n' + subject

        elif self.mode == 'mock-bayestar':
            subject = '*MOCK BAYESTAR TEST* \n\n' + subject
            
        elif self.mode == 'observation':
            subject = '*<!channel>*\n\n' + subject
        elif self.mode == 'dummy':
            subject=""
        if text.__contains__("has been running nonstop"):
            subject=""

        text = subject + text
        
        slack_msg = {'text': text,'subject':subject}

        return requests.post(self.token,data=json.dumps(slack_msg))
