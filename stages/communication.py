import smtplib
import os
import requests
import numpy as np
import yaml
from abc import ABC, abstractmethod
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from yaml.loader import SafeLoader



class Bot(ABC):

    @abstractmethod
    def from_config_file(self):
        pass 

    @abstractmethod
    def send_message(self):
        pass


class Email(Bot):

    def __init__(self, mode: str, email_login: str, email_auth: str) -> None:

        self.mode = mode
        self._email_login = email_login
        self.__email_auth = email_auth
        
    @classmethod
    def from_config_file(cls, config_path: str) -> 'Email':
        """
        A class method that creates an instance of Email
        
        Parameters:
        -----------
            config_path (str):
                A string to set config path location as yaml file.
        
        Returns:
        --------
            Email: an instance of Email.
        """

        with open(config_path) as f:
            email_config = yaml.load(f, Loader=SafeLoader)
        
        email = cls(email_config['mode'],
                    email_config['email_login'],
                    email_config['email_auth'])
        
        return email
    

    def send_message(self,
                     subject: str,
                     text: str,
                     email_list: str = os.path.join(os.environ['ROOT_DIR'], 'DESGW_O4_People.TXT'),
                     emergency: bool = False) -> None:
        """
        Send emails to a list of emails with a given subject and text.

        Parameters:
        -----------
        subject (str):
            Email Subject.
        text (str):
            Email text content.
        email_list (str):
            Path to an txt file containing an valid email list.
        emergency (bool):
            If True, will sent Emergency email to a emergency email list.
            Default values is False.

        Returns:
        --------
            None
        """

        if self.mode == 'test':
            return
        
        email_args = {
            'from_addr': 'DESGW Team',
            'smtpserver': '%s:%s' % ('smtp.gmail.com', 587),
            'gmail_login': self._email_login,
            'gmail_password': self.__email_auth,
        }

        from_addr = email_args['from_addr']
        smtpserver = email_args['smtpserver']
        gmail_login = email_args['gmail_login']
        gmail_password = email_args['gmail_password']

        print('Preparing email')


        people = np.genfromtxt(email_list,
                                dtype=[('name','S50'),('email','S50'),('phone','S50')],
                                delimiter=",",
                                skip_header=1)
        emails=[]

        test_emails = [
            'andsouzasanttos@gmail.com',
            'norafs@umich.edu'
        ]
        
        emergency_emails = [
            'andsouzasanttos@gmail.com',
            'norafs@umich.edu',
            'annis@fnal.gov',
            'mssantos@umich.edu'
        ]

        if not emergency:

            if self.mode == 'test':
                emails = test_emails
                subject = 'OFFLINE TESTING ' + subject
            else:
                for email in people['email']:
                    email = email.decode('utf-8')
                    emails.append(email.replace('\t','').replace(' ',''))

            for email in emails:
                msg = MIMEMultipart('alternative')
                msg['Subject'] = subject
                msg['From'] = from_addr
                msg['To'] = email
                payload = MIMEText(text)
                msg.attach(payload)

                with smtplib.SMTP(smtpserver) as server:
                    server.starttls()
                    server.login(gmail_login, gmail_password)
                    server.sendmail(from_addr=from_addr,
                                    to_addrs=[email],
                                    msg=msg.as_string())
                    
                    print('Send email success: {0}'.format(email))      

            phone_numbers = []
            for phone in people['phone']:
                    phone = phone.decode('utf-8')
                    phone = phone.replace('\t','').replace(' ','')
                    if phone != 'None':
                        phone_numbers.append(phone)

            for phone in phone_numbers:
                msg = MIMEMultipart('alternative')
                msg['Subject'] = subject
                msg['To'] = phone
                payload = MIMEText(text)
                msg.attach(payload)

                with smtplib.SMTP(smtpserver) as server:
                    server.starttls()
                    server.login(gmail_login, gmail_password)
                    server.sendmail(from_addr, phone, msg.as_string())
                    print('Send email success: {0}'.format(phone))

        if emergency:

            for email in emergency_emails:
                msg = MIMEMultipart('alternative')
                msg['Subject'] = 'LISTENER BROKE! RUN ALARMS!!'
                msg['From'] = from_addr
                msg['To'] = email
                payload = MIMEText('Check listener.py inside des machines.')
                msg.attach(payload)

                with smtplib.SMTP(smtpserver) as server:
                    server.starttls()
                    server.login(gmail_login, gmail_password)
                    server.sendmail(from_addr,[email],msg=msg.as_string())
                    print('Send email success: {0}'.format(email)) 

class Slack(Bot):

    def __init__(self, mode: str, ifttt_link: str) -> None:

        self.mode = mode
        self._link = ifttt_link

    @classmethod
    def from_config_file(cls, config_path: str) -> 'Slack':
        """
        A class method that creates an instance of Slack
        
        Parameters:
        -----------
            config_path (str):
                A string to set config path location as yaml file.
        
        Returns:
        --------
            Email:
                an instance of Email.
        """

        with open(config_path) as f:
            config = yaml.load(f, Loader=SafeLoader)
        

        slack = cls(config['mode'],
                    config['ifttt_link'])

        return slack
    
    def send_message(self, subject: str, text: str) -> None:
        """send message to slack workspace trough IFTTT Bot.
        
        Parameters:
        -----------
        subject (str):
            Slack Subject.
        text (str):
            Slack text content.

        Returns:
        --------
            None.
        """

        if self.mode == 'test':
            subject = 'OFFLINE TESTING ' + subject

        if self.mode == 'observation':
            subject = '@channel ' + subject

        post = {}
        post["value1"] = subject
        post["value2"] = text
        requests.post(self._link, data=post)
