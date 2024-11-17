import smtplib
import os
import numpy as np
import yaml
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from yaml.loader import SafeLoader

class EmailBot():

    def __init__(self, mode) -> None:

        ROOT = os.path.abspath(__file__)
        ROOT = os.path.dirname(ROOT)
        ROOT = os.path.dirname(ROOT)

        self.ROOT = ROOT
        self.CONFIG = os.path.join(ROOT, 'configs','communications.yaml')
        with open(self.CONFIG) as f:
            email_config = yaml.load(f, Loader=SafeLoader)

        self.mode = mode
        self._email_login = email_config['email_login']
        self._email_auth = email_config['email_auth']

    def send_email(self, subject, text, emergency: bool = False) -> None:
        """
        Send emails to a list of emails with a given subject and text.
        """

        email_args = {
            'from_addr': 'DESGW Team',
            'smtpserver': '%s:%s' % ('smtp.gmail.com', 587),
            'gmail_login': self._email_login,
            'gmail_password': self._email_auth,
        }

        from_addr = email_args['from_addr']
        smtpserver = email_args['smtpserver']
        gmail_login = email_args['gmail_login']
        gmail_password = email_args['gmail_password']

        print('Preparing email')

        people = np.genfromtxt(f"{os.path.join(self.ROOT,'DESGW_O4_People.TXT')}",
                                dtype=[('name','S50'),('email','S50'),('phone','S50')],
                                delimiter=",",
                                skip_header=1)
        emails=[]

        test_emails = [
            'seanmacb@umich.edu',
            'imcmahon@umich.edu'
        ]
        
        if not emergency:

            if self.mode =='test':
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
                    resp = server.sendmail(from_addr=from_addr,
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
                    resp = server.sendmail(from_addr, phone, msg.as_string())
                    print('Send email success: {0}'.format(phone))

        if emergency:
        
            emergency_list = {
                'names': [
                    'Andre Santos',
                    'Nora Sherman',
                    'Jim Annis',
                    'Marcelle Soares-Santos',
                    'Sean MacBride',
                    'Isaac McMahon'
                ],
                'emails': [
                    'andsouzasanttos@gmail.com',
                    'norafs@umich.edu',
                    'annis@fnal.gov',
                    'mssantos@umich.edu',
                    'sean.p.macbride@gmail.com',
                    'imcmahon@umich.edu'
                ]
            }

            for email in emergency_list['emails']:
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
