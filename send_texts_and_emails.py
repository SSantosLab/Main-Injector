import smtplib
import requests
import logging as log
import numpy as np
import re

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

FORMAT = '%(asctime)s %(message)s'
log.basicConfig(format=FORMAT)

def send_text(*args, **kwargs) -> None:
    pass

def postImagetoSLACK(impath,official=False):
    if not official:
        return
    post = {}
    post["value1"] = open(impath,'rb').read()
    post["value2"] = '' 
    r = requests.post('https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh',data=post)

def postToSLACK(subject,text,official=False,atchannel=False):

    if not official:
        return
    if atchannel:
        subject = '@channel ' + subject
        

    post = {}
    post["value1"] = subject
    post["value2"] = text
    requests.post("https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh", data=post)

def send(
    subject: str,
    text: str,
    official: bool = False,
 ) -> None:
    
    """
    Send emails to a list of emails with a given subject and text.
    if Official, the emails says that it's a OFICIAL LIGO ALERT.
    """

    email_args = {
        'from_addr': 'DESGW Team (TEST CASE)',
        'smtpserver': '%s:%s' % ('smtp.gmail.com', 587),
        'gmail_login': 'step.daily.report@gmail.com',
        'gmail_password': 'uxiywgbflgpqysah',
    }

    from_addr = email_args['from_addr']
    smtpserver = email_args['smtpserver']
    gmail_login = email_args['gmail_login']
    gmail_password = email_args['gmail_password']

    print('Preparing email')

    people = np.genfromtxt('DESGW_O4_People.TXT',
                            dtype=[('name','S50'),('email','S50'),('phone','S50')],
                            delimiter=",",
                            skip_header=1)
    emails=[]

    test_emails = ['andsouzasanttos@gmail.com',
                   'annis@fnal.gov',
                   'norafs@umich.edu']
    
    if not official:
        emails = test_emails
    
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
    
def send_emergencial_email():
    email_args = {
        'from_addr': 'DESGW Team (TEST CASE)',
        'smtpserver': '%s:%s' % ('smtp.gmail.com', 587),
        'gmail_login': 'step.daily.report@gmail.com',
        'gmail_password': 'uxiywgbflgpqysah',
    }

    from_addr = email_args['from_addr']
    smtpserver = email_args['smtpserver']
    gmail_login = email_args['gmail_login']
    gmail_password = email_args['gmail_password']

    emergency_list = {
        'names': ['Andre Santos', 'Nora Sherman', 'Jim Annis'],
        'emails': ['andsouzasanttos@gmail.com',
                    'norafs@umich.edu',
                    'annis@fnal.gov']
    }

    for email in emergency_list['emails']:
        msg = MIMEMultipart('alternative')
        msg['Subject'] = 'CHECKEVENT_KAFKA BROKE! RUN ALARMS!!'
        msg['From'] = from_addr
        msg['To'] = email
        payload = MIMEText('Check checkevent_kafka.py inside des machines.')
        msg.attach(payload)

        with smtplib.SMTP(smtpserver) as server:
            server.starttls()
            server.login(gmail_login, gmail_password)
            resp = server.sendmail(from_addr,[email],msg=msg.as_string())
            print('Send email success: {0}'.format(email)) 
    

    # try:
    #     # people = np.genfromtxt('/home/s1/desgw/PHONE_AND_EMAIL_LIST.TXT', dtype=[('name','S50'),('email','S50'),('phone','S50')], delimiter=",",skip_header=1)
    #     people = {
    #         'email': 'physics.programming@gmail.com'
    #     }
    # except:
    #     report = {}
    #     report["value1"] = '@channel WARNING, THE PHONE AND EMAIL LIST CRASHED... probably a formatting issue \n'
    #     report["value2"] = 'please fix asap.'
    #     requests.post("https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh", data=report)

    # print(people['name'])

    # me = 'automated-desGW@fnal.gov'
    # you = ['physics.programming@gmail.com']
    # if official:
    #     for email in people['email']:
    #         #if not p[0] == '#':
    #         you.append(email.replace('\t','').replace(' ',''))
    #else:
    #    you = ['djbrout@gmail.com']

    # if official:
    #     for phone in people['phone']:
    #         if not phone == 'none':
    #             #if not p[0] == '#':
    #             you.append(phone.replace('\t','').replace(' ',''))

    #     report = {}
    #     if atchannel:
    #         subject = '@channel '+subject
    #     report["value1"] = subject
    #     report["value2"] = text
    #     requests.post("https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh", data=report)            

    #else:
    #    you.append('2153008763@mms.att.net')

    # for y in you:
    #     try:
    #         msg['Subject'] = subject
    #         msg['From'] = me

    #         #msg['From'] = me
    #         msg['To'] = y

    #         s = smtplib.SMTP('localhost')
    #         s.sendmail(me, y, msg.as_string())
    #         s.quit()
    #     except:
    #         print('could not send alert! to ', y)

        

#send('This is just a test #6 DO NOT REPLY VIA EMAIL','LOOK AT SLACK PAGE, cross your fingers, should be instant...')

if __name__ == '__main__':
    send('Test', 'this is a test')
