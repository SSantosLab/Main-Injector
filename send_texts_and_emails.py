import smtplib
import requests
import logging as log
import numpy as np

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

FORMAT = '%(asctime)s %(message)s'
log.basicConfig(format=FORMAT)

def postImagetoSLACK(impath,official=False):
    if not official:
        return
    post = {}
    post["value1"] = open(impath,'rb').read()
    post["value2"] = '' 
    r = requests.post('https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh',data=post)

def postToSLACK(subject,text,official=False,atchannel=False):
    print(official)
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
        'subject': 'STEP Report'
    }

    from_addr = email_args['from_addr']
    smtpserver = email_args['smtpserver']
    gmail_login = email_args['gmail_login']
    gmail_password = email_args['gmail_password']

    print('Preparing email')

    msg = MIMEMultipart('alternative')
    msg['From'] = from_addr
    

    try: 
        people = np.genfromtxt('PHONE_AND_EMAIL_LIST.txt')
    except:
        people = {
            'names': ['Andre Santos'],
            'emails': ['andsouzasanttos@gmail.com']

        }
    
    with smtplib.SMTP(smtpserver) as server:
        server.starttls()
        server.login(gmail_login, gmail_password)

        for y in people.get('emails'):
            try:
                msg['Subject'] = subject
                msg['To'] = y
                payload = MIMEText(text)
                msg.attach(payload)
            except:
                log.info('Could not sent email to', y)

        resp = server.sendmail(from_addr, y, msg.as_string())
        print('Send email success: {0}'.format(y))


    # msg = MIMEText(text)
    

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