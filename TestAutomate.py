import subprocess ,os
from threading import Thread
import make_recycler_config
import argparse
import numpy as np
import smtplib
from email.mime.text import MIMEText

############# 8/14/18 update #############
### This scirpt is still a work in    ###
### progress. To run as is, make sure ###
### to run in a the directory that    ### 
### Main-Injector, gw_workflow, ,and  ###
### Post-Processing live.             ###
#########################################

############# Send emails to appropriate people when things fail #############

def send_email(error, where):
    
    text = 'There was an error with ATC pipeline during %s, with message % ' % (where, error)
    msg = MIMEText(text)
    # me == the sender's email address
    # you == the recipient's email address
    me = 'alyssag94@brandeis.edu'
    you = 'alyssag94@brandeis.edu'
    msg['Subject'] = 'ATC pipeline error'
    msg['From'] = me            
    msg['To'] = you
    s = smtplib.SMTP('localhost')
    s.sendmail(me, [you], msg.as_string())
    print('An email was sent to %s' % you)
    s.quit()


############# Make yaml file for recycler #############

#information needed to make the .yaml config file for recycler 
parser = argparse.ArgumentParser()
parser.add_argument('--camera', choices=['decam', 'hsc'], default='decam', help="what camera did we use, default=decam")
parser.add_argument('--res', type=str, choices=[64, 128, 256], default=128,
                    help="what resolution do you want the map, default=128") 
parser.add_argument('--debug', type=str, choices=[True,False], default=False, help="turn debugging on/off")
parser.add_argument('--propid', default='2017B-0110', help='proposal id for this run')
args=parser.parse_args()

#makeYaml takes (camera, res, propid, sendEmail=(default False), sendTexts=(default False), debug=(default False))
yamlName= make_recycler_config.makeYaml(camera=args.camera, res=args.res, propid=args.propid, debug=args.debug)

#this is a hack to make sure the true/false statements are capitalized. 
os.system("sed -i -e 's/false/False/g' "+yamlName) 
os.system("sed -i -e 's/true/True/g' "+yamlName)

#need this to live in the production directory 
os.system(str("mv ")+yamlName+str(" Main-Injector/"))


############# Use to update current environment because subprocess is shit #############

def source(script, update=1):
    pipe = subprocess.Popen(". %s > /dev/null; env" % script, stdout=subprocess.PIPE, shell=True)
    data = pipe.communicate()[0]
    env={}
    for line in data.splitlines():
        #print(line)
        splt1=line.split('=',1)
        if splt1[0] == 'BASH_FUNC_setup()':
            splt1[1] = '() {  eval `$EUPS_DIR/bin/eups_setup           "$@"`  \n}'
            #print(splt1[1])
        if splt1[0] == 'BASH_FUNC_unsetup()':
            splt1[1] = '() {  eval `$EUPS_DIR/bin/eups_setup --unsetup "$@"`  \n}'
            #print(splt1[1])
        if splt1[0] == 'BASH_FUNC_module()':
            splt1[1]='() {  eval `/usr/bin/modulecmd bash $*`   \n}'
            #print(splt1[1])
        if splt1[0] =='}':
            continue

        #print(splt1)
        env[splt1[0]]=splt1[1]
        
    if update:
        os.environ.update(env)

    return env


############# Main Injector #############

#not sure what kind of trigger will come from listener but something needs to be triggered here to start
#Add the make_recycler_config.py here or at the end of main injector? 

mainoutfile= open('test_main.out', 'w')
mainerrfile = open('test_main.err','w')

source('Main-Injector/SOURCEME')
print("Environment successfully set up for Main Injector")


start_main = subprocess.Popen(['python', '/data/des41.a/data/desgw/alyssa_test/Main-Injector/recycler.py'], 
                        stderr=subprocess.PIPE, stdout=subprocess.PIPE, cwd='Main-Injector/')

main_out, main_err = start_main.communicate()

mainoutfile.write(main_out)
mainerrfile.write(main_err)
mainoutfile.close()
mainerrfile.close()

rc = start_main.returncode
#if rc != 0:
    
#    send_email(error, where)

print('The return code for recycler is '+str(rc))
print('')
print('')
print("Finished recycler for main injector. Visit website for skymaps")
print('')
print('')
print("Moving on to image processing ...")
print('')
print('')


############# Image Processing ################ 

### The bash script setup_img_proc.sh, should have everything needed for img proc
### need to add the make_dagrc.py to make the dagmaker.rc file


################# create new season number #######################
### Y6 will start with 600 (417 for mock)   ###
import easyaccess
import fitsio

query = 'SELECT max(SEASON) from marcelle.SNFAKEIMG where SEASON < 800;'
connection=easyaccess.connect('destest')
connection.query_and_save(query,'testfile.fits')
data = fitsio.read('testfile.fits')
print(data[0][0])

newseason = (int(data[0][0]/100) + 1)*100
print("the season number for this event is "+str(newseason))
print('')
print('')
print('')
############################################################

### make config file for dagmaker ###
###       THIS NEEDS HELP         ###

#import make_dagrc
#make_dagrc.makeDagRC(season=newseason)
#dagrc_name= make_dagrc.makeDagRC(seasonval=417) #for mock run
#os.system('mv dagmaker.rc gw_workflow/dagmaker.rc') 

############################################################

source('Post-Processing/diffimg_setup.sh')
### run DAGMaker for all new exposures, ie the exposures in curatedExposure.list
os.system('. ./seasonCycler.sh') #find new exposures and create list

source('gw_workflow/setup_img_proc.sh')
print("Environment successfully set up for Image processing.")
print('')
print('')

explist = np.genfromtxt('new_curated.list', delimiter=' ', usecols=0) #this will just be curatedExposure.list for production

imgprocout = open('test_imgproc.out', 'w')
imgprocerr = open('test_imgproc.err', 'w')
for i in explist:
    EXPNUM = int(i)
    check = os.path.isdir('gw_workflow/mytemp_'+str(EXPNUM))
    print('does '+str(EXPNUM)+' mytemp directory exist? '+str(check))
    print('')
    check = False #only for test runs
    if check == False:   
        img1 = subprocess.Popen(['bash','-c', '/data/des41.a/data/desgw/alyssa_test/gw_workflow/DAGMaker.sh '+str(EXPNUM)] ,stdout = subprocess.PIPE, stderr=subprocess.PIPE, cwd='gw_workflow/') 
        
        im1out, im1err = img1.communicate()
        imgprocout.write(im1out)
        imgprocerr.write(im1err)
        
        rc1 = img1.returncode
        print('The return code for DAGMaker is '+str(rc1))
        print('Finished ./DAGMaker for exposure '+str(EXPNUM))
        print('')

        img2 = subprocess.Popen(['jobsub_submit_dag','--role=DESGW', '-G', 'des', 
                                 'file://desgw_pipeline_'+str(EXPNUM)+'.dag'], 
                                stdout = subprocess.PIPE, stderr=subprocess.PIPE, cwd='gw_workflow/')

        im2out, im2err = img2.communicate()
        imgprocout.write(im2out)
        imgprocerr.write(im2err)
        
        rc2 = img2.returncode
        print('The return code for this jobsub is '+str(rc2))
        print('Finished jobsub_submit_dag for exposure '+str(EXPNUM))
        print('')
        print('Look at test_imgproc.out for the jobid')
    else:
        print('Already processed exposure number '+str(EXPNUM))
imgprocout.close()
imgprocerr.close()

print("Finished image processing! Moving on to post processing...")
print('')
print('')

############# Run Post Processing #############
############# we have a script that will make the postproc_seasonnumber.ini 
############# just need to know where to look for the ligo id 

postprocout = open('test_postproc.out', 'w')
postprocerr = open('test_postproc.err', 'w')
source('Post-Processing/diffimg_setup.sh')
#print("Environment successfully set up for post processing.")

start_pp = ['bash','-c', '/data/des41.a/data/desgw/alyssa_test/Post-Processing/seasonCycler.sh']
#start_pp = ['id','.', '/data/des41.a/data/desgw/alyssa_test/Post-Processing/seasonCycler.sh']
postproc = subprocess.Popen(start_pp, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='Post-Processing/')

postproc_out, postproc_err = postproc.communicate()
rc3 = postproc.returncode
print('the return code for post processing is '+str(rc3))
print('')
postprocout.write(postproc_out)
postprocerr.write(postproc_err)

postprocout.close()
postprocerr.close()
print("Finished Post-Processing! Visit website for more information")





