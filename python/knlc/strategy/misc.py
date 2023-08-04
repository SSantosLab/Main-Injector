#import pyfits #deprecated lib --> new library is astropy.io.fits
#from sltools.catalog import ascii_data as asc;
#from sltools.image import imcp,segobjs;
from __future__ import print_function
import numpy
import os
#import time
#from sltools.image import add_noise
#from sltools.image import convolve
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import math
from math import sqrt
import datetime
import string

def is_X_running():
    """
    Function to check if X windows is running, adapted from:
    http://stackoverflow.com/questions/1027894/detect-if-x11-is-available-python

    Input:
     - NONE
    Output:
     - BOOLEAN : True if X window is running, False otherwise
    """
    import os
    if os.system("which xset") != 0:
        return False
    from subprocess import Popen, PIPE
    p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
    p.communicate()
    return p.returncode == 0

def get_string_file(string_test,search_dir,copy_dir=''):

    files=os.listdir(search_dir)
    string_test=string.rstrip(string_test,"\n")
    string_test=string.lstrip(string_test,"\n")
    matching = [s for s in files if string_test in s]
    if copy_dir!='':
        for i in xrange(0,len(matching)):
            os.system("cp "+matching+" "+copy_dir)

    return matching





def matplotlib_import():
    import matplotlib
    if not is_X_running():
       
       matplotlib.use('agg')

    return 0

def gauss(x, *p):
    A, mu, sigma = p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))

def gauss_norm(x, p):
    pi=3.141592653589793
    mu=p[0] 
    sigma= p[1]
    A=1/(sigma*sqrt(2*pi))
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))




def linear(x,*p):
    A, B = p
    return A*numpy.array(x)+B

def scale(x,*p):
    A = p
    return A*numpy.array(x)



def verify_kwarg(param_name,default_value,kwargs): 
    if param_name in kwargs.keys():
        param=kwargs[param_name]
    else:
        param=default_value
    return param

def quadrance(a, b):
    return math.sqrt(a**2+b**2)

def define_workdir(f):
    #f = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)
    return f
