from __future__ import print_function

#from io import open

import matplotlib
matplotlib.use('agg')
from matplotlib import rc
import matplotlib.patches as patches
#from sltools.image import sextractor as se
#import pyfits as fits #deprecated lib --> new library is astropy.io.fits
from astropy.io import fits
from astropy.io.fits import Column
#from sltools.catalog import ascii_data as asc;
#from sltools.image import imcp,segobjs;
import numpy
import numpy as np
import os
#import time
import matplotlib.pyplot as plt
import math
import datetime
from misc import verify_kwarg, quadrance
#from cbomcode.image.image import image_hist_plot
from astropy.coordinates import SkyCoord,Angle
import astropy.coordinates as coord
from astropy import units as u
import string
from math import sqrt
from cycler import cycler
import collections
import scipy.stats #import norm
from scipy.interpolate import interpn
import json
#import mocpy
from mocpy import MOC
from mocpy import WCS
import astropy_healpix as ah
import healpy as hp


def baye2moc(bay_file,prob_level=0.9):

    # check it out https://emfollow.docs.ligo.org/userguide/tutorial/skymaps.html

    #works with bayestar.multiorder.fits
    uniq,probdensity,data,nested=openbayeMulti(bay_file)
    if nested==True:
        return [],nested
    #with fits.open(bay_file) as hdul:
    #    hdul.info()
    #    hdul[1].columns

    #    data = hdul[1].data

    #uniq=data['UNIQ']
    #probdensity=data['PROBDENSITY']


    level, ipix = ah.uniq_to_level_ipix(uniq)
    area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)

    prob = probdensity * area
    try:
        moc_out=MOC.from_valued_healpix_cells(uniq, prob, cumul_to=prob_level)
    except:  
        print ('converting uniq to np.int64')
        moc_out=MOC.from_valued_healpix_cells(uniq.astype('int64'), prob, cumul_to=prob_level)
        
    return moc_out,nested

def check_inside_MOC(ras,decs,moc_data):
    return moc_data.contains(ra=ras*u.degree,dec=decs*u.degree)



def desifiles2vertices(path_name):
    
    filenames=os.listdir(path_name)

    for i in range(0,len(filenames)):
        name=filenames[i].split("/")[-1]

    return 0



def squarecut_catalogs(filenames=None,ralim=[148.0,152.0],declim=[0.2,4.2], ra_=None,dec_=None, colnames=['ra','dec']):
    mask=[]
    if (type(filenames)==type(None)) and (type(ra_)!=type(None)) and (type(dec_)!=type(None)) :
        filenames=[1]
        input_array=True 
    elif type(filenames)!=type(None):
        input_array=False
    else:
        print(' Error: Need valid file name or ra and dec arrays')
        return 0
    for i in range(0,len(filenames)):
        if input_array==False:
            cat_=open_fits_catalog(filenames[i])
            ra_=fits_catalog_col(cat_,col_name=colnames[0])
            dec_= fits_catalog_col(cat_,col_name=colnames[1])
        for j in range(0,len(ra_)):
            if (ra_[j]>=ralim[0] and ra_[j]<ralim[1]) and (dec_[j]>declim[0] and dec_[j]<declim[1]):
                mask.append(True)
            else:
                mask.append(False)
    return mask   
    

def check_legacy_files(filenames,ralim,declim, check_moc=None):
    """
    brick names (<brick>) have the format <AAAa>c<BBB> where A, a and B are digits and c is either the letter m or p (e.g. 1126p222). The names are derived from    the (RA, Dec) center of the brick. The first four digits are int(RA×10)
, followed by p to denote positive Dec or m to denote negative Dec ("plus"/"minus"), followed by three digits of int(Dec×10)
. For example the case 1126p222 corresponds to (RA, Dec) = (112.6, +22.2).
<brickmin> and <brickmax> denote the corners of a rectangle in (RA, Dec). Explicitly, <brickmin> has the format <AAA>c<BBB> where <AAA> denotes three digits of the minimum int(RA)
 in degrees, <BBB> denotes three digits of the minimum int(Dec)
 in degrees, and c uses the p/m ("plus"/"minus") format outlined in the previous bullet point. The convention is similar for <brickmax> and the maximum RA and Dec. For example 000m010-010m005 would correspond to a survey region limited by 0∘≤RA<10∘
 and −10≤Dec<−5.
    sub-directories are listed by the RA of the brick center, and sub-directory names (<AAA>) correspond to RA. For example 002 corresponds to brick centers between an RA of 2 and an RA of 3.
    <filter> denotes the g, r, or z
    band, using the corresponding letter.
    Note that it is not possible to go from a brick name back to an exact (RA, Dec) center (the bricks are not on 0.1° grid lines). The exact brick center for a given brick name can be derived from columns in the survey-bricks.fits.gz file (i.e. brickname, ra, dec).`
    """
    #sweep-350p000-360p005.fits
    roi_files=[]
    for i in range(0,len(filenames)):
        filename=filenames[i].split("/")[-1]
        filepath='/'.join(filename[:-1]) 
        array_name=filename.split("-")
        begin=array_name[1]
        end=array_name[2]
        end_aux=end.replace("p"," ") 
        begin_aux=begin.replace("p"," ")
        if begin_aux==begin:
            ra_begin=float(begin.split("m")[0])
            dec_begin=float(("-")+begin.split("m")[1])
        else:
            ra_begin=float(begin.split("p")[0])
            dec_begin=float(begin.split("p")[1])
        if end_aux==end:
            ra_end=float(end.split("m")[0])
            dec_end=float(("-")+(end.split("m")[1]).rstrip(".fits"))
        else:
            ra_end=float(end.split("p")[0])
            dec_end=float((end.split("p")[1]).rstrip(".fits"))

        if type(check_moc)==type(None):
            dx = min(ra_end, ralim[1]) - max(ra_begin, ralim[0])
            dy = min(dec_end, declim[1]) - max(dec_begin, declim[0])
        else:
            inter_size,inter_moc,moc_sq=check_intersection_square(ras=[ra_begin,ra_end,ra_end,ra_begin],decs=[dec_begin,dec_begin,dec_end,dec_end],moc=check_moc)
            if inter_size>0:
                dx=1
                dy=1
            else:
                dx=0
                dy=0
          
        if (dx>0) and (dy>0):  
            #print ("intersection found")
            #print(ralim)
            #print (declim)
            #print([ra_begin,ra_end])
            #print([dec_begin,dec_end]) 
            roi_files.append(filenames[i])
        

    return roi_files


def get_legacy_photz(filenames=None,moc=None,phz_path="",phot_err_lim=0.1,magr_lim=22.87, return_all_mags=False,cat_in=None,cat_pz=None, select_spec=False, correct_phot=True,fluxes=False, sel_ra=None,sel_dec=None, get_dist_selradec=False,remove_stars=True, color_cuts=True, replace_zphot=True, error_spec=0.0005, star_column=False,color_column=False):
    """
    TYPE !=PSF and (z_phot_u68-z_phot_l68)/2 < 0.1 and mag_i<21.
    magab= 22.5-2.5log10(flux)
    limiting magnitude was taken from https://iopscience.iop.org/article/10.3847/1538-3881/ab089d/pdf (table 4)
    """
    if type(filenames)==type(None) and type(cat_in)==type(None):
        print("Error: you need to define a catalog, exting function get_legacy_photo_z")
        return 0
    elif type(filenames)==type(None) and type(cat_in)!=type(None):
        filenames=[1]
    else:
        pass
    #if type(phz_path)==type(None):
    #    phz_path=
    _ids=[]
    _ras=[]
    _decs=[]
    _photz=[]
    _magnitude_r=[]
    _photz_err=[]

    mags_r=[]
    mags_g=[]
    mags_z=[]
    mags_w1=[]
    mags_w2=[]
    mags_w3=[]
    mags_w4=[]

    magserr_r=[]
    magserr_g=[]
    magserr_z=[]
    magserr_w1=[]
    magserr_w2=[]
    magserr_w3=[]
    magserr_w4=[]
    _photz_min=[]
    _photz_max=[]
    _zspec=[]
    _training=[]
    _mask_star=[]
    _mask_color=[] 

    stars_removed=0
    stars_removed_wzspec=0
    color_removed=0
    sepsamin=[]    
    for i in range(0,len(filenames)):#len(filenames)
        if type(cat_in)==type(None):
            cat_,names=open_fits_catalog(filenames[i])
            print('working on ',filenames[i],'. This is the ',i,'th file of', len(filenames))
        else:
            cat_=cat_in
        ra_psf=fits_catalog_col(cat_,col_name='RA')
        #cat_=select_subset_fits_catalog(cat_,field_name='TYPE',field_value='PSF',condition="!=")
        ra_=fits_catalog_col(cat_,col_name='RA')
        Type_=fits_catalog_col(cat_,col_name='TYPE')
        #print(len(ra_psf)-len(ra_),' objects were removed out of ',len(ra_psf))
        dec_=fits_catalog_col(cat_,col_name='DEC')
        flux_=fits_catalog_col(cat_,col_name='FLUX_R')
        if return_all_mags==True:
            flux_r=fits_catalog_col(cat_,col_name='FLUX_R')
            flux_g=fits_catalog_col(cat_,col_name='FLUX_G')
            flux_z=fits_catalog_col(cat_,col_name='FLUX_Z')
            flux_w1=fits_catalog_col(cat_,col_name='FLUX_W1')
            flux_w2=fits_catalog_col(cat_,col_name='FLUX_W2')
            flux_w3=fits_catalog_col(cat_,col_name='FLUX_W3')
            flux_w4=fits_catalog_col(cat_,col_name='FLUX_W4')

            flux_ivar_r=fits_catalog_col(cat_,col_name='FLUX_IVAR_R')
            flux_ivar_g=fits_catalog_col(cat_,col_name='FLUX_IVAR_G')
            flux_ivar_z=fits_catalog_col(cat_,col_name='FLUX_IVAR_Z')
            flux_ivar_w1=fits_catalog_col(cat_,col_name='FLUX_IVAR_W1')
            flux_ivar_w2=fits_catalog_col(cat_,col_name='FLUX_IVAR_W2')
            flux_ivar_w3=fits_catalog_col(cat_,col_name='FLUX_IVAR_W3')
            flux_ivar_w4=fits_catalog_col(cat_,col_name='FLUX_IVAR_W4')
            gaia_g=fits_catalog_col(cat_,col_name='GAIA_PHOT_G_MEAN_MAG')
            if correct_phot==True:
                mwtrans_r=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_R')
                mwtrans_g=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_G')
                mwtrans_z=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_Z')
                mwtrans_w1=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_W1')
                mwtrans_w2=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_W2')
                mwtrans_w3=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_W3')
                mwtrans_w4=fits_catalog_col(cat_,col_name='MW_TRANSMISSION_W4')
                #flux_rraw=flux_r.copy()
                flux_r=np.divide(flux_r,mwtrans_r)#fits_catalog_col(cat_,col_name='FLUX_R')
                flux_g=np.divide(flux_g,mwtrans_g)
                flux_z=np.divide(flux_z,mwtrans_z)
                flux_w1=np.divide(flux_w1,mwtrans_w1)
                flux_w2=np.divide(flux_w2,mwtrans_w2)
                flux_w3=np.divide(flux_w3,mwtrans_w3)
                flux_w4=np.divide(flux_w4,mwtrans_w4)

            if fluxes==False:
                magerr_r=1.086/np.array(flux_r*np.sqrt(flux_ivar_r)) #flux * sqrt(IVAR) = S/N 
                #https://www.sdss.org/dr14/manga/manga-tutorials/manga-faq/#WhydoyououtputIVAR(inversevariance)insteadoferrors?
                magerr_g=1.086/np.array(flux_g*np.sqrt(flux_ivar_g))
                magerr_z=1.086/np.array(flux_z*np.sqrt(flux_ivar_z))
                magerr_w1=1.086/np.array(flux_w1*np.sqrt(flux_ivar_w1)) 
                magerr_w2=1.086/np.array(flux_w2*np.sqrt(flux_ivar_w2))
                magerr_w3=1.086/np.array(flux_w3*np.sqrt(flux_ivar_w3))
                magerr_w4=1.086/np.array(flux_w4*np.sqrt(flux_ivar_w4))

                mag_r=22.5-(2.5*np.log10(flux_r))
                #mag_rraw=22.5-(2.5*np.log10(flux_rraw))
                mag_g=22.5-(2.5*np.log10(flux_g))
                mag_z=22.5-(2.5*np.log10(flux_z))
                mag_w1=22.5-(2.5*np.log10(flux_w1))
                mag_w2=22.5-(2.5*np.log10(flux_w2))
                mag_w3=22.5-(2.5*np.log10(flux_w3))
                mag_w4=22.5-(2.5*np.log10(flux_w4))
            else:
                magerr_r=np.array(flux_ivar_r) 
                magerr_g=np.array(flux_ivar_g)
                magerr_z=np.array(flux_ivar_z)
                magerr_w1=np.array(flux_ivar_w1) 
                magerr_w2=np.array(flux_ivar_w2)
                magerr_w3=np.array(flux_ivar_w3)
                magerr_w4=np.array(flux_ivar_w4)

                mag_r=np.array(flux_r)
                mag_g=np.array(flux_g)
                mag_z=np.array(flux_z)
                mag_w1=np.array(flux_w1)
                mag_w2=np.array(flux_w2)
                mag_w3=np.array(flux_w3)
                mag_w4=np.array(flux_w4)

        if type(moc)!=type(None):
            mask=check_inside_MOC(ras=ra_,decs=dec_,moc_data=moc)
        if (type(sel_ra)!=type(None)):
            #print (sel_ra[i])
            #print (sel_dec[i])
            #print(filenames[i])
            #print (ra_[:20])
            mask= np.where((ra_ > sel_ra[i][0]) & (ra_ < sel_ra[i][1]) & (dec_>sel_dec[i][0]) & (dec_<sel_dec[i][1]))
        #print ('len of ra')
        #print(len(ra_))
        if type(moc)!=type(None) or (type(sel_ra)!=type(None)):
            
            ra_=ra_[mask] 
            dec_=dec_[mask]
            flux_=flux_[mask]
            Type_=Type_[mask] 
            if return_all_mags==True:
                mag_r=mag_r[mask]
                mag_g=mag_g[mask]
                mag_z=mag_z[mask]
                mag_w1=mag_w1[mask]
                mag_w2=mag_w2[mask]
                mag_w3=mag_w3[mask]
                mag_w4=mag_w4[mask]
                magerr_r=magerr_r[mask]
                magerr_g=magerr_g[mask]
                magerr_z=magerr_z[mask]
                magerr_w1=magerr_w1[mask]
                magerr_w2=magerr_w2[mask]
                magerr_w3=magerr_w3[mask]
                magerr_w4=magerr_w4[mask]
                gaia_g=gaia_g[mask]
        #print ('len of ra after mask')
        #print(len(ra_))

 




        if len(ra_)>0:
            if type(cat_in)==type(None):
                filenamepz=filenames[i].split("/")[-1]
                #print(type(phz_path))
                #print(type("example str"))
                #print(type(phz_path)!=type("example str"))
                #print(isarraylike(phz_path))
                if type(phz_path)!=type("example str"):
                    #print ("in the wrong place") 
                    filenamepz=phz_path[i]+filenamepz.split(".")[0]+"-pz.fits" 
                else:
                    #print ("in the right place")
                    filenamepz=phz_path+filenamepz.split(".")[0]+"-pz.fits"
                #print (phz_path)
                #print (filenamepz)
                cat_pz,names=open_fits_catalog(filenamepz)
            phz_=fits_catalog_col(cat_pz,col_name='z_phot_median')
            phz_low=fits_catalog_col(cat_pz,col_name='z_phot_l68')
            phz_upper=fits_catalog_col(cat_pz,col_name='z_phot_u68')
            phz_err= (phz_upper.astype('float')-phz_low.astype('float'))/2.0

            specz_=fits_catalog_col(cat_pz,col_name='z_spec')
            training_=fits_catalog_col(cat_pz,col_name='training')

            if (type(moc)!=type(None)) or (type(sel_ra)!=type(None)):
                phz_=phz_[mask]
                phz_low=phz_low[mask]
                phz_upper=phz_upper[mask]
                phz_err=phz_err[mask]
                training_= training_[mask]
                specz_=specz_[mask]
            #remove PSF-like objects
            #print(Type_)
            #print('Remove PSF-LIKE')
            # Select objects with spectroscopy only



            phz_=phz_[Type_.astype('str') !='PSF']#select_subset_arr(Type_.astype('str'),val='PSF',condition="!=",into=phz_)
            ra_=ra_[Type_.astype('str') !='PSF']#select_subset_arr(Type_.astype('str'),val='PSF',condition="!=",into=ra_)
            dec_=dec_[Type_.astype('str') !='PSF']#select_subset_arr(Type_.astype('str'),val='PSF',condition="!=",into=dec_)
            phz_err=phz_err[Type_.astype('str') !='PSF']#select_subset_arr(Type_.astype('str'),val='PSF',condition="!=",into=phz_err)
            flux_=flux_[Type_.astype('str') !='PSF']#select_subset_arr(Type_.astype('str'),val='PSF',condition="!=",into=flux_)
            #print(len(ra_psf)-len(ra_),' PSF-like objects in legacy survey were removed out of ',len(ra_psf) )
            training_=training_[Type_.astype('str') !='PSF']
            specz_=specz_[Type_.astype('str') !='PSF']


            if return_all_mags==True:
                mag_r=mag_r[Type_.astype('str') !='PSF']
                mag_g=mag_g[Type_.astype('str') !='PSF']
                mag_z=mag_z[Type_.astype('str') !='PSF']
                mag_w1=mag_w1[Type_.astype('str') !='PSF']
                mag_w2=mag_w2[Type_.astype('str') !='PSF']
                mag_w3=mag_w3[Type_.astype('str') !='PSF']
                mag_w4=mag_w4[Type_.astype('str') !='PSF']

                magerr_r=magerr_r[Type_.astype('str') !='PSF']
                magerr_g=magerr_g[Type_.astype('str') !='PSF']
                magerr_z=magerr_z[Type_.astype('str') !='PSF']
                magerr_w1=magerr_w1[Type_.astype('str') !='PSF']
                magerr_w2=magerr_w2[Type_.astype('str') !='PSF']
                magerr_w3=magerr_w3[Type_.astype('str') !='PSF']
                magerr_w4=magerr_w4[Type_.astype('str') !='PSF']
                phz_low=phz_low[Type_.astype('str') !='PSF']
                phz_upper=phz_upper[Type_.astype('str') !='PSF']
                gaia_g=gaia_g[Type_.astype('str') !='PSF'] 

            if select_spec==True:
                phz_=select_subset_arr(specz_,val=0.0,condition=">",into=phz_)
                ra_=select_subset_arr(specz_,val=0.0,condition=">",into=ra_)
                dec_=select_subset_arr(specz_,val=0.0,condition=">",into=dec_)
                phz_err= select_subset_arr(specz_,val=0.0,condition=">",into=phz_err)
                flux_=select_subset_arr(specz_,val=0.0,condition=">",into=flux_)
                if return_all_mags==True:
                    mag_r=select_subset_arr(specz_,val=0.0,condition=">",into=mag_r)
                    mag_g=select_subset_arr(specz_,val=0.0,condition=">",into=mag_g)
                    mag_z=select_subset_arr(specz_,val=0.0,condition=">",into=mag_z)
                    mag_w1=select_subset_arr(specz_,val=0.0,condition=">",into=mag_w1)
                    mag_w2=select_subset_arr(specz_,val=0.0,condition=">",into=mag_w2)
                    mag_w3=select_subset_arr(specz_,val=0.0,condition=">",into=mag_w3)
                    mag_w4=select_subset_arr(specz_,val=0.0,condition=">",into=mag_w4)

                    magerr_r=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_r)
                    magerr_g=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_g)
                    magerr_z=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_z)
                    magerr_w1=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_w1)
                    magerr_w2=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_w2)
                    magerr_w3=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_w3)
                    magerr_w4=select_subset_arr(specz_,val=0.0,condition=">",into=magerr_w4)

                    phz_low=select_subset_arr(specz_,val=0.0,condition=">",into=phz_low)
                    phz_upper=select_subset_arr(specz_,val=0.0,condition=">",into=phz_upper)
                    gaia_g=select_subset_arr(specz_,val=0.0,condition=">",into=gaia_g)
                training_=select_subset_arr(specz_,val=0.0,condition=">",into=training_)
                specz_=select_subset_arr(specz_,val=0.0,condition=">",into=specz_)







            #remove negative fluxes
            #print('Remove Negative Fluxes')
            phz_=select_subset_arr(flux_,val=0.0,condition=">",into=phz_)
            ra_=select_subset_arr(flux_,val=0,condition=">",into=ra_)
            dec_=select_subset_arr(flux_,val=0,condition=">",into=dec_)
            phz_err=select_subset_arr(flux_,val=0,condition=">",into=phz_err)

            training_=select_subset_arr(flux_,val=0,condition=">",into=training_)
            specz_=select_subset_arr(flux_,val=0,condition=">",into=specz_)

            if return_all_mags==True:
                mag_r=select_subset_arr(flux_,val=0,condition=">",into=mag_r)
                mag_g=select_subset_arr(flux_,val=0,condition=">",into=mag_g)
                mag_z=select_subset_arr(flux_,val=0,condition=">",into=mag_z)
                mag_w1=select_subset_arr(flux_,val=0,condition=">",into=mag_w1)
                mag_w2=select_subset_arr(flux_,val=0,condition=">",into=mag_w2)
                mag_w3=select_subset_arr(flux_,val=0,condition=">",into=mag_w3)
                mag_w4=select_subset_arr(flux_,val=0,condition=">",into=mag_w4)

                magerr_r=select_subset_arr(flux_,val=0,condition=">",into=magerr_r)
                magerr_g=select_subset_arr(flux_,val=0,condition=">",into=magerr_g)
                magerr_z=select_subset_arr(flux_,val=0,condition=">",into=magerr_z)
                magerr_w1=select_subset_arr(flux_,val=0,condition=">",into=magerr_w1)
                magerr_w2=select_subset_arr(flux_,val=0,condition=">",into=magerr_w2)
                magerr_w3=select_subset_arr(flux_,val=0,condition=">",into=magerr_w3)
                magerr_w4=select_subset_arr(flux_,val=0,condition=">",into=magerr_w4)

                phz_low=select_subset_arr(flux_,val=0,condition=">",into=phz_low)
                phz_upper=select_subset_arr(flux_,val=0,condition=">",into=phz_upper)
                gaia_g=select_subset_arr(flux_,val=0,condition=">",into=gaia_g)

            flux_=select_subset_arr(flux_,val=0,condition=">",into=flux_)
            mag_=22.5-(2.5*np.log10(flux_))

            #remove_stars=True

            if remove_stars==True:
                #https://arxiv.org/pdf/2007.14950.pdf
                test_star=gaia_g-mag_
                mask_notingaia=gaia_g==0
                
                mask_gaia=np.logical_and(test_star>0.61,np.logical_not(mask_notingaia))

                mask_star=np.logical_or(mask_notingaia,mask_gaia)
                stars_removed_spec=len(mag_)-np.sum(mask_star)
                stars_removed_wzspec=stars_removed_wzspec+stars_removed_spec
                mask_star=np.logical_or(mask_star,specz_>0.0) # avoid removing confirmed specz         

                stars_removed_count=len(mag_)-np.sum(mask_star)
                stars_removed=stars_removed+stars_removed_count
                
                #print("Not Stars ="+str(np.sum(mask_star))+" out of "+str(len(mag_))) 
                #print("Not in gaia ="+str(np.sum(mask_notingaia)))
                #print("in gaia but satisfy the cut ="+str(np.sum(mask_gaia))) 
                #print("Stars removed ="+str(stars_removed))
                #print("gag - r ="+str(test_star[:50]))

                #print("gag ="+str(gaia_g[:50]))
                #print("r ="+str(mag_[:50]))
                #print()
                if star_column==False:
                    phz_=phz_[mask_star]
                    ra_=ra_[mask_star]
                    dec_=dec_[mask_star]
                    phz_err=phz_err[mask_star]
                    flux_=flux_[mask_star]
                    mag_=mag_[mask_star]
                    if return_all_mags==True:
                        mag_r=mag_r[mask_star]
                        mag_g=mag_g[mask_star]
                        mag_z=mag_z[mask_star]
                        mag_w1=mag_w1[mask_star]
                        mag_w2=mag_w2[mask_star]
                        mag_w3=mag_w3[mask_star]
                        mag_w4=mag_w4[mask_star]

                        magerr_r=magerr_r[mask_star]
                        magerr_g=magerr_g[mask_star]
                        magerr_z=magerr_z[mask_star]
                        magerr_w1=magerr_w1[mask_star]
                        magerr_w2=magerr_w2[mask_star]
                        magerr_w3=magerr_w3[mask_star]
                        magerr_w4=magerr_w4[mask_star]

                        phz_low=phz_low[mask_star]
                        phz_upper=phz_upper[mask_star]

                    training_=training_[mask_star]
                    specz_=specz_[mask_star]
           
            if replace_zphot:
                count_z=0
                nonspec_z=0
                print("Replacing z phot by zspec")#for i in range(0,len(ra_)):#len(ra_)
    
                for m in range(0,len(ra_)):
                    if  specz_[m] > 0.0: #zspec[i]
                        phz_err[m]=error_spec
                        phz_[m]=specz_[m]
                        count_z=count_z+1
                    else:
                        nonspec_z=nonspec_z+1
    

                print ("A total of ", count_z," was replaced for zpec info of ", nonspec_z)#"+str(gaia_g[:50]))



            if color_cuts==True:
                #https://arxiv.org/pdf/2007.14950.pdf

                if fluxes==True:
                    aux_r=22.5-(2.5*np.log10(mag_r))
                    aux_g=22.5-(2.5*np.log10(mag_g))
                    aux_z=22.5-(2.5*np.log10(mag_z)) 
                    test_color01=aux_g-aux_r
                    test_color02=aux_r-aux_z


                else: 
                    test_color01=mag_g-mag_r

                    test_color02=mag_r-mag_z

                
                mask_color1= np.logical_and(test_color01 < 4,test_color01 > -1)
                mask_color2= np.logical_and(test_color02 < 4,test_color02 > -1)
                mas_color_full=np.logical_and(mask_color1,mask_color2)
                mas_color_full=np.logical_or(mas_color_full,specz_>0.0) # avoid removing confirmed specz                 

                color_removed_count=len(ra_)-np.sum(mas_color_full)
                color_removed=color_removed+color_removed_count

                #print("Not Stars ="+str(np.sum(mask_star))+" out of"+str(len(mag_))) 
                #print("Not in gaia ="+str(np.sum(mask_notingaia)))
                #print("color_removed ="+str(color_removed))

                #print("color_removed g-r 50 ="+str(test_color01[:50]))
                #print("color_removed r-z 50 ="+str(test_color02[:50]))
                #print("color_removed g 50 ="+str(mag_g[:50]))
                #print("color_removed r 50 ="+str(mag_r[:50]))

                #print()
                if color_column==False:
                    phz_=phz_[mas_color_full]
                    ra_=ra_[mas_color_full]
                    dec_=dec_[mas_color_full]
                    phz_err=phz_err[mas_color_full]
                    flux_=flux_[mas_color_full]
                    mag_=mag_[mas_color_full]
                    if return_all_mags==True:
                        mag_r=mag_r[mas_color_full]
                        mag_g=mag_g[mas_color_full]
                        mag_z=mag_z[mas_color_full]
                        mag_w1=mag_w1[mas_color_full]
                        mag_w2=mag_w2[mas_color_full]
                        mag_w3=mag_w3[mas_color_full]
                        mag_w4=mag_w4[mas_color_full]

                        magerr_r=magerr_r[mas_color_full]
                        magerr_g=magerr_g[mas_color_full]
                        magerr_z=magerr_z[mas_color_full]
                        magerr_w1=magerr_w1[mas_color_full]
                        magerr_w2=magerr_w2[mas_color_full]
                        magerr_w3=magerr_w3[mas_color_full]
                        magerr_w4=magerr_w4[mas_color_full]

                        phz_low=phz_low[mas_color_full]
                        phz_upper=phz_upper[mas_color_full]

                    training_=training_[mas_color_full]
                    specz_=specz_[mas_color_full]
                #else:
                #    mask_star=mask_star[] 





            # ============= starting cuts
            #### 01 =======remove objects with phz<=0
            #print('Remove phz < 0 ')
            phz_sel=select_subset_arr(phz_,val=0.0,condition=">",into=phz_)
            ra_=select_subset_arr(phz_,val=0,condition=">",into=ra_)
            dec_=select_subset_arr(phz_,val=0,condition=">",into=dec_)
            phz_err= select_subset_arr(phz_,val=0,condition=">",into=phz_err)
            mag_=select_subset_arr(phz_,val=0,condition=">",into=mag_)

            training_=select_subset_arr(phz_,val=0,condition=">",into=training_)
            specz_=select_subset_arr(phz_,val=0,condition=">",into=specz_)
            if star_column==True:
                mask_star=select_subset_arr(phz_,val=0,condition=">",into=mask_star)
            if color_column==True:
                mas_color_full=select_subset_arr(phz_,val=0,condition=">",into=mas_color_full)
            if return_all_mags==True:
                mag_r=select_subset_arr(phz_,val=0,condition=">",into=mag_r)
                mag_g=select_subset_arr(phz_,val=0,condition=">",into=mag_g)
                mag_z=select_subset_arr(phz_,val=0,condition=">",into=mag_z)
                mag_w1=select_subset_arr(phz_,val=0,condition=">",into=mag_w1)
                mag_w2=select_subset_arr(phz_,val=0,condition=">",into=mag_w2)
                mag_w3=select_subset_arr(phz_,val=0,condition=">",into=mag_w3)
                mag_w4=select_subset_arr(phz_,val=0,condition=">",into=mag_w4)

                magerr_r=select_subset_arr(phz_,val=0,condition=">",into=magerr_r)
                magerr_g=select_subset_arr(phz_,val=0,condition=">",into=magerr_g)
                magerr_z=select_subset_arr(phz_,val=0,condition=">",into=magerr_z)
                magerr_w1=select_subset_arr(phz_,val=0,condition=">",into=magerr_w1)
                magerr_w2=select_subset_arr(phz_,val=0,condition=">",into=magerr_w2)
                magerr_w3=select_subset_arr(phz_,val=0,condition=">",into=magerr_w3)
                magerr_w4=select_subset_arr(phz_,val=0,condition=">",into=magerr_w4)

                phz_low=select_subset_arr(phz_,val=0,condition=">",into=phz_low)
                phz_upper=select_subset_arr(phz_,val=0,condition=">",into=phz_upper)




            #print('Mag Lim Cut')
            ### 02 === Magnitude limit cut
            ra_=select_subset_arr(mag_,val=magr_lim,condition="<",into=ra_)
            dec_=select_subset_arr(mag_,val=magr_lim,condition="<",into=dec_)
            phz_sel=select_subset_arr(mag_,val=magr_lim,condition="<",into=phz_sel)
            phz_err=select_subset_arr(mag_,val=magr_lim,condition="<",into=phz_err)

            training_=select_subset_arr(mag_,val=magr_lim,condition="<",into=training_)
            specz_=select_subset_arr(mag_,val=magr_lim,condition="<",into=specz_)
            if star_column==True:
                mask_star=select_subset_arr(mag_,val=magr_lim,condition="<",into=mask_star)
            if color_column==True:
                mas_color_full=select_subset_arr(mag_,val=magr_lim,condition="<",into=mas_color_full)
            if return_all_mags==True:
                mag_r=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_r)
                mag_g=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_g)
                mag_z=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_z)
                mag_w1=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_w1)
                mag_w2=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_w2)
                mag_w3=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_w3)
                mag_w4=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_w4)

                magerr_r=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_r)
                magerr_g=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_g)
                magerr_z=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_z)
                magerr_w1=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_w1)
                magerr_w2=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_w2)
                magerr_w3=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_w3)
                magerr_w4=select_subset_arr(mag_,val=magr_lim,condition="<",into=magerr_w4)

                phz_low=select_subset_arr(mag_,val=magr_lim,condition="<",into=phz_low)
                phz_upper=select_subset_arr(mag_,val=magr_lim,condition="<",into=phz_upper)


            mag_=select_subset_arr(mag_,val=magr_lim,condition="<",into=mag_)




            ### 03 === photo-z error cut
            #print('phz error cut')
            ra_=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=ra_)
            dec_=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=dec_)
            phz_sel=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=phz_sel)
            mag_=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_)

            training_=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=training_)
            specz_=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=specz_)
            if star_column==True:
                mask_star=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mask_star)
            if color_column==True:
                mas_color_full=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mas_color_full)


            if return_all_mags==True:
                mag_r=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_r)
                mag_g=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_g)
                mag_z=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_z)
                mag_w1=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_w1)
                mag_w2=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_w2)
                mag_w3=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_w3)
                mag_w4=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=mag_w4)

                magerr_r=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_r)
                magerr_g=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_g)
                magerr_z=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_z)
                magerr_w1=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_w1)
                magerr_w2=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_w2)
                magerr_w3=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_w3)
                magerr_w4=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=magerr_w4)
                phz_low=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=phz_low)
                phz_upper=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=phz_upper)

            phz_err=select_subset_arr(phz_err,val=phot_err_lim,condition="<",into=phz_err)

            # ========== End of cuts
 
            if ((type(sel_ra)!=type(None)) and (get_dist_selradec==True) and (len(ra_)>0)):
                central_ra=float(sel_ra[i][1]+sel_ra[i][0])/2.0
                central_dec=float(sel_dec[i][1]+sel_dec[i][0])/2.0 
                #ref = SkyCoord(ra=x_ref*u.degree, dec=y_ref*u.degree) # frame='icrs'  
                #idx, d2d, d3d = c.match_to_catalog_sky(ref)
                #print("the maximum, and median of idx")
                #print(max(idx))
                #print(numpy.median(numpy.array(idx)))
                #print(idx)
                #d2d=d2d*u.degree
                #print ("separation")
                #print (sel_ra[i])
                #print (sel_dec[i])
                #print(ra_[:10])
                #print(dec_[:10])
                c1 = SkyCoord(ra=central_ra*u.degree, dec=central_dec*u.degree)#, frame='icrs')
                c2 = SkyCoord(ra_*u.degree, dec_*u.degree )#frame='icrs'
                sep = c1.separation(c2)
                sepamin=sep.arcmin
                #print (sepamin)



            if len(ra_)>0:
                _ras=np.concatenate((_ras,ra_))
                _decs=np.concatenate((_decs,dec_))
                _photz=np.concatenate((_photz,phz_sel))
                _photz_err=np.concatenate((_photz_err,phz_err))
                _magnitude_r=np.concatenate((_magnitude_r,mag_))

                _zspec=np.concatenate((_zspec,specz_))
                _training=np.concatenate((_training,training_))
                if star_column==True:
                    _mask_star=np.concatenate((_mask_star,mask_star))
                if color_column==True: 
                    _mask_color=np.concatenate((_mask_color,mas_color_full))

                if ((type(sel_ra)!=type(None)) and (get_dist_selradec==True)):
                    sepsamin=np.concatenate((sepsamin,sepamin)) 
                if return_all_mags==True:
                    mags_r=np.concatenate((mags_r,mag_r))
                    mags_g=np.concatenate((mags_g,mag_g))
                    mags_z=np.concatenate((mags_z,mag_z))
                    mags_w1=np.concatenate((mags_w1,mag_w1))
                    mags_w2=np.concatenate((mags_w2,mag_w2))
                    mags_w3=np.concatenate((mags_w3,mag_w3))
                    mags_w4=np.concatenate((mags_w4,mag_w4))

                    magserr_r=np.concatenate((magserr_r,magerr_r))
                    magserr_g=np.concatenate((magserr_g,magerr_g))
                    magserr_z=np.concatenate((magserr_z,magerr_z))
                    magserr_w1=np.concatenate((magserr_w1,magerr_w1))
                    magserr_w2=np.concatenate((magserr_w2,magerr_w2))
                    magserr_w3=np.concatenate((magserr_w3,magerr_w3))
                    magserr_w4=np.concatenate((magserr_w4,magerr_w4))
                    _photz_min=np.concatenate((_photz_min,phz_low))
                    _photz_max=np.concatenate((_photz_max,phz_upper))



                ra_=[]
                dec_=[]
                phz_=[]
                phz_sel=[]
                phz_low=[]
                phz_upper=[]
                phz_err=[]
                mag_=[]
                zspec_=[]
                training_=[] 

                if return_all_mags==True:
                    mag_r=[]
                    mag_g=[]
                    mag_z=[]
                    mag_w1=[]
                    mag_w2=[]
                    mag_w3=[]
                    mag_w4=[]

                    magerr_r=[]
                    magerr_g=[]
                    magerr_z=[]
                    magerr_w1=[]
                    magerr_w2=[]
                    magerr_w3=[]
                    magerr_w4=[]

                    phz_low=[]
                    phz_upper=[]
    if star_column==True:
        stars_removed=np.logical_not(_mask_star)
    if color_column==True:
        color_removed=np.logical_not(_mask_color)    
    #print("Stars removed ="+str(stars_removed))
    #print("Stars removed before considering z spec ="+str(stars_removed_wzspec))
    if return_all_mags==True and get_dist_selradec==False:
        return _ras,_decs,mags_r,mags_g,mags_z,mags_w1,mags_w2,mags_w3,mags_w4,magserr_r, \
               magserr_g,magserr_z,magserr_w1,magserr_w2,magserr_w3,magserr_w4, \
              _photz,_photz_err,_photz_min,_photz_max, _zspec, _training,stars_removed, color_removed 
    elif return_all_mags==True and get_dist_selradec==True:
        return _ras,_decs,mags_r,mags_g,mags_z,mags_w1,mags_w2,mags_w3,mags_w4,magserr_r, \
               magserr_g,magserr_z,magserr_w1,magserr_w2,magserr_w3,magserr_w4, \
              _photz,_photz_err,_photz_min,_photz_max, _zspec, _training,stars_removed, color_removed,sepsamin 
    else:
        return _ras,_decs,_magnitude_r,_photz,_photz_err, _zspec, _training
                

def check_intersection_square(ras,decs,moc):
    vertices = np.array([[ras[0],  decs[0]        ],
                    [ras[1] ,  decs[1]],
                    [ras[2], decs[2]],
                    [ras[3], decs[3]]])
    skycoord = SkyCoord(vertices, unit="deg", frame="icrs")
    # A point that we say it belongs to the inside of the MOC
    #inside = SkyCoord(ra=10, dec=5, unit="deg", frame="icrs")
    moc_sq = MOC.from_polygon_skycoord(skycoord, max_depth=9)
    inter = moc_sq.intersection(moc)
    return inter.sky_fraction,inter,moc_sq



def openbayeMulti(bay_file):
    with fits.open(bay_file) as hdul:
        hdul.info()
        hdul[1].columns
        hdr = hdul[0].header
        data_b = hdul[1].data
    #print (data_b.names)
    #print(data_b.columns)
    #print (hdul[1].columns)
    #print (hdul[0].columns)
    #try:
    try:
        uni=data_b['UNIQ']
        probd=data_b['PROBDENSITY']
        nested=False
    except:
        uni=np.array([])
        probd=np.array([])
        nested=True        
    #except:
    #print ('File in NESTED ordering')
    #print (hdr)   
    #print (uni)
    #print (type(uni[0]))
    return uni,probd,data_b,nested

def hour2degree(ras,decs):

    c = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle, u.deg))

    return c.ra.degree,c.dec.degree


def open2df_cat(file_2df):
    """
    serial    I6      Database serial number (=SEQNUM)
    spectra   I1      Number of spectra obtained
    name      A10     2dFGRS name (=NAME)
    UKST      A3      UKST plate (=IFIELD)
    ra        A11     R.A. (B1950) 
    dec       A11     Dec. (B1950)
    ra2000    A11     R.A. (J2000)
    dec2000   A11     Dec. (J2000)
    BJG       F6.3    Final bj magnitude without extinction correction 
    BJSEL     F6.3    Final bj magnitude with extinction correction
    BJG_OLD   F6.3    Original bj magnitude without extinction correction
    BJSELOLD  F6.3    Original bj magnitude with extinction correction
    GALEXT    F5.3    Galactic extinction value
    SB_BJ     F6.3    SuperCosmos bj magnitude without extinction correction
    SR_R      F6.3    SuperCosmos R magnitude without extinction correction 
    z         F9.6    Best redshift (observed)
    z_helio   F9.6    Best redshift (heliocentric)
    obsrun    A5      Observation run of best spectrum
    quality   I1      Redshift quality parameter for best spectrum
                         (quality=1-5; reliable redshifts have quality>=3)
    abemma    I1      Redshift type (abs=1, emi=2, man=3)
    Z_ABS     F9.6    Cross-correlation redshift from best spectrum
    KBESTR    I1      Cross-correlation template from best spectrum
    R_CRCOR   F5.3    Cross-correlation R value from best spectrum
    Z_EMI     F9.6    Emission redshift from best spectrum
    NMBEST    I2      Number of emission lines for Z_EMI from best spectrum
    SNR       F6.2    Median S/N per pixel from best spectrum
    ETA_TYPE  F10.6   Eta spectral type parameter from best spectrum (-99.9 if none)
    """
    _id,ra1,ra2,ra3,dec1,dec2,dec3,z,z_helio,rank=open_ascii_cat(file_2df,skip_rows=1,unpack=True,delimiter=" ",usecols=(0,10,11,12,13,14,15,23,24,26))
    #print(len(out))
    #cat,names=open_fits_catalog(file_2df)
    #cat=select_subset_fits_catalog(cat,field_name='Column28',field_value=2,condition=">")
    #print (fits_catalog_col(cat,col_name='Column12')[:5])
    #print (fits_catalog_col(cat,col_name='Column13')[:5])
    #print (fits_catalog_col(cat,col_name='Column14')[:5])
    #print (type(fits_catalog_col(cat,col_name='Column14')[0]))
    # fits_catalog_col(cat,col_name='Column13').astype('str')
    ras=np.char.array(ra1)+' '+np.char.array(ra2)+' '+np.char.array(ra3)
    #print (ras[0])
    #print (decs[0])

    decs=np.char.array(dec1)+' '+np.char.array(dec2)+' '+np.char.array(dec3)
    #>>> c = SkyCoord('00 42 30 +41 12 00', unit=(u.hourangle, u.deg))
    #z_helio=fits_catalog_col(cat,col_name='Column26')
    c = SkyCoord(ra=ras, dec=decs, unit=(u.hourangle, u.deg))
    return _id,c.ra.degree,c.dec.degree,np.array(z_helio).astype('float'),np.array(rank).astype('int')






def plotbayestar(bay_file,ra=None,dec=None, title='Bayestar',out_name='Bayestar.png', fig=None,ax=None, source_colors=['navy'], show_legend=False, plot_center=[320,20], FOV=105, multi=True, obj_marker='.', obj_marker_size=0.5, dense_plot=False):

    if multi==True:
        uniq,probdensity,data,nested=openbayeMulti(bay_file)
    
        level, ipix = ah.uniq_to_level_ipix(uniq)
        area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)

        prob = probdensity * area

        cumul_to = [0.9, 0.5]#np.linspace(0.5, 0.9, 5)[::-1]
        colors = ['indianred', 'indianred']#, 'indianred', 'indianred', 'indianred']
        alphas=[0.4,0.9]#[0.3,0.7]#0.4,0.9
        try:
            mocs = [MOC.from_valued_healpix_cells(uniq, prob, cumul_to=c) for c in cumul_to]
        except:
            print ('converting uniq to np.int64')
            mocs = [MOC.from_valued_healpix_cells(uniq.astype('int64'), prob, cumul_to=c) for c in cumul_to]
    #else:
        #mocs=[MOC.from_fits_image(bay_file)]

 
    if type(fig)==type(None) and multi==True:
        fig = plt.figure(111, figsize=(15, 10))
        # Define a astropy WCS easily
        with WCS(fig, 
            fov=FOV * u.deg,
            center=SkyCoord(plot_center[0], plot_center[1], unit='deg', frame='icrs'),
            coordsys="icrs",
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
            for (moc, c, col,al) in zip(mocs, cumul_to, colors,alphas):
                print("Creating maps")
                moc.fill(ax=ax, wcs=wcs, alpha=al, linewidth=0, fill=True, color=col, label='confidence probability ' + str(round(c*100)) + '%')
                moc.border(ax=ax, wcs=wcs, alpha=al, color=col)

    else: #moc_multi==True:
        fig = plt.figure(111, figsize=(15, 10)) 
        with WCS(fig, 
            fov=FOV * u.deg,
            center=SkyCoord(plot_center[0], plot_center[1], unit='deg', frame='icrs'),
            coordsys="icrs",
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            
            import ligo.skymap.plot

            ax = plt.axes(projection='astro aitoff')#astro #rotate=Angle(0, u.degree) #radius='50 deg' #'astro zoom' center=SkyCoord(plot_center[0], plot_center[1], unit='deg', frame='icrs')
            ax.grid()
            ax.contour_hpx(bay_file ) #cmap='cylon' # # 'astro zoom' #levels=[0.1, 0.3, 0.5]
        
    #else:
    #    fig = plt.figure(111, figsize=(15, 10))
    #    probs = hp.read_map(fitsFile)
    #    nside = hp.pixelfunc.get_nside(probs)

        # Choose color map and set background to white
    #    cmap = cm.YlOrRd
    #    cmap.set_under("w")
    #    probs = hp.pixelfunc.ud_grade(probs, 64) #reduce nside to make it faster
    #    probs = probs/np.sum(probs)
    #    pixels = np.arange(probs.size)
    #    sample_points = np.array(hp.pix2ang(nside,pixels)).T

        ### plot 90% containment contour of GW PDF
    #    cumul_to = np.linspace(0.5, 0.9, 5)[::-1]
    #    theta_contour, phi_contour = compute_contours(cumul_to,probs)
    #    hp.projplot(theta_contour[0],phi_contour[0],linewidth=1,c='k',label='GW (90% C.L)')
    #    for i in range(1,len(theta_contour)):
    #        hp.projplot(theta_contour[i],phi_contour[i],linewidth=1,c='k')
        
        # Plot GW skymap in Mollweide projection
        #hp.mollview(probs,cbar=True,unit=r'Probability',min=0,max=3e-5,rot=180,cmap=cmap)
        #hp.graticule() # Set grid lines




    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    #for (moc, c, col) in zip(mocs, cumul_to, colors):
    #    moc.fill(ax=ax, wcs=wcs, alpha=0.4, linewidth=0, fill=True, color=col, label='confidence probability ' + str(round(c*100)) + '%')
    #    moc.border(ax=ax, wcs=wcs, alpha=0.4, color=col)
    if type(ra)!=type(None):
        from astropy.wcs.utils import skycoord_to_pixel
        if len(source_colors)==1: 

            sources = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) #frame='icrs'
            #x, y = skycoord_to_pixel(sources, wcs)
            #ax.plot(x,y,marker='.', c="k")
            x = sources.ra.wrap_at(180 * u.deg).radian
            y = sources.dec.radian

            ax.scatter(x,y, s=obj_marker_size, marker=obj_marker, c=source_colors[0])
    
        else:
            for i in range(0,len(source_colors)):
                sources = SkyCoord(ra=ra[i]*u.degree, dec=dec[i]*u.degree) #frame='icrs'
                #from astropy.wcs.utils import skycoord_to_pixel
                #x, y = skycoord_to_pixel(sources, wcs)
                x = sources.ra.wrap_at(180 * u.deg).radian
                y = sources.dec.radian
                #ax.plot(x,y,marker='.', c="k")
                if dense_plot==False:
                    ax.scatter(x,y, s=obj_marker_size, marker=obj_marker, c=source_colors[i])
                else:
                    ax=density_scatter(x, y, ax = ax,fig=fig, bins = 50, alpha = 0.5, s=5,marker='o', save=False)         
        

    dec_lim=ax.get_ylim()
    ra_lim=ax.get_xlim()
    print('ra and dec limits')
    print(dec_lim)
    print(ra_lim) 
    #ax.set_ylim(dec_lim[0],)
    ax.set_xlim()
    if show_legend==True:
        ax.legend()
    ax.set_xlabel('ra')
    #plt.xlabel('ra')
    ax.set_ylabel('dec')
    ax.set_title(title)
    ax.grid(b=True,color="black", linestyle="-")
    plt.savefig(out_name)
    return ax,fig


def saveMOC():
    return 0   

def readMOC():
    return 0



def addstr2listleft(mylist,mystring):
    return [ mystring+s for s in mylist]

def open_fits_catalog(fits_file): 
    hdu_list=fits.open(fits_file, ignore_missing_end=True)
    #print hdu_list
    hdu = hdu_list[1]    # table extensions can't be the first extension, so there's a dummy image extension at 0
    #print hdu.header
    cat_table = hdu.data
    cols=hdu.columns
    return cat_table, cols

def add_col_fits_catalog(cols,newcol,col_name,col_format="E", out_name="",overwrite=False, list_of_cols=False):
    if list_of_cols==False:
        new_cols = fits.ColDefs([Column(name=col_name, format=col_format,
                 array=newcol)])
    else:
        new_cols = fits.ColDefs([fits.Column(name=col_name[i], format=col_format[i],
                 array=newcol[i]) for i in range(0,len(col_name))])
    hdu = fits.BinTableHDU.from_columns(cols + new_cols)
    if out_name!="":
        hdu.writeto(out_name,overwrite=overwrite)
    return hdu

def append_fits_tables(hdu1,hdu2):
    nrows1 = hdu1.data.shape[0]
    nrows2 = hdu2.data.shape[0]
    nrows = nrows1 + nrows2
    hdu = fits.BinTableHDU.from_columns(hdu1.columns, nrows=nrows)
    for colname in hdu1.columns.names:
            hdu.data[colname][nrows1:] = hdu2.data[colname]
    return hdu

def get_posterior_json(Catalog_object,varname=None,posteriortype='validation', return_ids=True):

    cats = list(Catalog_object['catalogs'])
    
    """ Get indexes (in case any) """
    indexes = {'all': None}
    if 'metadata' in Catalog_object:
        if 'indexes' in Catalog_object['metadata']:
            indexes = Catalog_object['metadata']['indexes']
        
    """ Get real/estimated values """
    actuals = {ix: \
               {c: np.array(Catalog_object['catalogs'][c]['channels']['actuals']['data']) \
                    if indexes[ix] is None else \
                    np.array(Catalog_object['catalogs'][c]['channels']['actuals']['data'])[indexes[ix]] \
                    for c in cats} \
                for ix in indexes}
    posteriors = {ix: \
               {c: np.array(Catalog_object['catalogs'][c]['channels']['posteriors']['data']) \
                    if indexes[ix] is None else \
                    np.array(Catalog_object['catalogs'][c]['channels']['posteriors']['data'])[:,indexes[ix]] \
                    for c in cats} \
                for ix in indexes}
    if return_ids==True:
        _ids = {ix: \
                   {c: np.array(Catalog_object['catalogs'][c]['channels']['sample_id']['data']) \
                        if indexes[ix] is None else \
                        np.array(Catalog_object['catalogs'][c]['channels']['sample_id']['data'])[indexes[ix]] \
                        for c in cats} \
                    for ix in indexes}
    

    #_ids = {ix: \
    #           {c: np.array(Catalog_object['catalogs'][c]['channels']['sample_id']['data']) \
    #                if indexes[ix] is None else \
    #                np.array(Catalog_object['catalogs'][c]['channels']['sample_id']['data'])[indexes[ix]] \
    #                for c in cats} \
    #            for ix in indexes}


    if type(varname)==type(None):
        return actuals,posteriors
    else:
        pdf=posteriors[posteriortype][varname]
        truth=actuals[posteriortype][varname]
        trutht=np.transpose(truth)
        pdft=np.transpose(pdf)
        fake_id_len=len(trutht[0])
        #fake_id=range(1,fake_id_len)
        if return_ids==True:
            _sample_ids=_ids[posteriortype][varname]
            _ids_corr=np.transpose(_sample_ids)[0]
        else:
            _ids_corr=range(1,len(trutht[0])+1)
        return _ids_corr,trutht[0],pdft[0] 
def jsoncat2deeparray(_ids,truths,pdfs):
    cat=[list([_ids[i],truths[i]])+(pdfs[i].tolist()) for i in range(0,len(_ids))] 
    return np.array(cat)


#>>> import numpy as np
#>>> from cbomcode.tools import fits_cat as fc
#>>> a=fc.load_jsoncatalog(filename='/share/storage1/SLcosmology/SLregression/POLIOUTPUT_inputs_images_outputs_z_lens-z_source-einstein_radius-sigma_velocity_heteroscedastic_bayesian_valsplit_20_architecture_inception_outunits_128-64/catalog_lenspop_+aleatoric_linadj.json')
# b=fc.get_posterior_json(a)
def load_jsoncatalog(filename):
    #Data = json.load(open(filename))
    #Data =json.loads(open(filename).read().decode('utf-8'))
    #Data = json.loads(open(filename, encoding='utf-8').read())
    txt = open(filename, encoding='utf-8').read()
    txt = txt.replace('nan','-1')
    Data = json.loads(txt)
    # transform data into np array (json parser loads them as lists)
    for c in Data['catalogs']:
        for ch in Data['catalogs'][c]['channels']:
            try:
                float(Data['catalogs'][c]['channels'][ch]['data'][0][0])
                xx = np.array(Data['catalogs'][c]['channels'][ch]['data'])
                nanid = np.where(xx == 'nan')
                xx[nanid] = 0
                Data['catalogs'][c]['channels'][ch]['data'] = xx.astype(np.float)
                Data['catalogs'][c]['channels'][ch]['data'][nanid] = np.nan
            except:
                print('Leaving channel {} as imported'.format(ch))
            
    return Data

def json_recursivity(dct, level, convert = True):
        
    # number of tabs:
    tab = ''.join(['\t']*level)
        
    txt = ''
    
    # loop through all entries in dct
    if isinstance(dct,dict):
        
        txt = tab
        # we shall put this entry between brackets
        tab_b = ''.join(['\t']*(level+1))
        # tabulation one level further
        tab_e = ''.join(['\t']*(level+2))
        
        txt += '\n' + tab_e + '{\n'
        for e in dct:
            print('{}{} - {}'.format(tab_e,level+2,e))
            
            # add entry
            #txt += tab_b + '{\n'
            txt += '{}"{}":'.format(tab_e, e)
            
            # recursivity on dictionary
            txt += _recursivity(dct[e], level + 1)
            
            # remove last '\n' and close entry
            txt = txt[:-1] + ',\n'
            #txt += tab_e + '},'
        
        # discard last comma and place line jump
        txt = txt[:-2]
        txt += '\n'
        
        # close bracket
        txt += tab_e + '}\n'
            
    elif isinstance(dct,str) or isinstance(dct,int) or isinstance(dct,float):
        tab_b = ''.join(['\t']*(level+2))
        # simply a comment or a string value of some sort
        txt += '"{}"\n'.format(dct)
        print('{}{} - "{}"'.format(tab_b,level+2,dct))
    elif isinstance(dct,np.ndarray) or isinstance(dct,list) or isinstance(dct,tuple):
        tab_b = ''.join(['\t']*(level+2))
        dct = np.array(dct)
        print('{}{} - array with shape {}'.format(tab_b,level+2,dct.shape))
        if convert:
            data_tmp = np.array2string(dct,separator=',',threshold=np.inf).replace('. ','.0').replace('\n','')
            txt += '{}\n'.format(data_tmp).replace('.]','.0]')
        
    
    return txt

def recarray2fits_catalog(new_data): 
    hdu_new = fits.BinTableHDU(data=new_data)
    cols=hdu_new.columns
    return hdu_new,cols

def add_col_fits_cat_missing(cat_ftb,id_array,newcol,col_name,missing=-99.0,col_format="E", out_name="",overwrite=False,id_col="QUICK_OBJECT_ID"):

    _ids=cat_ftb.field(id_col)
    addcol=np.ones(len(_ids))*missing
    print("addcol len")
    print(len(addcol))
    for i in range(0,len(id_array)):
        _idaddcol=np.where(_ids==id_array[i])
        #print("matched_ids")
        #print(_idaddcol)
        addcol[_idaddcol]=newcol[i]
    hdu_aux,cols_ftb=recarray2fits_catalog(cat_ftb)
    hdu=add_col_fits_catalog(cols=cols_ftb,newcol=addcol,col_name=col_name,col_format=col_format, out_name=out_name,overwrite=False)
    
    #    for i in range(0,len(id_array)):
    return hdu
      

        




def fits_catalog_col(cat,col_name):
    return cat.field(str(col_name))

def selectfromlist(target_ids,_ids,_ids_feature):
    target_features=[]
    for i in range(0,len(target_ids)):
        mask=_ids ==target_ids[i]
        target_features.append(_ids_feature[mask][0])
    return np.array(target_features)


def fits2ascii(fits_file, **kwargs):
    delimiter=verify_kwarg("delimiter",",",kwargs)
    cat,names=open_fits_catalog(fits_file)
    ascii_file=fits_file.rstrip(".fits")
    ascii_file=ascii_file.lstrip(".fit")
    cat_ascii=open(ascii_file+".txt", "w")
    
    for i in range(0,len(names)):
        cat_ascii.write("#"+str(names[i])+"\n")
    for i in range(0,len(cat[names[0]])):
        for j in range(0,len(names)):
            line=str(cat[names[j]][i]).replace("\n", "")
            if j!=(len(names)-1):
                
                cat_ascii.write(line+delimiter)
            else:
                cat_ascii.write(line)
            
        if i!=(len(cat[names[0]])-1):
            cat_ascii.write("\n") 
        
    cat_ascii.close()
    return 0

def add_columns(cat,data,col_names,col_types):
    cols=[]
    for i in range(0,len(col_names)):
        cols.append(fits.Column(name=col_names[i], format=col_types[i], array=data[i]))
    new_cols = fits.ColDefs(cols)
    orig_cols=cat.columns
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    return hdu.data,hdu.columns

def estimate_psf_auto (cat,pix_scale,**kwargs):
    CS=verify_kwarg("CLASS_STAR",0.95,kwargs)
    MAG=verify_kwarg("mag",'MAG_AUTO',kwargs)
    cat=select_subset_fits_catalog(cat,field_name='FLUX_RADIUS',field_value=0,condition=">")

    cat=select_subset_fits_catalog(cat,field_name=MAG,field_value=90,condition="abs<")
    if ('MAG_ERR' in kwargs.keys()):
        cat=select_subset_fits_catalog(cat,field_name="MAGERR_AUTO",field_value=kwargs['MAG_ERR'],condition="<")
    if ('FLAGS' in kwargs.keys()):
        cat=select_subset_fits_catalog(cat,field_name='FLAGS',field_value=kwargs['FLAGS'],condition="==")
    cat=select_subset_fits_catalog(cat,field_name='CLASS_STAR',field_value=CS,condition=">")
    flux_r=fits_catalog_col(cat,"FLUX_RADIUS")
    mean_r=flux_r.mean()
    err_r=flux_r.std()
    median_r=numpy.median(flux_r)
    mean_r_deg=mean_r*pix_scale
    err_r_deg=pix_scale*err_r
    psf_sex_pix=numpy.median(fits_catalog_col(cat,"FWHM_IMAGE"))
    return mean_r,err_r,median_r, mean_r_deg,err_r_deg, psf_sex_pix

def find_nearest(x,value):
    x=numpy.array(x)
    x=numpy.nan_to_num(x)
    idx = (numpy.abs(x-value)).argmin()
    return idx

def filter_outlayers(x,sigmas=1,**kwargs):
    #if 'into' in kwargs.keys():
    #    y=kwargs['into']
    #else:
    #    y=numpy.array(x).copy()
    std_x=numpy.std(x)
    mean_x=numpy.mean(x)
    y=select_subset_arr(x,val=[mean_x-(std_x*sigmas),mean_x+(std_x*sigmas)],condition='interval',**kwargs)
    #y=select_subset_arr(x,val=[],condition='>',**kwargs)
    return y


def filter_columns(x,y,sigmas):
    """
    Filter outlyers and return x,y with same size. Y is the test column.
    """
    x=filter_outlayers(y,sigmas, into=x)
    y=filter_outlayers(y,sigmas)
    x = x[numpy.logical_not(numpy.isnan(y))]
    y = y[numpy.logical_not(numpy.isnan(y))]
    x = x[numpy.logical_not(numpy.isinf(y))]
    y = y[numpy.logical_not(numpy.isinf(y))]
    return x,y


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def filter_bined_data(bin_center,bin_vals,bin_valsx,sigmas=2,bin_valserr=None):
    """
    Filter data bined by make_bin and return bin_center,bin_vals 
    with same size suited for a scatter plot.
    """
    colx=[]
    colx_data=[]
    coly=[]
    colyerr=[]
    for i in range(0,len(bin_center)):

        bin_vals_aux=numpy.array(bin_vals[i]).copy()
        bin_vals_auxerr=numpy.array(bin_vals[i]).copy()
        bin_vals[i]=filter_outlayers(bin_vals[i],sigmas=sigmas)
        bin_valsx[i]=filter_outlayers(bin_vals_aux,sigmas=sigmas,into=bin_valsx[i])


        if bin_valserr!=None:
            bin_valserr[i]=filter_outlayers(bin_vals_auxerr,sigmas=sigmas,into=bin_valserr[i])

    for i in range(0,len(bin_center)):
        colx=colx+list(numpy.ones(len(bin_vals[i]))*bin_center[i])
        colx_data=colx_data+list(bin_valsx[i])
        coly=coly+list(bin_vals[i])
        if bin_valserr!=None:
            colyerr=colyerr+list(bin_valserr[i])
        else:
            colyerr=None 
    if bin_valserr!=None:
        return colx,coly,colx_data,colyerr
    else:
        return colx,coly,colx_data        


def estimate_psf (cat,pix_scale,lim,**kwargs):
    """
    be aware that flux_radius is the half of the diameter, so it is necessary to multiplicate by factor 2
    prior to this function
    e.g.
    cat.field('FLUX_RADIUS')[:]=cat.field('FLUX_RADIUS')[:]*2
    """
    CS=verify_kwarg("CLASS_STAR",0.0,kwargs)
    MAG=verify_kwarg("mag",'MAG_AUTO',kwargs)
 
    cat=select_subset_fits_catalog(cat,field_name='FLUX_RADIUS',field_value=0,condition=">")
    #cat=select_subset_fits_catalog(cat,field_name='FLAGS',field_value=0,condition="==")
    cat=select_subset_fits_catalog(cat,field_name=MAG,field_value=90,condition="abs<")
    cat=select_subset_fits_catalog(cat,field_name='CLASS_STAR',field_value=CS,condition=">")

    if ('mag_err' in kwargs.keys()):
        cat=select_subset_fits_catalog(cat,field_name="MAGERR_AUTO",field_value=kwargs['mag_err'],condition="<")
    if ('FLAGS' in kwargs.keys()):
        cat=select_subset_fits_catalog(cat,field_name='FLAGS',field_value=kwargs['FLAGS'],condition="==")
    cat=select_subset_fits_catalog(cat,field_name="FLUX_RADIUS",field_value=lim[0],condition=">")
    cat=select_subset_fits_catalog(cat,field_name="FLUX_RADIUS",field_value=lim[1],condition="<")
    
    cat=select_subset_fits_catalog(cat,field_name=MAG,field_value=lim[2],condition=">")
    cat=select_subset_fits_catalog(cat,field_name=MAG,field_value=lim[3],condition="<")


    flux_r=fits_catalog_col(cat,"FLUX_RADIUS")
    mean_r=flux_r.mean()
    err_r=flux_r.std()
    median_r=numpy.median(flux_r)
    mean_r_deg=mean_r*pix_scale
    err_r_deg=pix_scale*err_r

    return mean_r,err_r,median_r, mean_r_deg,err_r_deg


 





def math_cat_sky(cat,cat_ref,sep=0.00027, unit='degree', out='cat',cat_type=['fits','fits'],col_cat=[None,None],col_ref=[None,None],radec_name=["RA","DEC"],radec_nameref=["RA","DEC"],ref_unit='degree'):
    if unit=='arcsec':
        sep=sep*0.00027
    if cat_type[0]=='fits':
        x=fits_catalog_col(cat,radec_name[0])
        y=fits_catalog_col(cat,radec_name[1])
    elif cat_type[0]=='array':
        x=cat[0]
        y=cat[1] 
    else:
        cat_tp=np.transpose(cat)
        x=cat_tp[col_cat[0]]
        y=cat_tp[col_cat[1]]
    if cat_type[1]=='fits':
        x_ref=fits_catalog_col(cat_ref,radec_nameref[0])
        y_ref=fits_catalog_col(cat_ref,radec_nameref[1])
    elif cat_type[1]=='array':
        x_red=cat_ref[0]
        y_ref=cat_ref[1]
    else:
        cat_ref_tp=np.transpose(cat_ref)
        y_ref=np.array(cat_ref_tp[col_ref[1]]).astype('float')
        x_ref=np.array(cat_ref_tp[col_ref[0]]).astype('float')
    if ref_unit=='hour':
        x_ref,y_ref=hour2degree(x_ref,y_ref)

    
    c = SkyCoord(ra=x*u.degree, dec=y*u.degree) #frame='icrs'
    #print(type(x_ref[0]))
    #print(type(x[0]))  
    ref = SkyCoord(ra=x_ref*u.degree, dec=y_ref*u.degree) # frame='icrs'  
    idx, d2d, d3d = c.match_to_catalog_sky(ref)
    #print("the maximum, and median of idx")
    #print(max(idx))
    #print(numpy.median(numpy.array(idx)))
    #print(idx)
    d2d=d2d*u.degree
    match_ref=idx[numpy.where(numpy.array(d2d)<=sep)]
    match=numpy.where(numpy.array(d2d)<=sep)[0]
    
    if out=='cat':
        return cat[match],cat_ref[match_ref]
    else:
        return match,match_ref

    


def clear_unmatched(cat,matched,obj_img):

    return 0  


def isfits_catalog_col(col_test,col_names):
    if col_test in col_names.names:
        return True
    else:
        return False

from astropy.table import Table

def save_fits_catalog(file_name,data):
    hdu = fits.BinTableHDU(data=data)
    hdu.writeto(file_name, clobber=True)
    return 0

def create_fits_catalog_from_array(file_name,data_array,col_names):
    t = Table(data_array, names=col_names)
    t.write(file_name, format='fits', overwrite=True)
    return 0

def create_assoc_file_from_fits_cat(assoc_file,cat):
    fil = open(assoc_file, "w")
    x=fits_catalog_col(cat,col_name=X_IMAGE)
    y=fits_catalog_col(cat,col_name=Y_IMAGE)
    x,y=open_ascii_cat(cat_file,usecols=(3,4), unpack=True)
    for i in range(0,len(x)-1):
        fil.write(str(x[i])+" "+str(y[i])+"\n")
    fil.write(str(x[-1])+" "+str(y[-1]))        
    fil.close()
    return 0

def sex_assoc():
    """
    Warning: Assoc will work only if assoc variables are set as output in the .param file
    """



    return 0

def create_reg(x,y,file_name,**kwargs):
    color=verify_kwarg("color","green",kwargs)
    regunit=verify_kwarg("unit","pixel",kwargs)
    asc = open(file_name, "w")
    asc.write("# Region file format: DS9 version 4.1 \n")
    asc.write("global color="+color+" dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    if (regunit!='pixels') or (regunit!='pixel'):
        asc.write(" fk5 \n")
        r=verify_kwarg("radius",str(3.7200001),kwargs)
    else:
        asc.write(" image \n")  
        r=verify_kwarg("radius",str(20),kwargs)
    for i in range(0,len(x)-1):
        asc.write("circle("+str(x[i])+","+str(y[i])+","+str(r)+")\n")
    asc.write("circle("+str(x[-1])+","+str(y[-1])+","+str(r)+")\n")
    asc.close()
    return 0


def plot_kwargs(args,**kwargs):

    bin=verify_kwarg("bin",False,args)
    pos=verify_kwarg("pos",[0,0],args)
    save=verify_kwarg("save",True,args)
    work_dir=verify_kwarg("work_dir",str(os.getcwd()),args)
    time=datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d--%H%M%S-%f')
    if ('x_label' in kwargs.keys()):
        x_label=kwargs['x_label']
        if ('y_label' in kwargs.keys()):
            y_label=kwargs['y_label']
            plot_name=verify_kwarg("plot_name",str(work_dir+"/"+time+str(x_label)+"X"+str(y_label)+".png" ),args)
        else:
            plot_name=verify_kwarg("plot_name",str(work_dir+"/"+time+str(x_label)+".png" ),args)
    else:
        plot_name=verify_kwarg("plot_name",str(work_dir+"/"+time+".png" ),args)
    return pos,save,plot_name,bin

def fits_catalog_hist(cat,col, **kwargs):
    x_label=verify_kwarg("xlabel",str(col),kwargs)
    kwargs['xlabel']=x_label
    kwargs['pos'],kwargs['save'],kwargs['plot_name'],kwargs['bin']=plot_kwargs(kwargs, x_label=x_label)
    ran=verify_kwarg('range',None,kwargs)
    #xlim=verify_kwarg('xlim',[-40,80],kwargs)
    norm=verify_kwarg('normed',False,kwargs)
    bins=verify_kwarg('bins','auto',kwargs)
    col_dat=fits_catalog_col(cat,col)

    if ('col_index' in kwargs.keys()):
        col_dat=numpy.array(zip(*col_dat)[kwargs['col_index']])
    
    hist,bin_edges=cat_hist_plot(col_dat,plt_type="bar",**kwargs)
    return hist,bin_edges

def cat_hist_plot(image,**kwargs):
    ran=verify_kwarg('rang',None,kwargs)
    #xlim=verify_kwarg('xlim',[-40,80],kwargs)
    normed=verify_kwarg('normed',False,kwargs)
    bins=verify_kwarg('bins','auto',kwargs)
    legend_names=verify_kwarg('legend_names',[],kwargs)
    fill=verify_kwarg('fill',False,kwargs)
    marker=verify_kwarg('marker','b.',kwargs)
    hatch=verify_kwarg('hatch',None,kwargs)
    bar_scale=verify_kwarg('bar_scale',1.0,kwargs)
    shift=verify_kwarg('shift',0.0,kwargs)
    plt_type=verify_kwarg('plt_type',"line",kwargs)
    save=verify_kwarg('save',True,kwargs)
    fig=verify_kwarg('figure',None,kwargs)
    ax1=verify_kwarg('axis',None,kwargs)
    bins_corr=verify_kwarg('bins_corr',False,kwargs)
    line_sty=verify_kwarg('line_sty',None,kwargs)
    line_w=verify_kwarg('line_w',1,kwargs)
    alpha_bar=verify_kwarg('alpha',0.8,kwargs)
    time=datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d--%H%M%S')
    out_name=verify_kwarg('plot_name',time+'histogram.png',kwargs)
    fit_gauss=verify_kwarg('fit_gauss',False,kwargs)
    hist_weights=verify_kwarg('hist_weights',None,kwargs)
    addtext=''
    colors=['indianred','navy','g','b','r','y','m','k','c','b','g','r']
    if len(legend_names)==0:
        image=[image]
        marker=[marker]
        legend_names=['']
        alpha_bar=[alpha_bar]
        if not(isarraylike(hatch)):
            hatch=[hatch]
        if type(line_sty)==type(None):
            line_sty=['-']
        if type(hist_weights)==type(None):
            hist_weights=[None]
    else:
        if type(hatch)==type(None):
            hatch=[None for i in range(0,len(legend_names))]
        if type(line_sty)==type(None):
            line_sty=['-' for i in range(0,len(legend_names))]
        if not(isarraylike(alpha_bar)): 
            alpha_bar_new=[alpha_bar for i in range(0,len(legend_names))]
            alpha_bar=alpha_bar_new
        if not(isarraylike(marker)): 
            marker_new=[marker for i in range(0,len(legend_names))]
            marker=marker_new
        if not(isarraylike(hist_weights)): 
            weights_new=[None for i in range(0,len(legend_names))]
            hist_weights=weights_new

    if fig==None:
        fig = plt.figure()
    tex=verify_kwarg('tex',False,kwargs)
    if (tex == True):
        rc('text', usetex=True)
    else:
        rc('text', usetex=False)
    if ax1==None:
        ax1 = fig.add_subplot(111)
    if plt_type=="bar":
        #print "second bin edges"
        #print bin_edges
        for i in range(0,len(legend_names)):
            width_bar=[]

            if len(bins)>0 and (type(bins)!=type('string')):
               bin_edges=bins
            else:
               hist,bin_edges=cat_hist(image[i],bins=bins,normed=normed, bins_corr=bins_corr, weights=hist_weights[i])
            #print (bin_edges)
            bin_edges=np.array(bin_edges).astype('float')
            for k in range(0,len(bin_edges)-1):
                width_bar.append((bin_edges[k+1]-bin_edges[k])/2)
            width_bar=np.array(width_bar)
            print("width bar")
            print(len(width_bar))
            print (len(bin_edges))
            hist,bin_edges=cat_hist(image[i],rang=(min(image[i])-width_bar[0],max(image[i])+width_bar[-1]),bins=bins,normed=normed, bins_corr=bins_corr,weights=hist_weights[i])
            if fit_gauss==True:
                mean_,std_=scipy.stats.norm.fit(image[i])
                addtext=addtext+legend_names[i]
                addtext=addtext+"\n $\mu$: "+str(round(mean_,3))
                addtext=addtext+"\n $\sigma$: "+str(round(std_,3))+"\n "

            print ("this is width_bar")
            #print type(width_bar)
            print (width_bar)
            #for i in range(0,len(bin_edges)):
            #    bin_edges[i]=bin_edges[i]-(2*width_bar)
            print ("third bin edges")
            print (bin_edges[:-1])
            print (len(hist))
            print (hatch,alpha_bar,colors,fill)
            print (2*bar_scale*width_bar)
            print (len(2*bar_scale*width_bar))
            ax1.bar(bin_edges[:-1],hist,width=2*bar_scale*width_bar[0],alpha=alpha_bar[i],hatch=hatch[i], ec=colors[i], fill=fill, color=colors[i])

    elif plt_type=="step": 
        for i in range(0,len(legend_names)):
            hist,bin_edges=cat_hist(image[i],rang=ran,bins=bins,normed=normed, bins_corr=bins_corr,weights=hist_weights[i])    
            width_bar=(bin_edges[1]-bin_edges[0])/2
            bin_edges=bin_edges+width_bar+(i*shift)
            print("bin edges")
            print(bin_edges)
            #print(type([bin_edges[-2]]))
            bin_edges[-1]=bin_edges[-2]+0.0001
            hist=np.concatenate((hist,[0.0])) 
            print(line_sty)
            print("bin edges")
            print(bin_edges)
            ax1.step(bin_edges,hist,alpha=alpha_bar[i], color=colors[i],linestyle=line_sty[i],lw=line_w)#ls=line_sty
            #ax1.plot(bin_edges[:-1],hist, marker[i], markersize=5, color=colors[i]) #ls=line_sty #lw=0.5

    else: 
        for i in range(0,len(legend_names)):
            hist,bin_edges=cat_hist(image[i],rang=ran,bins=bins,normed=normed, bins_corr=bins_corr, weights=hist_weights[i])    
            width_bar=(bin_edges[1]-bin_edges[0])/2
            bin_edges=bin_edges+width_bar+(i*shift)
            #print("bin edges")
            #print(bin_edges)
            #ax1.step(bin_edges[:-1],hist,alpha=alpha_bar[i], color=colors[i])
            ax1.plot(bin_edges[:-1],hist, marker[i], markersize=5, color=colors[i]) #ls=line_sty #lw=0.5
    if ('ylim' in kwargs.keys()):
        ylim=kwargs['ylim']
        plt.ylim(ylim[0],ylim[1])
    if ('xlim' in kwargs.keys()):
        xlim=kwargs['xlim']
        plt.xlim(xlim[0],xlim[1])
    if ('xlabel' in kwargs.keys()):
        xlbl=kwargs['xlabel']
        plt.xlabel(xlbl, fontsize=20, labelpad=-5)
    if ('ylabel' in kwargs.keys()):
        ylbl=kwargs['ylabel']
        plt.ylabel(ylbl, rotation=0,fontsize=20, labelpad=15)
    if ('hide_y_ticks' in kwargs.keys()):
        if kwargs['hide_y_ticks']==True:
            ax1.axes.yaxis.set_ticks([])
    #plt.rcParams['xtick.labelsize']=50
    plt.xticks(fontsize=14)



    if ('scale' in kwargs.keys()) and ('mu' in kwargs.keys()) and ('sigma' in kwargs.keys()) and ('xlim' in kwargs.keys()):
        xm = numpy.linspace(xlim[0], xlim[1], 10000)  # 100 evenly spaced points
        popt=[kwargs['scale'],kwargs['mu'], kwargs['sigma']]
        ax1.plot(xm, gauss(xm, *popt),"g-",label="fit")

    
    full_text=''
    if ('text' in kwargs.keys()):
        if addtext!='':
           full_text=kwargs['text']+"\n "+addtext
        else:
           full_text=kwargs['text']
    else:
        full_text=addtext
    if full_text!='':
        pos=verify_kwarg("pos",[0.7,0.2],kwargs)
        plt.figtext(pos[0], pos[1], full_text)

    
    if save==True:    
        plt.savefig(out_name)
    return fig,ax1

def cat_hist(image,rang=None,bins='auto',**kwargs):
    bins_corr=verify_kwarg('bins_corr',False,kwargs)
    normed=verify_kwarg('normed',False,kwargs)
    weights=verify_kwarg('weights',None,kwargs)

    #print "hist info"
    #print type(image)
    #print len(image)
    #print bins
    
    hist, bin_edges = numpy.histogram(numpy.array(image),bins=bins,range=rang,density=normed,weights=weights)
    #print "first bin edges"
    #print bin_edges
    interval=(bin_edges[1]-bin_edges[0])/2




    if bins_corr==True:

        width_bar=[]
        for k in range(0,len(bin_edges)-1):
            width_bar.append((bin_edges[k+1]-bin_edges[k])/2)
        width_bar.append(width_bar[-1])
        width_bar=np.array(width_bar)
        for i in range(0,len(bin_edges)):
            bin_edges[i]=bin_edges[i]+width_bar[i]
    
    return hist, bin_edges










def fits_catalog_scatter(cat,col1,col2, **kwargs):
    x_label=verify_kwarg("xlabel",str(col1),kwargs)
    y_label=verify_kwarg("ylabel",str(col2),kwargs)
    line=verify_kwarg("line",False,kwargs)
    pos,save,plot_name,bin=plot_kwargs(kwargs,x_label=x_label,y_label=y_label)
      
    col1_dat=fits_catalog_col(cat,col1)
    col2_dat=fits_catalog_col(cat,col2)

    if ('col1_index' in kwargs.keys()):
        col1_dat=zip(*col1_dat)[kwargs['col1_index']]
    if ('col2_index' in kwargs.keys()):
        col2_dat=zip(*col2_dat)[kwargs['col2_index']]
    tex=verify_kwarg('tex',False,kwargs)
    if (tex == True):
        rc('text', usetex=True)

    fig = plt.figure()
    #plt.rc('text', usetex=True)
    ax1 = fig.add_subplot(111)
    if bin==True:
        bin_edges=verify_kwarg("bin_edges",numpy.arange(numpy.amin(col1_dat),numpy.amax(col1_dat),(numpy.amax(col1_dat)-numpy.amin(col1_dat))/len(col1_dat)),kwargs)
        col1_dat,col2_dat,col_err=make_bin(col1_dat,col2_dat,bin_edges)
        ax1.errorbar(col1_dat,col2_dat,yerr=col_err)
    else:
        ax1.scatter(col1_dat,col2_dat, s=1.0)
    if line==True:
        x_line=numpy.arange(0,numpy.amax(col1_dat)+1,0.5)
        y_line=numpy.arange(0,numpy.amax(col1_dat)+1,0.5)
        ax1.plot(x_line,y_line,'r-') 

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if ('xlim' in kwargs.keys()):
        plt.xlim(kwargs['xlim'][0],kwargs['xlim'][1])
    if ('ylim' in kwargs.keys()):
        plt.ylim(kwargs['ylim'][0],kwargs['ylim'][1])
    if ('box' in kwargs.keys()):
        box=kwargs['box']
        ax1.add_patch( patches.Rectangle(box[0],box[1],box[2],fill=False,edgecolor="red"))
    if ('text' in kwargs.keys()):
        plt.figtext(pos[0], pos[1], kwargs['text'])
    if ('title' in kwargs.keys()):
        plt.title(kwargs['title'])
    if save==True:
        plt.savefig(plot_name)
    axes = plt.gca()
    return fig, ax1


def start_plot(**kwargs):

    fig = plt.figure()

    tex=verify_kwarg('tex',False,kwargs)
    axis_eq=verify_kwarg('axis_eq',False,kwargs)
    if (tex == True):
        rc('text', usetex=True)
        #matplotlib.rcParams['text.latex.unicode']=False

    #plt.rc('text', usetex=True)
    axes = fig.add_subplot(111)

    x_label=verify_kwarg("xlabel",'x',kwargs)
    y_label=verify_kwarg("ylabel",'y',kwargs)
    if ('ylim' in kwargs.keys()):
        if type(kwargs['ylim'][0])==list:
            plt.ylim(kwargs['ylim'][0][0],kwargs['ylim'][0][1])
        else:
            plt.ylim(kwargs['ylim'][0],kwargs['ylim'][1])
    if type(x_label)==list:
        plt.xlabel(x_label[0])
    else:
        plt.xlabel(x_label)
    if type(y_label)==list:
        plt.ylabel(y_label[0])
    else:
        plt.ylabel(y_label)
    if axis_eq==True:
        plt.axis('equal')

    return fig,axes

def start_subplots(num_sub,**kwargs):
    x_label=verify_kwarg("xlabel",'x',kwargs)
    equal=verify_kwarg("axis1_eq",False,kwargs)
    
 
    if (type(x_label) is str):
        x_label=[x_label]
    y_label=verify_kwarg("ylabel",'y',kwargs)
    if (type(y_label) is str):
        y_label=[y_label]
    sx=verify_kwarg("sharex",True,kwargs)
    #number=verify_kwarg("num_sub",2,kwargs)
    fig, ax = plt.subplots(num_sub, sharex=sx)
    tex=verify_kwarg('tex',False,kwargs)
    if (tex == True):
        rc('text', usetex=True)
        #matplotlib.rcParams['text.latex.unicode']=False


    if sx==True:
        ax[-1].set_xlabel(x_label[0])
    elif len(x_label) >= num_sub:
        for i in range(0,num_sub):
            ax[i].set_xlabel(x_label[i])
    else:
        for i in range(0,num_sub):
            ax[i].set_xlabel('x')

    if len(y_label) >= num_sub:
        for i in range(0,num_sub):
            ax[i].set_ylabel(y_label[i])
    else:
        for i in range(0,num_sub):
            ax[i].set_ylabel('y')
    if 'title' in kwargs.keys():
        ax[0].set_title(kwargs['title'])


    if  equal==True:
        ax[0].axis('equal')
    if ('ylim' in kwargs.keys()):
        if (type(kwargs['ylim'][0]) is int) or (type(kwargs['ylim'][0]) is float):
            for i in range(0,num_sub):
                #ax[i].set_ylabel('y')
                ax[i].set_ylim(kwargs['ylim'][0],kwargs['ylim'][1])
        else:
            for i in range(0,num_sub):
                #ax[i].set_ylabel('y')
                ax[i].set_ylim(kwargs['ylim'][i][0],kwargs['ylim'][i][1])
    if ('xlim' in kwargs.keys()):
        if (type(kwargs['xlim'][0]) is int) or (type(kwargs['xlim'][0]) is float):
            for i in range(0,num_sub):
                #ax[i].set_ylabel('y')
                ax[i].set_xlim(kwargs['xlim'][0],kwargs['xlim'][1])
        else:
            for i in range(0,num_sub):
                #ax[i].set_ylabel('y')
                ax[i].set_xlim(kwargs['xlim'][i][0],kwargs['xlim'][i][1])


    return fig,ax

def start_regression_plot():

    return 0

def measure_uplow_std_limits(ybin, bin_stdevs, nsigma=1):
    lo = [ybin[i] - nsigma*bin_stdevs[i] for i in range(0,len(ybin))] 
    hi = [ybin[i] + nsigma*bin_stdevs[i] for i in range(0,len(ybin))]
        
    return lo, hi


def plot_summary_results(
        medians_1,medians_2,errors_1=None,
        errors_2=None,refmedian=None,refmedian_2=None,
        referror=None,referror_2=None,
        model_names=None, left_margin=0.2,
        plot_name='summary_plot.png',
        label1='$\\Omega_M$',label2='$w$',
        xmin=None, xmax=None,
        ymin=None, ymax=None,
        figsize=None,line_after=None):
    #(6, 8)
    if type(figsize)==type(None):
        hgt=float(len(medians_1))+0.5
        wth=8
        figsize=(wth,hgt)


    custom_cycler = (cycler(color=['r', 'g', 'b', 'y','c', 'm', 'y', 'k']))

    plt.rc('axes', prop_cycle=custom_cycler)

    fig = plt.figure(figsize=figsize, edgecolor='black', constrained_layout=True)         
    




    # set up grid of panels in figure

    grid = plt.GridSpec(1, 2, hspace=0.0, wspace=0.0, 
                        left=left_margin, right=0.9, bottom=0.3, top=1.0)#height_ratios=[4], width_ratios=[1,1]
    # main axis: predicted statistic vs. true statistic
    ax_var1 = fig.add_subplot(grid[0,0])
    
    # hist axis [x]: histogram of data on x axis on top
    ax_var2 = fig.add_subplot(grid[0,1], sharey=ax_var1)
    rang=np.flip(np.arange(0,float(len(model_names))/4,0.25))
    #ax_var1.set_yticks(rang, model_names)   
    plt.yticks(rang, model_names)
    if type(referror)!=type(None):
        ax_var1.fill_betweenx(rang,refmedian-referror[0],refmedian+referror[1],alpha=0.3, color='#bcc4d1')
    if type(referror_2)!=type(None):
        ax_var2.fill_betweenx(rang,refmedian_2-referror_2[0],refmedian_2+referror_2[1],alpha=0.3, color='#bcc4d1')

    plt.setp(ax_var2.get_yticklabels(), visible=False)
    if type(errors_1)!=type(None):
        lo,hi=np.transpose(errors_1)
        for i in range(0,len(medians_1)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var1.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')
            ax_var1.errorbar([medians_1[i]], [rang[i]], xerr=[[lo[i]],[hi[i]]], marker='o', markersize=5, linestyle='')
    else:
        for i in range(0,len(medians_1)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var1.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')
            ax_var1.plot([medians_1[i]], [rang[i]], marker='o', markersize=5, linestyle='')

    #print (lo2,hi2)
    #print('essa e o vetor original')
    #print(errors_2)
    if type(errors_2)!=type(None):
        lo2,hi2=np.transpose(errors_2)
        for i in range(0,len(medians_2)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var2.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')

            ax_var2.errorbar([medians_2[i]], [rang[i]], xerr=[[lo2[i]],[hi2[i]]], marker='o', markersize=5, linestyle='')
    else:
        for i in range(0,len(medians_2)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var2.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')
            ax_var2.plot([medians_2[i]], [rang[i]], marker='o', markersize=5, linestyle='')
    
    if type(refmedian)!=type(None):
        ax_var1.axvline(x=refmedian, ls='dashed', lw=0.5)
    if type(refmedian_2)!=type(None):
        ax_var2.axvline(x=refmedian_2, ls='dashed', lw=0.5)
    ax_var1.set_ylim((-0.1, float(len(model_names))/4)) 
    ax_var1.set_xlabel(label1)  
    ax_var2.set_xlabel(label2)  
    #high_error_var1,low_error_var1=np.transpose(real_data)
    #yerr_lo=real_me - real_lo
    #yerr_hi= real_hi - real_me

    #yerr_lor=realr_me - realr_lo
    #yerr_hir= realr_hi - realr_me

    #res_me=np.divide(real_me-real_tr,real_tr)
    #ax_main.errorbar(real_tr, real_me, yerr=[yerr_lo,yerr_hi], fmt='ko', markersize=5,label='DES')

    plt.savefig(plot_name)


    return 0

def plot_summary_results_h(
        medians_1,medians_2,errors_1=None,
        errors_2=None,refmedian=None,refmedian_2=None,
        referror=None,referror_2=None,
        model_names=None, left_margin=0.2,
        plot_name='summary_plot.png',
        label1='$\\Omega_M$',label2='$w$',
        xmin=None, xmax=None,
        ymin=None, ymax=None,
        figsize=(8,6),line_after=None):
    #(6, 8)
    if type(figsize)==type(None):
        hgt=float(len(medians_1))+0.5
        wth=8
        figsize=(hth,wgt)


    custom_cycler = (cycler(color=['r', 'g', 'b', 'y','c', 'm', 'y', 'k']))

    plt.rc('axes', prop_cycle=custom_cycler)

    fig = plt.figure(figsize=figsize, edgecolor='black', constrained_layout=True)         
    




    # set up grid of panels in figure

    grid = plt.GridSpec(2, 1, hspace=0.0, wspace=0.0, 
                        left=left_margin, right=0.9, bottom=0.3, top=1.0)#height_ratios=[4], width_ratios=[1,1]
    # main axis: predicted statistic vs. true statistic
    ax_var1 = fig.add_subplot(grid[0,0])
    
    # hist axis [x]: histogram of data on x axis on top
    ax_var2 = fig.add_subplot(grid[1,0], sharex=ax_var1)
    rang=np.flip(np.arange(0,float(len(model_names))/4,0.25))
    #ax_var1.set_yticks(rang, model_names)   
    plt.xticks(rang, model_names)
    if type(referror)!=type(None):
        ax_var1.fill_between(rang,refmedian-referror[0],refmedian+referror[1],alpha=0.3, color='#bcc4d1')
    if type(referror_2)!=type(None):
        ax_var2.fill_between(rang,refmedian_2-referror_2[0],refmedian_2+referror_2[1],alpha=0.3, color='#bcc4d1')

    plt.setp(ax_var1.get_xticklabels(), visible=False)
    if type(errors_1)!=type(None):
        lo,hi=np.transpose(errors_1)
        for i in range(0,len(medians_1)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var1.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')
            ax_var1.errorbar([medians_1[i]], [rang[i]], xerr=[[lo[i]],[hi[i]]], marker='o', markersize=5, linestyle='')
    else:
        for i in range(0,len(medians_1)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var1.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')
            ax_var1.plot([medians_1[i]], [rang[i]], marker='o', markersize=5, linestyle='')

    #print (lo2,hi2)
    #print('essa e o vetor original')
    #print(errors_2)
    if type(errors_2)!=type(None):
        lo2,hi2=np.transpose(errors_2)
        for i in range(0,len(medians_2)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var2.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')

            ax_var2.errorbar([medians_2[i]], [rang[i]], xerr=[[lo2[i]],[hi2[i]]], marker='o', markersize=5, linestyle='')
    else:
        for i in range(0,len(medians_2)):
            if type(line_after)!=type(None):
                if line_after==(i-1):
                    shift= float(rang[i-1]-rang[i])/2
                    ax_var2.axhline(y=rang[i]+shift, ls='-', lw=2.0, color='black')
            ax_var2.plot([medians_2[i]], [rang[i]], marker='o', markersize=5, linestyle='')
    
    if type(refmedian)!=type(None):
        ax_var1.axvline(x=refmedian, ls='dashed', lw=0.5)
    if type(refmedian_2)!=type(None):
        ax_var2.axvline(x=refmedian_2, ls='dashed', lw=0.5)
    ax_var1.set_ylim((-0.1, float(len(model_names))/4)) 
    ax_var1.set_xlabel(label1)  
    ax_var2.set_xlabel(label2)  
    #high_error_var1,low_error_var1=np.transpose(real_data)
    #yerr_lo=real_me - real_lo
    #yerr_hi= real_hi - real_me

    #yerr_lor=realr_me - realr_lo
    #yerr_hir= realr_hi - realr_me

    #res_me=np.divide(real_me-real_tr,real_tr)
    #ax_main.errorbar(real_tr, real_me, yerr=[yerr_lo,yerr_hi], fmt='ko', markersize=5,label='DES')

    plt.savefig(plot_name)


    return 0

def isarraylike(obj):
    if isinstance(obj, str):
        return False
    return isinstance(obj, collections.Sequence)



def plot_sum_single(data,
        medians,errors,
        plot_name='box_plot.png',limits=[-0.15,0.15],
        label1='$\\delta f/ f$',label2='$\%$',
        xmin=None, xmax=None,
        ymin=-0.5, ymax=+0.5,
        ticks=None,shift_in=None,
        figsize=(6, 4),line_after=None, axis=None, xaxis=None, xlabel=None, shift=0, alpha=1.0, hatch='//', color='IndianRed', marker='o', width = 0.35):

    

    #fig = plt.figure()
    
    #ax1 = fig.add_subplot(111)


    #bp4=ax1.boxplot(data, patch_artist=True, notch =False,  showfliers=False)
    
    #ax1.set_xticklabels(ticks,
    #                rotation=45, fontsize=8)

    #plt.savefig('bp1.png')
    #plt.clf()


    # set up grid of panels in figure
    custom_cycler = (cycler(color=['r', 'g', 'b', 'y','c', 'm', 'y', 'k']))
    color_shift=['r', 'g', 'b', 'y','c', 'm', 'y', 'k']

    plt.rc('axes', prop_cycle=custom_cycler)
    #fig = plt.figure(figsize=figsize) 
    if type(axis)==type(None):
        fig = plt.figure(figsize=figsize, edgecolor='black')         
    
        grid = plt.GridSpec(2, 1, hspace=0.0, wspace=0.0, 
                        left=0.15, right=0.9, bottom=0.1, top=0.9)#height_ratios=[4], width_ratios=[1,1]
    
        #grid = plt.GridSpec(2, 1, hspace=0.0, wspace=0.0, 
        #                    left=0.15, right=0.9, bottom=0.1, top=0.9)

        # main axis: predicted statistic vs. true statistic
        ax_var1 = fig.add_subplot(grid[0,0])

  
        # hist axis [x]: histogram of data on x axis on top
        ax_var2 = fig.add_subplot(grid[1, 0], sharex=ax_var1)

    else:
        fig=axis[1]
        ax_var1=axis[0][0]
        ax_var2=axis[0][1]
 


    ax_var1.axhline(y=0.2, ls='dashed', color='magenta', lw=0.5)
    ax_var1.axhline(y=0.1, ls='dashed', color='red', lw=0.5)
    ax_var1.axhline(y=0, ls='dashed', color='black', lw=0.5)
    ax_var1.axhline(y=-0.1, ls='dashed', color='red', lw=0.5)
    ax_var1.axhline(y=-0.2, ls='dashed', color='magenta', lw=0.5)
    if type(xaxis)==type(None):
        xaxis=range(len(medians))
    if type(medians)!=type(None):
        for i in range(0,len(medians)):
            print('Errors',errors[i])
            lo,hi=np.transpose(errors[i])
            print('Errors2',lo,hi)
            ax_var1.errorbar([xaxis[i]+shift],[medians[i]], yerr=[[lo],[hi]], marker=marker,markersize=5, linestyle='', color=color)

    #rang=np.arange(max(xparameter)+1,max(xparameter)+float(len(ext_names))+1,1)#np.flip(
    #rang=rang.tolist()
    #xaxis=np.array(xaxis+rang)
    #xaxisticks=xaxisticks+np.array(ext_names).tolist()
    #ax_var1.set_yticks(rang, model_names)
    #ext_names=xaxis+list(ext_names)   
    ax_var1.set_ylim((ymin, ymax))
    catastrophic=[]
    for i in range(0,len(data)):
      total_pts=len(data[i])
      subset=select_subset_arr(data[i],val=limits[0],condition=">")
      subset=select_subset_arr(subset,val=limits[1],condition="<")
      outl=float(total_pts-len(subset))/float(total_pts)*100
      catastrophic.append(outl) 
     
    ax_var2.bar(np.array(xaxis)+shift, catastrophic, width, color=color, alpha=alpha, hatch=hatch)
    #ax_var2.bar(ind, womenMeans, width,
    #         bottom=menMeans, yerr=womenStd)

    ax_var1.set_ylabel(label1, fontsize=12)  
    ax_var2.set_ylabel(label2, fontsize=12)
    if type(xlabel)!=type(None):
        ax_var2.set_xlabel(xlabel, fontsize=12)


    if type(ticks)!=type(None):
        plt.xticks(xaxis, ticks) #rotation=45
    if plot_name!='':
        plt.savefig(plot_name)
    return [ax_var1,ax_var2],fig


#    ax2 = fig.add_subplot(111)
#    bp2=ax2.boxplot(data, patch_artist=True, conf_intervals=errors,usermedians=medians, notch =False,  showfliers=False)
#    plt.savefig('bp2.png')
#    plt.cla()
#    ax3 = fig.add_subplot(111)
#    bp3=ax3.boxplot(data, patch_artist=True, conf_intervals=errors, notch =False,  showfliers=False)
#    plt.savefig('bp3.png')
#    plt.cla()
#    ax4 = fig.add_subplot(111)
#    bp4=ax4.boxplot(data, patch_artist=True, notch =False,  showfliers=False)
#    plt.savefig('bp4.png')


def plot_summary_asfunction(
        medians_1,medians_2,xparameter,errors_1=None,
        errors_2=None,refmedian=None,refmedian_2=None,
        referror=None,referror_2=None,
        medians_ext1=None,medians_ext2=None,errors_ext1=None,
        errors_ext2=None,
        ext_names=None, left_margin=0.2,
        plot_name='summaryasfunc_plot.png',
        label1='$\\Omega_M$',label2='$w$',
        xmin=None, xmax=None,
        ymin=None, ymax=None,
        horizontal=False,shift_in=None,
        figsize=None,colors=None,markers=None,ext_color='k',ext_marker="o",tex=True,xlabel='',mks=7):
    #(6, 8)
    #if type(figsize)==type(None):
    #    hgt=float(len(medians_1))+0.5
    #    wth=8
    #    figsize=(wth,hgt)

    medians_1=np.array(medians_1)
    medians_2=np.array(medians_2)
    xaxis=np.array(xparameter)
    
    custom_cycler = (cycler(color=['r', 'g', 'b', 'y','c', 'm', 'y', 'k']))
    color_shift=['r', 'g', 'b', 'y','c', 'm', 'y', 'k']
    marker_shifts=[".",",","o","v","^","<",">","1","2","3", "4","8","s","p","P","*"]


    plt.rc('axes', prop_cycle=custom_cycler)

    fig = plt.figure(figsize=(8, 6), edgecolor='black', constrained_layout=True)         
    
    


    # set up grid of panels in figure


    grid = plt.GridSpec(2, 1, hspace=0.0, wspace=0.0, 
                        left=0.15, right=0.9, bottom=0.1, top=0.9)#height_ratios=[4], width_ratios=[1,1]
    
    #grid = plt.GridSpec(2, 1, hspace=0.0, wspace=0.0, 
    #                    left=0.15, right=0.9, bottom=0.1, top=0.9)


    # main axis: predicted statistic vs. true statistic
    ax_var1 = fig.add_subplot(grid[0,0])
    
    # hist axis [x]: histogram of data on x axis on top
    ax_var2 = fig.add_subplot(grid[1, 0], sharex=ax_var1)
    xaxisticks=(np.array(xaxis).astype(str)).tolist()
    if tex==True:
       xaxisticks=list(['$'+s+'$' for s in xaxisticks]) 
    xaxis=xaxis.tolist()
    
    if type(ext_names)!=type(None):
        rang=np.arange(max(xparameter)+1,max(xparameter)+float(len(ext_names))+1,1)#np.flip(
        rang=rang.tolist()
        xaxis=np.array(xaxis+rang)
        xaxisticks=xaxisticks+np.array(ext_names).tolist()
        #ax_var1.set_yticks(rang, model_names)
        #ext_names=xaxis+list(ext_names)   
        plt.xticks(xaxis, xaxisticks)
    xaxis=xaxis.tolist()
    xaxis_aux=xaxis[:]
    ax_var1.axvline(max(xparameter)+0.5, ls='-', lw=2.0, color='black', zorder=10) 
    ax_var2.axvline(max(xparameter)+0.5, ls='-', lw=2.0, color='black', zorder=10)       



    if type(errors_1)!=type(None):
        for i in range(0,len(medians_1)):
            if not(isarraylike(medians_1[i])):
            
                lo,hi=np.transpose(errors_1[i])
                multiple=xaxis.count(xaxis[i])
                multiple_aux=xaxis_aux.count(xaxis[i])
                un=np.unique(xaxis)
                if multiple>1:
                    if len(un)>0 and type(shift_in)==type(None): 
                        shift=float(un[1]-un[0])/(2*multiple)
                    elif type(shift_in)!=type(None):
                        shift=shift_in
                    else:
                        shift=0 
                    col_ind=(multiple-multiple_aux)
                    shift_index_sc=(multiple-multiple_aux)-(float(multiple)/2)
                    xaxis_aux.remove(xaxis[i])
                    print('1- shift index: ',shift_index_sc,'shift: ',shift)
                else:
                    col_ind=0
                    shift=0
                    shift_index_sc=1
                    #print('xaxis ',xaxis)
                #print('multiple: ',multiple,'xaxis1: ',xaxis[1],'xaxis2: ',xaxis[0])
                #print('xaxis ',xaxis)
                #print('x axis with SHIFT')
                #print(xaxis[i]+(shift_index_sc*shift))
                if type(colors)!=type(None):
                    cl=colors[i]
                else:
                    cl=color_shift[col_ind]

                if type(markers)!=type(None):
                    mk=markers[i]
                else:
                    mk=markers_shift[col_ind]

                ax_var1.errorbar([xaxis[i]+(shift_index_sc*shift)],[medians_1[i]], yerr=[[lo],[hi]], marker=mk, color=cl,markersize=mks, linestyle='', zorder=10)
                #ax_var1.errorbar([xaxis[i]],[medians_1[i]], xerr=[[lo],[hi]], marker='o', markersize=5, linestyle='')
            else:
                shift=float(xaxis[0]-xaxis[1])/(2*len(medians_1[i]))
                shift_index=range(len(medians_1[i])) - (float(len(medians_1[i]))/2.) 
                for j in range(0,len(medians_1[i])):
                    lo,hi=np.transpose(errors_1[i])
                    ax_var1.errorbar([xaxis[i]+(shift_index[j]*shift)],[medians_1[i,j]], yerr=[[lo[j]],[hi[j]]], marker='o', markersize=mks, linestyle='', zorder=10)
    else:
        for i in range(0,len(medians_1)):
            if not(isarraylike(medians_1[i])):
                ax_var1.plot([xaxis[i]],[medians_1[i]], marker='o', markersize=mks, linestyle='', zorder=10)
            else:
                shift=float(xaxis[0]-xaxis[1])/(2*len(medians_1[i]))
                shift_index=range(len(medians_1[i])) - (float(len(medians_1[i]))/2.) 
                for j in range(0,len(medians_1[i])):
                    ax_var1.plot([xaxis[i]+(shift_index[j]*shift)],[medians_1[i,j]], marker='o', markersize=mks, linestyle='', zorder=10)



    xaxis_aux=xaxis[:]


    if type(errors_2)!=type(None):
        
        for i in range(0,len(medians_2)):
            if not(isarraylike(medians_2[i])):
                 
                lo,hi=np.transpose(errors_2[i])
                multiple=xaxis.count(xaxis[i])
                multiple_aux=xaxis_aux.count(xaxis[i])
                if multiple>1:
                    if len(un)>0 and type(shift_in)==type(None): 
                        shift=float(un[1]-un[0])/(2*multiple)
                    elif type(shift_in)!=type(None):
                        shift=shift_in
                    else:
                        shift=0
                    col_ind=(multiple-multiple_aux) 
                    shift_index_sc=(multiple-multiple_aux)-(float(multiple)/2)
                    xaxis_aux.remove(xaxis[i])
                    print('2 - shift index: ',shift_index_sc,'shift: ',shift)
                else:
                    shift=0
                    shift_index_sc=1
                    col_ind=0
                if type(colors)!=type(None):
                    cl=colors[i]
                else:
                    cl=color_shift[col_ind]
                if type(markers)!=type(None):
                    mk=markers[i]
                else:
                    mk=markers_shift[col_ind]

                ax_var2.errorbar([xaxis[i]+(shift_index_sc*shift)],[medians_2[i]], yerr=[[lo],[hi]],color=cl, marker=mk, markersize=mks, linestyle='', zorder=10)
                #ax_var1.errorbar([xaxis[i]],[medians_1[i]], xerr=[[lo],[hi]], marker='o', markersize=5, linestyle='')
            else:
                shift=float(xaxis[0]-xaxis[1])/(2*len(medians_2[i]))
                shift_index=range(len(medians_2[i])) - (float(len(medians_2[i]))/2.) 
                for j in range(0,len(medians_2[i])):
                    lo,hi=np.transpose(errors_2[i])
                    ax_var2.errorbar([xaxis[i]+(shift_index[j]*shift)],[medians_2[i,j]], yerr=[[lo[j]],[hi[j]]], marker='o', markersize=mks, linestyle='', zorder=10)
    else:
        for i in range(0,len(medians_2)):
            if not(isarraylike(medians_2[i])):
                ax_var2.plot([xaxis[i]],[medians_2[i]], marker='o', markersize=mks, linestyle='', zorder=10)
            else:
                shift=float(xaxis[0]-xaxis[1])/(2*len(medians_2[i]))
                shift_index=range(len(medians_2[i])) - (float(len(medians_2[i]))/2.) 
                for j in range(0,len(medians_2[i])):
                    ax_var2.plot([xaxis[i]+(shift_index[j]*shift)],[medians_2[i,j]], marker='o', markersize=mks, linestyle='', zorder=10)

    #plotting the external probes

    if type(errors_ext1)!=type(None):
        for i in range(0,len(medians_ext1)):
            lo,hi=np.transpose(errors_ext1[i])

            if type(ext_color)!=type(None):
                ax_var1.errorbar([rang[i]],[medians_ext1[i]], yerr=[[lo],[hi]],color=ext_color, marker=ext_marker, markersize=mks, linestyle='', zorder=10)
            else:
                ax_var1.errorbar([rang[i]],[medians_ext1[i]], yerr=[[lo],[hi]], marker=ext_marker, markersize=mks, linestyle='', zorder=10)
            #ax_var1.errorbar([xaxis[i]],[medians_1[i]], xerr=[[lo],[hi]], marker='o', markersize=5, linestyle='')
    else:
        for i in range(0,len(medians_2)):
            if type(ext_color)!=type(None):
                ax_var1.plot([rang[i]],[medians_ext1[i]], marker=ext_marker, markersize=mks, linestyle='',color=ext_color)
            else:
                ax_var1.plot([rang[i]],[medians_ext1[i]], marker=ext_marker, markersize=mks, linestyle='')


    if type(errors_ext2)!=type(None):
        for i in range(0,len(medians_ext2)):
            lo,hi=np.transpose(errors_ext2[i])
            if type(ext_color)!=type(None):
                ax_var2.errorbar([rang[i]],[medians_ext2[i]], yerr=[[lo],[hi]],color=ext_color, marker=ext_marker, markersize=mks, linestyle='', zorder=10)
            else:
                ax_var2.errorbar([rang[i]],[medians_ext2[i]], yerr=[[lo],[hi]], marker=ext_marker, markersize=mks, linestyle='', zorder=10)
            #ax_var1.errorbar([xaxis[i]],[medians_1[i]], xerr=[[lo],[hi]], marker='o', markersize=5, linestyle='')
    else:
        for i in range(0,len(medians_2)):
            if type(ext_color)!=type(None):
                ax_var2.plot([rang[i]],[medians_ext2[i]], marker=ext_marker, markersize=mks, linestyle='', color=ext_color)
            else:
                ax_var2.plot([rang[i]],[medians_ext2[i]], marker=ext_marker, markersize=mks, linestyle='')
            




    #print (lo2,hi2)
    #print('essa e o vetor original')
    #print(errors_2)





    plt.setp(ax_var2.get_xticklabels(), visible=True)
    plt.setp(ax_var2.get_yticklabels(), visible=True)
    plt.setp(ax_var1.get_xticklabels(), visible=True)
    plt.setp(ax_var1.get_yticklabels(), visible=True)   
    ax_var1.axhline(y=refmedian, ls='dashed', lw=0.5)
    ax_var2.axhline(y=refmedian_2, ls='dashed', lw=0.5)
    ax_var1.set_ylim((-0.05, 0.6)) 
    ax_var2.set_ylim((-4.0, 0.1)) 
    ax_var1.set_ylabel(label1)  
    ax_var2.set_ylabel(label2)
    ax_var2.set_xlabel(xlabel)    

    xmin1, xmax1 = ax_var1.get_xlim()
    if type(referror)!=type(None):
        ax_var1.fill_between(x=list([xmin1])+xaxis+list([xmax1]),y1=refmedian-referror[0],y2=refmedian+referror[1],alpha=0.3, color='#A9A9A9', zorder=1)  #bcc4d1
    xmin2, xmax2 = ax_var2.get_xlim()
    if type(referror_2)!=type(None):
        ax_var2.fill_between(x=list([xmin2])+xaxis+list([xmax2]),y1=refmedian_2-referror_2[0],y2=refmedian_2+referror_2[1],alpha=0.3, color='#A9A9A9', zorder=1)
    #xmin1, xmax1 = ax_var1.get_xlim()
    #xmin2, xmax2 = ax_var2.get_xlim()

    #ax_var1.set_xlim((xmax1, xmin1)) 
    ax_var2.set_xlim((xmax1, xmin1)) 


    #high_error_var1,low_error_var1=np.transpose(real_data)
    #yerr_lo=real_me - real_lo
    #yerr_hi= real_hi - real_me

    #yerr_lor=realr_me - realr_lo
    #yerr_hir= realr_hi - realr_me

    #res_me=np.divide(real_me-real_tr,real_tr)
    #ax_main.errorbar(real_tr, real_me, yerr=[yerr_lo,yerr_hi], fmt='ko', markersize=5,label='DES')

    plt.savefig(plot_name)


    return 0


def density_scatter( x , y, ax = None,fig=None, sort = True, bins = 20, alpha = 0.05, s=0.8,marker='o', out_name="testedense.png", save=True, color_bar=True, **kwargs )   :
    """
    Scatter plot colored by 2d histogram
    """
    from scipy.stats import gaussian_kde

    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins)

    #xy = np.vstack([x,y])
    #z = gaussian_kde(xy)(xy)

    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )
    #z = interpn( ( 0.05*(x_e[1:] + x_e[:-1]) , 0.05*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )
    #print("HIST 2D")
    #print( 0.05*(x_e[1:] + x_e[:-1])) 
    #print( 0.05*(y_e[1:]+y_e[:-1]))
    #print(x_e)
    #print(y_e)
    #print (data[0][0])
    #print (data[14][14])
    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    
    cax=ax.scatter( x, y, c=z,s=s,alpha=alpha )
    if color_bar==True: 
        fig.colorbar(cax,ax=ax)

    if save==True:
        plt.savefig("testedense.png")

    #ax.plot( x, y, c=z, **kwargs )
    return ax
def radec_density(ra,dec, ax = None,fig=None,sort = True, bins = 20, alpha = 0.05, s=0.8,marker='o', out_name="testedense.png", save=True, color_bar=True ):
    if ax is None :
        fig , ax = plt.subplots(projection="mollweide")
    ra = coord.Angle(ra*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(dec*u.degree)
    print(ax)
    print(fig)
    ax=density_scatter( x=ra.radian , y=dec.radian, ax = ax,fig=fig, sort = sort, bins = bins, alpha = alpha, s=s,marker=marker, save=False, color_bar=color_bar)
    #ax.scatter(ra.radian, dec.radian)
    ax.grid(True)
    if save==True:
        plt.savefig("testedense.png")

    #ax.plot( x, y, c=z, **kwargs )
    return ax



def plot_deviation_results(
        results,results_residuals,
        truth_values,individual_truth_values,
        individual_medians,plot_std=False,
        bin_edges=10,real_data=None,
        plot_dropout_err_sigmas=[True,True,True], plot_name='regression.png',
        xlabel='$\\theta_{E}$',ylabel='$\\theta_{E}$',
        xmin=None, xmax=None,
        ymin=None, ymax=None,
        color_std='magenta', color_err='green',
        color_pop='blue',plot_label=None,
        fill_label_std=None, fill_label_err=None,
        figsize=(6, 8),alpha_hist=0.4, 
        pop_error=None,pop_error_residuals=None,
        plot_real_residuals=False,plot_hist=True, color_median='gold', fontsize=14, fig=None, ax_diff=None, ax_hist_x=None, save=True, ylimits=[-0.8,0.8], ref_lines=True,is_rmse=False,plot_dropout_err_alpha=[0.5,0.4,0.3], multiple_entry=False, multiple_true_vals=False):
    

    multi_colors=['indianred','navy','lime','gold','purple','goldenrod','teal','b','olivedrab','r','c','m','dimgray','palevioletred','y']
    multi_ls=['solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1))]
    #'solid',dashdot,dotted,densely dashdot , densely dot

    if multiple_entry==True:
        xs_bin=[]
        #ys_bin=[]
        #ys_bin_std=[]
        #ys_bin_percentile1sig=[]
        #ys_bin_percentile2sig=[]
        #ys_bin_percentile3sig=[]
        #imeds=[]
        ysr_bin=[]
        ysr_bin_std=[]
        ysr_bin_percentile1sig=[]
        ysr_bin_percentile2sig=[]
        ysr_bin_percentile3sig=[]
        #yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig

    if multiple_entry==True:
        pop_error=None
        pop_error_residuals=None
        for i in range(0,len(results)):
            #x_bin,y_bin,y_bin_std,y_bin_percentile1sig,y_bin_percentile2sig,y_bin_percentile3sig=results[i]
            x_bin,yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig=results_residuals[i]
            xs_bin.append(x_bin)
            #ys_bin.append(y_bin)
            #ys_bin_std.append(y_bin_std)
            #ys_bin_percentile1sig.append(y_bin_percentile1sig)
            #ys_bin_percentile2sig.append(y_bin_percentile2sig)
            #ys_bin_percentile3sig.append(y_bin_percentile3sig)
            
            ysr_bin.append(yr_bin)
            ysr_bin_std.append(yr_bin_std)
            ysr_bin_percentile1sig.append(yr_bin_percentile1sig)
            ysr_bin_percentile2sig.append(yr_bin_percentile2sig)
            ysr_bin_percentile3sig.append(yr_bin_percentile3sig)
    else:
        #x_bin,y_bin,y_bin_std,y_bin_percentile1sig,y_bin_percentile2sig,y_bin_percentile3sig=results
        x_bin,yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig=results_residuals 









    #pop_err=pop_error
    pop_err_r=pop_error_residuals
    #x_bin,y_bin,y_bin_std,y_bin_percentile1sig,y_bin_percentile2sig,y_bin_percentile3sig=results
    #x_bin,yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig=results_residuals


    if type(xmin)==type(None):
        if type(bin_edges)!=type(10):
            xmin=min(bin_edges)
            xmax=max(bin_edges)
        else:
            xmin=min(x_bin)
            xmax=max(x_bin)



    if multiple_entry==True:
        imeds=individual_medians
    else:
        imed=individual_medians

    if (multiple_entry==True):
        itruths=individual_truth_values
    else:
        itruth=individual_truth_values



    #imed=individual_medians
    #itruth=individual_truth_values 

    #ty_bin_percentile1sig=np.transpose(y_bin_percentile1sig)
    #ty_bin_percentile2sig=np.transpose(y_bin_percentile2sig)
    #ty_bin_percentile3sig=np.transpose(y_bin_percentile3sig)

    #tyr_bin_percentile1sig=np.transpose(yr_bin_percentile1sig)
    #tyr_bin_percentile2sig=np.transpose(yr_bin_percentile2sig)
    #tyr_bin_percentile3sig=np.transpose(yr_bin_percentile3sig)
    

    if multiple_entry==True:
        #tys_bin_percentile1sig=[]
        #tys_bin_percentile2sig=[]
        #tys_bin_percentile3sig=[]

        tysr_bin_percentile1sig=[]
        tysr_bin_percentile2sig=[]
        tysr_bin_percentile3sig=[]

        for i in range(0,len(xs_bin)):
            #tys_bin_percentile1sig.append(np.transpose(ys_bin_percentile1sig[i]))
            #tys_bin_percentile2sig.append(np.transpose(ys_bin_percentile2sig[i]))
            #tys_bin_percentile3sig.append(np.transpose(ys_bin_percentile3sig[i]))

            tysr_bin_percentile1sig.append(np.transpose(ysr_bin_percentile1sig[i]))
            tysr_bin_percentile2sig.append(np.transpose(ysr_bin_percentile2sig[i]))
            tysr_bin_percentile3sig.append(np.transpose(ysr_bin_percentile3sig[i]))

    else:
        #ty_bin_percentile1sig=np.transpose(y_bin_percentile1sig)
        #ty_bin_percentile2sig=np.transpose(y_bin_percentile2sig)
        #ty_bin_percentile3sig=np.transpose(y_bin_percentile3sig)

        tyr_bin_percentile1sig=np.transpose(yr_bin_percentile1sig)
        tyr_bin_percentile2sig=np.transpose(yr_bin_percentile2sig)
        tyr_bin_percentile3sig=np.transpose(yr_bin_percentile3sig)







    if type(fig)==type(None):
       fig = plt.figure(figsize=figsize)         
    

    # set up grid of panels in figure
    grid = plt.GridSpec(2, 1, hspace=0.0, wspace=0.0, 
                        left=0.15, right=0.9, bottom=0.1, top=0.9)#,
    #                    height_ratios=[1, 4, 1], width_ratios=[4,1])
         
    # main axis: predicted statistic vs. true statistic
    #ax_main = fig.add_subplot(grid[1,0])
    
    # res axis: residual of statistic vs. true stratistic
    if type(ax_diff)==type(None):
        ax_diff = fig.add_subplot(grid[1,0])#, sharex=ax_main)

    # hist axis [x]: histogram of data on x axis on top
    if plot_hist==True:
        if type(ax_hist_x)==type(None):
            ax_hist_x = fig.add_subplot(grid[0, 0],sharex=ax_diff)
    else:
        ax_hist_x=None
    # hist axis [y]: histogram of data on y axis on right
    #ax_hist_y = fig.add_subplot(grid[1, 1], sharey=ax_main)
    



    # calculate upper and lower bounds [for systematic error bars]
    #if plot_std==True:
    #    lo_err, up_err = measure_uplow_std_limits(y_bin, y_bin_std, nsigma=1)#(ys_bin, ys_err_bin)

    # calculate upper and lower bounds on diff [for systematic error bars]
    if (plot_std==True) and (multiple_entry==False):
        lo_err_diff, up_err_diff = measure_uplow_std_limits(yr_bin, yr_bin_std, nsigma=1)
        

    # ================================
    # Plot data [sim]
    # ================================
    # plot: dropout medians [main axis]
    #ax_main.plot(x_bin,y_bin, color=color_err, alpha=0.9, label=plot_label)
    
    # plot: dropouts percentile errors [main axis]
    #if plot_dropout_err==True:
    #    ax_main.fill_between(x_bin, ty_bin_percentile1sig[0], ty_bin_percentile1sig[1], color=color_err, alpha=0.5, label=fill_label_err)
    #    ax_main.fill_between(x_bin, ty_bin_percentile2sig[0], ty_bin_percentile2sig[1], color=color_err, alpha=0.4)
    #    ax_main.fill_between(x_bin, ty_bin_percentile3sig[0], ty_bin_percentile3sig[1], color=color_err, alpha=0.3)
    
    # plot: one-to-one line [main axis]
    #ax_main.plot([xmin, xmax], [ymin, ymax], color='black', ls='dashed', alpha=0.9)
    
    # plot: errors from standard deviation [main axis]
    #if plot_std==True:
    #    ax_main.fill_between(x_bin, lo_err, up_err, color=color_std, alpha=0.4, label=fill_label_std)
    
    # plot: deltas [diff axis]
    #ax_diff.plot(x_bin,yr_bin, color=color_median, alpha=0.9)
    

    if multiple_entry==True:
       for i in range(0,len(xs_bin)):
            ax_diff.plot(xs_bin[i],ysr_bin[i], color=multi_colors[i], alpha=1.0, ls=multi_ls[i])
    else:
            ax_diff.plot(x_bin,yr_bin, color=color_median, alpha=0.9)
    


    # plot: statistical errors [diff axis]

    if multiple_entry==True:
        for i in range(0,len(xs_bin)):
            if (plot_dropout_err_sigmas[0]==True):
                ax_diff.fill_between(xs_bin[i], tysr_bin_percentile1sig[i][0], tysr_bin_percentile1sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[0])
            if (plot_dropout_err_sigmas[1]==True):    
                ax_diff.fill_between(xs_bin[i], tysr_bin_percentile2sig[i][0], tysr_bin_percentile2sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[1])
            if (plot_dropout_err_sigmas[2]==True):    
                ax_diff.fill_between(xs_bin[i], tysr_bin_percentile3sig[i][0], tysr_bin_percentile3sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[2])
    else:

        if (plot_dropout_err_sigmas[0]==True):
            ax_diff.fill_between(x_bin, tyr_bin_percentile1sig[0], tyr_bin_percentile1sig[1], color=color_err, alpha=plot_dropout_err_alpha[0], label=fill_label_err)
        if (plot_dropout_err_sigmas[1]==True):    
            ax_diff.fill_between(x_bin, tyr_bin_percentile2sig[0], tyr_bin_percentile2sig[1], color=color_err, alpha=plot_dropout_err_alpha[1])
        if (plot_dropout_err_sigmas[2]==True):    
            ax_diff.fill_between(x_bin, tyr_bin_percentile3sig[0], tyr_bin_percentile3sig[1], color=color_err, alpha=plot_dropout_err_alpha[2])




    #if plot_dropout_err_sigmas[0]==True:
    #    ax_diff.fill_between(x_bin, tyr_bin_percentile1sig[0], tyr_bin_percentile1sig[1], color=color_err, alpha=0.5, label=fill_label_err)
    #if plot_dropout_err_sigmas[1]==True:        
    #    ax_diff.fill_between(x_bin, tyr_bin_percentile2sig[0], tyr_bin_percentile2sig[1], color=color_err, alpha=0.4)
    #if plot_dropout_err_sigmas[2]==True:    
    #    ax_diff.fill_between(x_bin, tyr_bin_percentile3sig[0], tyr_bin_percentile3sig[1], color=color_err, alpha=0.3)
    
    # plot: systematic errors [main axis]
    if plot_std==True and (multiple_entry==False):
        ax_diff.fill_between(x_bin, lo_err_diff, up_err_diff, color=color_err, alpha=0.9, label=fill_label_std)
    
    # plot: errors from medians
    if type(pop_error)!=type(None):
        #xp_bin,yp_bin,yp_bin_std,yp_bin_percentile1sig,yp_bin_percentile2sig,yp_bin_percentile3sig=pop_err
        xp_bin,ypr_bin,ypr_bin_std,ypr_bin_percentile1sig,ypr_bin_percentile2sig,ypr_bin_percentile3sig=results_residuals=pop_err_r
        #typ_bin_percentile1sig=np.transpose(yp_bin_percentile1sig)
        typr_bin_percentile1sig=np.transpose(ypr_bin_percentile1sig)
        #ax_main.fill_between(xp_bin, typ_bin_percentile1sig[0], typ_bin_percentile1sig[1], color=color_pop, alpha=0.4)
        ax_diff.fill_between(xp_bin, typr_bin_percentile1sig[0], typr_bin_percentile1sig[1], color=color_pop, alpha=0.4)
    # plot: zero line [res axis]
    xmin_line, xmax_line = ax_diff.get_xlim()
    if ref_lines==True:
        ax_diff.axhline(y=0.2, ls='dashed', color='magenta', lw=0.5)
        ax_diff.axhline(y=0.1, ls='dashed', color='red', lw=0.5)
        ax_diff.axhline(y=0, ls='dashed', color='black', lw=0.5)
        ax_diff.axhline(y=-0.1, ls='dashed', color='red', lw=0.5)
        ax_diff.axhline(y=-0.2, ls='dashed', color='magenta', lw=0.5)

    if plot_hist==True:
    # ================================
    # Plot data [simulation]
    # ================================
    # plot: histogram on upper x axis

        if (multiple_entry==True) and (multiple_true_vals==True):
            for i in range(0,len(xs_bin)):
                ax_hist_x.hist(itruths[i], bin_edges, 
                   color=multi_colors[i], alpha=1.0, 
                   histtype='step', stacked=True)

        elif (multiple_entry==True) and (multiple_true_vals==False):
            ax_hist_x.hist(itruths[0], bin_edges, 
                   color=color_err, alpha=1.0, 
                   histtype='step')

        else:
            ax_hist_x.hist(itruth, bin_edges, 
                   color=color_err, alpha=alpha_hist, 
                   histtype='stepfilled')






        #ax_hist_x.hist(itruth, bin_edges, 
        #           color=color_err, alpha=alpha_hist, 
        #           histtype='stepfilled')
        #ax_hist_x.grid(True,alpha=0.2)

        plt.setp(ax_hist_x.get_yticklabels(), visible=False)
        plt.setp(ax_hist_x.get_xticklabels(), visible=False)
    #ax_diff.set_xlim((xmin, xmax))   
    #ax_diff.set_ylim((xmin, xmax))
    #ax_diff.set_ylabel(ylabel)   
    #ax_diff.grid(True, alpha=0.2) 

    #ax_hist_y.grid(True,alpha=0.2)
    #ax_main.legend()
    #ax_main.set_aspect('equal')
    
    # set plot options [residual axis]
    plt.setp(ax_diff.get_xticklabels(), visible=True,fontsize=fontsize)
    plt.setp(ax_diff.get_yticklabels(), visible=True,fontsize=fontsize)
    ax_diff.set_ylim((ylimits[0],ylimits[1]))
    ax_diff.set_xlabel(xlabel, fontsize=fontsize)
    ylabel_plain=ylabel.replace('$','')
    if is_rmse==False:
        ax_diff.set_ylabel('$\\frac{\\delta '+ylabel_plain+'}{'+ylabel_plain+'}$', rotation=0, fontsize=fontsize)
    else:
        ax_diff.set_ylabel('$RMSE_{'+ylabel_plain+'}$', rotation=90, fontsize=fontsize)
    ax_diff.grid(True, alpha=0.2)

    if save==True:
        plt.savefig(plot_name)
        return None,None,None
    else:
        return fig,ax_diff, ax_hist_x


def plot_regression_results(
        results,results_residuals,
        truth_values,individual_truth_values,
        individual_medians,plot_std=False,
        bin_edges=10,real_data=None,
        plot_dropout_err_sigmas=[True,True,True], plot_name='regression.png',
        xlabel='$\\theta_{E}$',ylabel='$\\theta_{E}$',
        xmin=None, xmax=None,
        ymin=None, ymax=None, resylim=(-0.8,0.8),
        color_std='magenta', color_err='green',color_median='gold',
        color_pop='blue',dashedlines=[0.0,0.1,0.2],plot_label=None,
        fill_label_std=None, fill_label_err=None,
        figsize=(6, 8),alpha_hist=0.4, 
        pop_error=None,pop_error_residuals=None,ylimits=None,
        plot_real_residuals=False, bin_edges_y=None, fontsize=14,
        res_axis=True, plot_dropout_err_alpha=[0.5,0.4,0.3],
        median_ls='-', real_label=None, multiple_entry=False, multiple_true_vals=False, vline=None):
    """ 
    Funtion to create Plots of true values vs predicted values.
    It also plots residuals and the histogram of true and 
    predicted.

    INPUT:
        results: <list> - a list with the format [bin centers, 
        bin medians, bin standard deviations,bin 1 sigma 
        intervals, i.e., [low,high] limit ,bin 2 sigma intervals 
        ,bin 3 sigma intervals].

        results_residuals: <list> - a list with the format [bin 
        centers, bin median residuals, bin standard deviation
        residuals,bin 1 sigma interval of residuals,bin 2 sigma 
        intervals of residuals, bin 3 sigma intervals of 
        residuals].

        truth_values: <array> - same as np.unique
        (individual_truth_values).

        individual_truth_value: <array> - The non bined truth 
        values for histogram.

        individual_medians: <array> - The non bined predicted 
        medians for histogram.

        plot_std <bool> - If True will plot a standard deviation
        in results and results_residuals

        bin_edges: <array> - Bin edges, used to define a scale 
        only, optional.

        real_data: data points, instead of binned curves to add
        in the plot. if plot_real_residuals=True it must follow 
        the format [true_valules,median,lowlim,highlim,
        median_residual,lowlim_residual,highlim_residual] if 
        plot_real_residuals=False the format is [true_valules,
        median,lowlim,highlim]

        plot_dropout_err: <bool> -  if True will plot the 
        1-sigma, 2-sigma and 3-sigma error in results_catalog 
        and results_residuals.
 
        plot_name: <str> : output plot name.


    OUTPUT:
       Saves a plot with a a path defined in the input plot_name.



    """

    multi_colors=['indianred','navy','g','gold','purple','goldenrod','teal','b','olivedrab','r','c','m','dimgray','palevioletred','y']
    multi_colorsm=['indianred','navy','lime','gold','purple','goldenrod','teal','b','olivedrab','r','c','m','dimgray','palevioletred','y']
    multi_ls=['solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1)),'solid','dashdot','dotted',(0, (3, 1, 1, 1)),(0, (1, 1))]
    #'solid',dashdot,dotted,densely dashdot , densely dot

    if multiple_entry==True:
        xs_bin=[]
        ys_bin=[]
        ys_bin_std=[]
        ys_bin_percentile1sig=[]
        ys_bin_percentile2sig=[]
        ys_bin_percentile3sig=[]
        #imeds=[]
        ysr_bin=[]
        ysr_bin_std=[]
        ysr_bin_percentile1sig=[]
        ysr_bin_percentile2sig=[]
        ysr_bin_percentile3sig=[]
        #yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig

    if multiple_entry==True:
        pop_error=None
        pop_error_residuals=None
        for i in range(0,len(results)):
            x_bin,y_bin,y_bin_std,y_bin_percentile1sig,y_bin_percentile2sig,y_bin_percentile3sig=results[i]
            x_bin,yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig=results_residuals[i]
            xs_bin.append(x_bin)
            ys_bin.append(y_bin)
            ys_bin_std.append(y_bin_std)
            ys_bin_percentile1sig.append(y_bin_percentile1sig)
            ys_bin_percentile2sig.append(y_bin_percentile2sig)
            ys_bin_percentile3sig.append(y_bin_percentile3sig)
            
            ysr_bin.append(yr_bin)
            ysr_bin_std.append(yr_bin_std)
            ysr_bin_percentile1sig.append(yr_bin_percentile1sig)
            ysr_bin_percentile2sig.append(yr_bin_percentile2sig)
            ysr_bin_percentile3sig.append(yr_bin_percentile3sig)
    else:
        x_bin,y_bin,y_bin_std,y_bin_percentile1sig,y_bin_percentile2sig,y_bin_percentile3sig=results
        x_bin,yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig=results_residuals 


        

    pop_err=pop_error
    pop_err_r=pop_error_residuals
    #x_bin,y_bin,y_bin_std,y_bin_percentile1sig,y_bin_percentile2sig,y_bin_percentile3sig=results
    #x_bin,yr_bin,yr_bin_std,yr_bin_percentile1sig,yr_bin_percentile2sig,yr_bin_percentile3sig=results_residuals


    if type(xmin)==type(None):
        if type(bin_edges)!=type(10):
            xmin=min(bin_edges)
            xmax=max(bin_edges)
        else:
            xmin=min(x_bin)
            xmax=max(x_bin)
    if type(bin_edges_y)==type(None):
        bin_edges_y=bin_edges

   
    if multiple_entry==True:
        imeds=individual_medians
    else:
        imed=individual_medians

    if (multiple_entry==True):
        itruths=individual_truth_values
    else:
        itruth=individual_truth_values

     



    if multiple_entry==True:
        tys_bin_percentile1sig=[]
        tys_bin_percentile2sig=[]
        tys_bin_percentile3sig=[]

        tysr_bin_percentile1sig=[]
        tysr_bin_percentile2sig=[]
        tysr_bin_percentile3sig=[]

        for i in range(0,len(xs_bin)):
            tys_bin_percentile1sig.append(np.transpose(ys_bin_percentile1sig[i]))
            tys_bin_percentile2sig.append(np.transpose(ys_bin_percentile2sig[i]))
            tys_bin_percentile3sig.append(np.transpose(ys_bin_percentile3sig[i]))

            tysr_bin_percentile1sig.append(np.transpose(ysr_bin_percentile1sig[i]))
            tysr_bin_percentile2sig.append(np.transpose(ysr_bin_percentile2sig[i]))
            tysr_bin_percentile3sig.append(np.transpose(ysr_bin_percentile3sig[i]))

    else:
        ty_bin_percentile1sig=np.transpose(y_bin_percentile1sig)
        ty_bin_percentile2sig=np.transpose(y_bin_percentile2sig)
        ty_bin_percentile3sig=np.transpose(y_bin_percentile3sig)

        tyr_bin_percentile1sig=np.transpose(yr_bin_percentile1sig)
        tyr_bin_percentile2sig=np.transpose(yr_bin_percentile2sig)
        tyr_bin_percentile3sig=np.transpose(yr_bin_percentile3sig)




    
    # ================================
    # Initialize figure canvas
    # ================================
    if xlabel==ylabel:
        ymax = xmax
        ymin = xmin
    # set up figure
    fig = plt.figure(figsize=figsize)         
    
    # set up grid of panels in figure
    if res_axis==True:
        grid = plt.GridSpec(3, 2, hspace=0.0, wspace=0.0, 
                        left=0.1, right=1.0, bottom=0.1, top=1.0,
                        height_ratios=[1, 4, 1], width_ratios=[4,1])
    else:
        grid = plt.GridSpec(2, 2, hspace=0.0, wspace=0.0, 
                        left=0.1, right=1.0, bottom=0.1, top=1.0,
                        height_ratios=[1, 4], width_ratios=[4,1])
         
    # main axis: predicted statistic vs. true statistic
    ax_main = fig.add_subplot(grid[1,0])
    
    # hist axis [x]: histogram of data on x axis on top
    ax_hist_x = fig.add_subplot(grid[0, 0], sharex=ax_main)
    
    # hist axis [y]: histogram of data on y axis on right
    ax_hist_y = fig.add_subplot(grid[1, 1], sharey=ax_main)
    
    # res axis: residual of statistic vs. true stratistic
    if res_axis==True:
        ax_diff = fig.add_subplot(grid[2,0], sharex=ax_main)
    

 
    
    # calculate upper and lower bounds [for systematic error bars]
    if (plot_std==True) and (multiple_entry==False):
        lo_err, up_err = measure_uplow_std_limits(y_bin, y_bin_std, nsigma=1)#(ys_bin, ys_err_bin)

    # calculate upper and lower bounds on diff [for systematic error bars]
    if (plot_std==True) and (multiple_entry==False):
        lo_err_diff, up_err_diff = measure_uplow_std_limits(yr_bin, yr_bin_std, nsigma=1)
        

    # ================================
    # Plot data [sim]
    # ================================
    # plot: dropout medians [main axis]

    if multiple_entry==True:
        for i in range(0,len(xs_bin)):
            ax_main.plot(xs_bin[i],ys_bin[i], color=multi_colorsm[i], alpha=1.0, ls=multi_ls[i])
    else:
        ax_main.plot(x_bin,y_bin, color=color_median, alpha=0.9, label=plot_label, ls=median_ls)
    
    # plot: dropouts percentile errors -  systematic errors from Neural Network [main axis]
    if multiple_entry==True:
        for i in range(0,len(xs_bin)):
            if plot_dropout_err_sigmas[0]==True:
                ax_main.fill_between(xs_bin[i], tys_bin_percentile1sig[i][0], tys_bin_percentile1sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[0])
            if plot_dropout_err_sigmas[1]==True:    
                ax_main.fill_between(xs_bin[i], tys_bin_percentile2sig[i][0], tys_bin_percentile2sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[1])
            if plot_dropout_err_sigmas[2]==True:
                ax_main.fill_between(xs_bin[i], ty_bin_percentile3sig[i][0], ty_bin_percentile3sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[2])


    else:
        if plot_dropout_err_sigmas[0]==True:
            ax_main.fill_between(x_bin, ty_bin_percentile1sig[0], ty_bin_percentile1sig[1], color=color_err, alpha=plot_dropout_err_alpha[0], label=fill_label_err)
        if plot_dropout_err_sigmas[1]==True:    
            ax_main.fill_between(x_bin, ty_bin_percentile2sig[0], ty_bin_percentile2sig[1], color=color_err, alpha=plot_dropout_err_alpha[1])
        if plot_dropout_err_sigmas[2]==True:
            ax_main.fill_between(x_bin, ty_bin_percentile3sig[0], ty_bin_percentile3sig[1], color=color_err, alpha=plot_dropout_err_alpha[2])
    
    # plot: one-to-one line [main axis]
    ax_main.plot([xmin, xmax], [ymin, ymax], color='black', ls='dashed', alpha=0.9, lw=0.5)
    



    # plot: errors from standard deviation [main axis]
    if (plot_std==True) and (multiple_entry==False):
        ax_main.fill_between(x_bin, lo_err, up_err, color=color_std, alpha=0.4, label=fill_label_std)

    if type(vline)!=type(None):
        ax_main.axvline(x=vline, ls='dashed', color='red', lw=1.0)

    
    # plot: deltas [diff axis]
    if res_axis==True:
        if multiple_entry==True:
            for i in range(0,len(xs_bin)):
                ax_diff.plot(xs_bin[i],ysr_bin[i], color=multi_colors[i], alpha=1.0, ls=multi_ls[i])
        else:
            ax_diff.plot(x_bin,yr_bin, color=color_median, alpha=0.9)
    
    # plot: dropouts percentile errors - systematic errors from Neural Network [diff axis]



    if multiple_entry==True:
        for i in range(0,len(xs_bin)):
            if (plot_dropout_err_sigmas[0]==True) and (res_axis==True):
                ax_diff.fill_between(xs_bin[i], tysr_bin_percentile1sig[i][0], tysr_bin_percentile1sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[0])
            if (plot_dropout_err_sigmas[1]==True) and (res_axis==True):    
                ax_diff.fill_between(xs_bin[i], tysr_bin_percentile2sig[i][0], tysr_bin_percentile2sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[1])
            if (plot_dropout_err_sigmas[2]==True) and (res_axis==True):    
                ax_diff.fill_between(xs_bin[i], tysr_bin_percentile3sig[i][0], tysr_bin_percentile3sig[i][1], color=multi_colors[i], alpha=plot_dropout_err_alpha[2])
    else:

        if (plot_dropout_err_sigmas[0]==True) and (res_axis==True):
            ax_diff.fill_between(x_bin, tyr_bin_percentile1sig[0], tyr_bin_percentile1sig[1], color=color_err, alpha=plot_dropout_err_alpha[0], label=fill_label_err)
        if (plot_dropout_err_sigmas[1]==True) and (res_axis==True):    
            ax_diff.fill_between(x_bin, tyr_bin_percentile2sig[0], tyr_bin_percentile2sig[1], color=color_err, alpha=plot_dropout_err_alpha[1])
        if (plot_dropout_err_sigmas[2]==True) and (res_axis==True):    
            ax_diff.fill_between(x_bin, tyr_bin_percentile3sig[0], tyr_bin_percentile3sig[1], color=color_err, alpha=plot_dropout_err_alpha[2])






    # plot: systematic errors [main axis]
    if (plot_std==True) and (res_axis==True) and (multiple_entry==False):
        ax_diff.fill_between(x_bin, lo_err_diff, up_err_diff, color=color_err, alpha=0.9, label=fill_label_std)
    
    # plot: errors from medians
    if type(pop_error)!=type(None):
        xp_bin,yp_bin,yp_bin_std,yp_bin_percentile1sig,yp_bin_percentile2sig,yp_bin_percentile3sig=pop_err
        xp_bin,ypr_bin,ypr_bin_std,ypr_bin_percentile1sig,ypr_bin_percentile2sig,ypr_bin_percentile3sig=results_residuals=pop_err_r
        typ_bin_percentile1sig=np.transpose(yp_bin_percentile1sig)
        typr_bin_percentile1sig=np.transpose(ypr_bin_percentile1sig)
        ax_main.fill_between(xp_bin, typ_bin_percentile1sig[0], typ_bin_percentile1sig[1], color=color_pop, alpha=0.4)
        if res_axis==True:  
            ax_diff.fill_between(xp_bin, typr_bin_percentile1sig[0], typr_bin_percentile1sig[1], color=color_pop, alpha=0.4)
    if res_axis==True:
        # plot: zero line [res axis]
        xmin_line, xmax_line = ax_diff.get_xlim()
        #dashedlines
        if multiple_entry==True:
            ax_diff.axhline(y=dashedlines[2], ls='dashed', color='black', lw=0.5)
            ax_diff.axhline(y=dashedlines[1], ls='dashed', color='black', lw=0.5)
            ax_diff.axhline(y=dashedlines[0], ls='dashed', color='black', lw=0.5)
            ax_diff.axhline(y=-1*dashedlines[1], ls='dashed', color='black', lw=0.5)
            ax_diff.axhline(y=-1*dashedlines[2], ls='dashed', color='black', lw=0.5)
        else:
            ax_diff.axhline(y=dashedlines[2], ls='dashed', color='darkviolet', lw=0.5)
            ax_diff.axhline(y=dashedlines[1], ls='dashed', color='peru', lw=0.5)
            ax_diff.axhline(y=dashedlines[0], ls='dashed', color='black', lw=0.5)
            ax_diff.axhline(y=-1*dashedlines[1], ls='dashed', color='peru', lw=0.5)
            ax_diff.axhline(y=-1*dashedlines[2], ls='dashed', color='darkviolet', lw=0.5)

    #ax_diff.plot([xmin, xmax], np.zeros(2), color='black', ls='dashed', alpha=0.4) 

    
    # ================================
    # Plot data [simulation]
    # ================================
    # plot: histogram on upper x axis

    if (multiple_entry==True) and (multiple_true_vals==True):
        for i in range(0,len(xs_bin)):
            ax_hist_x.hist(itruths[i], bin_edges, 
                   color=multi_colors[i], alpha=1.0, 
                   histtype='step', stacked=True)

    elif (multiple_entry==True) and (multiple_true_vals==False):
        ax_hist_x.hist(itruths[0], bin_edges, 
                   color=color_err, alpha=1.0, 
                   histtype='step')

    else:
        ax_hist_x.hist(itruth, bin_edges, 
                   color=color_err, alpha=alpha_hist, 
                   histtype='stepfilled')
        #hist_val,edges=np.histogram(itruth,bin_edges)
    #print('x histogram')
    #print(hist_val)
    # plot: histogram on right y axis
    #hist_val,edges=np.histogram(imed,bin_edges_y)
    #print('y histogram')
    #print(hist_val)

    if multiple_entry==True:
        for i in range(0,len(xs_bin)):
            blank = ax_hist_y.hist(imeds[i], bin_edges_y, 
                   color=multi_colors[i], alpha=1.0,
                   histtype='step', stacked=True,
                   orientation='horizontal')

    else:
        blank = ax_hist_y.hist(imed, bin_edges_y, 
                   color=color_err, alpha=alpha_hist,
                   histtype='stepfilled',
                   orientation='horizontal')
        
    # plot: sky data main_axis
    if type(real_data)!=type(None):
        print('real data')
        print(np.array(real_data).shape)
        if plot_real_residuals==False:
            _id_point,real_tr,real_me,real_lo,real_hi,real_odds,real_frac_err,real_frac_dev,real_half_interval=np.transpose(real_data)
            #real_tr,real_me,real_lo,real_hi=np.transpose(real_data)
        else:
            _id_point,real_tr,real_me,real_lo,real_hi,realr_me,realr_lo,realr_hi,real_odds,real_frac_err,real_frac_dev,real_half_interval=np.transpose(real_data)
            #real_tr,real_me,real_lo,real_hi,realr_me,realr_lo,realr_hi=np.transpose(real_data)
        yerr_lo=real_me - real_lo
        yerr_hi= real_hi - real_me
        if plot_real_residuals==True:
            yerr_lor=realr_me - realr_lo
            yerr_hir= realr_hi - realr_me

        #res_me=np.divide(real_me-real_tr,real_tr)
        ax_main.errorbar(real_tr, real_me, yerr=[yerr_lo,yerr_hi], fmt='ko', markersize=5,label=real_label)
        if (res_axis)==True and (plot_real_residuals==True):
            ax_diff.errorbar(real_tr, realr_me, yerr=[yerr_lor,yerr_hir], fmt='ko', markersize=5,label=real_label)


    #if xr is not None and yr is not None:
    #    if yr_err is not None:
    #        ax_main.errorbar(xr, yr, yerr=yr_err, color='black', fmt='o', label='DES')
    #    if yr_diff is not None:
    #        ax_diff.scatter(xr, yr_diff, color='black')
    #    if yr_diff is not None and yr_err is not None:
    #        ax_diff.errorbar(xr, yr_diff, yerr=yr_err, color='black', fmt='o')

    # ================================
    # Set plot options
    # ================================
    # set plot options [main axis]
    ax_main.set_xlabel(ylabel) 
    if type(ylimits)!=type(None):
        ax_main.set_ylim((ylimits[0],ylimits[1]))
  
    if res_axis==True:
        plt.setp(ax_main.get_xticklabels(), visible=False)
    else:
        plt.setp(ax_main.get_xticklabels(), visible=True, fontsize=fontsize)
    plt.setp(ax_main.get_yticklabels(), visible=True, fontsize=fontsize)
    if res_axis==True:
        plt.setp(ax_diff.get_yticklabels(), visible=False)
        plt.setp(ax_diff.get_xticklabels(), visible=True,fontsize=fontsize)
    plt.setp(ax_hist_x.get_xticklabels(), visible=False)
    plt.setp(ax_hist_x.get_yticklabels(), visible=False)
    plt.setp(ax_hist_y.get_yticklabels(), visible=False)
    plt.setp(ax_hist_y.get_xticklabels(), visible=False)
    #ax_main.set_ylabel(ylabel, fontsize=fontsize)
    #ax_main.autoscale(enable=True,axis='both')
    #ax_main.axis('equal')
    ax_main.set_xlim((xmin, xmax))   
    ax_main.set_ylim((ymin, ymax))
    ax_main.set_ylabel(ylabel, fontsize=fontsize)   
    ax_main.grid(True, alpha=0.2) 
    ax_hist_x.grid(True,alpha=0.2)
    ax_hist_y.grid(True,alpha=0.2)
    if type(real_label)!=type(None):
        ax_main.legend()
    #ax_main.set_aspect('equal')

    # set plot options [residual axis]
    if res_axis==True:
        ax_diff.set_ylim(resylim)
        ax_diff.set_xlabel(xlabel, fontsize=fontsize)
        ylabel_plain=ylabel.replace('$','')
        ax_diff.set_ylabel('$\\frac{\\delta '+ylabel_plain+'}{'+ylabel_plain+'}$',fontsize=fontsize)
        ax_diff.grid(True, alpha=0.2)
    else:
        ax_main.set_xlabel(xlabel, fontsize=fontsize)
    plt.savefig(plot_name)
    return 0


def get_mode(x,y):
    return x[y==max(y)]



def make_bin(col1_dat,col2_dat,bin_edges,return_vals=False,method="mean", remove_nan=True,mode_bin=[]):
    """
    col1= array like of the x values
    col2= array like of the y values
    Funtion that bins in x of a sample of (x,y) pairs.
    USe the return_vals if you want a list of y values in each bin

    """
    print (type(col2_dat))
    print (np.array(col2_dat).shape)
    #col1_dat=fits_catalog_col(cat,col1)
    #col2_dat=fits_catalog_col(cat,col2)
    #xy=zip(col1_dat,col2_dat)
    bin_vals = [[] for x in range(0,len(bin_edges)-1)]
    bin_valsx = [[] for x in range(0,len(bin_edges)-1)]
    bin_median=numpy.zeros(len(bin_edges)-1)
    bin_mean=numpy.zeros(len(bin_edges)-1)
    bin_std=numpy.zeros(len(bin_edges)-1)
    bin_center=numpy.zeros(len(bin_edges)-1)

    bin_percentile1sig=[[] for x in range(0,len(bin_edges)-1)]
    bin_percentile2sig=[[] for x in range(0,len(bin_edges)-1)]
    bin_percentile3sig=[[] for x in range(0,len(bin_edges)-1)]

    for i in range(0,len(col1_dat)):
        for j in range(0,len(bin_edges)-1):
            if bin_edges[j] <= col1_dat[i] <= bin_edges[j+1]:
                #print "data size"
                #try:
                #    print col1_dat.size
                #    print col2_dat.size
                #except:
                #    print len(col1_dat)
                #    print len(col2_dat)
                bin_vals[j].append(col2_dat[i])
                bin_valsx[j].append(col1_dat[i])            
    
    for i in range(0,len(bin_edges)-1):
        bin_center[i]= (bin_edges[i]+bin_edges[i+1])/2
    

    


    for i in range(0,len(bin_edges)-1):
        try:
            bin_mean[i]=numpy.mean(bin_vals[i])
        except:
            print(type(bin_vals[i]))
            print("Error!")

    
        bin_median[i]=numpy.median(bin_vals[i])
        bin_std[i]=numpy.std(bin_vals[i])
        if method=="percentile":
            if len(bin_vals[i])!=0:
                

                #median=np.percentile(measurements,50.0)
                lowlim=np.percentile(bin_vals[i],15.87)
                highlim=np.percentile(bin_vals[i],84.13)
                lowlim2sig=np.percentile(bin_vals[i],2.25)
                highlim2sig=np.percentile(bin_vals[i],97.6)
                lowlim3sig=np.percentile(bin_vals[i],0.15)
                highlim3sig=np.percentile(bin_vals[i],99.7)
            else:
                lowlim=0
                highlim=0
                lowlim2sig=0
                highlim2sig=0
                lowlim3sig=0
                highlim3sig=0
            bin_percentile1sig[i]=[lowlim,highlim]
            bin_percentile2sig[i]=[lowlim2sig,highlim2sig]
            bin_percentile3sig[i]=[lowlim3sig,highlim3sig]
            
        if bin_std[i]<0.001:
            print("zero bin std data:")
            print(bin_vals[i])
    if remove_nan==True:
        #ids2keep = np.where(~(np.isnan(bin_center) | np.isnan(bin_median) ))[0]
        ids2keep = (np.logical_not(np.isnan(bin_median)) )#| not(np.isnan(bin_median)) )

        #print (ids2keep)
        bin_center = bin_center[ids2keep]
        bin_std=bin_std[ids2keep]
        if method=='median':
            bin_median = bin_median[ids2keep]
         
        elif method=='percentile':
            bin_median = bin_median[ids2keep]
            bin_percentile1sig = np.array(bin_percentile1sig)[ids2keep]
            bin_percentile2sig = np.array(bin_percentile2sig)[ids2keep]
            bin_percentile3sig = np.array(bin_percentile3sig)[ids2keep]
        else:
            bin_mean = bin_mean[ids2keep]
    if return_vals==False:
        if method=='median':
            return bin_center,bin_median,bin_std
        elif method=='percentile':

            return bin_center,bin_median,bin_std,bin_percentile1sig,bin_percentile2sig,bin_percentile3sig       
        else:
            return bin_center,bin_mean,bin_std
    else:
        # if you want the a list of values inside each bin,
        # so you can perform your own bin metrics.
        return bin_center,bin_valsx,bin_vals




def mag_lim_plot(cat,**kwargs):
    
    MAG=verify_kwarg("MAG",'MAG_ISO',kwargs)
    loc=verify_kwarg("location",'(0,1)',kwargs)
    tex=verify_kwarg('tex',False,kwargs)
    

    try:
        MAG_AUX=MAG.split("_")
        MAG_ERR=MAG_AUX[0]+"ERR_"+MAG_AUX[1]
    except:
        MAG_ERR=''
    MAG_ERR=verify_kwarg("MAG",MAG_ERR,kwargs)
    
    mag_lim=verify_kwarg("mag_lim",0.1979,kwargs)
    bin=verify_kwarg("bin",False,kwargs)
    plot=verify_kwarg("plot_name",'maglim.png',kwargs)
    save_plot=verify_kwarg("save",True,kwargs)
    mag_err_aux=MAG.split("_")[-1]
    MAG_ERR="MAGERR"+"_"+mag_err_aux
    
    if ('bin' in kwargs.keys()) and ('bin_edges' in kwargs.keys()):  
        fig1,axes=fits_catalog_scatter(cat,col1=MAG,col2=MAG_ERR,plot_name=plot,save=False,bin=kwargs['bin'],bin_edges=kwargs['bin_edges'],tex=tex)
        x=kwargs['bin_edges']
    else:
        fig1,axes=fits_catalog_scatter(cat,col1=MAG,col2=MAG_ERR,plot_name=plot,save=False, tex=tex)
        x=fits_catalog_col(cat,MAG)
    
    if ('text' in kwargs.keys()):
        axes.text(loc[0], loc[1], kwargs['text'] , va='top', fontsize=9, transform=axes.transAxes)

    if ('plotmag' in kwargs.keys()):
        plotmag=kwargs['plotmag']
        newy=numpy.arange(0,mag_lim,0.001)
        newx=numpy.ones(len(newy))*plotmag
        axes.plot(newx,newy,"r-")
        xmax=plotmag
    else:
        xmax=max(x)
    
    xmin=min(x)
    newx=numpy.arange(int(round(xmin,0)),xmax,0.001)
    mag_lim_axis=numpy.ones(newx.shape)
    mag_lim_axis=mag_lim_axis*mag_lim
    axes.plot(newx,mag_lim_axis,"r-")

    if save_plot==True:
        plt.savefig(plot)
    return fig1, axes





def ascii_catalog_scatter(cat,col1,col2, **kwargs):
    x_label=verify_kwarg("xlabel","x",kwargs)
    y_label=verify_kwarg("ylabel","y",kwargs)
    col1_dat=select_column_ascii(cat,[int(col1)])
    col2_dat=select_column_ascii(cat,[int(col2)])
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(col1_dat,col2_dat)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if ('xlim' in kwargs.keys()):
        plt.xlim(options['xlim'][0],options['xlim'][1])
    if ('ylim' in kwargs.keys()):
        plt.ylim(options['ylim'][0],options['ylim'][1])
    time=datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d--%H%M%S-%f')
    plt.savefig(os.getcwd()+"/"+time+str(x_label)+"X"+str(y_label)+"txt.pdf")
    return fig

def list_radec(coord_file,**kwargs):
    comments=verify_kwarg("comments","#",kwargs)
    skiprows=verify_kwarg("skip_rows",0,kwargs)
    usecols=verify_kwarg("usecols",None,kwargs)
    unpack=verify_kwarg("unpack",False,kwargs)
    cat=open_ascii_cat(file_name=coord_file,comments=comments,skip_rows=skiprows,usecols=usecols,unpack=unpack)
    ra_pos=verify_kwarg("ra_col",1,kwargs)
    dec_pos=verify_kwarg("dec_col",2,kwargs)
    out_id=verify_kwarg("out_id",False,kwargs)
    id_pos=verify_kwarg("id_col",0,kwargs)
    ras=select_column_ascii(cat,ra_pos)
    ras=[float(ras[i]) for i in range(0,len(ras))]
    decs=select_column_ascii(cat,dec_pos)
    decs=[float(decs[i]) for i in range(0,len(ras))]
    if out_id==True:
        _ids=select_column_ascii(cat,id_pos)
        return _ids,ras,decs
    else:
        return ras,decs


def removelastcolumns(file_name,out_name,col_keep,**kwargs):
    comments=verify_kwarg("comments","#",kwargs)
    skiprows=verify_kwarg("skip_rows",0,kwargs)
    usecols=verify_kwarg("usecols",None,kwargs)
    cat=open_ascii_cat(file_name=file_name,comments=comments,skip_rows=skiprows,usecols=usecols,unpack=False, delimiter=',')
    new_cat=[]
    #print cat[0]
    cat=np.array(cat).astype(type('str'))
    #print cat[0]
    new_cat=cat[:,:col_keep].copy()
    #for i in range(0,len(cat)):
    #    new_cat.append(cat[i][:col_keep])
    save_ascii_cat(new_cat,out_name,header=None,identifier="#",delimiter=',') 
    
def open_ascii_cat(file_name,**kwargs):
    comments=verify_kwarg("comments","#",kwargs)
    skiprows=verify_kwarg("skip_rows",0,kwargs)
    sk_last=verify_kwarg("skip_footer",0,kwargs)
    #delimiter=verify_kwarg("delimiter"," , ",kwargs)#delimiter=" , "delimiter 
    usecols=verify_kwarg("usecols",None,kwargs) #Which columns to read, with 0 being the first. For example, usecols = (1,4,5) will extract the 2nd, 5th and 6th columns. The default, None, results in all columns being read.
    unpack=verify_kwarg("unpack",False,kwargs) 
    vartype=verify_kwarg("vartype",type('str'),kwargs) 
    if 'delimiter' in kwargs.keys():
        #print('I am skipping ',skiprows ,' rows')
        data=numpy.loadtxt(file_name,dtype=vartype,comments=comments,delimiter=kwargs['delimiter'],skiprows=skiprows,usecols=usecols,unpack=unpack)
    else:
        try:
            data=numpy.loadtxt(file_name,dtype=vartype,comments=comments,skiprows=skiprows,usecols=usecols,unpack=unpack)
        except:
            header=get_header_ascii(file_name,identifier=comments)
            m_v=0
            f_v=0
            cols=range(0,len(header))
            print("WARNING: table format is wrong.")
            data=numpy.genfromtxt(file_name,dtype=vartype,comments=comments,skip_header=skiprows, skip_footer=sk_last, missing_values=m_v, filling_values=f_v, usecols=cols, unpack=unpack)
    return data
def get_lenstool_rms(file_name,source=False,**kwargs):
    """
    Possible kwargs and default values
    ('plt_type',"line")
    ('rang',None,)
    #xlim=verify_kwarg('xlim',[-40,80],kwargs)
    ('normed',False)
    """
    skiprows=verify_kwarg("skip_rows",2,kwargs)
    skipfooter=verify_kwarg("skip_footer",4,kwargs)
    chi2=get_lenstool_chi2(file_name)
    data=open_ascii_cat(file_name,skip_rows=skiprows,unpack=True, skip_footer=skipfooter)
    #for i in range(1,6):
    #    print("line",len(data)-i,data[-i])
    #data=data[:-5,:]
    #print("line",len(data)-1,data[-1])
    #data[i, field]
    if source==True:
        ind=8
    else:
        ind=9
    tdata=numpy.transpose(data)
    tdata=select_subset_ascii(tdata,field=3,field_value=2,condition="<")
    data=numpy.transpose(tdata)
    rms=np.divide(np.sum(np.square(data[ind].astype("float"))),float(len(tdata[ind])))
    rms=sqrt(rms)
    alphabet=list(string.ascii_lowercase)
    _ids=data[1]
    num_img=len(_ids)
    for i in range(0,len(alphabet)):
        _ids=[s.replace(alphabet[i] , '') for s in _ids]
    num_fam=len(np.unique(_ids))
    return rms,num_fam,num_img,chi2
#def get_lenstool_nspec

def get_lenstool_chi2(file_name):
    fo = open(file_name, "r")
    lines = fo.readlines()
    chiline=lines[-2]
    chiline=chiline.replace("chitot","")
    chiline=chiline.replace("\n","")
    chi2tot=float(chiline)
    return chi2tot
    #chitot           92.67

def get_lenstool_dof(num_fam,num_img,mcmcfile="Bayes.dat", model="wcdm"):
     
    header=get_header_ascii(mcmcfile,identifier='#')
    free_params=len(header)-3 #remove Nsample ln(Lhood) Chi2 
    if model=="wcdm":
        free_params=free_params-1 #also account that the fact that omega_m+omega_l=1
    dof=(2*num_img) -free_params - (2*num_fam)
    return dof
def get_nz_optimized(mcmcfile="Bayes.dat"):
    header=get_header_ascii(mcmcfile,identifier='#')
    nz=0
    for i in range(0,len(header)):
        #header_aux=header[i]
        header_aux=header[i].replace("Redshift","")
        if header_aux !=header[i]:
            nz=nz+1
    return nz 

def save_ascii_cat(data,file_name,header=None,identifier="#",**kwargs):
    delimiter=verify_kwarg("delimiter"," ",kwargs)
    #method read/write append rwa was created to avoid \n in the end of the initial file.
    #cat,names=open_fits_catalog(fits_file)
    #ascii_file=string.rstrip(fits_file,".fits")
    #ascii_file=string.lstrip(ascii_file,".fit")
    #file_name=string.rstrip(file_name,".txt")
    
    method=verify_kwarg("method","w",kwargs)
    if method=='rwa':
        method='r'
    
    if method=='r':
       cat_ascii = open(file_name, 'r') 
       cat_data= cat_ascii.readlines()
       cat_ascii.close()
       cat_ascii=open(file_name, 'w')
       for i in range(0,len(cat_data)):
           cat_data[i]=cat_data[i].lstrip("\n")
           cat_ascii.write(cat_data[i].rstrip("\n")+'\n')
       #cat_ascii.write(cat_data[i]
    else:
        cat_ascii=open(file_name, method)   
    if header!=None:
        for i in range(0,len(header)):
            header[i]=header[i].rstrip("\n")
            header[i]=header[i].lstrip(identifier)
            cat_ascii.write(identifier+header[i]+"\n")
    

    for i in range(0,len(data)):
        for j in range(0,len(data[0])):
            line=str(data[i][j]).replace("\n",'')
            if j !=(len(data[0])-1):
                cat_ascii.write(line+delimiter)
            else:
                cat_ascii.write(line)
        if i!= (len(data)-1):
            cat_ascii.write("\n")
    cat_ascii.close()
    return 0






def string2array(a,delimiter=" "):
    a=a.lstrip("[")
    a=a.lstrip("(")
    a=a.rstrip("]")
    a=a.rstrip(")")
    a=a.split(delimiter)
    a=numpy.array(a)
    return a
def get_header_ascii(file_name,identifier='#'):

    f = open(file_name, 'r')
    lines=f.readlines()
    header=[]
    for i in range(0,len(lines)):
        lines_aux=lines[i].lstrip(identifier)
        if lines_aux!=lines[i]:
            lines_hdr=lines_aux.rstrip('\n')
            lines_hdr=lines_hdr.rstrip(' ')
            header.append(lines_hdr)
    return header
        
def add_column_ascii(data,column,header=None,col_name='', identifier='#',place=None):
    if place==None:
        out=numpy.append(data, column, axis=1)
    else:
        column=[[x] for x in column]
        out=numpy.insert(data,[place],column,axis=1)
    if (header!=None) and (col_name!=''):
        if place==None:
            header.append(identifier+col_name)
        else:
            header.insert(place,identifier+col_name)
        return out,header
    else:
        return out
   

def npy2array(file_name):

    data = np.load(file_name)
    # write to fits file
    #hdu = fits.PrimaryHDU(data)
    #hdu.writeto(myfile_fits)

    # read from fits file and show shape (test)
    #hdu = fits.open(myfile_fits)
    #scidata = hdu[0].data
    #print scidata.shape

    ### result: (20000, 3, 64, 64)

    return data









def select_column_ascii(data,column=0,header=None,col_name='',identifier='#'):
    if (header!=None) and (col_name!=''):
        for i in range(0,len(header)):
            test_line=header[i].lstrip(identifier)
            test_line=test_line.rstrip("\n")
            if test_line==col_name:
               column=i
    return data[:, column]

def remove_column_ascii(data,columns):
    IN_columns = [i for i in range(numpy.shape(data)[1]) if i not in columns]
    return data[:,IN_columns]

def arcslenstool2reg(lenstool_file,reg_file,delimiter=" "):
    ra,dec=open_ascii_cat(lenstool_file,usecols=(1,2),unpack=True,delimiter=delimiter)
    reg=open(reg_file, "w")
    reg.write("# Region file format: DS9 version 4.1 \n")
    reg.write("global color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n")
    reg.write("fk5 \n")
    for i in range(0,len(ra)-1):
        reg.write('circle('+str(ra[i])+','+str(dec[i])+',2.6") \n')
    reg.write('circle('+str(ra[len(ra)-1])+','+str(dec[len(dec)-1])+',2.6")')
    reg.close()

    return 0

def lenstool_images_distribution(lenstool_file, other_sources=[],shift_bins=0.08, plot_name="sources_dist.png",images=False,**kwargs):
    """
    Possible kwargs and default values
    ('plt_type',"line")
    ('rang',None,)
    #xlim=verify_kwarg('xlim',[-40,80],kwargs)
    ('normed',False)
    """
    kwargs['xlabel']='z'
    
    kwargs['plt_type']=verify_kwarg('plt_type',"bar",kwargs)
    kwargs['alpha']=verify_kwarg('alpha',0.5,kwargs)
    hatches=verify_kwarg('hatches',[None for i in range(len(lenstool_file)+1)],kwargs)
    #ran=verify_kwarg('rang',None,kwargs)
    #xlim=verify_kwarg('xlim',[-40,80],kwargs)
    #norm=verify_kwarg('normed',False,kwargs)
    #bins=verify_kwarg('bins','auto',kwargs)
    #marker=verify_kwarg('marker','b.',kwargs)
    #shift=verify_kwarg('shift',0.0,kwargs)
    #plt_type=verify_kwarg('plt_type',"line",kwargs)

    fig=verify_kwarg('figure',None,kwargs)
    ax1=verify_kwarg('axis',None,kwargs)
    fmt=[]
    markers=['.','^','s','o','D','v','<','*','s','p','*']
    colors=['b','g','r','c','m','k','y','c','b','g','r']
    for j in range(0,len(colors)):
        fmt.append(colors[j]+markers[j])

    if type(lenstool_file) is str:
        lenstool_file=[lenstool_file]
    interval=int(-float(len(lenstool_file))/2)
    if images==True: 
        fig,axes=start_subplots(2,**kwargs)
    else:
        fig,axes=start_plot(**kwargs) 
        axes=[axes]  
    for i in range(0,len(lenstool_file)):
        _id, ra, dec, a, b, theta, z, mag=open_ascii_cat(lenstool_file[i],unpack=True)
        z=numpy.array(z).astype('float')
        z_sources=numpy.unique(z)
        print("This are the lenstool sources")
        print(z_sources) 
        #if i==6:
        #    print "this is z"
        #    print z
        #    print lenstool_file[i]
        #==== kwargs for the cat_hist_plot function
        kwargs['save']=False
        kwargs['bins']=verify_kwarg('bins',numpy.arange(0,8,0.5),kwargs)
        kwargs['figure']=fig
        kwargs['axis']=axes[0]
        kwargs['shift']=(interval+i)*shift_bins
        kwargs['marker']=fmt[i]
        fig,axes[0]=cat_hist_plot(z_sources,hatch=hatches[i],**kwargs)

        if images==True:
            kwargs['axis']=axes[1]
            fig,axes[1]=cat_hist_plot(z,**kwargs)
            fig.subplots_adjust(hspace=0)
            axes[1].set_ylabel('Images')
        axes[0].set_ylabel('Sources')
        axes[0].set_xlim(0,7)

    if len(other_sources)>0:
        print("this are the other sources")
        print(other_sources)
        kwargs['alpha']= kwargs['alpha']*0.8
        kwargs['save']=False
        kwargs['bins']=verify_kwarg('bins',numpy.arange(0,8,0.5),kwargs)
        kwargs['figure']=fig
        kwargs['axis']=axes[0]
        kwargs['shift']=(interval+(i+1))*shift_bins
        kwargs['marker']=fmt[i+1]
        fig,axes[0]=cat_hist_plot(other_sources,hatch=hatches[-1],**kwargs)


    plt.savefig(plot_name)
    
    return 0


def reg2arcslenstool(reg_file,lenstool_file,out_file):
    out=open(out_file, "w")
    out_list=[]
    reg=open(reg_file, "r")
    _id, ra, dec, a, b, theta, z, mag=open_ascii_cat(lenstool_file,unpack=True)
    lines=reg.readlines()
    out_list.append("#REFERENCE 0 \n")
    raerr=0.0002
    decerr=0.000002
    for i in range(3,len(lines)):
        line=lines[i].lstrip("circle(")
        line=line.rstrip(")")
        rareg,decreg,radius=line.split(",")
        print(rareg)
        print(decreg)
        print(dec)
        for j  in range(0,len(ra)):
            if (float(rareg)<(float(ra[j])+raerr)) and (float(rareg)>(float(ra[j])-raerr)) and (float(decreg)<(float(dec[j])+decerr)) and (float(decreg)>(float(dec[j])-decerr)):
                out_list.append(str(_id[j])+" "+str(ra[j])+" "+str(dec[j])+" "+str(a[j])+" "+str(b[j])+" "+str(theta[j])+" "+str(z[j])+" "+str(mag[j])+"\n")            
    out_list[-1]=out_list[-1].rstrip("\n")
    for i in range(0,len(out_list)):
        out.write(out_list[i])
    out.close()

    return 0 

def select_subset_ascii(data,field,field_value,condition):

 #FIXME I think only 'in' is not working 
    if (condition=="=" or condition=="=="):
        line_ind = [i for i in range(numpy.shape(data)[0]) if float(data[i, field])==float(field_value)]
        return  data[line_ind,:]
    elif condition==">":
        line_ind = [i for i in range(numpy.shape(data)[0]) if float(data[i, field])>float(field_value)]
        return  data[line_ind,:]
    elif condition=="<":
        #col=[data[i, field] for i in range(numpy.shape(data)[0])]
        #print(col)
        line_ind = [i for i in range(numpy.shape(data)[0]) if float(data[i, field])<float(field_value)]
        return  data[line_ind,:]
    elif condition=="interval":
        min_val=min(field_value)
        max_val=max(field_value)
        line_ind = [i for i in range(numpy.shape(data)[0]) if float(data[i, field])>float(min_val)]
        extractedData=data[line_ind,:]
        line_ind2 = [i for i in range(numpy.shape(extractedData)[0]) if float(data[i, field])<float(max_val)]
        return extractedData[line_ind2,:]
    elif condition=="in":
        line_ind = [i for i in range(numpy.shape(data)[0]) if float(data[i, field]) in field_value]
        return  data[line_ind,:]
    else:
        print("no valid condition entry. Condition must be a string: '>','<','==' or 'interval'.")
        return 0

def select_subset_arr(x,val,condition, return_mask=False, **kwargs):
    x=numpy.array(x)
    if 'test_column' in kwargs.keys():
        x_test=numpy.array(zip(*x)[kwargs['test_column']])
    else:
        x_test=numpy.array(x).copy()
    if 'into' in kwargs.keys():
        x_out=numpy.array(kwargs['into']).copy()
    else:
        x_out=numpy.array(x).copy()
        
    if (condition=="=" or condition=="=="):
        return  x_out[x_test == val]
    elif condition==">":
        return  x_out[x_test > val]
    elif condition=="<":
        return  x_out[x_test < val]
    elif condition=="interval":
        min_val=min(val)
        max_val=max(val)
        x_aux=(x_out[x_test < max_val]).copy()
        x_aux_test=(x_test[x_test < max_val]).copy()
        #if 'test_column' in kwargs:
        #    x_aux_test=numpy.array(zip(*x_aux)[kwargs['test_column']])
        #else:
        #    x_aux_test=x
        return x_aux[x_aux_test > min_val]
        
           
    elif condition=="in": #fix me
        return numpy.array([x_out[i] for i in range(0,len(val)) if val[i] in x_test]) #val e o subsample
    elif condition=="abs>":
        return x_out[abs(x_test) > abs(val)]
    elif condition=="abs<":
        return x_out[abs(x_test) < abs(val)]
    else:
        print("no valid condition entry. Condition must be a string: '>','<','==' or 'interval'.")
        return 0

def select_id_fits_catalog(cat,_id):
    subset = cat[cat.field('ID') == str(_id)]
    return subset

def match_cat(matchfield,cat_ref,cat_test):
    col_ref=fits_catalog_col(cat_ref,col_name=matchfield)
    col_test=fits_catalog_col(cat_test,col_name=matchfield)
    cat_matched=cat_ref

    for i in range(0,len(col_ref)):
        #print cat_test[matchfield]
        ind=np.where(col_test==col_ref[i])#cat_test_aux=cat_test[matchfield == col_ref[i]]
        #print ind[0]
        if len(ind[0])==0:
            cat_matched=select_subset_fits_catalog(cat_matched,field_name=matchfield,field_value=col_ref[i],condition="!=")
            #print "eliminated   "+str(col_ref[i])

    return cat_matched

def select_subset_fits_catalog(cat,field_name,field_value,condition,**kwargs):

    if 'into' in kwargs.keys():
        cat_out=kwargs['into']#.copy()
    else:
        cat_out=cat        

    if 'col_index' in kwargs:
        test_field=numpy.array(zip(*cat.field(str(field_name)))[kwargs['col_index']])
    else:
        test_field=cat.field(str(field_name))

    if (condition=="id"):
        return  cat_out[field_value]
    if (condition=="=" or condition=="=="):
        return  cat_out[test_field == field_value]
    if (condition=="!="):
        return  cat_out[test_field != field_value]
    elif condition==">":
        return  cat_out[test_field > field_value]
    elif condition=="<":
        return  cat_out[test_field < field_value]
    elif condition=="interval":
        min_val=min(field_value)
        max_val=max(field_value)
        cat_aux=cat_out[test_field < max_val]
        cat_aux_test=cat[test_field < max_val]
        if 'col_index' in kwargs:
            test_field_aux=numpy.array(zip(*cat_aux_test.field(str(field_name)))[kwargs['col_index']])
        else:
            test_field_aux=cat_aux_test.field(str(field_name))

        return cat_aux[test_field_aux > min_val]
    elif condition=="in": #fix me
        return  cat_out[numpy.where(numpy.in1d(field_value, cat.field(str(field_name))))[0]]
    elif condition=="abs>":
        return cat_out[abs(test_field) > abs(field_value)]
    elif condition=="abs<":
        return cat_out[abs(test_field) < abs(field_value)]
    else:
        print("no valid condition entry. Condition must be a string: '>','<','==' or 'interval'.")
        return 0

