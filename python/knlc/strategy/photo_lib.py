
from __future__ import division
from __future__ import print_function
import matplotlib
matplotlib.use('agg')
import numpy as np
import fits_cat as fc
from fits_cat import make_bin
from matplotlib import rc
from misc import quadrance,define_workdir,verify_kwarg,linear,scale
#from cbomcode.image.image import *
#import pyfits as fits
from math import log10,exp,sqrt, atan, log
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit
import scipy
import numpy as np
from misc import gauss_norm, gauss, verify_kwarg
from sklearn.metrics import mean_squared_error
import matplotlib.cm as cm
import itertools
import random
import cosmology as co
import constants as const
from scipy.integrate import simps
import math
from astropy.io import fits
import seaborn as sns
from scipy.interpolate import interp1d
import h5py
import pandas as pd



def spec_phot_weights(cat_spec,cat_phot,colors=False, nei_num=20, id_col='QUICK_OBJECT_ID',bands=['MAG_G','MAG_R','MAG_I','MAG_Z']):

    n_phot=len(cat_phot)
    n_spec=len(cat_spec)    
    #dist = numpy.linalg.norm(mags_spec-mags_spec)

    neis=[]

    weights_spec=[]

    mags_spec=np.transpose([cat_spec.field(band) for band in bands])
    #for j in range(0,len(bands)):
    

    for i in range(0,n_spec):
        cat_cut=remove_dist(mags_spec[i], cat_spec,cut_lim=0.5, bands=bands)
        if ((i+1) % 50) ==0:
            print('Object number '+str(i+1)+' of '+str(n_spec))
            print ('Number of objects to be considered before cut '+str(len(cat_spec)))
            print ('Number of objects to be considered after cut '+str(len(cat_cut)))
        d2=0
        for j in range(0,len(bands)):
            xa=mags_spec[i][j]-cat_cut.field(bands[j])
            if j==0:
                d2=xa**2
            else:
                d2=d2+xa**2
        #mags_spec_cut=np.transpose([cat_cut.field(band) for band in bands])
        #d2=np.linalg.norm((mags_spec[i]-mags_spec_cut),axis=1)
        # with \approx 550 objects cur it magnitudes difference
        # higher than 0.5 the linalg performed 376 and 372 s while the one-for-loop 352 and 362
        #d_sum=sum(d2-d3)
        #if d_sum>0
        #    print ('The difference beetween two methods were '+str(i+1)+' of '+str(n_spec)


        #people = numpy.array(people)
        #ages = numpy.array(ages)
        inds = d2.argsort()
        cat_select = cat_cut[inds]
        cat_select=cat_select[:nei_num]
        maglims=[]
        for j in range(0,len(bands)):
            maglims.append(cat_select.field(bands[j])[-1])
        radius=np.abs(mags_spec[i]-maglims)
        
        cat_phot_cut=cat_phot.copy()
        for j in range(0,len(bands)):
            cat_phot_cut=fc.select_subset_fits_catalog(cat_phot_cut,field_name=bands[j],field_value=mags_spec[i][j]+radius[j],condition="<")
            cat_phot_cut=fc.select_subset_fits_catalog(cat_phot_cut,field_name=bands[j],field_value=mags_spec[i][j]-radius[j],condition=">")
            #FIXME: calculate the radius in the photometric catalog properly in more than one band.
        weights_spec.append(float(len(cat_phot_cut))/float(len(cat_select)))     
        
    weights_spec=(1/(n_phot))*np.array(weights_spec)
    return weights_spec      
         


def extended_class(spread_model,spread_model_error):
    extend_val=((spread_model + 3*spread_model_error ) > 0.005)*1
    extend_val+=((spread_model + spread_model_error ) > 0.003)*1
    extend_val+=((spread_model - spread_model_error ) > 0.003)*1
    extend_val[spread_model==-1]=-9
    extend_val[(spread_model==0) & (spread_model_error==0)]=-9
    extend_val[(spread_model==1) & (spread_model_error==1)]=-9
    return(extend_val)




def remove_dist(mags,cat,cut_lim,bands):

    for i in range(0,len(bands)): 
        cat_cut=fc.select_subset_fits_catalog(cat,field_name=bands[i],field_value=mags[i]+cut_lim,condition="<")
        cat_cut=fc.select_subset_fits_catalog(cat_cut,field_name=bands[i],field_value=mags[i]-cut_lim,condition=">")
    return cat_cut



def knbullaload(file_name,sigma=0.05,out_dir="",telescope_area=10.014, \
object_size=3,exposure=90,filter_ids=[293,294,295,296,297,298],filter_file= \
'/home/cleciobom/lib/eazy-photoz/filters/FILTER.RES.latest', \
filter_names=['u','g','r','i','z','Y'],filter_source='eazy', filter_unit='A', thoughput=0.525, \
atmtransmission='/home/cleciobom/lib/eazy-photoz-master2/cptrans_zm_23_10.txt', \
mirror_reflectivity='/home/cleciobom/lib/eazy-photoz-master2/alreflectivity.txt', \
nmirrors=2, sky_flux='skybg_50_10.datALL.txt' , zs=np.arange(0.01,0.5,0.02),plot=False, tex = False,kn_id=0, maglim=22.0, filtermaglim=2):
    """
    nphX_mejY_phiZ_TW.txt
    X: number of Monte Carlo packets
    Y: total ejecta mass
    Z: half-opening angle of the lanthanide-rich component
    W: Temperature at 1 day (T0)
    
    Fluxes are given in erg s-1 cm-2 A-1 at the distance of 10 pc
    """
    print ("working on... ",file_name)
  
    kn_file = open(file_name, 'r') 
    kn_data= kn_file.readlines()
    kn_file.close()
    file_name_var=(file_name.split("/")[-1]).split("_")
    #print(file_name_var)
    mc_p= float(file_name_var[0].lstrip("nph"))
    ejm=float(file_name_var[1].lstrip("mej"))
    phi=float(file_name_var[2].lstrip("phi"))
    temp=file_name_var[3].lstrip("T")
    temp=float(temp.rstrip(".txt"))
  
    #v_angle=float(kn_data[0])
    cosva=[]
    stepva=1/(float(kn_data[0])-1.0)
    print(kn_data[0])
    for i in range(0,int(kn_data[0])):
        cosva.append(i*stepva)
    nwave=float(kn_data[1])
    timebin=float(kn_data[2].split()[0])
    t_i=float(kn_data[2].split()[1])
    t_f=float(kn_data[2].split()[2])
    dt=(t_f-t_i)/timebin
    times=[]
    i=0
    print(timebin)
    for i in range(0,int(timebin)):
        times.append(t_i+((i+0.5)*dt))

    cosva_wav=[]
    #kn_id=0
    for m in range(0,len(zs)):
        for i in range(0,len(cosva)):
            kn_id=kn_id+1
            wave=[]
            fluxes=[]
            wave_txt=[]
            if i==0:
                ini=3
            else:
                ini=ini+j+1
            for j in range(0,int(nwave)):
                wave_txt.append(kn_data[ini+j].split()[0])
                fluxes.append(kn_data[ini+j].split()[1:])
            fluxes_tp=np.transpose(fluxes)
            mags_time=[]
            mags_time_err=[]
            mags_time_noise=[]
            for k in range(0,int(timebin)):
                
                wave=np.array(wave_txt).astype(np.double)
                flux_tp=np.array(fluxes_tp[k]).astype(np.double)
                #print('flux_tp') 
                #print(flux_tp)
                z_wave,z_flux,z_eff=refshift_tempflux2(wavelen=wave.copy(),flux=flux_tp.copy(),zin=2.37 * 10E-9,zout=zs[m])#refshift_templum(wavelen=wave,flux=flux_tp,zout=zs[m])
                #print('flux_z') 
                #print(z_flux)
                wavelen_err,flux_total_err,flux_total,mags,mags_err,mags_noise=add_spectro_err(flux=z_flux.copy(),wavelen=z_wave.copy(),telescope_area=telescope_area,object_size=object_size,exposure=exposure,filter_ids=filter_ids,filter_file=filter_file,filter_source=filter_source, filter_unit=filter_unit, thoughput=thoughput,atmtransmission=atmtransmission,mirror_reflectivity=mirror_reflectivity,nmirrors=nmirrors, sky_flux=sky_flux,noise=sigma)
                kilonova_id_str=str(kn_id)
                kilonova_id_str.zfill(5)
                KN_name='KN'+kilonova_id_str+'_'+str(mc_p)+'_'+str(ejm)+'_'+str(phi)+'_'+str(temp)+'_'+str(cosva[i])+'_'+str(times[k])+'_'+str(zs[m])
                if (mags[filtermaglim]<= maglim):
                    hf = h5py.File(out_dir+KN_name+".h5", 'w')
                    hf.create_dataset('flux_erg_s_a_pure', data=flux_tp)
                    hf.create_dataset('wavelen_pure', data=wave)
                    #print ("working on... ",file_name)
                    #print ('Kilonova'+kilonova_id_str)                
                    #print (len(wave))
                    #print (nwave)
                    #print (wave)
                    #print (len(flux_tp))
                    #print (flux_tp)

                    hf.create_dataset('flux_erg_s_a_z', data=z_flux)
                    hf.create_dataset('wavelen_z', data=z_wave)

                    hf.create_dataset('flux_erg_s_a_noise', data=flux_total_err)
                    hf.create_dataset('flux_erg_s_a_noiseless', data=flux_total)
                    hf.create_dataset('wavelen_noise', data=wavelen_err)
                    for l in range(0,len(filter_names)):
                        hf.create_dataset('mag_'+filter_names[l], data=mags[l], dtype=type(2.0))
                        hf.create_dataset('magerr_'+filter_names[l], data=mags_err[l], dtype=type(2.0))
                      

                    hf.create_dataset('mc_pac', data=mc_p)
                    hf.create_dataset('ejm', data=ejm)
                    hf.create_dataset('phi', data=phi)
                    hf.create_dataset('temperature', data=temp) 
                    hf.create_dataset('cosva', data=cosva[i])
                    hf.create_dataset('time', data=times[k])
                    hf.create_dataset('z', data=z_eff)

                    hf.create_dataset('noise_level', data=sigma)
                    hf.create_dataset('original_file_name', data=file_name)
                    #hf.create_dataset('exposure_s', data=exposure)
                    #hf.create_dataset('telescope_mirror_area_sqm', data=telescope_area)
                    #hf.create_dataset('telescope_nmirror', data=nmirrors)
                    #hf.create_dataset('detector_thoughput', data=thoughput)        

                    hf.close()
                mags_time.append(mags)
                mags_time_err.append(mags_err)
                mags_time_noise.append(mags_noise)
                if plot==True and (mags[filtermaglim]<= maglim) :
                    #model= h5py.File(out_dir+str(mc_p)+str(ejm)+str(phi)+str(temp)+str(cosva[i])+str(times[k])+".h5",'r')
                    df_model_pure = pd.DataFrame({'LAMBDA':wave, 'Flux':flux_tp})
                    df_model_z = pd.DataFrame({'LAMBDA':z_wave, 'Flux_z':z_flux})
                    df_model_sky = pd.DataFrame({'LAMBDA':wavelen_err, 'Flux_sky':flux_total})
                    df_model_noise = pd.DataFrame({'LAMBDA':wavelen_err, 'Flux_noise':flux_total_err})  
                    

                    norm = df_model_pure['Flux'].median()
                    df_model_pure['normFlux'] = df_model_pure['Flux']/norm
 
                    norm = df_model_z['Flux_z'].median()
                    df_model_z['normFlux_z'] = df_model_z['Flux_z']/norm

                    norm = df_model_sky['Flux_sky'].median()
                    df_model_sky['normFlux_sky'] = df_model_sky['Flux_sky']/norm

                    norm = df_model_noise['Flux_noise'].median()
                    df_model_noise['normFlux_noise'] = df_model_noise['Flux_noise']/norm  



                    if (tex == True):
                        rc('text', usetex=True)
                    else:
                        rc('text', usetex=False)
                    ax = df_model_pure.plot('LAMBDA','normFlux',grid=True)
                    df_model_z.plot('LAMBDA', 'normFlux_z', ax=ax)
                    df_model_sky.plot('LAMBDA', 'normFlux_sky', ax=ax)
                    df_model_noise.plot('LAMBDA', 'normFlux_noise', ax=ax)
                    title = "EJM="+str(ejm)+" Cos VA ="+str(cosva[i])+" t="+str(times[k])+" z="+str(zs[m])+" mag_"+filter_names[0]+"="+str(mags[0])
                    #plt.figtext(pos[0], pos[1], title)
                    #ax.set_ylim(0,10)
                    plt.title(title)
                    plt.savefig(out_dir+KN_name+"_spec.png")
            if plot==True:
                mags_time_err_tp=np.transpose(mags_time_err)
                mags_time_noise_tp=np.transpose(mags_time_noise)
                mags_time_tp=np.transpose(mags_time)
                for o in range(0,len(filter_names)):
                    fig_mag = plt.figure()
                    ax1 = fig_mag.add_subplot(111)
                    c = itertools.cycle(["r", "b", "g","u","k","c","m"])
                    #iter(cm.rainbow(np.linspace(0, 1, R)))
                    #,colors=next(c)
                    masknone=np.array(mags_time_err_tp[o]>-90)
                    print (masknone)
                    print (times)
                    print (mags_time_err_tp[o])
                    print (mags_time_noise_tp[o])
                    ax1.errorbar(np.array(times).astype('float')[masknone],mags_time_err_tp[o][masknone],yerr=mags_time_noise_tp[o][masknone], ls='--',markersize=2, alpha=0.5)
                    ax1.set_ylim(30, 15)  # decreasing time

                    plt.savefig(out_dir+KN_name+"_mag_"+filter_names[o]+".png")              

             #np.save(out_dir+str(mc_p)+str(ejm)+str(phi)+str(temp)+str(cosva[i])+str(times[k])+".h5",[wave,fluxes_tp[k]])
             #np.save(out_param_dir+str(mc_p)+str(ejm)+str(phi)+str(temp)+str(cosva[i])+str(times[k])+".npz",[mc_p,ejm,phi,temp,cosva[i],times[k]])
     
     #cosva_wav.append(wave)
     #cosva_flux.append(fluxes)  


    return kn_id







    #last=desc_data[-1]
            #print desc_data[:5]
            #print description_file
            #print desc_data[5:]
            #print desc_data[-1]
            #print desc_data[-1].split()
#            last_filter=len(desc_data)#int(desc_data[-1].split()[0])
#            cat_ascii=open(description_file, "w")




def magerr2sn(magerr):
    return 1.0875/magerr #n_counts
def sn2magerr(sn):
    return 1.0875/sn

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def absmag2loglum(absmag):
    return 0.4*(-48.6-absmag)

def loglum2absmag(loglum):
    return -2.5*loglum-48.6

def l_lim(absmag,mag,mag_lim):
    return (absmag2loglum(absmag)+0.4*(mag-mag_lim))

def get_l_min(absmags,mags,mag_lim,zs,z, kcorr=False):
    mag_low=np.percentile(mags,0.8)
    mags=np.array(mags)
    zs=np.array(zs)
    mask=mags>mag_low
    #print (mask)
    mags=mags[mask]
    absmags=absmags[mask]
    zs=zs[mask]
    l_maglim=l_lim(absmags,mags,mag_lim)
    
    mag_atmaglim=absmag2mag(loglum2absmag(l_maglim),zs, kcorr=kcorr)
    l_cut=np.percentile(l_maglim,0.95)
    if (loglum2absmag(l_cut) < -10) or (loglum2absmag(l_cut)>0.0) :
        #print ('absmag_cut')
        #print(loglum2absmag(l_cut))
        #print (l_cut)
        #print (mag_lim) 
        #print (np.transpose([absmags,loglum2absmag(l_maglim),mags]))
        abs_low_hist, abs_low_bin_edges = np.histogram(loglum2absmag(l_maglim),density=False)
        abs_low_bin_center=abs_low_bin_edges[:-1]+(float(abs_low_bin_edges[1]-abs_low_bin_edges[0])/2.0)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.bar(abs_low_bin_center,abs_low_hist,width=abs_low_bin_edges[1]-abs_low_bin_edges[0],alpha=0.8,hatch='\\', fill=False, ec='indianred')
        plt.savefig("absmag_hist_z_"+str(mag_lim)+"_"+str(loglum2absmag(l_cut))+"_"+str(z)+".png")
        plt.close('all')
        mag_l_maglim, mag_l_maglim_edges = np.histogram(mag_atmaglim,density=False)
        
        mag_l_maglim_center=mag_l_maglim_edges[:-1]+(float(mag_l_maglim_edges[1]-mag_l_maglim_edges[0])/2.0)
        fig = plt.figure()
        ax2 = fig.add_subplot(111)
        ax2.bar(mag_l_maglim_center,mag_l_maglim,width=mag_l_maglim_edges[1]-mag_l_maglim_edges[0],alpha=0.8,hatch='\\', fill=False, ec='indianred')
        plt.savefig("mag_at_lim_hist_z_"+str(mag_lim)+"_"+str(loglum2absmag(l_cut))+"_"+str(z)+".png")
    return loglum2absmag(l_cut)

def mag2absmag(mag,z,kcorr=False):
    if kcorr==False:
        return (mag -  lumn_dist(z))#dist_mod(z)) #
    else:
        return (mag -  dist_mod(z))#dist_mod(z)) #
def absmag2mag(absmag,z, kcorr=False):
    if kcorr==False:
        return absmag+ lumn_dist(z)# dist_mod(z) #
    else:
        return absmag+ dist_mod(z)# dist_mod(z) #



def dist_mod(z_vals,d=''):
    z = np.load(d+'dist_mod_kcorr.npz')['z']
    dm = np.load(d+'dist_mod_kcorr.npz')['dm']
    f_dm = interp1d(z, dm, kind='cubic', bounds_error=False, fill_value=0)
    return f_dm(z_vals)

def lumn_dist(z):
    H0 = 67.74
    h = H0/100.
    WV = 0.6911
    WR = 4.165E-5/(h*h)
    WM = 1. - WV - WR
    c = 299792.458
    n = 10000
    az = 1.0/(1.+z)
    DCMR = 0.
    DTT = 0.
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = np.sqrt((WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)
    DCMR = (1.-az)*DCMR/n
    DA = az*DCMR
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    # convert to distance modulus and return
    DL_mag = 5.*np.log10(DL_Mpc*1.e5) - 2.5*np.log10(1.+z)
    return DL_mag


def maglim_zbins(mags,zs,z_edges=np.arange(0.0,3.0,0.01),plot=False,mags_edges=np.arange(17,21.2,0.1),out_name=".png"):
    maglims=[]
    mags_z_median=[]
    zbin_center,zs_bined,mags_bined=fc.make_bin(zs,mags,bin_edges=z_edges,return_vals=True,method="mean", remove_nan=True)
    for i in range(0,len(zbin_center)):
        mags_z_median.append(np.median(mags_bined[i]))
        mag_hist, bin_edges_mag = np.histogram(mags_bined[i],bins=mags_edges,density=False)
        magbin_center=bin_edges_mag[:-1]+(float(mags_edges[1]-mags_edges[0])/2.0)
        #magbin_center,magbin_median,magbin_std,magbin_percentile1sig,magbin_percentile2sig,magbin_percentile3sig=fc.make_bin(mags_bined[i],mags_bined[i],bin_edges=mags_edges,return_vals=False,method='percentile', remove_nan=True)
        if len(mag_hist)>0:
            maglim=magbin_center[mag_hist==max(mag_hist)]
            print(maglim)
            if len(maglim)>1:
                maglim=np.array([maglim[0]])
            if maglim<18.0:
                print("Limiting magnitude="+str(maglim)+"in z="+str(zbin_center[i]))
                print(np.transpose([magbin_center,mag_hist]))
        else:
            maglim=[-99.0]
        maglims.append(maglim[0])
        if plot==True:
           fig = plt.figure()
           ax1 = fig.add_subplot(111)
           ax1.bar(magbin_center,mag_hist,width=mags_edges[1]-mags_edges[0],alpha=0.8,hatch='\\', fill=False, color='indianred')
           plt.savefig("mag_hist_z_"+str(zbin_center[i])+out_name)
           plt.close('all')
    if plot==True:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.plot(zbin_center,maglims,color='indianred')
        #ax2.plot(zbin_center,mags_z_median,alpha=0.8, color='navy')
        plt.savefig("maglim_z_"+out_name)

    return zbin_center,maglims,mags_bined,zs_bined    
       

def get_filter_centers_eazy(filter_file,filter_ids):
    centers=[]
    for i in range(0,len(filter_ids)):
        centers.append(get_mean_filter_eazy(filter_file,filter_ids[i],clear=0.01))
    return np.array(centers)


def open_template(temp_file,temp_type='SED',work_dir=''):

    if temp_type=='SED':
        wavelen,SED=fc.open_ascii_cat(temp_file, unpack=True, skip_rows=5, vartype=type(1.0))
        #fluxlambda=np.divide(np.power(10,SED.copy()),wavelen)
        fluxlambda=SED.copy()#np.divide(SED.copy(),wavelen)
        fluxlambda=fluxlambda*10E35
        nu,fluxnu=fluxlambda2fluxnu(wavelen.copy(),fluxlambda.copy())
        mag=fluxnu2mag(fluxnu.copy())

        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)

        
        ax0.plot(wavelen,fluxlambda,c='b')
        plt.savefig(work_dir+"flux_comp.png") 
    if temp_type=='sdss':
        f = fits.open(temp_file)
        tbdata = f[1].data
        fluxlambda=tbdata['model']*10E-17
        wavelen=np.power(10,tbdata['loglam'])
        nu,fluxnu=fluxlambda2fluxnu(wavelen.copy(),fluxlambda.copy())
        mag=fluxnu2mag(fluxnu.copy())

        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        
        ax0.plot(wavelen,fluxlambda,c='b')
        plt.savefig(work_dir+"flux_comp.png") 

    else:
        wavelen,mag=fc.open_ascii_cat(temp_file, unpack=True, skip_rows=5, vartype=type(1.0))
        fluxnu=mag2fluxnu(np.array(mag).astype('float'))
        fluxlambda=fluxnu2fluxlambda(wavelen,fluxnu)
        nu,fluxnu_rev=fluxlambda2fluxnu(wavelen.copy(),fluxlambda.copy())
        mag_rev=fluxnu2mag(fluxnu_rev.copy())
        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        ax0.plot(wavelen,mag_rev+0.5,c='g')
        ax0.plot(wavelen,mag,c='b')
        ax0.plot(wavelen,mag_rev-0.5,c='g')
        plt.gca().invert_yaxis()
        plt.savefig(work_dir+"mag_comp.png") 
    return mag.copy(),fluxlambda.copy(),fluxnu.copy(),wavelen.copy()


def find_nearestxy(x,y,target):
    #print "find near args"
    #a=np.array([3,4,5])
    #c=a-2.1
    #print c
    #print x
    #print target
    B=np.array(x)-target
    #print B
    idx = (np.abs(np.array(x).astype('float')-target)).argmin()
    #print "results nearest"
    #print (np.abs(np.array(x).astype('float')-target)).min()
    #print x[idx]
    #print y[idx]
    #print "#### the interval was"
    #if idx<(len(x)-1):
    #    interval= abs(x[idx]-x[idx+1])
    #else:
    #    interval= abs(x[idx]-x[idx-1])

    if target > x[idx]:
       if idx<(len(x)-1):
           interval= abs(x[idx]-x[idx+1])
       else:
           interval= abs(x[idx]-x[idx-1])

    else:
        if idx>0:
            interval= abs(x[idx]-x[idx-1])
        else:
            interval= abs(x[idx]-x[idx+1])

    #if idx>0:
    #    interval= abs(x[idx]-x[idx-1])
    #else:
    #    interval= abs(x[idx]-x[idx+1])

    #print interval
    #print abs(x[idx]-target)
    if interval > abs(x[idx]-target):
        flag=0
    else: 
        flag=1
        print("interval")
        print(idx)
        if idx==0:
           flag=2
        print(interval)
        print(x[idx])
        print(min(x))
        try:
            print(x[idx+1])
        except:
            pass
        print(x[idx-1])
        print(target)
        print(abs(x[idx]-target))
        print(y[idx])
    return y[idx],x[idx],flag

#ax1.set_ylim(-3E-3,+3E-3)

def fluxjan2mag(flux):
    """
    input in jankys
    """
    return -2.5*np.log10(np.array(flux).astype('float'))+8.90

def mag2fluxjan(mag):

    expoents=(-0.4)*(np.array(mag).astype('float')-8.9)
    return np.power(10,expoents)

def fluxnu2mag(flux):
    """
    input in erg s-1 cm-2 hz-1
    """
    return (-2.5*np.log10(np.array(flux).astype('float')))-48.60

def mag2fluxnu(mag):
    """
    input in mag AB, output in erg s-1 cm-2 hz-1
    """
    expoents=(-0.4)*(np.array(mag).astype('float')+48.60)
    return np.power(10,expoents)

def lambda2nu(wavelen):
    c=2997924580000000000.0
    nu=c/wavelen
    return nu

def fluxlambda2mag(wavelen,flux):
    return  (-2.5*np.log10( np.array(flux).astype('float') ))- (5.0 * np.log10( np.array(wavelen).astype('float') )) - 2.402


def fluxnu2magfilter(flux):
    #mags=
    return 0

def refshift_templum(wavelen,flux,zout):
    """
    WARNING: FLUX DENSITY MUST BE IN units of wavelenght
    """
    #for i in range(0,len(wavelen)):
    #    wavelen[i]+=zout*wavelen[i]
    print ("redshifting...")
    print (zout)
    print(wavelen)
    wavelen=(1+zout)*(np.array(wavelen).astype(type(1.0)))
    flux=np.array(flux).astype('double')
    DL=co.lum_dist(0, z_final=zout, omega0_m=0.3, omega0_k=0, omega0_l=0.7,w=-1,unit='km', model='lambdacdm')
    DL=DL*10E5
    z_plus_one=1.0+zout
    flux=(1.0/z_plus_one)*flux
    #print("the k correction was    ")
    flux=flux/((4*const.pi)*(DL**2))
    return wavelen,flux

def redshift(wavelen,z):
    wavelen=np.array(wavelen).astype(type(1.0))
    wavelen=(1+z)*wavelen
    return wavelen


def refshift_tempflux(wavelen,flux,zin,zout):
    """
    WARNING: FLUX DENSITY MUST BE IN units of wavelenght
    
    """
    #print "this are the wavelen axis before and after the redshift"
    #print wavelen#[:10]
    for i in range(0,len(wavelen)):
        wavelen[i]+=((zout-zin)/(1+zin))*wavelen[i]
    #print wavelen#[:10]
    flux=np.array(flux).astype('double')
    DLout=co.lum_dist(0, z_final=zout, omega0_m=0.3, omega0_k=0, omega0_l=0.7,w=-1,unit='km', model='lambdacdm')
    DLout=DLout*10E5

    DLin=co.lum_dist(0, z_final=zin, omega0_m=0.3, omega0_k=0, omega0_l=0.7,w=-1,unit='km', model='lambdacdm')
    DLin=DLin*10E5
    print("the cosmological distances")
    print(DLin)
    print(zin)
    print(DLout)
    print(zout)
    DLratio_inv= DLin/DLout
    z_plus_one_out=1.0+zout
    z_plus_one_in=1.0+zin
    flux=(z_plus_one_in/z_plus_one_out)*flux
    print("the k correction was    ")
    print(DLratio_inv)
    flux=flux*(DLratio_inv**2)
    return wavelen,flux


def refshift_tempflux2(wavelen,flux,zin,zout):
    """
    WARNING: FLUX DENSITY MUST BE IN units of wavelenght
    
    """
    zadd=zout-zin
    zeff=zin+zadd+(zin*zadd)
    #print "this are the wavelen axis before and after the redshift"
    #print wavelen#[:10]
    for i in range(0,len(wavelen)):
        wavelen[i]+=zadd*wavelen[i]
    #print wavelen#[:10]

    flux=np.array(flux).astype('double')
    DLeff=co.lum_dist(0, z_final=zeff, omega0_m=0.3, omega0_k=0, omega0_l=0.7,w=-1,unit='km', model='lambdacdm')
    DLeff=DLeff*10E5

    DLin=co.lum_dist(0, z_final=zin, omega0_m=0.3, omega0_k=0, omega0_l=0.7,w=-1,unit='km', model='lambdacdm')
    DLin=DLin*10E5
    #print "the cosmological distances"
    #print DLin
    #print zin
    #print DLeff
    #print zeff
    DLratio_inv= DLin/DLeff
    z_plus_one_eff=1.0+zeff
    z_plus_one_in=1.0+zin
    flux=(z_plus_one_in/z_plus_one_eff)*flux # https://arxiv.org/pdf/astro-ph/9905116.pdf
    #print "the k correction was    "
    #print DLratio_inv
    flux=flux*(DLratio_inv**2)
    return wavelen,flux,zeff




def fluxlambda2fluxjan(wavelen,flux):
    """
    fluxnu in jankys
    """
    wavelen=np.array(wavelen).astype('float')
    flux=3.34*10E4*(wavelen*wavelen)*np.array(flux).astype('float')
    c=2997924580000000000.0
    nu=c/wavelen
    
    return nu,flux

def fluxlambda2fluxnu(wavelen,flux):
    """
    output in erg s-1 cm-2 hz-1
    """
    wavelen=np.array(wavelen).astype('float')
    c=2997924580000000000.0
    #c2=c*c
    flux=(1.0/c)*(wavelen*wavelen)*np.array(flux).astype('float')
    
    nu=c/wavelen
    
    return nu,flux

def fluxnu2fluxlambda(wavelen,flux):
    """
    output in erg s-1 cm-2 hz-1
    """
    wavelen=np.array(wavelen).astype('float')
    c=2997924580000000000.0
    #c2=c*c
    flux=(1.0/(wavelen*wavelen))*(c)*np.array(flux).astype('float')
    
    #nu=c/wavelen
    
    return flux

def get_median_filter_eazy(filter_file,filter_id,clear=0.01):
    wavelen,out=get_filter_eazy(filter_file,filter_id)
    wavelen,out=fc.open_ascii_cat(filter_name+'.res', unpack=True, vartype=type('float'))
    wavelen=np.array(wavelen).astype('float')
    out=np.array(out).astype('float')
    mask=np.where(out>0.01)
    out_clear=out[mask]
    wavelen_clear=wavelen[mask]
    #print wavelen
    return np.median(wavelen_clear)

def get_mean_filter_eazy(filter_file,filter_id,clear=0.01):
    wavelen,out=get_filter_eazy(filter_file,filter_id)#fc.open_ascii_cat(filter_name+'.res', unpack=True, vartype=type('float'))
    wavelen=np.array(wavelen).astype('float')
    out=np.array(out).astype('float')
    mask=np.where(out>0.01)
    out_clear=out[mask]
    wavelen_clear=wavelen[mask]
    out_mean=np.mean(out_clear)
    out_mean,wavelen_mean,flag=find_nearestxy(wavelen_clear,out_clear,out_mean)
    #print wavelen
    return wavelen_mean

def get_std_filter_eazy(filter_file,filter_id,clear=0.01):

    #FIXME
    
    wavelen,out=get_filter_eazy(filter_file,filter_id)#fc.open_ascii_cat(filter_name+'.res', unpack=True, vartype=type('float'))
    wavelen=np.array(wavelen).astype('float')
    out=np.array(out).astype('float')
    mask=np.where(out>0.01)
    out_clear=out[mask]
    wavelen_clear=wavelen[mask]
    out_std=np.std(out_clear)
    out_mean=np.mean(out_clear)
    out_plus,wavelen_plus,flagplus=find_nearestxy(wavelen_clear,out_clear,out_mean+out_std)
    out_minus,wavelen_minus,flagminus=find_nearestxy(wavelen_clear,out_clear,out_mean-out_std)

    if flagplus==0 and flagminus==0:
        wavelen_std=float(wavelen_minus+wavelen_plus)/2.
    elif flagminus!=0:
        wavelen_std=wavelen_plus
    else:
        wavelen_std=wavelen_minus

    return wavelen_std


def clear_filter(wavs,mags,mags_err=[]):

    wavs=np.array(wavs).astype('float')
    mags=np.array(mags).astype('float')
    mags_err=np.array(mags_err).astype('float')
    mask=np.array(np.logical_and(mags>=-90, mags<=90))
    print("the number of elliminated is")
    print(len(mags[mask])-len(mags))
    print(mask.shape)
    print(wavs.shape)
    print(mags.shape)
    print(mags_err.shape)
    print(wavs[mask].shape)
    print(mags[mask].shape)
    print(mags_err[mask].shape)
    return wavs[mask],mags[mask],mags_err[mask]
    

def convolve_curve(x,y,curve_file, unit="a"):

   
   wavelen,response=fc.open_ascii_cat(curve_file, unpack=True, vartype=type(1.0))
   wavelen=np.array(wavelen).astype('float')
   response=np.array(response).astype('float')

   #print "conv limits"
   #print max(wavelen)
   #print min(wavelen)
   #print max(response)
   #print min(response)


   if unit=='mum':
       wavelen=wavelen*10000
   #print "conv limits 2"
   #print max(wavelen)
   #print min(wavelen)
   #print max(response)
   #print min(response)
   #y_int  = np.interp(wavelen,x,y,left=0.0,right=0.0)              #Interpolate to common wavelength axis
   response_int  = np.interp(x,wavelen,response,left=0.0,right=0.0)

   #print "conv limits 3"
   #print max(wavelen)
   #print min(wavelen)
   #print max(x)
   #print min(x)
   #print max(response_int)
   #print min(response_int)
   #print len(response_int)

   #print max(y)
   #print min(y)
   #print len(y)


   #filt_conv=np.convolve(y,response_int,mode='same')
   filt_conv=y*response_int
   #print "conv limits 4"
   #print max(filt_conv)
   #print min(filt_conv)
   #print len(filt_conv)


   mask=np.where(filt_conv!=0.0)#[0])
   #print mask
   
   filt_conv=filt_conv[mask]
   #print "conv limits 5"
   #print max(filt_conv)
   #print min(filt_conv)
   #print len(filt_conv)

   wavelen_conv=x[mask]
   #Iobject  = simps(flux_phot*filt_int*wavelen,wavelen)                     #Denominator
   #Ibackground= simps(sky_int*filt_int*wavelen,wavelen)
   return wavelen_conv,filt_conv


def eazy_generate_filters(filter_file,description_file,filter_name,R=100,lambda0=4230,limits=[3500,13500],norm=True,centered_in=None, filter_type='gauss',thoughput=0.525,atmtransmission='',mirror_reflectivity='',nmirrors=2, filter_file_full='',desc_file_full='',**kwargs):

    lambs=[]
    fwhms=[]

    filter_names=[]
    filter_ids=[]
    lambs_aux=[]
    std_aux=[]
    amp_aux=[]
    if type(centered_in)!=type(None):

        
        for j in range(0,len(centered_in)):
            
            print("####### Filter ",str(centered_in[j]))
            #print get_mean_filter_eazy(filter_file,centered_in[j])
            wavelen,response=get_filter_eazy(filter_file,centered_in[j])
            idx=np.where(response==max(response))
            print(min(wavelen))
            print(max(wavelen))
            print(response[idx[0]][0])
            print(wavelen[idx[0]][0])
            lambs_aux.append(wavelen[idx[0]][0])#get_mean_filter_eazy(filter_file,centered_in[j])
            std_resp=np.std(response)
            idxs1 = np.abs(np.array(response).astype('float')-(response[idx[0]]+std_resp)).argmin()
            idxs2 = np.abs(np.array(response).astype('float')-(response[idx[0]]-std_resp)).argmin()
            print((wavelen[idxs1]-wavelen[idxs2])/2) 
             
            std_aux.append(abs(wavelen[idxs1]-wavelen[idxs2])/2) #get_std_filter_eazy(filter_file,centered_in[j])
            amp_aux.append(thoughput)
    else:
     
        lamb=limits[0]
        fwhm=(lamb**2)/(R*lambda0)
        std=fwhm/(2*sqrt(2*log(2)))
        lamb=lamb+fwhm
        for k in range(0,R):
            lambs_aux.append(lamb)
            fwhm=(lamb**2)/(R*lambda0)
            std_aux.append(fwhm/(2*sqrt(2*log(2))))
            #x=np.arange(lamb-3*std,lamb+3*std,2.0)
            #print x
            amp_aux.append(thoughput)#0.3-(0.15*(float(i)/100.0))
            lamb=lamb+(fwhm )
    print(lambs_aux)
    print(amp_aux)  
    for i in range(0,len(lambs_aux)):

        
        #fwhm=(lamb**2)/(R*lambda0)
        #std=fwhm/(2*sqrt(2*log(2)))
        x=np.arange(lambs_aux[i]-3*std_aux[i],lambs_aux[i]+3*std_aux[i],2.0)
        #print x
        #amp=0.3-(0.15*(float(i)/100.0))
        if filter_type!='top':
            y=gauss(x,*[amp_aux[i],lambs_aux[i],std_aux[i]])
        else:
            y=np.zeros(x.size)
            mask=np.where(np.logical_and(x<=lambs_aux[i]+std_aux[i], x>=lambs_aux[i]-std_aux[i]))
            y[mask]=thoughput
        #y=gauss_norm(x,[lamb,std])

        mask=np.where(np.logical_and(x>=limits[0], x<=limits[1]))
        x_data=x[mask]#np.around(x[mask],5)
        y_data=y[mask]#np.around(y[mask],5)
        #print "filters limits"
        #print max(x_data)
        #print min(x_data)
        #print max(y_data)
        #print min(y_data)
        if (atmtransmission!='') and (len(x_data)>0):
            x_data,y_data=convolve_curve(x_data,y_data,curve_file=atmtransmission, unit="mum")

        #print "filters limits 2"
        #print max(x_data)
        #print min(x_data)
        #print max(y_data)
        #print min(y_data)


        if (mirror_reflectivity!='')  and (len(x_data)>0):
           for m in range(1,nmirrors+1):
               x_data,y_data=convolve_curve(x_data,y_data,curve_file=mirror_reflectivity, unit="mum")
        #print "filters limits 3"
        #print max(x_data)
        #print min(x_data)
        #print max(y_data)
        #print min(y_data)


        

        if len(x_data)>0:
            x_data_th=x_data.copy()
            y_data_th=y_data.copy()
            y_data=np.array(y_data)/max(y_data)


            lambs.append(lambs_aux[i])
            counter=range(1,len(x_data)+1)
            counter=np.array(counter).astype(int)
            filter_data=np.transpose([counter,x_data,y_data])
            filter_data_full=np.transpose([counter,x_data_th,y_data_th])

            #print "filters limits 4"
            #print max(x_data_th)
            #print min(x_data_th)
            #print max(y_data_th)
            #print min(y_data_th)



            #dt = np.dtype('int, float, float')
            filter_data=filter_data.tolist()#=np.array(filter_data,dt)
            filter_data_full=filter_data_full.tolist()
            for l in range(0,len(filter_data)):
                filter_data[l][0]=int(round(float(filter_data[l][0]),0))
                filter_data_full[l][0]=int(round(float(filter_data[l][0]),0))
                #print "###### filter_Data l 0 "
                #print type(filter_data[l][0])
                #print filter_data[l][0]
            if i==0:
                appendix='\n'
            else:
                appendix="\n"
            fc.save_ascii_cat(data=filter_data,file_name=filter_file,header=[appendix+str(len(x_data))+"  "+filter_name+str(i)],method="a", identifier="")
             
            if filter_file_full!='':
                fc.save_ascii_cat(data=filter_data_full,file_name=filter_file_full,header=[appendix+str(len(x_data_th))+"  "+filter_name+str(i)],method="a", identifier="")
            desc = open(description_file, 'r') 
            desc_data= desc.readlines()
            desc.close()
            #last=desc_data[-1]
            #print desc_data[:5]
            #print description_file
            #print desc_data[5:]
            #print desc_data[-1]
            #print desc_data[-1].split()
            last_filter=len(desc_data)#int(desc_data[-1].split()[0])
            cat_ascii=open(description_file, "w")
            for m in range(0,len(desc_data)):
                desc_data[m]=desc_data[m].rstrip(' ')
                cat_ascii.write(desc_data[m].rstrip('\n')+'\n')
            #if i!=0:
            #    cat_ascii.write("\n")
            cat_ascii.write(" "+str(last_filter+1)+"  "+filter_name+str(i))
            cat_ascii.close()
            #if desc_file_full!='':
               
            
            filter_ids.append(last_filter+1) 
            filter_names.append(filter_name+str(i))
            
        

    #print "novamente os lambs e  filter_id"       
    #print lambs
    #print filter_ids
    return np.array(lambs),np.array(fwhms), np.array(filter_names),np.array(filter_ids)

def plot_filters(filter_file,filter_ids,plot_name,transmission='',reflectivity='',limits=[3500,13500]):

    fig = plt.figure()

    ax1 = fig.add_subplot(111)

    
    c = itertools.cycle(["r", "b", "g","u","k","c","m"])
    #iter(cm.rainbow(np.linspace(0, 1, R)))

    #,colors=next(c)
    #print "#### This is the filter ids"
    #print filter_ids
    for i in range(0,len(filter_ids)):
            #print "filter id "+str(filter_ids[i])
            wavelen,out=get_filter_eazy(filter_file,filter_ids[i])#fc.open_ascii_cat(filters_name[i]+".res", unpack=True, vartype=type('float'))
            ax1.plot(np.array(wavelen).astype('float'),out)

  
    if transmission!='':
       wavelenT,responseT=fc.open_ascii_cat(transmission, unpack=True, vartype=type(1.0))
       wavelenT=np.array(wavelenT).astype('float')*10000
       ax1.plot(np.array(wavelenT).astype('float'),responseT, color='#82a67d', label='sky transmission', ls='--',markersize=2, alpha=0.5)
    if reflectivity!='':

       wavelenR,responseR=fc.open_ascii_cat(reflectivity, unpack=True, vartype=type(1.0))
       wavelenR=np.array(wavelenR).astype('float')*10000
       ax1.plot(np.array(wavelenR).astype('float'),responseR, color='#af6d52', label='Mirror reflectivity', ls=':',markersize=2, alpha=0.5)
    ax1.legend(loc=1, borderaxespad=0., fontsize=8)
    ax1.set_ylim(0,1.05)
    ax1.set_xlim(limits[0],limits[1])
    plt.savefig(plot_name)
    
    return 0

def open_phot_results(cat_file,origin='eazy',cols=(0,20,21,22,24,1)):

    if origin=='bpz' or origin=='BPZ':
        _ids,z,zmax,zmin,odds,zspec=fc.open_ascii_cat(cat_file.rstrip("_bpz.cat")+"_bpz.cat", usecols=(0,1,2,3,5,9), unpack=True, skip_rows=0, vartype='float')
    elif origin=='ascii':
         _ids,z,zmax,zmin,odds,zspec=fc.open_ascii_cat(cat_file, delimiter=",",unpack=True, usecols=cols, vartype='float')
    else:
        #eazy
        _ids,zspec,z,odds,zmin,zmax=fc.open_ascii_cat(cat_file, usecols=(0,1,5,8,9,10), unpack=True, skip_rows=0, vartype='float')
    return _ids,z,zmax,zmin,odds,zspec


def get_sigmaNMAD(delta_z,zspec):

    #delta_z=z-zspec
    median_delta_z=np.median(delta_z)
    #print "median of delta_Z"
    #print median_delta_z
    argument=(delta_z-median_delta_z)/(1.0+zspec)
    #print "median of argument"
    #print argument
    sigmaz=1.48*np.median(np.absolute(argument))
    #print sigmaz
    return sigmaz 


    


def plot_redshifts(files=None,plot_name="photo_z.png", text="",legend="",origin='eazy',asfunctionof=[],labelasfunction="x",sigmas=1, dense=True, catastrophic_lim=0.15, shift=0.03, golden_sel_mag=22.5, golden_sel_mag_min=17.5,plot_gold=False,select_gold_mag=[],golden_sel_odds=0.9,zgold=1.0,plot_z_uncertainties=False,cols=(0,20,21,22,24,1), weights=None, gridsz=(50,50),denseclean=3,delimiter=" ", mag_cat=None,mag_col=None, uncertainties_lim=0.045,zbins=[np.arange(0,9,0.01),60],in_ids=None,in_z=None,in_zmax=None,in_zmin=None,in_odds=None,in_zspec=None, out_format='.png',reduced_stats=False, text_pos=[0.65, 0.23], **kwargs):
    # select_gold_mag=[] shall be a list with the format [catalog_file_path,number_of_mag_column] it assumes that the first column is the _id. _ids in diferent files must be unique # ids must be integer
    #if type(bpz_file)==type('str'):

    if type(files)==type(None):
        files_loop=['array_entry']
        
    else:
        files_loop=files 
        

    rc('text', usetex=True)
    text_cp=text
    if type(legend)==type('string'):
        legend=[legend for x in range(0,len(files_loop))]
    
    colors=['g','b','r','y','m','k','c','b','g','r']

    markers=['o','s','^','s','o','D','v','<','>','s','p','*']
    n=int(float(len(files_loop))/2)



    deltaz_reduced_full=np.array([])
    z_spec_full=np.array([])
    z_full=np.array([])
    _ids_full=np.array([])
    odds_full=np.array([])
    zmax_full=np.array([])
    zmin_full=np.array([])
    mag_full=np.array([])
    for i in range(0,len(files_loop)):
        if type(files)!=type(None):
            _ids,z,zmax,zmin,odds,zspec=open_phot_results(files[i],origin=origin)
        else:
            _ids=in_ids.copy()
            z=in_z.copy()
            zmax=in_zmax.copy()
            zmin=in_zmin.copy()
            odds=in_odds.copy()
            zspec=in_zspec.copy()
        zmax=zmax-z
        zmin=z-zmin
        #zerr=(zmax-zmin)*0.5
        #zerr=np.array(zerr)
        #zerr_reduced=zerr/(1.0+zspec)
        deltaz_reduced=(z-zspec)/(1.0+zspec)
        deltaz_reduced_full=np.concatenate((deltaz_reduced_full,deltaz_reduced))
        z_spec_full=np.concatenate((z_spec_full,zspec))
        z_full=np.concatenate((z_full,z))
        zmin_full=np.concatenate((zmin_full,zmin))
        zmax_full=np.concatenate((zmax_full,zmax))
        _ids_full=np.concatenate((_ids_full,_ids))
        odds_full=np.concatenate((odds_full,odds))
        if type(mag_cat)!=type(None):
            #print(mag_cat[i])
            #print(mag_col)
            mag_=fc.select_column_ascii(mag_cat[i],column=mag_col,header=None,col_name='',identifier='#')
            mag_=np.array(mag_).astype(type(1.0))
            mag_full=np.concatenate((mag_full,mag_))   
    uncertainties_full=(zmax_full-zmin_full)*0.5




    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    if dense==False:
        ax1.errorbar(z_spec_full+((i-n)*shift),z_full,yerr=[zmin_full,zmax_full], fmt=colors[i]+markers[i],markersize=1,label=legend[i])
    else:
        print ("histogram")
        print (zbins)
        ed,h=np.histogram(z_full, bins=zbins[0])
        print (ed)
        print (h)
        
        ax1=fc.density_scatter(np.array(z_spec_full), np.array(z_full), ax = ax1, fig=fig, bins = [zbins[0],zbins[0]], alpha = 0.5, s=8,marker='o', save=False) 
    #print z_spec_full        
        #ax1.plot(z_spec_full+((i-n)*shift),z_full,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
    
    ax1.set_ylabel("$z_{phot}$")
    ax1.set_xlabel("$z_{spec}$")
    if len(legend)>1:
        ax1.legend(loc=1, borderaxespad=0., fontsize=8)

    if ('xlim' in kwargs.keys()):
        xlim=kwargs['xlim']
        ax1.set_xlim(xlim[0],xlim[1])
        ax1.set_ylim(xlim[0],xlim[1])
    else:
        ax1.set_ylim(min(z_spec_full)-0.1,max(z_spec_full)+0.4)
    x_line=np.arange(min(z_spec_full)-0.1,max(z_spec_full)+0.4,0.05)
    ax1.plot(x_line,x_line,'c--')
    plt.savefig(plot_name)
    plot_nameerr=plot_name.split(".png")[0]


    plt.clf()
    #ax7 = fig.add_subplot(111)
    #ax7= sns.kdeplot(z_spec_full, z_full, shade=True,cmap="rainbow", ax=ax7, cbar=True, shade_lowest=False, n_levels=30)

    #ax7.set_ylabel("$z_{phot}$")
    #ax7.set_xlabel("$z_{spec}$")
    #if ('xlim' in kwargs.keys()):
    #    xlim=kwargs['xlim']
    #    ax7.set_xlim(xlim[0],xlim[1])
    #    ax7.set_ylim(xlim[0],xlim[1])
    #else:
    #    ax7.set_ylim(min(z_spec_full)-0.1,max(z_spec_full)+0.4)
    #ax7.plot(x_line,x_line,'c--')

    #plt.savefig(plot_nameerr+"kde.png")


    if len(files_loop)>1:
        plt.clf()
        ax2 = fig.add_subplot(111)
        #golden_lim_plot=np.ones(2)*golden_lim
        for i in range(0,len(files_loop)):
            if type(files)!=type(None):
                _ids,z,zmax,zmin,odds,zspec=open_phot_results(files[i],origin=origin)
            else:
                _ids=in_ids.copy()
                z=in_z.copy()
                zmax=in_zmax.copy()
                zmin=in_zmin.copy()
                odds=in_odds.copy()
                zspec=in_zspec.copy()

            #_ids,z,zmax,zmin,odds,zspec=open_phot_results(files[i],origin=origin)
            #zmax=zmax-z
            #zmin=z-zmin
            #zerr=(zmax-zmin)*0.5
            #zerr=np.array(zerr)
            deltaz_reduced=(z-zspec)/(1.0+zspec)
            ax2.plot(zspec+((i-n)*shift),deltaz_reduced,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
            #ax2.plot(zspec+((i-n)*shift),deltaz_reduced,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
        ax2.axhline(y=catastrophic_lim, ls='dashed', color='red', lw=0.5)
        ax2.axhline(y=catastrophic_lim*-1, ls='dashed', color='red', lw=0.5)
        if len(legend)>1:
            ax2.legend(loc=1, borderaxespad=0., fontsize=8)
        ax2.set_ylabel("$\Delta z/(1+z_{\rm spec})$")
        ax2.set_xlabel("$z_{spec}$")
        ax2.set_ylim(-2,+2)

        #print zspec
    

        ax2.text(0.6, 0.1, text, va='top', fontsize=8, transform=ax1.transAxes)
        plt.savefig(plot_nameerr+"ERR"+out_format)

    #if dense==True:
    plt.clf()
        
    mask1=deltaz_reduced_full<abs(catastrophic_lim) #and  deltaz_reduced_full<0.15
    deltaz_reduced_cat=deltaz_reduced_full[mask1]
    mask2=deltaz_reduced_cat>(-1*abs(catastrophic_lim))
    deltaz_reduced_cat=deltaz_reduced_cat[mask2]
    cat100=float(len(deltaz_reduced_cat))/float(len(deltaz_reduced_full))
    sigmaz_=get_sigmaNMAD(deltaz_reduced_full,z_spec_full)
    bias=get_photo_z_bias(deltaz_reduced_full,z_spec_full)

    #if type(weights)==type(None):
    #    delta_z_full=z_full-z_spec_full
    #else:
    #    delta_z_full=weights*(z_full-z_spec_full)
    #bin_center_z,bin_vals_z,bin_std=fc.make_bin(z_spec_full,delta_z_full,bin_edges=np.arange(0.0,max(z_spec_full)+0.1,0.1), return_vals=False)
    mean_bias,median_bias,sigma68=get_zstatistic(z_full,z_spec_full,weights=weights)
    text="$\sigma_{NMAD}$: $"+str(round(sigmaz_,3))+"$"
    if reduced_stats==False:
        text+="\n $\sigma_{68}$: $"+str(round(sigma68,3))+"$"
        text+="\n Catastrofic: $"+str(round(1.0-round(cat100,3),4))+"$"
        #text+="\n Bias: $"+str(round(bias,4))+"$"
        text+="\n Mean Bias: $"+str(round(mean_bias,3))+"$"
   
    text+="\n Median Bias: $"+str(round(median_bias,4))+"$"
    if reduced_stats==False:
        text+="\n Number of objects: $"+str(len(z_spec_full))+"$"


        #fc.cat_hist_plot(np.array(z_spec_full),plot_name=plot_nameerr+'zsepcfull-histV3.png', xlabel='$z_{spec}$',rang=(0,0.8))
        #plt.clf()
    ax5 = fig.add_subplot(111)
    
    if dense==True:
        ax5=fc.density_scatter(np.array(z_spec_full), np.array(deltaz_reduced_full), ax = ax5,fig=fig, bins = zbins, alpha = 0.5, s=5,marker='o', save=False)         
        #ax5.plot(z_spec_full,deltaz_reduced_full,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
    else:
        ax5.plot(np.array(z_spec_full), np.array(deltaz_reduced_full),colors[0]+markers[0], markersize=1.0)
        

    ax5.axhline(y=2*sigma68, ls='dashed', color='red', lw=0.5)
    ax5.axhline(y=2*sigma68*-1, ls='dashed', color='red', lw=0.5)
    
    #ax5.axhline(y=catastrophic_lim, ls='dashed', color='red', lw=0.5)
    #ax5.axhline(y=catastrophic_lim*-1, ls='dashed', color='red', lw=0.5)

    ax5.set_ylabel("$\Delta z/(1+z_{\rm spec})$", fontsize=20, labelpad=1)
    ax5.set_xlabel("$z_{spec}$", fontsize=20, labelpad=-1)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if len(legend)>1:
        ax5.legend(loc=1, borderaxespad=0., fontsize=8)
    if ('xlim' in kwargs.keys()):
        xlim=kwargs['xlim']
        ax5.set_xlim(xlim[0],xlim[1])
    if ('ylim' in kwargs.keys()):
        ylim=kwargs['ylim']
    else: 
        ylim=[-1,1]
    ax5.set_ylim(ylim[0],ylim[1])

    ax5.text(text_pos[0], text_pos[1], text, va="top", fontsize=12, transform=ax5.transAxes)
    plt.savefig(plot_nameerr+"ERRTotal"+out_format)

    plt.clf()
    if dense==True:      
        ax7 = fig.add_subplot(111)
        ax7.axhline(y=0.8, ls='dashed', color='red', lw=0.5)
        ax7.axhline(y=0.9, ls='dashed', color='red', lw=0.5)
        ax7=fc.density_scatter(np.array(z_spec_full), np.array(odds_full), ax = ax7,fig=fig, bins = zbins, alpha = 0.5, s=5,marker='o', save=False)         
        #ax5.plot(z_spec_full,deltaz_reduced_full,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
        ax7.set_ylim(0.6,+1)
        ax7.set_ylabel("Odds")
        ax7.set_xlabel("$z_{spec}$")
        if len(legend)>1:
            ax7.legend(loc=1, borderaxespad=0., fontsize=8)
        if ('xlim' in kwargs.keys()):
            xlim=kwargs['xlim']
            ax7.set_xlim(xlim[0],xlim[1])
        ax7.text(0.65, 0.23, text, va="top", fontsize=8, transform=ax5.transAxes)
        plt.savefig(plot_nameerr+"Odds"+out_format)



        plt.clf()

        ax8 = fig.add_subplot(111)
        ax8.axhline(y=0.05, ls='dashed', color='red', lw=0.5)
        ax8.axhline(y=0.05, ls='dashed', color='red', lw=0.5)
        ax8=fc.density_scatter(np.array(z_spec_full), np.array(uncertainties_full), ax = ax8,fig=fig, bins = zbins, alpha = 0.5, s=5,marker='o', save=False)         
        #ax5.plot(z_spec_full,deltaz_reduced_full,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
        ax8.set_ylim(-0.1,+0.1)
        ax8.set_ylabel("$\sigma z$")
        ax8.set_xlabel("$z_{spec}$")
        if len(legend)>1:
            ax8.legend(loc=1, borderaxespad=0., fontsize=8)
        if ('xlim' in kwargs.keys()):
            xlim=kwargs['xlim']
            ax8.set_xlim(xlim[0],xlim[1])
        ax8.text(0.65, 0.23, text, va="top", fontsize=8, transform=ax5.transAxes)
        plt.savefig(plot_nameerr+"dense_uncertainties"+out_format)
        #plt.clf()
        #fc.cat_hist_plot(z_spec_full,plot_name=plot_nameerr+'zsepcfull-hist.png', xlabel='$z_{spec}$',rang=(0,0.8))
        #plt.clf()
        #fc.cat_hist_plot(np.array(z_spec_full),plot_name=plot_nameerr+'zsepcfull-histV2.png', xlabel='$z_{spec}$',rang=(0,0.8))

        #ax4 = fig.add_subplot(111)
        #ax4= sns.kdeplot(z_spec_full, deltaz_reduced_full, shade=True,cmap="rainbow", ax=ax4, cbar=True, shade_lowest=False, n_levels=30)
        #ax4.axhline(y=catastrophic_lim, ls='dashed', color='red', lw=0.5)
        #ax4.axhline(y=catastrophic_lim*-1, ls='dashed', color='red', lw=0.5)
        
        
        

         
        #ax4.set_ylim(-1,+1)
        #ax4.set_ylabel("$\Delta z/(1+z)$")
        #ax4.set_xlabel("$z_{spec}$")
        #ax4.text(0.7, 0.23, text, va='top', fontsize=8, transform=ax4.transAxes)
        #if ('xlim' in kwargs.keys()):
        #    xlim=kwargs['xlim']
        #    ax4.set_xlim(xlim[0],xlim[1])
        #plt.savefig(plot_nameerr+"ERRdense.png")


    if plot_gold==True:
        plt.clf()
        data_cat=fc.open_ascii_cat(select_gold_mag[0], delimiter=delimiter,unpack=False)
        #print data_cat.shape
        magi_cat=fc.select_column_ascii(data_cat,column=select_gold_mag[1],header=None,col_name='',identifier='#')
        #print magi_cat.shape
        _idfromcat=fc.select_column_ascii(data_cat,column=0,header=None,col_name='',identifier='#')
        #print _idfromcat.shape 
        #print _ids_full.shape
        #print _ids_full[:10]
        #print _idfromcat[:10]
        print(magi_cat[0])
        print(data_cat[0])
        magi_cat=np.array(magi_cat).astype(type(1.0))
        _ids_full=np.array(_ids_full).astype(type(1))
        _idfromcat=np.array(_idfromcat).astype(type(1))
        print ('Here comes the problem ',len(magi_cat),len(mag_full), len(_ids_full),len(_idfromcat))
        print (magi_cat[:5])
        print (mag_full[:5])
        print ('Ids ')
        print (_idfromcat[:5])
        print (_ids_full[:5]) 
        mag_i=fc.selectfromlist(target_ids=_ids_full,_ids=_idfromcat,_ids_feature=magi_cat)
        print ('After the match ',len(mag_i),len(mag_full), len(_ids_full),len(_idfromcat))

        print (mag_i[:5])
        print (mag_full[:5])
        print ('Ids ')
        print (_idfromcat[:5])
        print (_ids_full[:5])

        mask_gold=(mag_i>golden_sel_mag_min) & (odds_full>golden_sel_odds) & (z_spec_full < zgold) &  (mag_i < golden_sel_mag) & (uncertainties_full<uncertainties_lim)
        
        numzspec_before=np.where(z_spec_full<0.11)


        z_spec_gold=z_spec_full[mask_gold]

        numzspec_after=np.where(z_spec_full<0.11)
        numzspec_gold=np.where(z_spec_gold<0.11)
        print ('n of zspec below 0.11 ',numzspec_before[0][:5],numzspec_after[0][:5], numzspec_gold[0][:5])
        print ('n of zspec below 0.11 ',len(numzspec_before[0]),len(numzspec_after[0]), len(numzspec_gold[0]))
        z_full_gold=z_full[mask_gold]
        deltaz_reduced_gold=deltaz_reduced_full[mask_gold]
        uncertainties_full_gold=uncertainties_full[mask_gold]
        frac_total=float(len(deltaz_reduced_gold))/float(len(deltaz_reduced_full))         
        if type(mag_cat)!=type(None):
            mag_full_gold=mag_full[mask_gold] 
        mask1=deltaz_reduced_gold<abs(catastrophic_lim) #and  deltaz_reduced_full<0.15
        deltaz_reduced_gold_cat=deltaz_reduced_gold[mask1]
        mask2=deltaz_reduced_gold_cat>(-1*abs(catastrophic_lim))
        deltaz_reduced_gold_cat=deltaz_reduced_gold_cat[mask2]
        cat100=float(len(deltaz_reduced_gold_cat))/float(len(deltaz_reduced_gold))
   


#        text="$\sigma_{NMAD}$: $"+str(round(sigmaz_,3))+"$"
#        text+="\n $\sigma_{68}$: $"+str(round(sigma68,3))+"$"
#        text+="\n Catastrofic: $"+str(round(1.0-round(cat100,3),4))+"$"
#        text+="\n Bias: $"+str(round(bias,3))+"$"
        




        sigmaz_=get_sigmaNMAD(deltaz_reduced_gold,z_spec_gold)
        text_gold="$\sigma_{NMAD}$: "+str(round(sigmaz_,3))
        bias=get_photo_z_bias(deltaz_reduced_gold,z_spec_gold)
        mean_bias,median_bias,sigma68=get_zstatistic(z_full_gold,z_spec_gold,weights=weights)
        text_gold+="\n $\sigma_{68}$: $"+str(round(sigma68,3))+"$"
        text_gold+="\n Catastrofic: $"+str(round(1.0-round(cat100,3),4))+"$"
        text_gold+="\n Bias: "+str(round(bias,3))
        text_gold+="\n Mean Bias: $"+str(round(mean_bias,3))+"$"
        text_gold+="\n Median Bias: $"+str(round(median_bias,3))+"$" 
        text_gold+="\n fraction of objects: "+str(round(frac_total,3))

        ax6 = fig.add_subplot(111)
        ax6.axhline(y=2*sigma68, ls='dashed', color='red', lw=0.5)
        ax6.axhline(y=2*sigma68*-1, ls='dashed', color='red', lw=0.5)
        #ax6.axhline(y=catastrophic_lim, ls='dashed', color='red', lw=0.5)
        #ax6.axhline(y=catastrophic_lim*-1, ls='dashed', color='red', lw=0.5)
        ax6=fc.density_scatter(np.array(z_spec_gold), np.array(deltaz_reduced_gold), ax = ax6,fig=fig, bins = zbins, alpha = 0.5, s=5,marker='o', save=False)         
        #ax6.plot(z_spec_gold,deltaz_reduced_gold,colors[i]+markers[i], markersize=1.0,  label=legend[i])#ls=":"ls=":"
        ax6.set_ylim(-1,+1)

        ax6.set_ylabel("$\Delta z/(1+z)$")
        ax6.set_xlabel("$z_{spec}$")
        if len(legend)>1:
            ax6.legend(loc=1, borderaxespad=0., fontsize=8)
        if ('xlim' in kwargs.keys()):
            xlim=kwargs['xlim']
            ax6.set_xlim(xlim[0],xlim[1])
        ax6.text(0.65, 0.23, text_gold, va='top', fontsize=8, transform=ax6.transAxes)
        plt.savefig(plot_nameerr+"ERRTotal_gold"+out_format) #plot_name=work_dir+'L-level-hist.png'
        #plt.clf()
        #fc.cat_hist_plot(mag_i,plot_name=plot_nameerr+'magi-hist.png', xlabel='mag_i')
        #plt.clf()   
        #fc.cat_hist_plot(z_spec_gold,plot_name=plot_nameerr+'zsepcgold-hist.png', xlabel='$z_{spec}$', rang=(0,0.8))
    #ax1.text(0.6, 0.1, text, va='top', fontsize=8, transform=ax1.transAxes)
    #ax2.plot(zspec,zerr_reduced,colors[i]+markers[i], markersize=5.0, ls=":", label=legend[i])#ls=":"
    if len(asfunctionof)>0:
        plt.clf()
        ax3 = fig.add_subplot(111)
        for i in range(0,len(files_loop)):
            if type(files)!=type(None):
                _ids,z,zmax,zmin,odds,zspec=open_phot_results(files[i],origin=origin)
            else:
                _ids=in_id.copy()
                z=in_z.copy()
                zmax=in_zmax.copy()
                zmin=in_zmin.copy()
                odds=in_odds.copy()
                zspec=in_zspec.copy()
            #_ids,z,zmax,zmin,zspec=open_phot_results(files[i],origin=origin)
            #zmax=zmax-z
            #zmin=z-zmin
            zerr=(zmax-zmin)*0.5
            zerr=np.absolute(zerr)
            #print(zerr)
            #print(np.array(asfunctionof))
            if dense==True:
                ax3=fc.density_scatter(np.array(asfunctionof), np.array(zerr), ax = ax3,fig=fig, bins = zbins, alpha = 0.5, s=5,marker='o', save=False)         
            else:
                ax3.plot(np.array(asfunctionof)+((i-n)*shift),zerr,c=colors[i], markersize=5.0, ls=":", label=legend[i])
            ax3.set_ylabel("$\sigma z$")
            ax3.set_xlabel(labelasfunction)
        plt.savefig(plot_nameerr+labelasfunction+"vsERR"+out_format)
        plt.clf()
        ax3 = fig.add_subplot(111)
        for i in range(0,len(files_loop)):
            if type(files)!=type(None):
                _ids,z,zmax,zmin,odds,zspec=open_phot_results(files[i],origin=origin)
            else:
                _ids=in_id.copy()
                z=in_z.copy()
                zmax=in_zmax.copy()
                zmin=in_zmin.copy()
                odds=in_odds.copy()
                zspec=in_zspec.copy()


            #_ids,z,zmax,zmin,zspec=open_phot_results(files[i],origin=origin)
            #zmax=zmax-z
            #zmin=z-zmin
            zerr=(zmax-zmin)*0.5
            zerr=np.absolute(zerr)
            if dense==True:
                ax3=fc.density_scatter(np.array(asfunctionof), np.array(z), ax = ax3,fig=fig, bins = zbins, alpha = 0.5, s=5,marker='o', save=False)         
            else:
                ax3.plot(np.array(asfunctionof)+((i-n)*shift),z,c=colors[i], markersize=5.0, ls=":", label=legend[i])
            ax3.set_ylabel("$z$")
            ax3.set_xlabel(labelasfunction)
            plt.savefig(plot_nameerr+labelasfunction+"vsZ"+out_format)
    if plot_gold==True:
        if type(mag_cat)!=type(None): 
            return z_spec_full,z_full,deltaz_reduced_full,uncertainties_full, z_spec_gold,z_full_gold, deltaz_reduced_gold, uncertainties_full_gold,mag_full,mag_full_gold
        else:
            return z_spec_full,z_full,deltaz_reduced_full,uncertainties_full, z_spec_gold,z_full_gold, deltaz_reduced_gold,uncertainties_full
    else:
        if type(mag_cat)!=type(None): 
            return z_spec_full,z_full,deltaz_reduced_full,uncertainties_full,mag_full
        else:
            return z_spec_full,z_full,deltaz_reduced_full,uncertainties_full


def get_zstatistic(z,zspec,weights=None):
    if type(weights)==type(None):
        deltaz=z-zspec
        mean_bias=np.mean(deltaz)
    else:
        deltaz=weights*(z-zspec)
        mean_bias=np.average(deltaz,weights=weights)
    median_bias=np.percentile(deltaz,50)
    p16=np.percentile(deltaz,15.85)
    p84=np.percentile(deltaz,84.05)
    sigma68=0.5*(p84-p16)
    return mean_bias,median_bias,sigma68

def get_photo_z_bias(delta_z,zspec):

    return np.median(delta_z/(1.0+zspec))

def modest_class(catalog_,mtype='gal',sm_key='SPREAD_MODEL_I',smerr_key='SPREADERR_MODEL_I',mag_key='MAG_AUTO_I',wsm_key='WAVG_SPREAD_MODEL_I'):
    """
    modest_class selection criteria to separate star and galaxies. 
    The definitions can be seen in table 6 of https://arxiv.org/pdf/1708.01531.pdf
    mtype=gal represents high confidance galaxies (class 1) and mtype !=gal represents everything that is not a star (not class 0,2,3)
    """
    #  table 6


    test_field=np.array(catalog_.field(str(sm_key)))+((5/3)*np.array(catalog_.field(str(smerr_key))))

    #print(len(test_field))
    print(' tamanho do catalogo',len(catalog_))
    #print(len(test_field > 0.002))
    #print(test_field > 0.002)
    if mtype=='highpure':
        cat_out=catalog_[test_field > 0.005]
    if mtype=='highcomplete':
        cat_out=catalog_[test_field > 0.002]
    if mtype=='highpure2':
        cat_out=catalog_[test_field > 0.01]
    #test_field=np.array(cat_out.field(str(sm_key)))+((5/3)*np.array(cat_out.field(str(smerr_key))))
    #test_field=np.abs(test_field)
    #cat_out=cat_out[test_field > 0.002]
    if mtype=='highgal':
        cat_out=catalog_[test_field > 0.005]
        #test_field2=not((cat_out.field(str(wsm_key)) <0.002) and (cat_out.field(str(mag_key))<21.5))
        test_field_aux=np.abs(cat_out.field(str(wsm_key)))
        cat_out=cat_out[test_field_aux>0.002]
        test_field_aux_2=cat_out.field(str(mag_key))
        cat_out=cat_out[test_field_aux_2>21.5]
    
    return cat_out


def crazy_colorcut(cat,magi_key='MAG_AUTO_I',magr_key='MAG_AUTO_R',magg_key='MAG_AUTO_G',magz_key='MAG_AUTO_Z'):#wsm_key='WAVG_SPREAD_MODEL_I')
    # throw away g-r,r-i,i-z > 4 or <-1

    test_field1=np.array(cat.field(str(magg_key)))-(np.array(cat.field(str(magr_key))))
    #test_field2=np.array(cat.field(str(magr_key)))-(np.array(cat.field(str(magi_key))))
    #test_field3=np.array(cat.field(str(magi_key)))-(np.array(cat.field(str(magz_key))))
    cat_out=cat[test_field1 < 4]
    test_field1=np.array(cat_out.field(str(magg_key)))-(np.array(cat_out.field(str(magr_key))))
    cat_out=cat_out[test_field1>-1]
    
    test_field2=np.array(cat_out.field(str(magr_key)))-(np.array(cat_out.field(str(magi_key))))
    cat_out=cat_out[test_field2 < 4]
    test_field2=np.array(cat_out.field(str(magr_key)))-(np.array(cat_out.field(str(magi_key))))
    cat_out=cat_out[test_field2>-1]

    test_field3=np.array(cat_out.field(str(magi_key)))-(np.array(cat_out.field(str(magz_key))))
    cat_out=cat_out[test_field3 < 4]
    test_field3=np.array(cat_out.field(str(magi_key)))-(np.array(cat_out.field(str(magz_key))))
    cat_out=cat_out[test_field3>-1]


    return cat_out

def deltaz_hist(zphot,zspec,legend_names=[],out_name='deltaz_hist.png', bins=10,normed=False,alpha=None, fill=False, hatch=['//','-', '+', 'x', '\\', '*', 'o', 'O', '.'],fit_gauss=True):
    if type(alpha)==type(None):
        alpha=[1.0 for i in range(0,len(zphot))]
    deltaz_reduced=[]
    for i in range(0,len(zphot)):
        z=np.array(zphot[i])
        ztr=np.array(zspec[i])
        deltaz=(z-ztr)/(1.0+ztr)
        deltaz_reduced.append(deltaz)
    
    fc.cat_hist_plot(deltaz_reduced,plot_name=out_name, xlabel='$\delta_{z}/{1+z}$',bins=bins,legend_names=legend_names,normed=normed,plt_type="bar",alpha=alpha, hatch=hatch, fit_gauss=fit_gauss)


def fraction_cuts(num_objects,num_objects_full,norm_factor,cut_labels):



    fig = plt.figure(figsize=figsize)         
    grid = plt.GridSpec(3, 1, hspace=0.0, wspace=0.0, 
                        left=0.15, right=0.9, bottom=0.1, top=0.9)

    ax_cumulative = fig.add_subplot(grid[0,0])#, sharex=ax_main)

    
    ax_absolute = fig.add_subplot(grid[1, 0],sharex=ax_bias)


    ax_cumulative_norm_step = fig.add_subplot(grid[0,0])#, sharex=ax_main)

def zbin_plots(zphot,var_ref,zspec,bins,uncertainties=None,surveys=[],figsize=(6, 8),plot_name="zbins.png",outliers=True, xlabel='$z_{phot}$', print_bias=True, outliers_joint=True):
    
    colors=['navy','indianred','g','b','r','y','m','k','c','b','g','r']

    markers=['o','s','^','s','o','D','v','<','>','s','p','*']

    fig = plt.figure(figsize=figsize)
    if outliers_joint==False:         
        grid = plt.GridSpec(4, 1, hspace=0.0, wspace=0.0, 
                        left=0.15, right=0.9, bottom=0.1, top=0.9)
    else:
        grid = plt.GridSpec(3, 1, hspace=0.01, wspace=0.0, 
                        left=0.2, right=0.9, bottom=0.1, top=0.9)

    
    ax_bias = fig.add_subplot(grid[0,0])#, sharex=ax_main)

    
    ax_sigma = fig.add_subplot(grid[1, 0],sharex=ax_bias)



    if outliers==True:
        if outliers_joint==False:
            ax_out2 = fig.add_subplot(grid[2, 0],sharex=ax_bias)
            ax_out3 = fig.add_subplot(grid[3, 0],sharex=ax_bias)
        else:
            ax_out2 = fig.add_subplot(grid[2, 0],sharex=ax_bias)
            #ax_out3 = ax_out2

    if len(surveys)==0:
        zphot=[zphot]
        zspec=[zspec]
        var_ref=[var_ref]
        survey=[" "]
    if type(uncertainties)!=type(None):
        if len(surveys)==0:
            uncertainties=[uncertainties]
    medians_un=[]
    lun=[]
    hun=[]
    for i in range(0,len(zphot)):


        bin_center,var_valx,deltaz_bin=fc.make_bin(col1_dat=var_ref[i],col2_dat=zphot[i]-zspec[i],bin_edges=bins,return_vals=True,method="mean")

        if type(uncertainties)!=type(None):
            bin_un,var_valx_un,un_bin=fc.make_bin(col1_dat=var_ref[i],col2_dat=uncertainties[i],bin_edges=bins,return_vals=False,method="mean")
            median_un_bin=[]
            uper_bin=[]
            low_bin=[]
        #out2sig_bin=[]
        #out3sig_bin=[]


        median_bias_bin=[]
        sigma68_bin=[]
        out2sig_bin=[]
        out3sig_bin=[] 
        print("Analisando o Survey", surveys[i]) 
        for j in range(len(deltaz_bin)):
            if len(deltaz_bin[j])!=0:
                median_bias_bin.append(np.percentile(deltaz_bin[j],50))
                print("Objects in the bin centered in ", str(bin_center[j])," : ",str(len(deltaz_bin[j])) ) 
                p16=np.percentile(deltaz_bin[j],15.85)
                p84=np.percentile(deltaz_bin[j],84.05)
                sig=0.5*(p84-p16)
                sigma68_bin.append(sig)

                if type(uncertainties)!=type(None):
                    p16un=np.percentile(un_bin[j],15.85)
                    p84un=np.percentile(un_bin[j],84.05)
                    median_un_bin.append(np.percentile(deltaz_bin[j],50))
                    uper_bin.append(p84un)
                    low_bin.append(p16un)
            
                total=float(len(deltaz_bin[j]))
                #print("total de objetos no bin ",total ) 
                sig2_h=np.percentile(deltaz_bin[j],97.725)
                sig2_l=np.percentile(deltaz_bin[j],2.275)
                #print("2d Sigma Intervals") 
                #print(sig2_h)
                #print(sig2_l)
                delta_aux1=fc.select_subset_arr(deltaz_bin[j],sig2_h,">")
                delta_aux2=fc.select_subset_arr(deltaz_bin[j],sig2_l,"<")
                #print("NUmero de objetos no intervalo superior e inferior")
                #print(len(delta_aux1))
                #print (len(delta_aux2))
                out2=float(len(delta_aux1)+len(delta_aux2))/total
                #print("Fracao de outliers")
                #print(out2)
                delta_aux1=[]
                delta_aux2=[]
                sig3_l=np.percentile(deltaz_bin[j],0.135) 
                sig3_h=np.percentile(deltaz_bin[j],99.865)
                delta_aux1=fc.select_subset_arr(deltaz_bin[j],sig3_h,">")
                delta_aux2=fc.select_subset_arr(deltaz_bin[j],sig3_l,"<")
                out3=float(len(delta_aux1)+len(delta_aux2))/total
                delta_aux1=[]
                delta_aux2=[]
                out2sig_bin.append(out2)
                out3sig_bin.append(out3)
            else:
                print("Empty bin", bins[j],"  ",bins[j+1])
                median_bias_bin.append(np.nan) 
                sigma68_bin.append(np.nan)
                out2sig_bin.append(np.nan)
                out3sig_bin.append(np.nan)
        ax_sigma.plot(bin_center,sigma68_bin, color=colors[i],marker=markers[i], alpha=0.9,markersize=1.0,  label=surveys[i])
        if print_bias==True:
            print("Bias in ", surveys[i])
            print(np.transpose([np.array(bin_center).astype(str),np.array(median_bias_bin).astype(str)]))
        ax_bias.plot(bin_center,median_bias_bin, color=colors[i],marker=markers[i], alpha=0.9,markersize=1.0,  label=surveys[i])
        if outliers==True:
            ax_out2.plot(bin_center,out2sig_bin, color=colors[i],marker=markers[i], alpha=0.9,markersize=1.0,  label=surveys[i])
            print("Bias in ", surveys[i])
            print(np.transpose([np.array(bin_center).astype(str),np.array(out2sig_bin).astype(str),np.array(out3sig_bin).astype(str)]))
            if outliers_joint==True:
                out3color=colors[i+1]
                out3marker=markers[i+1]
                ax_out2.plot(bin_center,out3sig_bin, color=out3color,marker=out3marker, alpha=0.9,markersize=1.0,  label=surveys[i])
            else:        
                ax_out3.plot(bin_center,out3sig_bin, color=out3color,marker=out3marker, alpha=0.9,markersize=1.0,  label=surveys[i])
        if type(uncertainties)!=type(None):
            medians_un.append(median_un_bin)
            lun.append(uper_bin)
            hun.append(low_bin)
        ax_sigma.set_ylim(+0.001,0.030)
        ax_bias.set_ylim(-0.011,+0.0050)
        ax_bias.axhline(y=0, ls='dashed', color='black', lw=0.5)

        ax_bias.set_xlim(bins[0],bins[-1])
        ax_bias.tick_params(labelsize=14)
        #plt.setp(ax_bias.get_xticklabels(), visible=True, fontsize=14)
        #plt.setp(ax_sigma.get_xticklabels(), visible=True, fontsize=14)
        #plt.yticks()
        ax_bias.set_ylabel("$\overline{\Delta z}$", fontsize=14, rotation=0)
        ax_sigma.set_ylabel("$\sigma_{68}$", fontsize=14, rotation=0)

        plt.setp(ax_sigma.get_xticklabels(), visible=False)

        plt.setp(ax_bias.get_xticklabels(), visible=False)

        if outliers==True:
            if outliers_joint==False:
                ax_out3.set_ylabel("$out_{3\sigma}$", fontsize=16, rotation=0)
                ax_out2.set_ylabel("$out_{2\sigma}$", fontsize=16, rotation=0)
                ax_out3.set_xlabel(xlabel,fontsize=16)
                ax_out2.set_ylim(0.04,+0.06)
                ax_out3.set_ylim(0,+0.005)
                plt.setp(ax_out2.get_xticklabels(), visible=False)
            else:
                ax_out2.set_ylabel("$outliers{2\sigma}$", fontsize=16, rotation=0)
                ax_out2.set_xlabel(xlabel,fontsize=16)
                ax_out2.set_ylim(0.0001,+0.095)



        else:
            ax_sigma.set_xlabel(xlabel, fontsize=16)
            plt.setp(ax_sigma.get_xticklabels(), visible=True)


        
        if  len(surveys)>1:   
            ax_bias.legend(loc=1, borderaxespad=0., fontsize=8)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig(plot_name)
    if type(uncertainties)!=type(None):
        plt.clf()
        ax1 = fig.add_subplot(111)
        for k in range(0,len(uncertainties)):
            medians_un.append(median_un_bin)
            #lun.append(uper_bin)
            #hun.append(low_bin)
            ax1.errorbar(bin_un,medians_un[k],yerr=[lun[k],hun[k]], alpha=0.9,markersize=1.0,label=surveys[k])
            ax1.set_ylabel("$\sigma z$")
            ax1.set_xlim(bins[0],bins[-1])
            ax1.set_ylim(-0.2,0.2)
            if  len(surveys)>1:   
                ax1.legend(loc=1, borderaxespad=0., fontsize=8)
            plt.savefig(plot_name.rstrip(".png")+"uncertainties.png")

    return 0


def get_filter_eazy(filter_file,filter_id):
    desc = open(filter_file, 'r') 
    desc_data= desc.readlines()
    desc.close()
    count=1
    current_line=0
    wavelenF=[]
    response=[]
    for i in range(0,filter_id-1):
        filter_lines=int(desc_data[current_line].split()[0])
        current_line+=filter_lines+1
        count+=1
    #print "current line"
    #print current_line
    filter_lines=int(desc_data[current_line].split()[0])
    #current_line+=1
    filter_name=desc_data[current_line]
    #print "##### current filter"
    #print filter_name 
    for j in range(0,filter_lines):
        current_line+=1
        data_aux=desc_data[current_line].split()
        wavelenF.append(data_aux[1])
        response.append(data_aux[2]) 
    
    return np.array(wavelenF).astype('float'),np.array(response).astype('float')    
    
def get_filter_cols(filter_file,col,noatm=False):
    if noatm==False:
        wavelenF,filt=fc.open_ascii_cat(filter_file, unpack=True,usecols=(0,col), vartype=type(1.0))
    else:
        wavelenF,filt,atm=fc.open_ascii_cat(filter_file, unpack=True,usecols=(0,col,-1), vartype=type(1.0))
        filt=np.divide(filt,atm)
    return wavelenF,filt  

def template2magnitudes(template,filter_file,filter_ids=None,sigma_err=0.000, plot=True,plot_dir='', filter_source='eazy', filter_unit='A'):
    
    c_AAs     = 2.99792458e18                          # Speed of light in Angstrom/s
    wavelen,flux=fc.open_ascii_cat(template, unpack=True, vartype=type(1.0))
    mAB=[]
    if filter_source=='bpz' or filter_source=='BPZ' or type(filter_ids)==type(None):
        filter_ids=filter_file

    for i in range(0,len(filter_ids)):
        if filter_source=='eazy':
            wavelenF,filt=get_filter_eazy(filter_file,int(filter_ids[i]))#
        elif filter_source=='ctio-decam':
            wavelenF,filt=get_filter_cols(filter_file,int(filter_ids[i]))
        elif filter_source=='ctio-decam-noatm':
             wavelenF,filt=get_filter_cols(filter_file,int(filter_ids[i]),noatm=True)
        else:
            wavelenF,filt=fc.open_ascii_cat(filter_ids[i], unpack=True, vartype=type(1.0))      
        if filter_unit=='nm':
            wavelenF=wavelenF*10.0  
        wavelen=np.array(wavelen).astype(float)
        wavelenF=np.array(wavelenF).astype(float)
        filt=np.array(filt).astype(float)
        flux=np.array(flux).astype(float)#
        filt_int  = np.interp(wavelen,wavelenF,filt,left=0.0,right=0.0)              #Interpolate to common wavelength axis
        integrand=flux*filt_int*wavelen
        mask=integrand>0.0
        I1        = simps(flux[mask]*filt_int[mask]*wavelen[mask],wavelen[mask])
        mask2=filt_int>0.0
        I2        = simps(  filt_int[mask2]/wavelen[mask2],wavelen[mask2])                     #Numerator
        fnu       = (I1/I2) / c_AAs                          #Average flux density https://astronomy.stackexchange.com/questions/16286/how-can-i-convolve-a-template-spectrum-with-a-photometric-filter-response-spectr
        #mag_cal=-2.5*log10(fnu) - 48.6
        
        #I2a=simps(wavelen[mask2]*filt_int[mask2],wavelen[mask2])

        ##I3=simps(filt_int[mask2]/wavelen[mask2],wavelen[mask2]) 
        #mag_cal=-2.5*log10(fl_mean) - 2.5*log10((I2/I3)/c_AAs)-48.6

        #print "magnitudes comparison"
        #fl_mean=(I1/I2a)
        #print fl_mean
        #mag_cal2=-2.5*log10(fl_mean) - 2.5*log10((I2a/I3)/c_AAs)-48.6
        #print mag_cal
        #fl_mean=fl_mean / c_AAs
        #print fl_mean
        #mag_cal=-2.5*log10(fl_mean) - 2.5*log10((I2/I3)/c_AAs)-48.6
        #print mag_cal
        mag_cal=-2.5*log10(fnu)-48.6

        if sigma_err>0.0:
            mag_cal=random.gauss(mag_cal, sigma_err)
        #if math.isnan(mag_cal):
        #    mag_cal=-99
        mAB.append(mag_cal)            #AB magnitude


    return np.array(mAB).astype('float')


def add_spectro_err(flux,wavelen,telescope_area,object_size,exposure,filter_ids=[295],filter_file='/home/cleciobom/lib/eazy-photoz/filters/FILTER.RES.latest',filter_source='eazy', filter_unit='A', thoughput=0.525,atmtransmission='/home/cleciobom/lib/eazy-photoz-master/cptrans_zm_23_10.txt',mirror_reflectivity='/home/cleciobom/lib/eazy-photoz-master/alreflectivity.txt',nmirrors=2, sky_flux='skybg_50_10.datALL.txt',noise=0.05):

    if (atmtransmission!='') and (len(wavelen)>0):
        wavelen,flux=convolve_curve(wavelen,flux,curve_file=atmtransmission, unit="mum")

        #print "filters limits 2"
        #print max(x_data)
        #print min(x_data)
        #print max(y_data)
        #print min(y_data)


    #if (mirror_reflectivity!='')  and (len(wavelen)>0):
    #    for m in range(1,nmirrors+1):
    #        wavelen,flux=convolve_curve(wavelen,flux,curve_file=mirror_reflectivity, unit="mum")
    #print "filters limits 3"

    #Fluxes are given in erg s-1 cm-2 A-1 at the distance of 10 pc

    #here we start calculating the number of photons in order to get the error

    c_AAs     = 2.99792458e18                          # Speed of light in Angstrom/s
    hp=6.626069311e-27                                 # constant of planck in ergs/s
    pi=3.14159265359	
    area_obj=pi*(object_size*object_size)
    # calculating number of photons to estimate the error, not used in the moment
    #wavelen_sky,flux_sky=fc.open_ascii_cat(sky_flux, unpack=True, vartype=type(1.0))
    #wavelen_sky=np.array(wavelen_sky).astype(float)*10.0 # convert from nn to Angstrom
    #flux_sky=np.array(flux_sky).astype(float)*0.1 # convert from nn-1 to Angstrom-1
    #flux_sky=flux_sky*area_obj # consider the size of imaged object
    #nu=c_AAs/np.array(wavelen).astype(float)# convert wavelenght to frequency 
    
    
    #flux_phot=np.array(flux).astype(float)*10e4 #convert cm^2 to m^2
    print ('Flux')
    print (flux)
    #flux_phot=np.array(flux_phot).astype(float)/(hp*nu)

    #print type(flux_sky)
    #print type(telescope_area)
    #print type(exposure)
    #print type(flux_phot)
    #flux_sky=flux_sky*telescope_area*exposure*thoughput
    #flux_phot=flux_phot*telescope_area*exposure*thoughput


    #flux_sky_interp  = np.interp(wavelen,wavelen_sky,flux_sky,left=0.0,right=0.0)
        
    #flux_total_phot_err=flux_phot#+flux_sky_interp
    #flux_total_err=np.random.poisson(flux_total)

    #flux_total= flux_total*(hp*nu)
    #flux_total_err= (np.sqrt(flux_total_phot_err)*(hp*nu)*10e-4)/(telescope_area*exposure)

    mask_f=flux>0
    flux_total_err=max(flux[mask_f])-min(flux[mask_f])*noise
    flux_error_added=np.random.normal(flux,(flux-min(flux[mask_f]))*noise)#flux
    #mags_on_filters[k]=random.gauss(mags_on_filters[k], errors_on_filters[k])
    print ('Flux w error')
    print (flux_error_added)
    print ('Flux')
    print (flux)
    #print (hp*nu)

    #integrandobj=flux_phot*filt_int*wavelen
    #mask=integrandobj>0.0
    #Iobject  = simps(flux_phot[mask]*filt_int[mask],wavelen[mask])#*wavelen[mask]
    #print "Iobj antes"
        
    #print Iobject
    #print "fim"
    #integrandback=sky_int*filt_int*wavelen


    #if filter_source=='bpz' or filter_source=='BPZ' or type(ref_filter)==type(None):
    #    ref_filter=filter_file
    mags=[]
    mags_err=[]#flux
    mags_noise=[]
    for i in range(0,len(filter_ids)):
        if filter_source=='eazy':
            wavelenF,filt=get_filter_eazy(filter_file,int(filter_ids[i]))#
        elif filter_source=='ctio-decam':
            wavelenF,filt=get_filter_cols(filter_file,int(filter_ids[i]))
        elif filter_source=='ctio-decam-noatm':
            wavelenF,filt=get_filter_cols(filter_file,int(filter_ids[i]),noatm=True)
        else:
            wavelenF,filt=fc.open_ascii_cat(filter_ids[0], unpack=True, vartype=type(1.0))

        if filter_unit=='nm':
            wavelenF=wavelenF*10.0
    
        wavelen=np.array(wavelen).astype(float)
        wavelenF=np.array(wavelenF).astype(float)
        filt=np.array(filt).astype(float)
        flux=np.array(flux).astype(float)#
        filt_int  = np.interp(wavelen,wavelenF,filt,left=0.0,right=0.0)              #Interpolate to common wavelength axis
        integrand=flux*filt_int*wavelen
        mask=integrand>0.0
        #print('Non zero elements')
        #print(np.count_nonzero(mask))
        #print (flux)
        if np.count_nonzero(mask)!=0:
            I1        = simps(flux[mask]*filt_int[mask]*wavelen[mask],wavelen[mask])
            mask2=filt_int>0.0
            I2        = simps(  filt_int[mask2]/wavelen[mask2],wavelen[mask2])                     #Numerator
            fnu       = (I1/I2) / c_AAs                          #Average flux density https://astronomy.stackexchange.com/questions/16286/how-can-i-convolve-a-template-spectrum-with-a-photometric-filter-response-spectr
            if (fnu)>0:
                mag_cal=-2.5*log10(fnu)-48.6
                mag_cal_err=random.gauss(mag_cal,noise)
                
            else:
                mag_cal=-99
                mag_cal_err=-99

        else:
            mag_cal=-99
            mag_cal_err=-99
        mags.append(mag_cal)
        mags_err.append(mag_cal_err)#flux
        mags_noise.append(noise)
    #mag_cal=0 wavelen_err,flux_total_err,flux_total,mags,mags_err
    return wavelen,flux_error_added,flux,mags,mags_err,mags_noise 

 

def template2magerr(template,filter_file,skyflux,telescope_area,filter_ids=None,object_size=3,exposure=100, plot=True,plot_dir='', cut_sn=5, filter_source='eazy',filter_unit='A'):
    
    c_AAs     = 2.99792458e18                          # Speed of light in Angstrom/s
    hp=6.626069311e-27                                 # constant of planck in ergs/s
    pi=3.14159265359	
    area_obj=pi*(object_size*object_size)
    wavelen,flux=fc.open_ascii_cat(template, unpack=True, vartype=type(1.0))
    wavelen_sky,flux_sky=fc.open_ascii_cat(skyflux, unpack=True, vartype=type(1.0))
    #print wavelen_sky[:10]
    wavelen_sky=np.array(wavelen_sky).astype(float)*10.0 # convert from nn to Angstrom
    #print wavelen_sky[:10]
    flux_sky=np.array(flux_sky).astype(float)*0.1 # convert from nn-1 to Angstrom-1
    flux_sky=flux_sky*area_obj # consider the size of imaged object
    
    nu=c_AAs/np.array(wavelen).astype(float)
    flux=np.array(flux).astype(float)*10e4 #convert cm^2 to m^2
    flux_phot=np.array(flux).astype(float)/(hp*nu)

    #print type(flux_sky)
    #print type(telescope_area)
    #print type(exposure)
    #print type(flux_phot)
    flux_sky=flux_sky*telescope_area*exposure
    flux_phot=flux_phot*telescope_area*exposure
    

    mag_errs=[]

    if filter_source=='bpz' or filter_source=='BPZ' or type(filter_ids)==type(None):
        filter_ids=filter_file

    for i in range(0,len(filter_ids)):
        if filter_source=='eazy':
            wavelenF,filt=get_filter_eazy(filter_file,int(filter_ids[i]))#
        elif filter_source=='ctio-decam':
            wavelenF,filt=get_filter_cols(filter_file,int(filter_ids[i]))
        elif filter_source=='ctio-decam-noatm':
            wavelenF,filt=get_filter_cols(filter_file,int(filter_ids[i]),noatm=True)
        else:
            wavelenF,filt=fc.open_ascii_cat(filter_ids[i], unpack=True, vartype=type(1.0))

        if filter_unit=='nm':
            wavelenF=wavelenF*10.0 
   
        #wavelenF,filt=get_filter_eazy(filter_file,int(filter_ids[i]))#fc.open_ascii_cat(filt_curve[i]+".res", unpack=True, vartype=type('float'))
        wavelen=np.array(wavelen).astype(float)
        wavelenF=np.array(wavelenF).astype(float)
        filt=np.array(filt).astype(float)
        flux=np.array(flux).astype(float)
        filt_int  = np.interp(wavelen,wavelenF,filt,left=0.0,right=0.0)              #Interpolate to common wavelength axis
        sky_int  = np.interp(wavelen,wavelen_sky,flux_sky,left=0.0,right=0.0)
        
        integrandobj=flux_phot*filt_int*wavelen
        mask=integrandobj>0.0
        Iobject  = simps(flux_phot[mask]*filt_int[mask],wavelen[mask])#*wavelen[mask]
        #print "Iobj antes"
        
        #print Iobject
        #print "fim"
        integrandback=sky_int*filt_int*wavelen
        mask2=integrandback>0.0
        #print integrandback
        Ibackground= simps(sky_int[mask2]*filt_int[mask2],wavelen[mask2]) #*wavelen[mask2]
        #print Ibackground
        #print "fim2"
        back_counts=Iobject+Ibackground
        #if back_counts==0:
        #    print "back counts are 0 in template "+template+" \n filter:"+str(filter_ids[i])
        sn_counts=Iobject/sqrt(back_counts)
        mag_err=1.0875/sn_counts
        print("Error Analysis")
        print(Iobject)
        print(Ibackground)
        print("SN:",str(Iobject),"/","sqrt("+str(sqrt(back_counts)),")=",str(sn_counts))
        print("magerr=",str(mag_err))

        #sn_counts=Iobject/sqrt(back_counts)
        
        if math.isnan(mag_err):
            print("NAN magerr")
            print("sn_counts, Iobj, back")
            print(sn_counts)
            print(Iobject)
            print(back_counts)
            print(flux_phot*filt_int*wavelen)
            mag_err=-99
        if sn_counts <= cut_sn:
            #print "SN_cut"
            #print sn_counts
            mag_err=-99
        mag_errs.append(mag_err)
        #print "numerator"
        #print I1
        #I2        = simps(  filt_int/wavelen,wavelen)                     #Numerator
        #print "Denominador"
        #print I2

        #fnu       = (I1/I2) / c_AAs                          #Average flux density
        #mag_cal=-2.5*log10(fnu) - 48.6
        #if sigma_err>0.0:
        #    mag_cal=random.gauss(mag_cal, sigma_err)
        #if math.isnan(mag_cal):
        #    mag_cal=-99

        #mAB.append(mag_cal)            #AB magnitude


    return np.array(mag_errs).astype(float)





def template2magnitudes_naive(template,filt, sigma_err=0.002, plot=True,plot_dir=''):

    wavelen,flux=fc.open_ascii_cat(template, unpack=True, vartype='float')
    wavelen=np.array(wavelen).astype('float')
    flux=np.array(flux).astype('float')
    nu,fluxnu=fluxlambda2fluxnu(wavelen.copy(),flux.copy())
    mags=fluxnu2mag(fluxnu.copy()) # filter_unit='nm'
    
    
    mags=np.array(mags).astype('float')


    
   
    #print mags[:50]
    #print "this is the first 50 mags"
    #filt_center=filt[0]
    #print filt
    #print "this is filt"
    mags_filt=[]
    mags_filt_err=[]
    nearwavs=[]
    flag_count=0
    for i in range(0,len(filt)):
        #print "working on filter "+str(i)
        mag_filt,nearwav,flag=find_nearestxy(wavelen,mags,filt[i])
        #dist_near=(abs(nearwav-filt[i]))/abs(wavelen[0]-wavelen[1])
        #print " ##### this is mag"
        #print mag_filt
        
        mag_filt=random.gauss(mag_filt, sigma_err)
        print(" ##### this is mag after and the template")
        print(mag_filt)
        print(template)
        #print "this is distnear"
        #print str(dist_near)+'    '+str(nearwav)+'    '+str(mag_filt)+'    '+str(abs(wavelen[0]-wavelen[1]))+'    '+str(abs(wavelen[-1]-wavelen[-2]))
        if (mag_filt > 29.5) or (mag_filt < 5) or flag > 0:
           print("eliminei")
           print(flag)
           if flag==2:
               flag_count+=1
           mag_filt=99.0
        #print "the nearest is "+str(mag_filt)
        mags_filt.append(round(mag_filt,3))
        mags_filt_err.append(sigma_err*3)
        nearwavs.append(nearwav)
    

    if plot==True:
        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.plot(wavelen,mags)
    
        #print mags_filt[:20]
        print("Antes de filtrar")
        print(str(len(nearwavs)),'    ',str(len(mags_filt)),'    ',str(len(mags_filt_err)))
        #print nearwavs
        #print mags_filt
        #print mags_filt_err
        nearwavs_plot,mags_filt_plot,mags_filt_err_plot=clear_filter(nearwavs,mags_filt,mags_filt_err)
        print("depois de filtrar")
        print(str(len(nearwavs_plot)),'    ',str(len(mags_filt_plot)),'    ',str(len(mags_filt_err_plot)))
        print("elliminated due to out of range")
        print(flag_count)
        #print nearwavs_plot
        #print mags_filt_plot
        #print mags_filt_err_plot
        ax1.errorbar(nearwavs_plot,mags_filt_plot,yerr=mags_filt_err_plot,c='g', fmt='o')
        plt.gca().invert_yaxis()
        plot_name=template.split("/")[-1]+".png"
        plt.savefig(plot_dir+plot_name)


    return mags_filt,mags_filt_err,nearwavs

def build_eazy_in(file_name, flux_objs,flux_objs_err,filter_names,filterids, **kwargs):
    """
    WARNING: mags_objs is a list of mags_filt

    """
    file_name=file_name.rstrip('.cat')
    cat = open(file_name+'.cat', 'w')
    #col = open(file_name+'.columns', 'w')
    zs=verify_kwarg("zs",[],kwargs)
    zperr=verify_kwarg("zperr",np.ones(len(filter_names))*0.01,kwargs)
    m0=verify_kwarg("m0",None,kwargs)
    check=verify_kwarg("check",False,kwargs)
    check_dir=verify_kwarg("check_dir",'',kwargs)
    nearwavs=verify_kwarg("check_wavelen",[],kwargs)
    #print nearwavs
    #print type(nearwavs)
    zpoff=verify_kwarg("zpoff",np.zeros(len(filter_names)),kwargs)
    _ids=verify_kwarg("_ids",range(1,len(flux_objs)+1),kwargs)
    
    header='# id   '
    for i in range(0,len(filterids)):
        header+="F"+str(filterids[i])+"    E"+str(filterids[i])
        if i!=(len(filterids)-1):        
            header+="    "
        if i==(len(filterids)-1):
            if len(zs)>0:
                header+="    z_spec\n"
            else:
                header+="\n"
    cat.write(header)
    for j in range(0,len(flux_objs)):
        #print np.array(mags_objs).shape
        #print _ids
        cat.write(str(_ids[j])+'    ')
        for k in range(0,len(flux_objs[0])):
            #print(j,"    ids ",_ids[j], " ", k, "    mags ",str(flux_objs[j][k]))
            #print("    zspec", zs[j])
            #print "linha    "+str(k)
            cat.write(str(flux_objs[j][k])+"    "+str(flux_objs_err[j][k]))
            if k!=(len(flux_objs[0])-1):
                cat.write("    ")
            if k ==(len(flux_objs[0])-1) and len(zs)>0:
                cat.write("    "+str(zs[j]))   
            if k ==(len(flux_objs[0])-1) and j!=(len(flux_objs)-1):
                cat.write("\n")
    cat.close()
    """
    if (check==True) and (len(nearwavs)!=0):
        _data=fc.open_ascii_cat(file_name+'.cat', unpack=False, vartype='float')
        #print _data.shape

        try:
            temp_id=_data[0][0]
        except:
            _data=np.array([_data])
        #print "new data shape"
        #print _data.shape
   
        for l in range(0,len(_data)):
            mags_filt=[]
            mags_filt_err=[]
            temp_ids=[]
            temp_id=_data[l][0]
            #print "inside L"
            #print l
            for m in np.arange(1,len(_data[0])-1,2):
                #print "inside m"
                #print m
                #print _data[l][m]
                mags_filt.append(_data[l][m])
                mags_filt_err.append(_data[l][m+1])
                
            #if l==0:
            #    print mags_filt[:20]
            fig = plt.figure()

            ax1 = fig.add_subplot(111)
            #ax1.plot(wavelen,mags)
            #print np.array(nearwavs).size
            #print np.array(mags_filt_err)
            #print "Antes de filtrar II vez"
            #print nearwavs
            #print mags_filt
            #print mags_filt_err
            #print "lenghts I"
            #print str(len(nearwavs))+'    '+str(len(mags_filt))+'    '+str(len(mags_filt_err))
            #print "before the buildbpz cleaning"
            #print nearwavs.shape
            #print np.array(mags_filt).shape
            #print np.array(mags_filt_err).shape
            nearwavs_plot,mags_filt,mags_filt_err=clear_filter(nearwavs,mags_filt,mags_filt_err)
            #print "Depois de filtrar II vez"
            #print "lenghts II"
            #print str(len(nearwavs_plot))+'    '+str(len(mags_filt))+'    '+str(len(mags_filt_err))
            #nearwavs,mags_filt,mags_filt_err=clear_filter(nearwavs,mags_filt,mags_filt_err)
            #print nearwavs
            #print mags_filt
            #print mags_filt_err

            ax1.errorbar(np.array(nearwavs_plot).astype('float'),np.array(mags_filt).astype('float'),yerr=np.array(mags_filt_err).astype('float'),c='g', fmt='o')
            plt.gca().invert_yaxis()
            plt.savefig(check_dir+str(temp_id)+".png")
    """
    return 0











