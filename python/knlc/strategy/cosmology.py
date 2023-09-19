import matplotlib
matplotlib.use('agg')
import numpy 
import matplotlib.pyplot as plt
from . import cosmoconst as cosmo
from . import constants as const
from scipy import integrate
from math import sqrt
import random
import itertools as it
from .misc import verify_kwarg, quadrance
"""
 Just a toy model to exemplify the mcmc. Its double einstein ring system, with only \omega_r and \omega_m as free parameters in a point-like lens

"""

#Cosmology

def ang_dis(z_init, z_final, omega0_m, omega0_k, omega0_l,w=-1, unit='km', model='wcdm',aux1=0.0):
    """
    Compute the angular diameter distance

    Input:
     - z_init  float : initial red-shift
     - z_final float : final red-shift
    Output:
     - angular diameter distance in megaparsec
    """
    #megaparsec = 3.08568025E19 #[kilo-meters]
    #H0 = 100.0 * self.__h / 3.08568025E19
    #c_light = 299792.458/megaparsec
    #return c_light/H0/(1.0 + z_final)*self.int_flat(z_init, z_final)
    return (1.0/(1.0 + z_final))*com_dist(z_init, z_final, omega0_m, omega0_k, omega0_l,w=w, unit=unit, model=model,aux1=aux1) 

def com_dist(z_init, z_final, omega0_m, omega0_k, omega0_l,w=-1,unit='km', model='wcdm',aux1=0.0, omega0_g=False):
    """
    Compute the comoving distance

    Input:
     - z_init  float : initial red-shift
     - z_final float : final red-shift
    Output:
     - FIX ME

    args=


    """
    c1 = const.c/cosmo.H0_sec
    if unit=='m' or unit =='meters':
       c1=c1*10E3
    if unit=='cm' or unit =='centimeters':
       c1=c1*10E5

    return c1*int_fflat(z_init, z_final, omega0_m, omega0_k, omega0_l, w=w, model=model,aux1=aux1, omega0_g=omega0_g) 




def lum_dist(z_init, z_final, omega0_m, omega0_k, omega0_l,w=-1,unit='km', model='wcdm',**kwargs):
    """
    Compute the luminosity distance

    Input:
     - z_init  float : initial red-shift
     - z_final float : final red-shift
    Output:
     - FIX ME
    """
    one_plus_z_final=(1.0 + z_final)
    return (one_plus_z_final)*com_dist(z_init, z_final, omega0_m, omega0_k, omega0_l,w=w,unit=unit, model=model,**kwargs)

def int_fflat(z_init, z_final, omega0_m, omega0_k, omega0_l, w=-1,model='wcdm', aux1=0.0,omega0_g=False):
    """
    Compute the eq (3.4) in my msc thesis

    Input:
     - z_init  float : initial red-shift
     - z_final float : final red-shift

    Output:
     - 'I'
    """
    
    out = integrate.quad(fadm_hubble_param_inv, z_init, z_final, \
                args=(omega0_m, omega0_k, omega0_l,w,model,aux1,omega0_g))
    #print out
    if omega0_k <= 10E-8:
        return out[0]
    else:
        print('ERROR in int_flat: you showld use a flat cosmology')
        return 0.0   

def fadm_hubble_param_inv(z_in, omega0_m, omega0_k, omega0_l=0.7, w=-1.0, model='wcdm', aux1=0.0,omega0_g=False,**kwargs):
    """
    omega_g=2.469E-5
    """
    if omega0_g ==True:# kwargs.keys():
        h=cosmo.h # hubble reduced const
        N_eff=3.04 # Komatsu et al. 2011
        omega0_g=((2.469E-5)/(h*h))*(1+0.2271*N_eff) 
    else:
        omega0_g=0.0
    one_plus_z = 1.0 + z_in

    if model=='lambdacdm':
        out = sqrt(omega0_m*one_plus_z**3.0 + omega0_k*one_plus_z**2.0 + \
               omega0_l )
    elif model=='CLP':
        wa=aux1
        omega0_l=1-omega0_m - omega0_g
        if type(z_in)!=type(1.0):
           z_in=np.array(z_in).astype(float)
           one_plus_z=np.array(one_plus_z).astype(float) 
        arg_exp=(-3*wa*z_in)/one_plus_z
        f_z=one_plus_z**(3*(1+w+wa))*np.exp(arg_exp)
        out = sqrt(omega0_g*one_plus_z**4.0 + omega0_m*one_plus_z**3.0 + omega0_l*f_z)
    elif model=='IDE':
        omega0_l=1-omega0_m - omega0_g
        delta=aux1
        fz=one_plus_z**(3*(1+w))
        f2z=((delta*fz)+(3*w*one_plus_z**(3-delta))) /(delta+3*w)
        out = sqrt(omega0_g*one_plus_z**4.0 + omega0_m*f2z + omega0_l*f_z)
    else:
        omega0_l=1.0-abs(omega0_m - omega0_g)
        #print "cosmo params"
        #print omega0_l
        #print omega0_m
        #print omega0_g
        #print one_plus_z
        #print w
        #print "fim cosmo"
        out = sqrt(omega0_g*one_plus_z**4.0 + omega0_m*one_plus_z**3.0 + omega0_l*one_plus_z**(3*(1+w)))
       
 
    return 1.0/out


 


def sim_einstein_ring(n=500,zs_interval=[1.5,6.5],zl_interval=[0.1,1.0],omega0_m=0.3, omega0_l=0.7,error=0.1):
    """ Einstein Ring for a point mass lens"""
    result=[]
    zs=[]
    mass=[]
    zl=[]
    omega0_k=0
    for i in range(0,n):
        zl_aux=round(random.uniform(zl_interval[0],zl_interval[1]),2)
        zs_aux=round(random.uniform(zs_interval[0],zs_interval[1]),2)
        mass_exp= int(random.randint(13, 14))
        mass_aux=10**mass_exp
        mass_aux=mass_aux*const.solar_mass
        c_ms=const.c*10E3 # convert km/sec to m/sec
        ER=ang_dis(zl_aux, zs_aux,omega0_m,omega0_k,omega0_l) *4*const.g_newton*mass_aux

        ER=ER/(c_ms**2)
        ER=ER/(ang_dis(0, zs_aux,omega0_m,omega0_k,omega0_l)*ang_dis(0, zl_aux,omega0_m,omega0_k,omega0_l))
        ER=sqrt(ER)
        if (ER*3600) > 1.0:
            ER=ER+(random.uniform(error*-1,error)/3600.0)
            result.append(ER)
            zs.append(zs_aux)
            zl.append(zl_aux)
            mass.append(mass_aux)         


    return  numpy.array(result),numpy.array(zs),numpy.array(zl),numpy.array(mass)



def dist_ratio(cosmo_params,zl,zs,model='wcdm'):
    
    #print cosmo_params
    if (len(cosmo_params)==3) and (model=='wcdm' or model=='lambdacdm'):
        w,omega_l,omega_m =cosmo_params
    else:
        w,omega_l,omega_m,aux1,aux2 =cosmo_params
    zs1=zs[0]
    zs2=zs[1]
    omega_k=0
    if model=='IDE':
        dist_ratio=ang_dis(zl, zs1, omega_m,omega_k,omega_l, w, model=model, aux1=aux1)*ang_dis(0, zs2, omega_m,omega_k,omega_l, w, model=model, aux1=aux1)
        dist_ratio=dist_ratio/(ang_dis(0, zs1, omega_m,omega_k,omega_l, w,model=model, aux1=aux1)*ang_dis(zl, zs2, omega_m,omega_k,omega_l, w, model=model, aux1=aux1))
    elif model=='CLP':
        dist_ratio=ang_dis(zl, zs1, omega_m,omega_k,omega_l, w, model=model, aux1=aux1)*ang_dis(0, zs2, omega_m,omega_k,omega_l, w, model=model, aux1=aux1)
        dist_ratio=dist_ratio/(ang_dis(0, zs1, omega_m,omega_k,omega_l, w,model=model, aux1=aux1)*ang_dis(zl, zs2, omega_m,omega_k,omega_l, w, model=model, aux1=aux1))
    else:
        dist_ratio=ang_dis(zl, zs1, omega_m,omega_k,omega_l, w, model=model)*ang_dis(0, zs2, omega_m,omega_k,omega_l, w, model=model)
        dist_ratio=dist_ratio/(ang_dis(0, zs1, omega_m,omega_k,omega_l, w,model=model)*ang_dis(zl, zs2, omega_m,omega_k,omega_l, w, model=model))
    return dist_ratio


def dist_ratio_arr(zl,zs,w,omega_l,omega_m, **kwargs):
    fix_zs1=verify_kwarg("fix_zs1",-1,kwargs)
    model=verify_kwarg("model",'wcdm',kwargs)
    if fix_zs1==-1:
        sources_comb=numpy.array(list(it.combinations(zs,2)))
    else:
        zs=zs[numpy.where(zs!=fix_zs1)]
        sources_aux=[zs,[fix_zs1]]
        sources_comb=numpy.array(list(it.product(*sources_aux)))
    
    dist_ratio_mean=[]
    dist_ratio_std=[]
    dist_ratio_median=[]

    for i in xrange(0,len(sources_comb)):
        dist_ratio_aux=[w,omega_l,omega_m]#[sources_comb[i],[zl],w,omega_l,omega_m]
        dist_ratio_out=numpy.apply_along_axis(dist_ratio,0,dist_ratio_aux, zl=zl,zs=sources_comb[i])
        del dist_ratio_aux
        dist_ratio_mean.append(numpy.mean(dist_ratio_out))
        dist_ratio_std.append(numpy.std(dist_ratio_out))
        dist_ratio_median.append(numpy.median(dist_ratio_out))
        del dist_ratio_out
            
    return dist_ratio_mean, dist_ratio_std, dist_ratio_median, sources_comb
    

    
def dist_ratio_curve(zl,zs, w,omega_l,omega_m, **kwargs):
    aux1=verify_kwarg("aux1",0.0,kwargs)
    aux2=verify_kwarg("aux2",0.0,kwargs)

    zlim=verify_kwarg("zlims",[zl+0.1,5*zl],kwargs)
    curve_out=[]    

    z=numpy.arange(zlim[0],zlim[1],0.1)
    sources_comb=numpy.array(list(it.product(*[z,[zs]])))
    input_aux=[w,omega_l,omega_m,aux1,aux2]#[sources_comb,zl,w,omega_l,omega_m]
    #curve_in=numpy.array(list(it.product(*input_aux)))
    
    for i in xrange(0,len(sources_comb)):
        dist_ratio_out=dist_ratio(input_aux,zl=zl,zs=sources_comb[i])
        curve_out.append(dist_ratio_out)

    return z, numpy.array(curve_out)
###### CMB model fitting ##########################

def la_acustic_scale_cmb(z,omegaM,w,model='wcdm',aux1=0.0):
    la=const.pi*com_dist(0.0, z, omega0_m=omegaM, omega0_k=0.0, omega0_l=(1-omegaM),w=w, model=model,aux1=aux1, omega0_g=True)
    la=la/com_sound_horizon(z, omega0_m=omegaM, omega0_k=0.0, omega0_l=(1-omegaM),w=w, model='wcdm',h=0.7,aux1=aux1)
    return la


def R_shift_cmb(z,omegaM,w,model='wcdm',aux1=0.0,unit='km'):

    c1 = const.c/cosmo.H0_sec
    if unit=='m' or unit =='meters':
       c1=c1*10E3
    if unit=='cm' or unit =='centimeters':
       c1=c1*10E5

    R_cmb=(sqrt(omegaM)/c1)*com_dist(0.0, z, omega0_m=omegaM, omega0_k=0.0, omega0_l=(1-omegaM),w=w, model=model,aux1=aux1, omega0_g=True )
    return R_cmb

def z_dec(omegaM,h,ob):
    # small scale cosmological perturbation
    # an analytic approach Hu and Sugiyama 1996
    omegaB=ob#0.044939 #magana 2017
    obhh=ob*h*h#0.02228#omegaB*h*h
    ohh=omegaM*h*h
    g1=0.0783*(obhh**(-0.238))
    g1=g1/(1.0+39.5*(obhh**(0.763)))
    g2=0.560/(1.0+21.1*(obhh**(1.81)))
    z_out=1048.0*(1.0+0.00124*(obhh**(-0.738)))
    z_out=z_out*(1.0+g1*(ohh**(g2)))

    return z_out

##############BAO STUFF ################

def com_sound_horizon(z, omega0_m, omega0_k, omega0_l,w=-1,unit='km', model='wcdm',h=0.7,aux1=0.0):
    """
    Compute the comoving sound horizon distance

    Input:
     - z_init  float : initial red-shift
     - z_final float : final red-shift
    Output:
     - FIX ME

    args=


    """
    c1 = const.c/cosmo.H0_sec
    if unit=='m' or unit =='meters':
       c1=c1*10E3
    if unit=='cm' or unit =='centimeters':
       c1=c1*10E5

    return c1*rsound(z,omega0_m,w,aux1=aux1,h=0.7,obhh=0.02228)#int_soundflat(z, omega0_m, omega0_k, omega0_l, w=w, model=model,h=h,aux1=aux1)


def int_soundflat(z_final, omega0_m, omega0_k, omega0_l, w=-1,model='wcdm', h=0.7,aux1=0.0):
    """
    Compute the eq (3.4) in my msc thesis

    Input:
     - z_final float : final red-shift

    Output:
     - 'I'
    """
    
    out = integrate.quad(soundflat, z_final, 10E30, \
                args=(omega0_m, omega0_k, omega0_l,w,model,h,aux1))
    #print out
    if omega0_k <= 10E-8:
        return out[0]
    else:
        print('ERROR in int_flat: you showld use a flat cosmology')
        return 0.0   


def soundflat(z_in, omega0_m, omega0_k=0.0, omega0_l=0.7, w=-1, model='wcdm',h=0.7, aux1=0.0):
    tcmb=2.72548 # doi:10.1088/0004-637X/707/2/916
    obhh=0.02228#omegaB*h*h
    R=31500.0*obhh*((tcmb/2.7)**(-4))
    cz=1.0/sqrt(3*(1+(R/(1+z_in))))
    outsf=cz*fadm_hubble_param_inv(z_in, omega0_m=omega0_m, omega0_k=omega0_k, omega0_l=omega0_l, w=w, model=model, aux1=aux1, omega0_g=True)
    return outsf
#sound speed

def cs(a,h,ob):
    f1=(3.0*ob)/(4.0*2.469e-5*pow(h,-2.0))
    f2=3.0*(1.0+f1*a)
    return 1.0/(sqrt(f2))

#sound horizon
def rsound(z,omega0_m,w,aux1=0.0,h=0.7,obhh=0.02228,model='wcdm'):
    ob=(obhh)/(h*h) 
    li=1.0/(1.0+z)
    omega0_k=0.0
    #fadm_hubble_param_inv(z_in, omega0_m=omega0_m, omega0_k=omega0_k, omega0_l=omega0_l, w=w, model=model, aux1=aux1, omega0_g=True)
    
    return integrate.quad(lambda x:cs(x,h,ob)*(1/(x*x))*fadm_hubble_param_inv((1.0/x)-1.0,omega0_m=omega0_m,omega0_k=omega0_k, omega0_l=0.7, w=w, model=model, aux1=aux1, omega0_g=True),0,li,)[0]

####################### end of BAO STUFF

 
"""
def einstein_ring(zs=zs,zl=zl,mass=M,omega0_m=omega0_m_prior, omega0_l=omega0_l_prior):
    
    result=[]
    for i in range(0,len(zs)):
        omega0_k=0
        c_ms=const.c*10E3 # convert km/sec to m/sec
        ER=ang_dis(zl[i], zs[i],omega0_m,omega0_k,omega0_l) *4*const.g_newton*mass[i]
        ER=ER/(c_ms**2)
        ER=ER/(ang_dis(0, zs[i],omega0_m,omega0_k,omega0_l)*ang_dis(0, zl[i],omega0_m,omega0_k,omega0_l))
        result.append(sqrt(ER))
    return  numpy.array(result)


# (We define the error as 2 std)
y = pymc.Normal('y', mu=einstein_ring, tau=1. / (edata/2) ** 2, observed=True, value=ydata)
 
# package the full model in a dictionary
model1 = dict(omega0_m=omega0_m_prior, omega0_l=omega0_l_prior, y_model=einstein_ring, y=y)
 
# run the basic MCMC:
S = pymc.MCMC(model1)
S.sample(iter=1000000, burn=500000)
 
# extract the traces and plot the results
pymc_trace_unifo = [S.trace('omega0_m')[:],
              S.trace('omega0_l')[:]]
 
plot_MCMC_results(zl, ydata, pymc_trace_unifo)
print()
print("omega0_m mean= {:.2f}".format(pymc_trace_unifo[0].mean()))
print("omega0_l  mean= {:.2f}".format(pymc_trace_unifo[1].mean()))


"""
