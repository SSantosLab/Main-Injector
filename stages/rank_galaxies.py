from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp 
import numpy as np
from scipy.stats import norm
from scipy.special import gammaincinv
from scipy.special import gammaincc
import cosmolopy.magnitudes as mag
import pandas as pd
import os


def find_galaxy_list(map_path: str,
                     event_name: str,
                     galax : str,
                     out_path: str,
                     credzone: float = 0.99,
                     nsigmas_in_d: int = 3,
):
        
    #loading the map:
    skymap = hp.read_map(map_path, field=None)
    prob, distmu, distsigma, distnorm = skymap
   

    #map parameters:
    npix = len(prob)
    nside = hp.npix2nside(npix)

    #galaxy parameters(RA, DEC to theta, phi):
    galax = np.load(galax,allow_pickle=True)
    galax = (galax[np.where(galax[:,3]>0),:])[0] #no distance<0
    theta = 0.5 * np.pi - np.pi*(galax[:,2])/180
    phi = np.deg2rad(list(galax[:,1]))
    d = np.array(galax[:,3])

    #finding given percent probability zone(default is 99%):
    probcutoff = 1
    probsum = 0
    npix99 = 0

    sortedprob = np.sort(prob)
    while probsum<credzone:
        probsum = probsum+sortedprob[-1]
        probcutoff = sortedprob[-1]
        sortedprob = sortedprob[:-1]
        npix99 = npix99+1

    area = npix99 * hp.nside2pixarea(nside, degrees=True)

    #converting galaxy coordinates to map pixels:
    ipix = hp.ang2pix(nside, list(theta), list(phi))


    #calculating probability for galaxies by the localization map:
    p = prob[ipix]
    distp = (norm(distmu[ipix], distsigma[ipix]).pdf(d) * distnorm[ipix])

    #cuttoffs- 99% of probability by angles and 3sigma by distance:
    inddistance = np.where(np.abs(d-distmu[ipix])<nsigmas_in_d*distsigma[ipix])
    indcredzone = np.where(p>=probcutoff)

    #if no galaxies
    if (galax[np.intersect1d(indcredzone,inddistance)]).size == 0:
        while probsum < 0.99995:
            if sortedprob.size == 0:
                break
            probsum = probsum + sortedprob[-1]
            probcutoff = sortedprob[-1]
            sortedprob = sortedprob[:-1]
            npix99 = npix99 + 1
        inddistance = np.where(np.abs(d - distmu[ipix]) < 5 * distsigma[ipix])
        indcredzone = np.where(p >= probcutoff)

    #apply the cuts
    ipix = ipix[np.intersect1d(indcredzone, inddistance)]
    p = p[np.intersect1d(indcredzone, inddistance)]
    p = (p * (distp[np.intersect1d(indcredzone, inddistance)]))  ##d**2?


    galax = galax[np.intersect1d(indcredzone, inddistance)]
    if galax.size == 0:
        print ("no galaxies in field")
        print ("99.995% of probability is ", npix99*hp.nside2pixarea(nside, degrees=True), "deg^2")
        return


    # normalized luminosity to account for mass:
    luminosity = mag.L_nu_from_magAB(galax[:, 4] - 5 * np.log10(galax[:, 3] * (10 ** 5)))
    luminosityNorm = luminosity / np.sum(luminosity)

    
    #Create a df of the galaxies that I can return
    df = pd.DataFrame(data = galax, columns = ['ID', 'RA', 'Dec', 'Dist', 'imag'])
    score = list(p*luminosityNorm)
    scoresum = np.sum(score)
    df['Normalized Score'] = [i/scoresum for i in score]
    newdf = df.sort_values(by='Normalized Score', ascending=False).reset_index(drop=True)
    percent_loc = [(i / len(newdf)) * 100 for i in list(newdf.index)]
    newdf['Percent'] = percent_loc
    newdf.to_csv(out_path, index=False)
    
    return