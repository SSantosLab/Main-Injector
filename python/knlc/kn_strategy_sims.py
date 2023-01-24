# A module to print brightness estimates

import os
import sys
import math
import glob
import os.path

from math import log10
from argparse import ArgumentParser
from os.path import join, basename, dirname
import healpy as hp
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib

import astropy.units as u
import matplotlib.pyplot as plt

from astropy.cosmology import z_at_value
from scipy.stats import uniform
from scipy.stats import norm
from scipy.interpolate import interp1d
from astropy.cosmology import WMAP9 as cosmo

matplotlib.use("Agg")

np.set_printoptions(threshold=sys.maxsize)


def verify_kwarg(param_name, default_value, kwargs):
    if param_name in kwargs.keys():
        param = kwargs[param_name]
    else:
        param = default_value
    return param


def get_header_ascii(file_name, identifier='#'):

    f = open(file_name, 'r')
    lines = f.readlines()
    header = []
    for i in range(0, len(lines)):
        lines_aux = lines[i].lstrip(identifier)
        if lines_aux != lines[i]:
            lines_hdr = lines_aux.rstrip('\n')
            lines_hdr = lines_hdr.rstrip(' ')
            header.append(lines_hdr)
    return header


def open_ascii_cat(file_name, **kwargs):
    comments = verify_kwarg("comments", "#", kwargs)
    skiprows = verify_kwarg("skip_rows", 0, kwargs)
    sk_last = verify_kwarg("skip_footer", 0, kwargs)
    usecols = verify_kwarg("usecols", None, kwargs)
    unpack = verify_kwarg("unpack", False, kwargs)
    vartype = verify_kwarg("vartype", type('str'), kwargs)
    if 'delimiter' in kwargs.keys():
        data = np.loadtxt(file_name, dtype=vartype, comments=comments,
                          delimiter=kwargs['delimiter'], skiprows=skiprows,
                          usecols=usecols, unpack=unpack)
    else:
        try:
            data = np.loadtxt(file_name, dtype=vartype, comments=comments,
                              skiprows=skiprows, usecols=usecols, unpack=unpack)
        except:
            header = get_header_ascii(file_name, identifier=comments)
            m_v = 0
            f_v = 0
            cols = range(0, len(header))
            print("WARNING: table format is wrong.")
            data = np.genfromtxt(file_name, dtype=vartype, comments=comments,
                                 skip_header=skiprows, skip_footer=sk_last,
                                 missing_values=m_v, filling_values=f_v,
                                 usecols=cols, unpack=unpack)
    return data


tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]


markers20 = ["*", "D", ".", "v", "^", "<", ">", "1",
             "2", "3", "4", "8", "s", "p", "P", ",", "o"]


for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


def plt_style():
    plt.rcParams.update({
                        'lines.linewidth': 1.0,
                        'lines.linestyle': '-',
                        'lines.color': 'black',
                        'font.family': 'serif',
                        'font.weight': 'normal',
                        'font.size': 13.0,
                        'text.color': 'black',
                        'text.usetex': False,
                        'axes.edgecolor': 'black',
                        'axes.linewidth': 1.0,
                        'axes.grid': False,
                        'axes.titlesize': 'x-large',
                        'axes.labelsize': 'x-large',
                        'axes.labelweight': 'normal',
                        'axes.labelcolor': 'black',
                        'axes.formatter.limits': [-4, 4],
                        'xtick.major.size': 7,
                        'xtick.minor.size': 4,
                        'xtick.major.pad': 8,
                        'xtick.minor.pad': 8,
                        'xtick.labelsize': 'medium',
                        'xtick.minor.width': 1.0,
                        'xtick.major.width': 1.0,
                        'ytick.major.size': 7,
                        'ytick.minor.size': 4,
                        'ytick.major.pad': 8,
                        'ytick.minor.pad': 8,
                        'ytick.labelsize': 'medium',
                        'ytick.minor.width': 1.0,
                        'ytick.major.width': 1.0,
                        'legend.numpoints': 1,
                        'legend.shadow': False,
                        'legend.frameon': False})


# Handle command-line arguments


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)

    return math.sqrt(variance)


def get_distances(eventmap):

    pb, distmu, distsigma, distnorm = hp.read_map(eventmap,
                                                  field=range(4),
                                                  dtype=[np.float64,
                                                         np.float64,
                                                         np.float64,
                                                         np.float64])  # clecio
    return 0


class KNCalc(): 
    def __init__(self,
                 skymap_name,
                 distance,
                 distance_err,
                 time_delay,
                 loglan=None,
                 vk=None,
                 logmass=None,
                 info_file="knsed_info.txt",
                 kn_weight_type="gaussian",
                 plot=False,
                 use_map=None,
                 area_covered=0.9,
                 reduce_mapresolution=False,
                 set_mjd_correction=False,
                 m_exp_kncalc=False,
                 deep_coverage=0.5):

        self.event = skymap_name
        self.sw_mexp = False
        self.distance = distance
        self.distance_err = distance_err

        self.loglan = loglan
        self.vk = vk
        self.logmass = logmass

        self.set_mjd_correction = set_mjd_correction

        if float(time_delay) > 400.8:
            MSG = "Currently, only times delays less thant 400.8 hours " +\
                  "(16.7 days) post merger are supported."

            print(MSG)
            sys.exit()

        self.delta_mjd = round(float(time_delay) / 24.0, 1)
        delta_mjd_full = float(time_delay) / 24.0

        delta_mjd_later = round(self.delta_mjd+0.3, 1)

        if (delta_mjd_full != self.delta_mjd) and (self.set_mjd_correction):
            self.set_mjd_correction = True
            if delta_mjd_full > self.delta_mjd:
                delta_mjd_corr = self.delta_mjd+0.1
                if (delta_mjd_full-self.delta_mjd) >= 0.05:
                    delta_mjd_later = self.delta_mjd+0.4
            else:
                delta_mjd_corr = self.delta_mjd-0.1
            delta_mjd_corr = round(delta_mjd_corr, 1)
        else:
            self.set_mjd_correction = False

        # Set directory for lookup table
        knlc_dir = os.getenv("DESGW_DIR", "./")
        if knlc_dir != "./":
            knlc_dir = knlc_dir + "/knlc/"

        # Choose lookup table based on time_delay
        if delta_mjd_later < 2.3:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry.csv')
        elif delta_mjd_later < 4.7:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry_2.csv')
        elif delta_mjd_later < 7.1:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry_3.csv')
        elif delta_mjd_later < 9.5:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry_4.csv')
        elif delta_mjd_later < 11.9:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry_5.csv')
        elif delta_mjd_later < 14.3:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry_6.csv')
        elif delta_mjd_later < 16.7:
            df_later = pd.read_csv(knlc_dir+'data/grouped_photometry_7.csv')

        if self.delta_mjd < 2.3:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry.csv')
        elif self.delta_mjd < 4.7:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_2.csv')
        elif self.delta_mjd < 7.1:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_3.csv')
        elif self.delta_mjd < 9.5:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_4.csv')
        elif self.delta_mjd < 11.9:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_5.csv')
        elif self.delta_mjd < 14.3:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_6.csv')
        elif self.delta_mjd < 16.7:
            df = pd.read_csv(knlc_dir+'data/grouped_photometry_7.csv')

        if self.set_mjd_correction:
            if delta_mjd_corr < 2.3:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry.csv')
            elif delta_mjd_corr < 4.7:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry_2.csv')
            elif delta_mjd_corr < 7.1:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry_3.csv')
            elif delta_mjd_corr < 9.5:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry_4.csv')
            elif delta_mjd_corr < 11.9:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry_5.csv')
            elif delta_mjd_corr < 14.3:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry_6.csv')
            elif delta_mjd_corr < 16.7:
                df_corr = pd.read_csv(knlc_dir+'data/grouped_photometry_7.csv')

        df_later['ZMEAN'] = np.mean(df_later[['ZMIN', 'ZMAX']].values, axis=1)
        df['ZMEAN'] = np.mean(df[['ZMIN', 'ZMAX']].values, axis=1)

        mean_z = z_at_value(cosmo.luminosity_distance,
                            float(self.distance) * u.Mpc)
        template_df_mean = df[(df['ZMIN'].values < mean_z) &\
                              (df['ZMAX'].values > mean_z) &\
                              (df['DELTA_MJD'].values == self.delta_mjd)]
        template_df_mean = template_df_mean.copy().reset_index(drop=True)

        sed_filenames, kn_inds, vks, loglans, logmass_s = open_ascii_cat(
                                                                    info_file,
                                                                    unpack=True
                                                                )
        loglans = np.array(loglans).astype('float')
        kn_inds = np.array(kn_inds).astype('float')
        vks = np.array(vks).astype('float')
        logmass_s = np.array(logmass_s).astype('float')

        if template_df_mean.empty:
            checkzmax = df['ZMAX'].values > mean_z

            if not all(checkzmax):
                print(f"Object is too far away. mean_z is {mean_z}. Exiting.")
                self.Flag = 0
                return
            else:
                print("Something wrong with knlc photometry template library. Exiting")
                sys.exit()

        template_df_mean['WEIGHT'] = 1.0 / template_df_mean.shape[0]
        self.template_df_mean = template_df_mean

        # Full distance calculation
        template_df_full = df[df['DELTA_MJD'].values ==
                              self.delta_mjd].copy().reset_index(drop=True)
        template_df_later = df_later[df_later['DELTA_MJD'].values == delta_mjd_later].copy(
        ).reset_index(drop=True)

        if self.set_mjd_correction == True:

            template_df_corr = df_corr[df_corr['DELTA_MJD'].values == float(
                delta_mjd_corr)].copy().reset_index(drop=True)
            mag_g = template_df_full['MAG_g'].values
            mag_r = template_df_full['MAG_r'].values
            mag_i = template_df_full['MAG_i'].values
            mag_z = template_df_full['MAG_z'].values

            mag_g_corr = template_df_corr['MAG_g'].values
            mag_r_corr = template_df_corr['MAG_r'].values
            mag_i_corr = template_df_corr['MAG_i'].values
            mag_z_corr = template_df_corr['MAG_z'].values

            if delta_mjd_corr > self.delta_mjd:
                print("correcting magnitudes >")
                template_df_full['MAG_g'] = [
                    np.interp(delta_mjd_full,
                              [self.delta_mjd,delta_mjd_corr],
                              [mag_g[l], mag_g_corr[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_g))
                ]
                template_df_full['MAG_r'] = [
                    np.interp(delta_mjd_full,
                              [self.delta_mjd, delta_mjd_corr],
                              [mag_r[l], mag_r_corr[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_r))
                ]
                template_df_full['MAG_i'] = [
                    np.interp(delta_mjd_full,
                              [self.delta_mjd, delta_mjd_corr],
                              [mag_i[l], mag_i_corr[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_i))
                ]
                template_df_full['MAG_z'] = [
                    np.interp(delta_mjd_full,
                              [self.delta_mjd, delta_mjd_corr],
                              [mag_z[l], mag_z_corr[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_z))
                ]

            if delta_mjd_corr < self.delta_mjd:
                print("correcting magnitudes <")
                print(mag_g[0])
                template_df_full['MAG_g'] = [
                    np.interp(delta_mjd_full,
                              [delta_mjd_corr, self.delta_mjd],
                              [mag_g_corr[l], mag_g[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_g))
                ]
                template_df_full['MAG_r'] = [
                    np.interp(delta_mjd_full,
                              [delta_mjd_corr, self.delta_mjd],
                              [mag_r_corr[l], mag_r[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_r))
                    ]
                template_df_full['MAG_i'] = [
                    np.interp(delta_mjd_full,
                             [delta_mjd_corr, self.delta_mjd],
                             [mag_i_corr[l], mag_i[l]],
                             left=0,
                             right=0,
                             period=None) for l in range(0, len(mag_i))
                ]
                template_df_full['MAG_z'] = [
                    np.interp(delta_mjd_full,
                              [delta_mjd_corr, self.delta_mjd],
                              [mag_z_corr[l], mag_z[l]],
                              left=0,
                              right=0,
                              period=None) for l in range(0, len(mag_z))
                    ]

        loglans_full = []
        vks_full = []
        logmass_full = []
        for i in range(0, len(template_df_full['SIM_TEMPLATE_INDEX'])):
            _id = template_df_full['SIM_TEMPLATE_INDEX'][i]
            if _id == 27.5:
                _id = 99
                template_df_full['SIM_TEMPLATE_INDEX'][i] = 99
                template_df_later['SIM_TEMPLATE_INDEX'][i] = 99
            #print (_id)
            loglans_full.append(loglans[np.where(kn_inds == _id)[0]][0])
            vks_full.append(vks[np.where(kn_inds == _id)[0]][0])
            logmass_full.append(logmass_s[np.where(kn_inds == _id)[0]][0])

        loglans_full = np.array(loglans_full)
        vks_full = np.array(vks_full)
        logmass_full = np.array(logmass_full)
        mass_full = 10 ** logmass_full
        mass_ = 10 ** logmass_s
        preprocessed_weights = False

        if use_map:
            map_name = basename(self.event)
            path2map = dirname(self.event)

            if ".fits.gz" in map_name:
                map_name = map_name.rstrip(".fits.gz")
            if ".fits" in map_name:
                map_name = map_name.rstrip(".fits")

            lowres_path = join(path2map, "lowres")
            weights_path = join(path2map, "weights")
            use_map_info = join(lowres_path,f"{map_name}_mapinfo.npy")

            if not os.path.isdir(lowres_path):
                os.mkdir(lowres_path)
            
            if not os.path.isdir(weights_path):
                os.mkdir(weights_path)

            if m_exp_kncalc == False:
                use_map_weights = join(path2map, "weights"),

                use_map_weights_info = join(
                    path2map,
                    "weights",
                    f"{map_name}ac{area_covered}info.npy"
                )

                try:

                    weights_pre = np.load(use_map_weights)
                    use_map = ""
                    print("Using preprocessed weights for ", use_map_weights)
                    preprocessed_weights = True
                    area_deg_info, resolution = np.load(use_map_weights_info)
                    self.Resolution = resolution
                    self.area_deg = area_deg_info  # num_pix_covered*resolution

                except:
                    print(f"Failed to load preprocessed weights for {use_map_weights}")
            else:
                use_map_weights = join(path2map, "weights", f"{map_name}ac{area_covered}ad{deep_coverage}.npy")
                use_map_weights_info = join(path2map, "weights", f"{map_name}ac{area_covered}ad{deep_coverage}info.npy")

                try:
                    weights_pre, weights_pre_deep = np.load(use_map_weights)
                    area_deg_info, area_deg_deep_info, resolution = np.load(
                        use_map_weights_info
                    )
                    self.area_deg_deep = area_deg_deep_info  # num_pix_covered_deep*resolution
                    self.Resolution = resolution
                    self.area_deg = area_deg_info  # num_pix_covered*resolution
                    use_map = ""
                    print("Using preprocessed weights for ", use_map_weights)
                    preprocessed_weights = True
                except:
                    print("Failed to load preprocessed weights for ",
                          use_map_weights)

        if use_map != "":
            print("Using maps, not the distance")
            use_map_lowres = basename(self.event)
            path2map = dirname(self.event)
            use_map_lowres = join(path2map, "lowres", use_map_lowres)
            print("Tryng to open low resolution map")
            try:
                pb, distmu, distsigma = hp.read_map(use_map_lowres, field=range(
                    3), dtype=[np.float64, np.float64, np.float64])
                distmu_hr_average, distmu_std, distsigma_hr_average, distsigma_std = np.load(
                    use_map_info)
                reduce_mapresolution = False
                read_lowres_map = True
            except:
                print(
                    "Failed to open low resolution map, opening high resolution map", skymap_name)
                pb, distmu, distsigma, distnorm = hp.read_map(skymap_name, field=range(4), dtype=[
                                                              np.float64, np.float64, np.float64, np.float64])  # clecio dtype='numpy.float64'
                read_lowres_map = False

            pb_check = pb[np.logical_not(np.isinf(distmu))]
            distsigma_check = distsigma[np.logical_not(np.isinf(distmu))]
            distmu_check = distmu[np.logical_not(np.isinf(distmu))]

            pb_check = pb_check[np.logical_not(np.isinf(distsigma_check))]
            distmu_check = distmu_check[np.logical_not(
                np.isinf(distsigma_check))]
            distsigma_check = distsigma_check[np.logical_not(
                np.isinf(distsigma_check))]

            distmu_check_average = np.average(distmu_check, weights=pb_check)
            distsigma_check_average = np.average(
                distsigma_check, weights=pb_check)
            try:
                mean_z = z_at_value(cosmo.luminosity_distance, float(
                    distmu_check_average) * u.Mpc)
            except:
                print('Warning: Object too close with distance ',
                      str(distmu_check_average))
                self.Flag = 0
                return

            mean_z68 = z_at_value(cosmo.luminosity_distance, float(
                distmu_check_average+distsigma_check_average) * u.Mpc)
            z_max = max(template_df_full['ZMEAN'].values)
            luminosity_dist_max = max(cosmo.luminosity_distance(
                np.unique(template_df_full['ZMEAN'].values)))

            print(f"Mean redshift: {mean_z}")
            print(f"Mean z68: {mean_z68}")
            print(f"z max: {z_max}")
            if float(mean_z68) > float(z_max):
                print('Warning: Object too far away zmax= ',
                      str(z_max), ' z_event ', str(mean_z))
                self.Flag = 0
                return

            if np.isnan(pb).any():
                print('Warning: prob map contains nan')
                print('number of nans'+sum(np.isnan(pb))+' from '+len(pb))

            highres_sum = sum(pb)
            NSIDE = hp.npix2nside(pb.shape[0])
            resolution = (hp.nside2pixarea(NSIDE, degrees=True))

            if reduce_mapresolution == True:

                flag_inf_mu = False
                flag_inf_sigma = False
                if np.isinf(distmu).any():
                    flag_inf_mu = True
                    print('Warning: number of infs in distance array ' +
                          str(sum(np.isinf(distmu)))+' from '+str(len(distmu)))
                    print('prob of correspondend infs region ', sum(
                        pb[np.isinf(distmu)]), ' from ', sum(pb[np.logical_not(np.isinf(distmu))]))

                    pb_hr = pb[np.logical_not(np.isinf(distmu))]

                    distsigma_hr = distsigma[np.logical_not(np.isinf(distmu))]
                    distmu_hr = distmu[np.logical_not(np.isinf(distmu))]
                    distmu_hr_average = np.average(distmu_hr, weights=pb_hr)
                    distmu_std = weighted_avg_and_std(distmu_hr, weights=pb_hr)
                    distmu[np.isinf(distmu)] = 10000  # 10**25
                else:
                    distmu_hr_average = np.average(distmu, weights=pb)
                    distmu_std = weighted_avg_and_std(distmu, weights=pb)

                if np.isinf(distsigma).any():
                    flag_inf_sigma = True
                    print('Warning: number of infs in distance sigma array ' +
                          str(sum(np.isinf(distsigma)))+' from '+str(len(distsigma)))
                    print('prob of correspondend infs region ', sum(
                        pb[np.isinf(distsigma)]), ' from ', sum(pb[np.logical_not(np.isinf(distsigma))]))
                    pb_hr_sigma = pb[np.logical_not(np.isinf(distsigma))]
                    distsigma_hr = distsigma[np.logical_not(np.isinf(distmu))]
                    distsigma_hr_average = np.average(
                        distsigma_hr, weights=pb_hr_sigma)
                    distsigma_std = weighted_avg_and_std(
                        distsigma_hr, weights=pb_hr_sigma)
                    distsigma[np.isinf(distsigma)] = 10000

                else:
                    distsigma_hr_average = np.average(distsigma, weights=pb)
                    distsigma_std = weighted_avg_and_std(distsigma, weights=pb)
                
                np.save(use_map_info,
                        [
                            distmu_hr_average,
                            distmu_std,
                            distsigma_hr_average,
                            distsigma_std
                        ]
                )

                highres_sum = sum(pb)
                NSIDE = hp.npix2nside(pb.shape[0])
                res_high = (hp.nside2pixarea(NSIDE, degrees=True))

                target_res = 2.0
                final_nside_exp = math.ceil(
                    math.log((NSIDE*math.sqrt(res_high))/target_res, 2))  # int(NSIDE/8.0)
                final_nside = 2**final_nside_exp
                final_nside = int(final_nside)

                res_low = hp.nside2pixarea(final_nside, degrees=True)
                res_final = hp.nside2pixarea(final_nside, degrees=True)
                pb = hp.ud_grade(pb, final_nside, power=-2)
                lowres_sum = sum(pb)
                distsigma = hp.ud_grade(distsigma, final_nside)
                distmu = hp.ud_grade(distmu, final_nside)  # power=-2
                print('saving low resolution map')
                hp.write_map(use_map_lowres, m=[pb, distmu, distsigma], nest=False, dtype=None, fits_IDL=True,
                             coord=None, partial=False, column_names=None, column_units=None, extra_header=(), overwrite=False)

                if flag_inf_mu == True:
                    print('prob of reduced correspondend region due to infs ', sum(
                        pb[np.abs(distmu) > (distmu_hr_average+(3*distmu_std))]), ' from ', sum(pb))
                    pb = pb[np.abs(distmu) < (
                        distmu_hr_average+(3*distmu_std))]
                    distsigma = distsigma[np.abs(distmu) < (
                        distmu_hr_average+(3*distmu_std))]
                    distmu = distmu[np.abs(distmu) < (
                        distmu_hr_average+(3*distmu_std))]
                if flag_inf_sigma == True:
                    print('prob of reduced correspondend region due to sigma infs ', sum(pb[np.abs(
                        distsigma) > (distsigma_hr_average+(3*distsigma_std))]), ' from ', sum(pb))
                    pb = pb[np.abs(distsigma) < (
                        distsigma_hr_average+(3*distsigma_std))]
                    distmu = distmu[np.abs(distsigma) < (
                        distsigma_hr_average+(3*distsigma_std))]
                    distsigma = distsigma[np.abs(distsigma) < (
                        distsigma_hr_average+(3*distsigma_std))]  # (distmu_hr_average+(3*distmu_std))

                distmu_lowres_average = np.average(distmu, weights=pb)
                distsigma_lowres_average = np.average(distsigma, weights=pb)

                print("Reducing map resolution ")
                print("Previous resolution and new resolution (deg) ",
                      res_high, " ", res_low)
                print("Sums - high res map prob sum ", highres_sum,
                      " low res map pb sum ", lowres_sum)  # " ",sum(new_pb))
                print("average distance in high and low res map with its standard deviation",
                      distmu_hr_average, "+-", distmu_std, "and ", distmu_lowres_average)
                print("average distance sigma in high and low res map with its standard deviation",
                      distsigma_hr_average, "+-", distsigma_std, "and ", distsigma_lowres_average)
                resolution = res_low
            if read_lowres_map == True:
                print('prob of reduced correspondend region due to infs ', sum(
                    pb[np.abs(distmu) > (distmu_hr_average+(3*distmu_std))]), ' from ', sum(pb))
                pb = pb[np.abs(distmu) < (distmu_hr_average+(3*distmu_std))]
                distsigma = distsigma[np.abs(distmu) < (
                    distmu_hr_average+(3*distmu_std))]
                distmu = distmu[np.abs(distmu) < (
                    distmu_hr_average+(3*distmu_std))]
                print('prob of reduced correspondend region due to sigma infs ', sum(pb[np.abs(
                    distsigma) > (distsigma_hr_average+(3*distsigma_std))]), ' from ', sum(pb))
                pb = pb[np.abs(distsigma) < (
                    distsigma_hr_average+(3*distsigma_std))]
                distmu = distmu[np.abs(distsigma) < (
                    distsigma_hr_average+(3*distsigma_std))]
                distsigma = distsigma[np.abs(distsigma) < (
                    distsigma_hr_average+(3*distsigma_std))]  # (distmu_hr_average+(3*distmu_std))

            if np.isinf(distmu).any():
                print('Warning: number of infs in distance array ' +
                      str(sum(np.isinf(distmu)))+' from '+str(len(distmu)))
                print('prob of correspondend infs region ', sum(
                    pb[np.isinf(distmu)]), ' from ', sum(pb[np.logical_not(np.isinf(distmu))]))

                pb = pb[np.logical_not(np.isinf(distmu))]
                distsigma = distsigma[np.logical_not(np.isinf(distmu))]
                distmu = distmu[np.logical_not(np.isinf(distmu))]

            if np.isnan(distmu).any():
                print('Warning: distance map contains nan')
                print('number of nans '+str(sum(np.isnan(distmu))) +
                      ' from '+str(len(distmu)))
                print('prob of correspondend nan region ', sum(
                    pb[np.isnan(distmu)]), ' from ', sum(pb[np.logical_not(np.isnan(distmu))]))
                pb = pb[np.logical_not(np.isnan(distmu))]
                distsigma = distsigma[np.logical_not(np.isnan(distmu))]
                distmu = distmu[np.logical_not(np.isnan(distmu))]

            idx_sort = np.argsort(pb)
            idx_sort_up = list(reversed(idx_sort))

            sum_ = 0.
            id_c = 0
            sum_full = 0
            id_full = 0
            area_max = sum(pb)
            id_deep = 0

            while (sum_full < min(area_max-0.01, 0.98)) and (id_full < len(idx_sort_up)):
                this_idx = idx_sort_up[id_full]
                sum_full = sum_full+pb[this_idx]
                id_full = id_full+1
            total_area = id_full*resolution
            print("Total event area (deg)="+str(id_full*resolution)+" ")
            while (sum_ < area_covered) and (id_c < len(idx_sort_up)):
                this_idx = idx_sort_up[id_c]
                sum_ = sum_+pb[this_idx]

                if m_exp_kncalc == True:
                    if sum_ < deep_coverage:
                        id_deep = id_deep+1

                id_c = id_c+1

            idx_sort_full = idx_sort_up[:id_full]

            idx_sort_uncovered = idx_sort_up[id_c:]
            if m_exp_kncalc == True:
                print("number of deep coverage pixels and total covered are pixels=" +
                      str(id_deep)+" "+str(id_c))
            if (id_deep == 0) and (m_exp == True):
                print(
                    "===== Warning: number of deep coverage is 0, switching to single exposure mode")
                m_exp_kncalc = False
                self.sw_mexp = True
            else:
                self.sw_mexp = False
            if m_exp_kncalc == False:
                idx_sort_cut = idx_sort_up[:id_c]

            else:
                idx_sort_cut = idx_sort_up[id_deep:id_c]
                idx_sort_deep = idx_sort_up[:id_deep]
                pb_covered_deep = pb[idx_sort_deep]
                num_pix_covered_deep = len(pb_covered_deep)
                distmu_covered_deep = distmu[idx_sort_deep]
                distsigma_covered_deep = distsigma[idx_sort_deep]

            pb_covered = pb[idx_sort_cut]
            num_pix_covered = len(pb_covered)
            distmu_covered = distmu[idx_sort_cut]
            distsigma_covered = distsigma[idx_sort_cut]

            try:
                distmu_average = np.average(distmu_covered, weights=pb_covered)
                distsigma_average = np.average(
                    distsigma_covered, weights=pb_covered)
            except:
                self.Flag = 0
                return
            print("the distance average="+str(distmu_average))
            print("the distance sigma="+str(distsigma_average))
            if m_exp_kncalc == False:
                print("the area covered="+str(area_covered))
                print("the prob covered="+str(np.sum(pb_covered)))

            distmu_full = distmu[idx_sort_full]
            pb_full = pb[idx_sort_full]
            distsigma_full = distsigma[idx_sort_full]
            if (distmu_full < 0.0).any():
                print('Warning: prob of distmu <0.0 in the full region region ', sum(
                    pb_full[distmu_full < 0.0]), ' from ', sum(pb_full[distmu_full > 0.0]))
                distmu_full[distmu_full < 0.0] = 0.01
            if (distmu_covered < 0.0).any():
                print('Warning: prob of distmu<0.0 region ', sum(
                    pb_covered[distmu_covered < 0.0]), ' from ', sum(pb_covered[distmu_covered > 0.0]))
                distmu_covered[distmu_covered < 0.0] = 0.01
            if m_exp_kncalc == True:
                if (distmu_covered < 0.0).any():
                    print('Warning: prob of distmu<0.0 in the deep region ', sum(
                        pb_covered_deep[distsigma_covered_deep < 0.0]), ' from ', sum(pb_covered_deep[distsigma_covered_deep > 0.0]))
                    distsigma_covered_deep[distsigma_covered_deep < 0.0] = 0.01

            pb_vol_norm = np.sum(np.multiply(distmu_full, pb_full))

            pb_vol_covered_all = np.sum(
                np.multiply(distmu_covered, pb_covered))

            pb_vol_covered_all = pb_vol_covered_all/pb_vol_norm

            print("the prob volume covered (except deep region if any)=" +
                  str(pb_vol_covered_all))

            if np.isnan(pb_full).any():
                print('Warning: prob map full cut contains nan')
                print('number of nans'+sum(np.isnan(pb_full)) +
                      ' from '+len(pb_full))

            for k in range(0, len(pb_covered)):
                weights_pix = [norm.pdf(x.value, loc=float(distmu_covered[k]), scale=float(
                    distsigma_covered[k])) for x in cosmo.luminosity_distance(template_df_full['ZMEAN'].values)]
                weights_pix_norm = [norm.pdf(x.value, loc=float(distmu_covered[k]), scale=float(
                    distsigma_covered[k])) for x in cosmo.luminosity_distance(np.unique(template_df_full['ZMEAN'].values))]


                pb_vol_covered = (
                    pb_covered[k] * distmu_covered[k])/pb_vol_norm
                if np.sum(weights_pix_norm) == 0.0:
                    print("the weights pix sum is 0, skipping pixel")
                    continue
                weights_pix = (weights_pix / np.sum(weights_pix_norm))
                weights_pix = weights_pix*float(pb_vol_covered)
                if k == 0:
                    weights_pix_area = weights_pix
                else:
                    weights_pix_area = np.add(weights_pix_area, weights_pix)

            print("the weights sum="+str(np.sum(weights_pix_area)))
            template_df_full['WEIGHT'] = weights_pix_area
            if m_exp_kncalc == True:
                print("Getting the weights sum of Deep region with " +
                      str(len(pb_covered_deep))+" voxels")
                for k in range(0, len(pb_covered_deep)):
                    weights_pix_deep = [norm.pdf(x.value, loc=float(distmu_covered_deep[k]), scale=float(
                        distsigma_covered_deep[k])) for x in cosmo.luminosity_distance(template_df_full['ZMEAN'].values)]
                    weights_pix_norm_deep = [norm.pdf(x.value, loc=float(distmu_covered_deep[k]), scale=float(
                        distsigma_covered_deep[k])) for x in cosmo.luminosity_distance(np.unique(template_df_full['ZMEAN'].values))]

                    pb_vol_covered_deep = (
                        pb_covered_deep[k] * distmu_covered_deep[k])/pb_vol_norm
                    if np.sum(weights_pix_norm) == 0.0:
                        print("the weights pix sum is 0, skipping pixel")
                        continue
                    weights_pix_deep = (
                        weights_pix_deep / np.sum(weights_pix_norm_deep))
                    weights_pix_deep = weights_pix_deep * \
                        float(pb_vol_covered_deep)
                    if k == 0:
                        weights_pix_area_deep = weights_pix_deep
                    else:
                        weights_pix_area_deep = np.add(
                            weights_pix_area_deep, weights_pix_deep)

                print("the weights sum of Deep area=" +
                      str(np.sum(weights_pix_area_deep)))
                template_df_full['WEIGHT_deep'] = weights_pix_area_deep
                np.save(use_map_weights, [
                        weights_pix_area, weights_pix_area_deep])
                np.save(use_map_weights_info, [
                        num_pix_covered*resolution, num_pix_covered_deep*resolution, resolution])
            else:
                np.save(use_map_weights, weights_pix_area)
                np.save(use_map_weights_info, [
                        num_pix_covered*resolution, resolution])

            self.Resolution = resolution
            self.area_deg = num_pix_covered*resolution
            if m_exp_kncalc == True:
                self.area_deg_deep = num_pix_covered_deep*resolution
                # use_map_weights_info

        elif preprocessed_weights == True:

            template_df_full['WEIGHT'] = weights_pre
            if m_exp_kncalc == True:
                template_df_full['WEIGHT_deep'] = weights_pre_deep

        else:
            weights = [norm.pdf(x.value, loc=float(self.distance), scale=float(
                self.distance_err)) for x in cosmo.luminosity_distance(template_df_full['ZMEAN'].values)]
            weights_norm = [norm.pdf(x.value, loc=float(self.distance), scale=float(
                self.distance_err)) for x in cosmo.luminosity_distance(np.unique(template_df_full['ZMEAN'].values))]

            template_df_full['WEIGHT'] = (
                weights / np.sum(weights_norm))  # *area_covered
            print("the weights sum (no map)=" +
                  str(np.sum(template_df_full['WEIGHT'])))
            print("probability covered (no map)="+str(area_covered))

        template_df_later['WEIGHT'] = template_df_full['WEIGHT'].values
        if m_exp_kncalc == True:
            template_df_later['WEIGHT_deep'] = template_df_full['WEIGHT_deep'].values

        self.delta_mjd_later = delta_mjd_later
        self.template_df_later = template_df_later
        self.template_df_full = template_df_full
        self.Flag = 1
        self.mult_exp = m_exp_kncalc
        return


def telescope_time(exptime, area_deg, field_of_view_area=3.0, m_exp=False):
    if m_exp == False:
        num_exposures = area_deg/field_of_view_area
        return num_exposures*(30.0 + exptime)
    else:
        total_time = 0
        for i in range(0, len(exptime)):
            num_exposures = area_deg[i]/field_of_view_area
            total_time += num_exposures*(30.0 + exptime[i])
        return total_time


# Functions to calculate metrics of interest
def weighted_average(quantity, weights):
    print('weighted_average')
    print(quantity.shape)
    print(weights.shape)
    a = np.dot(quantity, weights) / np.sum(weights)
    return a


def weighted_average_multi(quantity, weights):

    for i in range(0, len(weights)):
        if i == 0:
            out_average = weights[i].copy()
        else:
            out_average = np.array(out_average) * weights[i]
    out_average = out_average/sum(out_average)

    print('combined weights')
    print(len(weights))
    # print(out_average)
    return weighted_average(quantity, out_average)


def get_all_mags(data, use_knmodel_weights=False, kn_type='red'):
    if use_knmodel_weights == True:
        g_m = weighted_average_multi(data['MAG_g'].values,
                                     [
                                      data['WEIGHT'].values,
                                      data[f'WEIGHT_loglan_{kn_type}'].values,
                                      data[f'WEIGHT_vk_{kn_type}'].values,
                                      data[f'WEIGHT_mass_{kn_type}'].values
                                     ])

        g_merr = weighted_average_multi(data['MAGERR_g'].values,
                                        [
                                          data['WEIGHT'].values,
                                          data[f'WEIGHT_loglan_{kn_type}'].values,
                                          data[f'WEIGHT_vk_{kn_type}'].values,
                                          data[f'WEIGHT_mass_{kn_type}'].values
                                        ])

        r_m = weighted_average_multi(data['MAG_r'].values,
                                     [
                                      data['WEIGHT'].values,
                                      data[f'WEIGHT_loglan_{kn_type}'].values,
                                      data[f'WEIGHT_vk_{kn_type}'].values,
                                      data[f'WEIGHT_mass_{kn_type}'].values
                                     ])

        r_merr = weighted_average_multi(data['MAGERR_r'].values,
                                        [
                                          data['WEIGHT'].values,
                                          data[f'WEIGHT_loglan_{kn_type}'].values,
                                          data[f'WEIGHT_vk_{kn_type}'].values,
                                          data[f'WEIGHT_mass_{kn_type}'].values
                                        ])

        i_m = weighted_average_multi(data['MAG_i'].values,
                                     [
                                       data['WEIGHT'].values,
                                       data[f'WEIGHT_loglan_{kn_type}'].values,
                                       data[f'WEIGHT_vk_{kn_type}'].values,
                                       data[f'WEIGHT_mass_{kn_type}'].values
                                     ])

        i_merr = weighted_average_multi(data['MAGERR_i'].values,
                                        [
                                          data['WEIGHT'].values,
                                          data[f'WEIGHT_loglan_{kn_type}'].values,
                                          data[f'WEIGHT_vk_{kn_type}'].values,
                                          data[f'WEIGHT_mass_{kn_type}'].values
                                        ])

        z_m = weighted_average_multi(data['MAG_z'].values,
                                     [
                                      data['WEIGHT'].values,
                                      data[f'WEIGHT_loglan_{kn_type}'].values,
                                      data['WEIGHT_vk_'+kn_type].values,
                                      data[f'WEIGHT_mass_{kn_type}'].values
                                     ])

        z_merr = weighted_average_multi(data['MAGERR_z'].values,
                                        [
                                          data['WEIGHT'].values,
                                          data[f'WEIGHT_loglan_{kn_type}'].values,
                                          data['WEIGHT_vk_'+kn_type].values,
                                          data[f'WEIGHT_mass_{kn_type}'].values
                                        ])

        out_dict = {'g_mag': g_m,
                    'r_mag': r_m,
                    'i_mag': i_m,
                    'z_mag': z_m,
                    'g_magerr': g_merr,
                    'r_magerr': r_merr,
                    'i_magerr': i_merr,
                    'z_magerr': z_merr}

    else:

        out_dict = {'g_mag': weighted_average(data['MAG_g'].values,
                                            data['WEIGHT'].values),
                    'r_mag': weighted_average(data['MAG_r'].values,
                                            data['WEIGHT'].values),
                    'i_mag': weighted_average(data['MAG_i'].values,
                                            data['WEIGHT'].values),
                    'z_mag': weighted_average(data['MAG_z'].values,
                                            data['WEIGHT'].values),
                    'g_magerr': weighted_average(data['MAGERR_g'].values,
                                                data['WEIGHT'].values),
                    'r_magerr': weighted_average(data['MAGERR_r'].values,
                                                data['WEIGHT'].values),
                    'i_magerr': weighted_average(data['MAGERR_i'].values,
                                                data['WEIGHT'].values),
                    'z_magerr': weighted_average(data['MAGERR_z'].values,
                                                data['WEIGHT'].values)}
    return out_dict


def gw170817(data):
    blue_df = data[data['SIM_TEMPLATE_INDEX'].values == 178.0]
    red_df = data[data['SIM_TEMPLATE_INDEX'].values == 224.0]
    print('get_all_mags call')
    blue_mags = get_all_mags(blue_df)
    red_mags = get_all_mags(red_df)

    return blue_mags, red_mags


def print_dict(d, title='', outfile=None):

    outdata = [title + '\n--------------------------------------------------------------\n',
               "\tg\t\tr\t\ti\t\tz\n",
               "--------------\t--------------\t--------------\t--------------\t\n",
               "%.2f +/- %.2f\t%.2f +/- %.2f\t%.2f +/- %.2f\t%.2f +/- %.2f\n\n" % (d['g_mag'], d['g_magerr'],
                                                                                   d['r_mag'], d['r_magerr'],
                                                                                   d['i_mag'], d['i_magerr'],
                                                                                   d['z_mag'], d['z_magerr'])]
    if outfile:
        stream = open(outfile + '.txt', 'a+')
        stream.writelines(outdata)
        stream.close()
    else:
        for row in outdata:
            print(row[:-1])

    return


def calc_mag_fractions(data,
                       use_knmodel_weights=False,
                       kn_type='red',
                       info_file="knsed_info.txt",
                       kn_weight_type="uniform",
                       kn_weight_sigma=1.0,
                       model_weights_path='',
                       m_exp=False):

    if use_knmodel_weights == True:
        print('===== Using KN model priors '+kn_weight_type +
              ' , considering gw170817 with '+kn_type+' component')

    percentile_levels = np.linspace(0.0, 100.0, 101)
    percentile_dict = {}
    percentile_dict_deep = {}
    prob_dict = {}

    try:
        model_weights = np.load(
            model_weights_path+'kn_weights_type'+kn_type+'prior'+kn_weight_type+'.npy')
        ids_total = np.load(
            model_weights_path+'kn_weights_ids_type'+kn_type+'prior'+kn_weight_type+'.npy')
        print('loading KN model weights ...')
    except:
        print('Warning: Building model weights from scratch. This is a slow process.')
        model_weights, ids_total = get_model_weights(
            kn_weight_type=kn_weight_type, kn_type=kn_type, info_file=info_file, kn_weight_sigma=kn_weight_sigma)
        np.save(model_weights_path+'kn_weights_type'+kn_type +
                'prior'+kn_weight_type+'.npy', model_weights)
        np.save(model_weights_path+'kn_weights_ids_type' +
                kn_type+'prior'+kn_weight_type+'.npy', ids_total)

    mags_range = np.arange(14, 28, 0.2)
    prob_at_maglim = {}
    prob_at_maglim_deep = {}
    for band in ['g', 'r', 'i', 'z']:
        prob_at_maglim['prob_'+band] = []
        data_test = {'MAG_'+band: np.array(data['MAG_'+band].values).astype(
            'float'), 'ids_': data['SIM_TEMPLATE_INDEX'].values, 'weights_z': data['WEIGHT']}
        if m_exp == True:
            prob_at_maglim_deep['prob_'+band] = []
            data_test_deep = {'MAG_'+band: np.array(data['MAG_'+band].values).astype(
                'float'), 'ids_': data['SIM_TEMPLATE_INDEX'].values, 'weights_z': data['WEIGHT_deep']}
        for j in range(0, len(mags_range)):
            # all_objs=data.copy()
            mask_det = (data_test['MAG_'+band] < mags_range[j])
            # print(mask_det)
            detected_objs = {'MAG_'+band: data_test['MAG_'+band][mask_det],
                             'ids_': data_test['ids_'][mask_det], 'weights_z': data_test['weights_z'][mask_det]}

            if m_exp == True:
                detected_objs_deep = {'MAG_'+band: data_test_deep['MAG_'+band][mask_det],
                                      'ids_': data_test_deep['ids_'][mask_det], 'weights_z': data_test_deep['weights_z'][mask_det]}

                if len(detected_objs_deep['MAG_'+band]) > 0:
                    prob_at_maglim_unweighted_deep = calc_prob_redshift(
                        pd.DataFrame(data=detected_objs_deep), band=band, quiet=True)
                    prob_at_maglim_weighted_deep = []
                    prob_at_maglim_weighted_deep = [prob_at_maglim_unweighted_deep['prob_'+band][i]*model_weights[np.where(
                        ids_total == prob_at_maglim_unweighted_deep['_ids'][i])] for i in range(0, len(prob_at_maglim_unweighted_deep['prob_'+band]))]
                    prob_at_maglim_deep['prob_'+band].append(
                        float(sum(prob_at_maglim_weighted_deep)))
                else:
                    prob_at_maglim_deep['prob_'+band].append(0.0)

            if len(detected_objs['MAG_'+band]) > 0:

                prob_at_maglim_unweighted = calc_prob_redshift(
                    pd.DataFrame(data=detected_objs), band=band, quiet=True)
                prob_total_unweighted = calc_prob_redshift(
                    pd.DataFrame(data=data_test), band=band, quiet=True)

                prob_at_maglim_weighted = []

                prob_at_maglim_weighted = [prob_at_maglim_unweighted['prob_'+band][i]*model_weights[np.where(
                    ids_total == prob_at_maglim_unweighted['_ids'][i])] for i in range(0, len(prob_at_maglim_unweighted['prob_'+band]))]
                prob_at_maglim['prob_' +
                               band].append(float(sum(prob_at_maglim_weighted)))
            else:
                prob_at_maglim['prob_'+band].append(0.0)
        prob_at_maglim['prob_' +
                       band] = np.array(prob_at_maglim['prob_'+band])*100.0
        if m_exp == True:
            prob_at_maglim_deep['prob_' +
                                band] = np.array(prob_at_maglim_deep['prob_'+band])*100.0
        percentile_dict['%s_cutoff' % band] = np.interp(
            percentile_levels, prob_at_maglim['prob_'+band], mags_range, left=0, right=100, period=None)
        if m_exp == True:
            percentile_dict_deep['%s_cutoff' % band] = np.interp(
                percentile_levels, prob_at_maglim_deep['prob_'+band], mags_range, left=0, right=100, period=None)

    if m_exp == True:
        return percentile_dict, percentile_dict_deep
    return percentile_dict


def get_model_weights(kn_weight_type="uniform", kn_type='red', info_file="knsed_info.txt", kn_weight_sigma=1.0):

    sed_filenames, kn_inds, vks, loglans, logmass_s = open_ascii_cat(
        info_file, unpack=True)  # usecols=(0,1)
    loglans = np.array(loglans).astype('float')
    kn_inds = np.array(kn_inds).astype('float')

    vks = np.array(vks).astype('float')
    logmass_s = np.array(logmass_s).astype('float')

    # https://science.sciencemag.org/content/358/6370/1583
    mass_red = 0.035

    loglan_red = -2.0

    vk_red = 0.15

    mass_blue = 0.025

    loglan_blue = -5.0

    vk_blue = 0.25

    weights_dict = {}  # template_df_full

    mass_ = 10 ** logmass_s

    if kn_weight_type == "gaussian":
        print("gaussian prior")
        mass_blue_err = 0.001*10
        loglan_blue_err = 1.0*10
        vk_blue_err = 0.01*10
        mass_red_err = 0.015*10
        loglan_red_err = 0.5*10
        vk_red_err = 0.03*10

        weights_loglan_red = [norm.pdf(x, loc=float(
            loglan_red), scale=kn_weight_sigma*float(loglan_red_err)) for x in loglans]
        weights_loglan_red_norm = [norm.pdf(x, loc=float(
            loglan_red), scale=kn_weight_sigma*float(loglan_red_err)) for x in np.unique(loglans)]
        weights_dict['WEIGHT_loglan_red'] = weights_loglan_red / \
            np.sum(weights_loglan_red_norm)

        weights_vks_red = [norm.pdf(x, loc=float(
            vk_red), scale=kn_weight_sigma*float(vk_red_err)) for x in vks]
        weights_vks_red_norm = [norm.pdf(x, loc=float(
            vk_red), scale=kn_weight_sigma*float(vk_red_err)) for x in np.unique(vks)]
        weights_dict['WEIGHT_vk_red'] = weights_vks_red / \
            np.sum(weights_vks_red_norm)

        weights_mass_red = [norm.pdf(x, loc=float(
            mass_red), scale=kn_weight_sigma*float(mass_red_err)) for x in mass_]
        weights_mass_red_norm = [norm.pdf(x, loc=float(
            mass_red), scale=kn_weight_sigma*float(mass_red_err)) for x in np.unique(mass_)]
        weights_dict['WEIGHT_mass_red'] = weights_mass_red / \
            np.sum(weights_mass_red_norm)

        weights_loglan_blue = [norm.pdf(x, loc=float(
            loglan_blue), scale=kn_weight_sigma*float(loglan_blue_err)) for x in loglans]
        weights_loglan_blue_norm = [norm.pdf(x, loc=float(
            loglan_blue), scale=kn_weight_sigma*float(loglan_blue_err)) for x in np.unique(loglans)]
        weights_dict['WEIGHT_loglan_blue'] = weights_loglan_blue / \
            np.sum(weights_loglan_blue_norm)

        weights_vks_blue = [norm.pdf(x, loc=float(
            vk_blue), scale=kn_weight_sigma*float(vk_blue_err)) for x in vks]
        weights_vks_blue_norm = [norm.pdf(x, loc=float(
            vk_blue), scale=kn_weight_sigma*float(vk_blue_err)) for x in np.unique(vks)]
        weights_dict['WEIGHT_vk_blue'] = weights_vks_blue / \
            np.sum(weights_vks_blue_norm)

        weights_mass_blue = [norm.pdf(x, loc=float(
            mass_blue), scale=kn_weight_sigma*float(mass_blue_err)) for x in mass_]
        weights_mass_blue_norm = [norm.pdf(x, loc=float(
            mass_blue), scale=kn_weight_sigma*float(mass_blue_err)) for x in np.unique(mass_)]
        weights_dict['WEIGHT_mass_blue'] = weights_mass_blue / \
            np.sum(weights_mass_blue_norm)
    if kn_weight_type == "gaussian_narrow":
        print("gaussian prior narrow")
        mass_blue_err = 0.001
        loglan_blue_err = 1.0
        vk_blue_err = 0.01
        mass_red_err = 0.015
        loglan_red_err = 0.5
        vk_red_err = 0.03

        weights_loglan_red = [norm.pdf(x, loc=float(
            loglan_red), scale=kn_weight_sigma*float(loglan_red_err)) for x in loglans]
        weights_loglan_red_norm = [norm.pdf(x, loc=float(
            loglan_red), scale=kn_weight_sigma*float(loglan_red_err)) for x in np.unique(loglans)]
        weights_dict['WEIGHT_loglan_red'] = weights_loglan_red / \
            np.sum(weights_loglan_red_norm)

        weights_vks_red = [norm.pdf(x, loc=float(
            vk_red), scale=kn_weight_sigma*float(vk_red_err)) for x in vks]
        weights_vks_red_norm = [norm.pdf(x, loc=float(
            vk_red), scale=kn_weight_sigma*float(vk_red_err)) for x in np.unique(vks)]
        weights_dict['WEIGHT_vk_red'] = weights_vks_red / \
            np.sum(weights_vks_red_norm)

        weights_mass_red = [norm.pdf(x, loc=float(
            mass_red), scale=kn_weight_sigma*float(mass_red_err)) for x in mass_]
        weights_mass_red_norm = [norm.pdf(x, loc=float(
            mass_red), scale=kn_weight_sigma*float(mass_red_err)) for x in np.unique(mass_)]
        weights_dict['WEIGHT_mass_red'] = weights_mass_red / \
            np.sum(weights_mass_red_norm)

        weights_loglan_blue = [norm.pdf(x, loc=float(
            loglan_blue), scale=kn_weight_sigma*float(loglan_blue_err)) for x in loglans]
        weights_loglan_blue_norm = [norm.pdf(x, loc=float(
            loglan_blue), scale=kn_weight_sigma*float(loglan_blue_err)) for x in np.unique(loglans)]
        weights_dict['WEIGHT_loglan_blue'] = weights_loglan_blue / \
            np.sum(weights_loglan_blue_norm)

        weights_vks_blue = [norm.pdf(x, loc=float(
            vk_blue), scale=kn_weight_sigma*float(vk_blue_err)) for x in vks]
        weights_vks_blue_norm = [norm.pdf(x, loc=float(
            vk_blue), scale=kn_weight_sigma*float(vk_blue_err)) for x in np.unique(vks)]
        weights_dict['WEIGHT_vk_blue'] = weights_vks_blue / \
            np.sum(weights_vks_blue_norm)

        weights_mass_blue = [norm.pdf(x, loc=float(
            mass_blue), scale=kn_weight_sigma*float(mass_blue_err)) for x in mass_]
        weights_mass_blue_norm = [norm.pdf(x, loc=float(
            mass_blue), scale=kn_weight_sigma*float(mass_blue_err)) for x in np.unique(mass_)]
        weights_dict['WEIGHT_mass_blue'] = weights_mass_blue / \
            np.sum(weights_mass_blue_norm)

    if kn_weight_type == "uniform":
        epsilon = 0.00000001
        mass_blue_err = 0.001*10
        loglan_blue_err = 1.0*10
        vk_blue_err = 0.01*10
        mass_red_err = 0.015*10
        loglan_red_err = 0.5*10
        vk_red_err = 0.03*10

        print("uniform prior")
        weights_loglan_red = [uniform.pdf(x, loc=float(
            loglan_red)-float(loglan_red_err), scale=2*float(loglan_red_err+epsilon)) for x in loglans]
        weights_loglan_red_norm = [uniform.pdf(x, loc=float(loglan_red)-float(
            loglan_red_err), scale=2*float(loglan_red_err+epsilon)) for x in np.unique(loglans)]
        weights_dict['WEIGHT_loglan_red'] = weights_loglan_red / \
            np.sum(weights_loglan_red_norm)

        vk_red_err_ext = 0.05000001
        weights_vks_red = [uniform.pdf(x, loc=float(
            vk_red)-float(vk_red_err_ext), scale=2*float(vk_red_err_ext)) for x in vks]
        weights_vks_red_norm = [uniform.pdf(x, loc=float(
            vk_red)-float(vk_red_err_ext), scale=2*float(vk_red_err_ext)) for x in np.unique(vks)]
        weights_dict['WEIGHT_vk_red'] = weights_vks_red / \
            np.sum(weights_vks_red_norm)

        weights_mass_red = [uniform.pdf(x, loc=float(
            mass_red)-float(mass_red_err), scale=2*float(mass_red_err+epsilon)) for x in mass_]
        weights_mass_red_norm = [uniform.pdf(x, loc=float(
            mass_red)-float(mass_red_err), scale=2*float(mass_red_err+epsilon)) for x in np.unique(mass_)]
        weights_dict['WEIGHT_mass_red'] = weights_mass_red / \
            np.sum(weights_mass_red_norm)

        weights_loglan_blue = [uniform.pdf(x, loc=float(
            loglan_blue)-float(loglan_blue_err), scale=2*float(loglan_blue_err+epsilon)) for x in loglans]
        weights_loglan_blue_norm = sum([uniform.pdf(x, loc=float(loglan_blue)-float(
            loglan_blue_err), scale=2*float(loglan_blue_err+epsilon)) for x in np.unique(loglans)])
        weights_dict['WEIGHT_loglan_blue'] = weights_loglan_blue / \
            weights_loglan_blue_norm

        vk_blue_err_ext = 0.1000001
        weights_vks_blue = [uniform.pdf(x, loc=float(
            vk_blue)-float(vk_blue_err_ext), scale=2*float(vk_blue_err_ext)) for x in vks]
        weights_vks_blue_norm = sum([uniform.pdf(x, loc=float(
            vk_blue)-float(vk_blue_err_ext), scale=2*float(vk_blue_err_ext)) for x in np.unique(vks)])
        weights_dict['WEIGHT_vk_blue'] = weights_vks_blue / \
            np.sum(weights_vks_blue_norm)  # np.sum(weights_vks_blue)

        weights_mass_blue = [uniform.pdf(x, loc=float(
            mass_blue)-float(mass_blue_err), scale=2*float(mass_blue_err+epsilon)) for x in mass_]
        weights_mass_blue_norm = [uniform.pdf(x, loc=float(mass_blue)-float(
            mass_blue_err), scale=2*float(mass_blue_err+epsilon)) for x in np.unique(mass_)]
        weights_dict['WEIGHT_mass_blue'] = weights_mass_blue / \
            np.sum(weights_mass_blue_norm)

        if np.isnan(weights_dict['WEIGHT_vk_blue']).any():
            print('WEIGHT_vk_blue NAN')
        if np.isnan(weights_dict['WEIGHT_mass_blue']).any():
            print('WEIGHT_mass_blue NAN')
        if np.isnan(weights_dict['WEIGHT_loglan_blue']).any():
            print('WEIGHT_loglan_blue NAN')

        if np.isnan(weights_dict['WEIGHT_vk_red']).any():
            print('WEIGHT_vk_red NAN')
            # print(vks)
        if np.isnan(weights_dict['WEIGHT_mass_red']).any():
            print('WEIGHT_mass_red NAN')
        if np.isnan(weights_dict['WEIGHT_loglan_red']).any():
            print('WEIGHT_loglan_red NAN')


    if kn_weight_type == "uniform_loglan":
        print("uniform cut in loglan")
        weights_loglan_blue = [uniform.pdf(
            x, loc=float(loglan), scale=0.1) for x in loglans]
        template_df_full['WEIGHT_loglan_blue'] = weights_loglan_blue / \
            np.sum(weights_loglan_blue)
        weights_loglan_red = [uniform.pdf(x, loc=float(
            loglan), scale=0.1) for x in loglans_full]
        weights_dict['WEIGHT_loglan_red'] = weights_loglan_red / \
            np.sum(weights_loglan_red)

        weights_vks_red = [uniform.pdf(x, loc=min(
            vks_full)-1, scale=max(vks_full)+1) for x in vks]
        weights_dict['WEIGHT_vk_red'] = weights_vks_red / \
            np.sum(weights_vks_red)

        weights_mass_red = [uniform.pdf(x, loc=min(
            mass_full)-1, scale=max(mass_full)+1) for x in mass_full]
        weights_dict['WEIGHT_mass_red'] = weights_mass_red / \
            np.sum(weights_mass_red)

        weights_vks_blue = [uniform.pdf(x, loc=min(
            vks_full)-1, scale=max(vks_full)+1) for x in vks]
        weights_dict['WEIGHT_vk_blue'] = weights_vks_blue / \
            np.sum(weights_vks_blue)

        weights_mass_blue = [uniform.pdf(x, loc=min(
            mass_full)-1, scale=max(mass_full)+1) for x in mass_]
        weights_dict['WEIGHT_mass_blue'] = weights_mass_blue / \
            np.sum(weights_mass_blue)
    if kn_type == 'red':
        model_weights = np.array(weights_dict['WEIGHT_vk_red'])*np.array(
            weights_dict['WEIGHT_mass_red'])*np.array(weights_dict['WEIGHT_loglan_red'])
    if kn_type == 'blue':
        model_weights = np.array(weights_dict['WEIGHT_vk_blue'])*np.array(
            weights_dict['WEIGHT_mass_blue'])*np.array(weights_dict['WEIGHT_loglan_blue'])
    if kn_type == 'blue_GW170817':
        model_weights = np.zeros(len(vks))
        model_weights[177] = 1.0
    if kn_type == 'red_GW170817':
        model_weights = np.zeros(len(vks))
        model_weights[223] = 1.0
    return model_weights, kn_inds


def calc_prob_redshift(data, band, quiet=True):

    data_out = pd.DataFrame(data=[get_weighted_redshift_prob_dict(
        data[data['ids_'].values == y], band=band, quiet=quiet) for y in np.unique(data['ids_'].values)])
    data_out['_ids'] = np.unique(data['ids_'].values)

    if quiet == False:
        print(band)
        print(data_out['prob_'+band].shape)
        print(np.sum(data_out['prob_'+band]))
        # print()
    return data_out


def get_weighted_redshift_prob_dict(data_obj, band, quiet=True):

    out_dict = {'prob_'+band: sum(data_obj['weights_z'].values)}
    if quiet == False:
        print('prob sum in z is='+str(sum(data_obj['weights_z'].values)))
    return out_dict


def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def mags_of_percentile(cutoff, percentile_dict):
    if cutoff < 1.0:
        cutoff *= 100

    index = int(round(cutoff))
    return {band: percentile_dict['%s_cutoff' % band][index] for band in ['g', 'r', 'i', 'z']}


def make_output_csv(cutoffs, percentile_dict, outfile=None, return_df=False, write_answer=False, flt='', fraction=90.0, datadir='./'):

    out_data = [mags_of_percentile(cutoff, percentile_dict)
                for cutoff in cutoffs]
    out_df = pd.DataFrame(out_data)
    out_df['PERCENTILE'] = cutoffs
    if outfile:
        out_df.to_csv(outfile + '.csv', index=False)

    if write_answer:
        if fraction < 1.0:
            fraction *= 100
        # get closest index to fraction
        closest_index = np.argmin(
            np.abs(float(fraction) - out_df['PERCENTILE'].values))
        stream = open(datadir+'answer_%s.txt' % flt, 'w+')
        #stream = open('answer_%s.txt' %flt, 'w+')
        stream.write('%.2f' % out_df[flt].values[closest_index])
        stream.close()

    if return_df:
        return out_df


def make_plot_days(time, fracs, plotname="detection_distance.png", _title="", day=True, maglim=False):
    color_dict = {'i': 'brown', 'g': 'green', 'r': 'red', 'z': 'dimgray'}
    plt.figure()
    bands = ['g', 'r', 'i', 'z']

    if maglim == False:
        time_levels = np.arange(time[0], time[-1], 5)
        for i in range(0, len(bands)):
            fracs_interp = interp1d(time, fracs[i])
            plt.plot(time_levels, fracs_interp(time_levels), lw=2,
                     label=bands[i], color=color_dict[bands[i]])
    else:

        for i in range(0, len(bands)):
            time_levels = np.arange(time[i][0], time[i][-1], 5)
            fracs_interp = interp1d(time, fracs[i])
            plt.plot(time_levels[i], fracs_interp(
                time_levels[i]), lw=2, label=bands[i], color=color_dict[bands[i]])

    if day:
        plt.xlabel("Hours after merger", fontsize=14)  # days after merger
    else:
        plt.xlabel("Distance (MPC)", fontsize=14)  # days after merger

    if maglim == True:
        plt.xlabel("lim magnitude", fontsize=14)  # days after merger
    plt.ylabel("percent (at lim magnitude) ", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.ylim(0, 100)
    plt.legend(fontsize=14, loc='upper right', frameon=False)
    plt.title(_title, fontsize=14)
    plt.savefig(plotname)
    plt.close()


def make_plot(percentile_dict, title='', outfile=None, fraction=None, teff_type='moony'):
    plt.figure()

    percentile_levels = np.arange(len(percentile_dict['g_cutoff']))

    color_dict = {'i': 'brown', 'g': 'green', 'r': 'red', 'z': 'dimgray'}

    interps = {}

    if teff_type == 'moony':
        teffs = {'g': 0.05, 'r': 0.15, 'i': 0.45, 'z': 0.6}
    else:
        teffs = {'g': 0.7, 'r': 0.8, 'i': 0.7, 'z': 0.6}

    fractions_at_maglim = []
    for band in ['g', 'r', 'i', 'z']:
        m0 = get_m0(band, teff=teffs[band])
        plt.plot(percentile_dict['%s_cutoff' % band],
                 percentile_levels, lw=2, label=band, color=color_dict[band])
        plt.axvline(x=m0, color=color_dict[band], lw=1, ls=':')
        interps[band] = interp1d(
            percentile_dict['%s_cutoff' % band], percentile_levels)
        if m0 < percentile_dict['%s_cutoff' % band][0]:
            fractions_at_maglim.append(float(0.0))
        elif m0 > percentile_dict['%s_cutoff' % band][-1]:
            fractions_at_maglim.append(float(99.0))

        else:
            fractions_at_maglim.append(float(interps[band](m0)))

        if band == 'z':
            plt.text(19.1, 40, '90s 10$\sigma$ limit',
                     verticalalignment='center', fontsize=12)

            plt.xlabel("magnitude", fontsize=14)
            plt.ylabel("percent (<magnitude) ", fontsize=14)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

    if fraction:
        if fraction < 1.0:
            fraction *= 100

    plt.legend(fontsize=14, loc='upper left', frameon=False)
    plt.ylim(0, 100)
    plt.xlim(18, 26)
    plt.title(title, fontsize=14)

    if outfile:
        plt.savefig(outfile)

    plt.close()

    return fractions_at_maglim


def get_percentile_at_exp_time(exptime, band, percentile_dict, teff_type='notmoony', return_maglim=False):

    if teff_type == 'moony':
        teffs = {'g': 0.05, 'r': 0.15, 'i': 0.45, 'z': 0.6}
    else:
        teffs = {'g': 0.7, 'r': 0.8, 'i': 0.7, 'z': 0.6}
    m0 = get_m0(band)

    mag = m0 + 1.25 * np.log10(exptime*teffs[band] / 90.0)

    df = make_output_csv(np.linspace(0.0, 100.0, 101),
                         percentile_dict, outfile=None, return_df=True)

    index_of_closest_mag = np.argmin(np.abs(mag - df[band].values))
    percentile = df.iloc[index_of_closest_mag]['PERCENTILE']
    if return_maglim == False:
        return percentile
    else:
        #print ('mag lim in '+str(band)+' is '+str(mag))
        return percentile, mag


def get_exptime(m0, mag):
    return 90.0 * 10 ** ((mag - m0) / 1.25)


def get_m0(band, teff=1.0):
    out_dict = {'g': 23.4+(1.25*log10(teff)), 'r': 23.1+(1.25*log10(teff)), 'i': 22.5+(
        1.25*log10(teff)), 'z': 21.8+(1.25*log10(teff)), 'Y': 20.3+(1.25*log10(teff))}
    return out_dict[band]


def make_exptime_plot(percentile_dict, title='', outfile=None, teff=1.0, exposures=[90, 100, 110, 120, 130, 140, 150, 160, 170, 180], teff_type='notmoony'):
    percentile_levels = np.arange(len(percentile_dict['g_cutoff']))

    color_dict = {'i': 'brown', 'g': 'green', 'r': 'red', 'z': 'dimgray'}

    if teff_type == 'moony':
        teffs = {'g': 0.05, 'r': 0.15, 'i': 0.45, 'z': 0.6}
    else:
        teffs = {'g': 0.7, 'r': 0.8, 'i': 0.7, 'z': 0.6}


    fig, ax1 = plt.subplots()

    for band in ['g', 'r', 'i', 'z']:
        m0 = get_m0(band, teff=teffs[band])
        exptimes = get_exptime(m0, percentile_dict['%s_cutoff' % band])
        ax1.plot(exptimes*teffs[band], percentile_levels,
                 lw=2, label=band, color=color_dict[band])

    ax1.axvline(x=90.0, ls='--', lw=0.5)

    percentiles_at_90_sec = {band: get_percentile_at_exp_time(
        90.0, band, percentile_dict, teff_type=teff_type) for band in ['g', 'r', 'i', 'z']}
    bands_exp = ['g', 'r', 'i', 'z']
    fraction_atexposure = []
    for i in range(0, len(bands_exp)):
        fraction_atexposure.append([get_percentile_at_exp_time(
            exposure_time, bands_exp[i], percentile_dict, teff_type=teff_type) for exposure_time in exposures])
    ax1.text(105, 25, "90 sec in g = %.1f %%" % (
        percentiles_at_90_sec['g']), horizontalalignment='left', verticalalignment='center', fontsize=11)
    ax1.text(105, 20, "90 sec in r = %.1f %%" % (
        percentiles_at_90_sec['r']), horizontalalignment='left', verticalalignment='center', fontsize=11)
    ax1.text(105, 15, "90 sec in i = %.1f %%" % (
        percentiles_at_90_sec['i']), horizontalalignment='left', verticalalignment='center', fontsize=11)
    ax1.text(105, 10, "90 sec in z = %.1f %%" % (
        percentiles_at_90_sec['z']), horizontalalignment='left', verticalalignment='center', fontsize=11)

    ax1.set_xlabel("exptime * teff", fontsize=14)
    ax1.set_ylabel("percent (< 10$\sigma$ limiting mag)", fontsize=14)
    ax1.set_xlim(30, 300)
    ax1.set_ylim(0, 100)

    ax1.set_title(title, fontsize=14)
    ax1.grid()
    ax1.legend(fontsize=14, loc=4)

    if outfile:
        plt.savefig(outfile)

    plt.close()

    return fraction_atexposure


def get_detection_later(p_, p_later, start_time, delay, start_later):
    delay_days = float(delay)/(60*60*24)
    start_day = (float(start_time)/24.0)
    return np.interp(start_day+delay_days, [start_day, start_later], [p_, p_later])


if __name__ == '__main__':

    parser = ArgumentParser(__doc__)

    parser.add_argument('--input',
                      help='Skymap Location.')

    parser.add_argument('--output', default=None,
                      help='Output location for csv file.')
                      
    args = parser.parse_args()
    area_deg_fix = 40.0

if 1 == 1:
    teff_kind = 'moony'  # 'notmoony'
    options = {
        'distance': 150,
        'distance_err': 50,
        'time_delay': 60.0,
        'fraction': 90,
        'magplot_file': 'kn_mag_plot_clecio_loglan5.png',
        'expplot_file': 'kn_exp_plot_clecio_test_loglan5.png',
        'report_file': 'kn_report_clecio_test',
        'filter': None
    }

    distances = [50, 75, 150, 200, 250, 300]
    loglan = -1.0
    kntype = 'blue'  # 'blue'
    kn_weights = True
    w_type = "gaussian_narrow"  # 'gaussian'#

    event_path = args.input
    out_dir = args.output
    sufix = '_'+teff_kind+'_'+kntype+'_'

    event_list = glob.glob(join(event_path, '*.fits.gz'))
    map_mode = 1  # 1 for map 0 for dist
    area_covered_deep = [0.8, 0.7, 0.5, 0.3, 0.7,
                         0.5, 0.3, 0.5, 0.4, 0.3]  # [0.8,0.7]#
    m_exp_run = True  # FIXME just to call attention to this variable.
    m_exp = m_exp_run

    if m_exp_run == True:
        area_covered = [0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7]
    else:
        area_covered = [0.9, 0.85, 0.8, 0.75, 0.7]

    if map_mode == 1:
        second_loop = area_covered
        second_loop_legend = "Region Coverage"
        areas_cov_deep = area_covered_deep
        distances = [250.0 for i in range(0, len(area_covered))]
    else:
        second_loop = distances
        eventmap = ""
        event_list = [""]
        second_loop_legend = "Distance (MPC)"
        m_exp = False
        m_exp_run = False

    if m_exp == False:
        area_covered_deep = [0.0 for i in range(0, len(second_loop))]

    plot_dist_name = 'detection_distance_loglan5.png'
    plot_td_name = 'detection_time_delay_loglan5.png'
    time_delays = [12.0, 24.0, 36.0, 48.0, 60.0, 72.0, 84.0, 96.0]

    filters_comb = ['gg', 'gr', 'gi', 'gz', 'rg', 'rr', 'ri',
                    'rz', 'ig', 'ir', 'ii', 'iz', 'zg', 'zr', 'zi', 'zz']

    if m_exp == False:
        exposure_times_calc = [60.0, 90.0, 120.0,
                               200.0, 300.0, 600.0, 1200.0, 2400.0, 3600.0]
    else:
        exposure_times_calc = [60.0, 90.0, 120.0, 200.0, 300.0,
                               300.0, 600.0, 600.0, 1200.0, 1200, 2400.0, 3600.0, 300.0]

    exposure_times_calc_deep = [90.0, 120.0, 200.0, 300.0, 600.0,
                                1200.0, 1200.0, 2400.0, 2400.0, 3600.0, 3600.0, 5400.0, 5400.0]
    maglims_exp = np.zeros(
        (len(['g', 'r', 'i', 'z']), len(exposure_times_calc)))
    maglims_exp_l = np.zeros(
        (len(['g', 'r', 'i', 'z']), len(exposure_times_calc)))
    fractions_all = np.zeros((len(['g', 'r', 'i', 'z']), len(
        exposure_times_calc), len(second_loop), len(time_delays)),)
    fractions_all_later = np.zeros((len(['g', 'r', 'i', 'z']), len(
        exposure_times_calc), len(second_loop), len(time_delays)),)

day_delays = np.array(time_delays)/24.0
hours_per_night = 8.0
time_per_night_sec = hours_per_night*60*60  # 21600
mjd_correction = True





for e in event_list:
    if map_mode == 1:

        print(" ============================== \n ==============================\n")
        print(f"starting the run for file {e}")
        print(" ============================== \n ==============================\n")
        plot_name_ = basename(e)
        plot_name_ = plot_name_.rstrip(".fits.gz")
        plot_name_ = plot_name_.rstrip(".fits")
        plot_name_ = out_dir+plot_name_+sufix
    else:
        plot_name_ = out_dir+"GW_sim"+sufix

    if os.path.isfile(plot_name_+"_allconfig.csv"):
        print('This event event strategy already exists. Skipping')
        continue

    warning_area_deep = []
    area_deg_arr = np.zeros((len(time_delays), len(second_loop)))
    area_deg_arr_deep = np.zeros((len(time_delays), len(second_loop)))
    for j in range(0, len(time_delays)):

        for i in range(0, len(second_loop)):  # time_delays
            options['distance'] = distances[i]
            options['time_delay'] = time_delays[j]
            if map_mode == 1:
                _map_file = e
            else:
                _map_file = ""  # e
            m_exp = m_exp_run
            kn_calc = KNCalc(
                e,
                float(options['distance']),
                float(options['distance_err']),
                float(options['time_delay']),
                float(loglan),
                kn_weight_type=w_type,
                plot=False,
                use_map=True,
                area_covered=second_loop[i],
                reduce_mapresolution=True,
                set_mjd_correction=mjd_correction,
                m_exp_kncalc=m_exp,
                deep_coverage=area_covered_deep[i]
            )

            if (kn_calc.sw_mexp == True) and (m_exp_run == True):
                warning_area_deep.append(area_covered_deep[i])

            if kn_calc.Flag == 0:
                print(' I will pass')
                continue
            m_exp = kn_calc.mult_exp

            if map_mode == 1:
                area_deg_arr[j][i] = kn_calc.area_deg
                if m_exp == True:
                    area_deg_arr_deep[j][i] = kn_calc.area_deg_deep

            else:
                area_deg_arr[j][i] = area_deg_fix

            if m_exp == False:
                percentile_dict = calc_mag_fractions(
                    kn_calc.template_df_full, use_knmodel_weights=kn_weights, kn_type=kntype,  kn_weight_type=w_type)
                percentile_dict_later = calc_mag_fractions(
                    kn_calc.template_df_later, use_knmodel_weights=kn_weights, kn_type=kntype,  kn_weight_type=w_type)
            else:
                percentile_dict, percentile_dict_deep = calc_mag_fractions(
                    kn_calc.template_df_full, use_knmodel_weights=kn_weights, kn_type=kntype,  kn_weight_type=w_type, m_exp=m_exp)
                percentile_dict_later, percentile_dict_later_deep = calc_mag_fractions(
                    kn_calc.template_df_later, use_knmodel_weights=kn_weights, kn_type=kntype,  kn_weight_type=w_type, m_exp=m_exp)

            cutoffs = mags_of_percentile(
                float(options['fraction']), percentile_dict)
            cutoff_dict = {'%s_mag' % k: v for k, v in cutoffs.items()}
            for band in ['g', 'r', 'i', 'z']:
                cutoff_dict['%s_magerr' % band] = 0.00

            plot_title = "%s +/- %s Mpc  -- %.2f Days After Merger" % (
                options['distance'], options['distance_err'], float(options['time_delay']) / 24.0)
            # if options['magplot_file']:
            fractions_maglim = make_plot(
                percentile_dict, title=plot_title, outfile=options['magplot_file'], fraction=options['fraction'], teff_type=teff_kind)
            fractions_maglim_later = make_plot(percentile_dict_later, title=plot_title,
                                               outfile=options['magplot_file'], fraction=options['fraction'], teff_type=teff_kind)
            if m_exp == True:
                fractions_maglim_deep = make_plot(
                    percentile_dict_deep, title=plot_title, outfile=options['magplot_file'], fraction=options['fraction'], teff_type=teff_kind)
                fractions_maglim_later_deep = make_plot(
                    percentile_dict_later_deep, title=plot_title, outfile=options['magplot_file'], fraction=options['fraction'], teff_type=teff_kind)

            for l in range(0, len(exposure_times_calc)):
                expt = exposure_times_calc[l]

                frac_at_exptg, maglimg = get_percentile_at_exp_time(
                    expt, band='g', percentile_dict=percentile_dict, teff_type=teff_kind, return_maglim=True)
                frac_at_exptr, maglimr = get_percentile_at_exp_time(
                    expt, band='r', percentile_dict=percentile_dict, teff_type=teff_kind, return_maglim=True)
                frac_at_expti, maglimi = get_percentile_at_exp_time(
                    expt, band='i', percentile_dict=percentile_dict, teff_type=teff_kind, return_maglim=True)
                frac_at_exptz, maglimz = get_percentile_at_exp_time(
                    expt, band='z', percentile_dict=percentile_dict, teff_type=teff_kind, return_maglim=True)

                frac_at_exptg_l, maglimg_l = get_percentile_at_exp_time(
                    expt, band='g', percentile_dict=percentile_dict_later, teff_type=teff_kind, return_maglim=True)
                frac_at_exptr_l, maglimr_l = get_percentile_at_exp_time(
                    expt, band='r', percentile_dict=percentile_dict_later, teff_type=teff_kind, return_maglim=True)
                frac_at_expti_l, maglimi_l = get_percentile_at_exp_time(
                    expt, band='i', percentile_dict=percentile_dict_later, teff_type=teff_kind, return_maglim=True)
                frac_at_exptz_l, maglimz_l = get_percentile_at_exp_time(
                    expt, band='z', percentile_dict=percentile_dict_later, teff_type=teff_kind, return_maglim=True)

                if m_exp == True:
                    expt_deep = exposure_times_calc_deep[l]
                    frac_at_exptg_deep, maglimg_deep = get_percentile_at_exp_time(
                        expt_deep, band='g', percentile_dict=percentile_dict_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_exptr_deep, maglimr_deep = get_percentile_at_exp_time(
                        expt_deep, band='r', percentile_dict=percentile_dict_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_expti_deep, maglimi_deep = get_percentile_at_exp_time(
                        expt_deep, band='i', percentile_dict=percentile_dict_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_exptz_deep, maglimz_deep = get_percentile_at_exp_time(
                        expt_deep, band='z', percentile_dict=percentile_dict_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_exptg_l_deep, maglimg_l_deep = get_percentile_at_exp_time(
                        expt_deep, band='g', percentile_dict=percentile_dict_later_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_exptr_l_deep, maglimr_l_deep = get_percentile_at_exp_time(
                        expt_deep, band='r', percentile_dict=percentile_dict_later_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_expti_l_deep, maglimi_l_deep = get_percentile_at_exp_time(
                        expt_deep, band='i', percentile_dict=percentile_dict_later_deep, teff_type=teff_kind, return_maglim=True)
                    frac_at_exptz_l_deep, maglimz_l_deep = get_percentile_at_exp_time(
                        expt_deep, band='z', percentile_dict=percentile_dict_later_deep, teff_type=teff_kind, return_maglim=True)

                maglims_exp[0][l] = round(maglimg, 1)
                maglims_exp[1][l] = round(maglimr, 1)
                maglims_exp[2][l] = round(maglimi, 1)
                maglims_exp[3][l] = round(maglimz, 1)

                maglims_exp_l[0][l] = round(maglimg_l, 1)
                maglims_exp_l[1][l] = round(maglimr_l, 1)
                maglims_exp_l[2][l] = round(maglimi_l, 1)
                maglims_exp_l[3][l] = round(maglimz_l, 1)

                if m_exp == False:

                    fractions_all[0][l][i][j] = frac_at_exptg
                    fractions_all[1][l][i][j] = frac_at_exptr
                    fractions_all[2][l][i][j] = frac_at_expti
                    fractions_all[3][l][i][j] = frac_at_exptz

                    fractions_all_later[0][l][i][j] = frac_at_exptg_l
                    fractions_all_later[1][l][i][j] = frac_at_exptr_l
                    fractions_all_later[2][l][i][j] = frac_at_expti_l
                    fractions_all_later[3][l][i][j] = frac_at_exptz_l
                else:
                    fractions_all[0][l][i][j] = frac_at_exptg + \
                        frac_at_exptg_deep
                    fractions_all[1][l][i][j] = frac_at_exptr + \
                        frac_at_exptr_deep
                    fractions_all[2][l][i][j] = frac_at_expti + \
                        frac_at_expti_deep
                    fractions_all[3][l][i][j] = frac_at_exptz + \
                        frac_at_exptz_deep

                    fractions_all_later[0][l][i][j] = frac_at_exptg_l + \
                        frac_at_exptg_l_deep
                    fractions_all_later[1][l][i][j] = frac_at_exptr_l + \
                        frac_at_exptr_l_deep
                    fractions_all_later[2][l][i][j] = frac_at_expti_l + \
                        frac_at_expti_l_deep
                    fractions_all_later[3][l][i][j] = frac_at_exptz_l + \
                        frac_at_exptz_l_deep

    if kn_calc.Flag == 0:
        print(' I will pass to the next event')
        continue

    if map_mode == 1:
        cadence_matrix = np.zeros((len(filters_comb), len(second_loop), len(
            exposure_times_calc)*len(exposure_times_calc)))  # ,len(time_delays)*len(time_delays)))
        cadence_matrix = cadence_matrix.tolist()
        bands_plot = ['g', 'r', 'i', 'z']
        for m in range(0, len(cadence_matrix)):
            for i in range(0, len(cadence_matrix[0])):
                for k in range(0, len(cadence_matrix[0][0])):
                    cadence_matrix[m][i][k] = []
        exp_comblegend = []
        day_delays_comb = []
        day_delays_comb_number = []
        exp_comb_number = []
        exp_comb_number_deep = []

        for j_1 in range(0, len(day_delays)):
            for j_2 in range(0, len(day_delays)):
                if j_1 <= j_2:
                    day_delays_comb.append(
                        str(day_delays[j_1])+'+'+str(day_delays[j_2]))
                    day_delays_comb_number.append(
                        [day_delays[j_1], day_delays[j_2]])
        for k_1 in range(0, len(exposure_times_calc)):
            for k_2 in range(0, len(exposure_times_calc)):

                exp_comb_number.append(
                    [exposure_times_calc[k_1], exposure_times_calc[k_2]])
                if (m_exp_run == True):  # or (kn_calc.sw_mexp):
                    exp_comb_number_deep.append(
                        [exposure_times_calc_deep[k_1], exposure_times_calc_deep[k_2]])  # FIXME 20220322
                    exp_comblegend.append(str(round(exposure_times_calc[k_1]/60.0, 1))+'D'+str(round(exposure_times_calc_deep[k_1]/60.0, 1))+'+'+str(
                        round(exposure_times_calc[k_2]/60.0, 1))+'D'+str(round(exposure_times_calc_deep[k_2]/60.0, 1)))
                else:
                    exp_comblegend.append(str(round(
                        exposure_times_calc[k_1]/60.0, 1))+'+'+str(round(exposure_times_calc[k_2]/60.0, 1)))

        m = 0
        k = 0
        j = 0
        telescope_time1_all = []
        prob1_all = []
        prob2_all = []
        telescope_time2_all = []
        area_deg_all = []
        area_deg_all_deep = []
        fov = 3.0
        for m_1 in range(0, len(bands_plot)):  # in range(0,4)
            for m_2 in range(0, len(bands_plot)):  # in range(0, 4)
                for i in range(0, len(second_loop)):  # in range(0, 10) if m_exp_run == True in range(0, 5) otherwise
                    k = 0
                    for k_1 in range(0, len(exposure_times_calc)):
                        for k_2 in range(0, len(exposure_times_calc)):
                            j = 0
                            for j_1 in range(0, len(time_delays)):
                                for j_2 in range(0, len(time_delays)):
                                    if j_1 <= j_2:
                                        if (m_exp_run == False) or (area_covered_deep[i] in warning_area_deep):
                                            if area_covered_deep[i] in warning_area_deep:
                                                print('area deep warning for ', str(
                                                    area_covered_deep[i]))
                                                area_deg_all_deep.append(0.0)
                                            telescope_time1 = telescope_time(
                                                exptime=exposure_times_calc[k_1], area_deg=area_deg_arr[j_1][i], field_of_view_area=fov)
                                            telescope_time2 = telescope_time(
                                                exptime=exposure_times_calc[k_2], area_deg=area_deg_arr[j_2][i], field_of_view_area=fov)
                                        else:
                                            # FIXME is K_1 20220322
                                            area_deg_all_deep.append(
                                                area_deg_arr_deep[j_1][i])
                                            telescope_time1 = telescope_time(exptime=[exposure_times_calc[k_1], exposure_times_calc_deep[k_1]], area_deg=[
                                                                             area_deg_arr[j_1][i], area_deg_arr_deep[j_1][i]], field_of_view_area=fov, m_exp=m_exp_run)
                                            telescope_time2 = telescope_time(exptime=[exposure_times_calc[k_2], exposure_times_calc_deep[k_2]], area_deg=[
                                                                             area_deg_arr[j_2][i], area_deg_arr_deep[j_2][i]], field_of_view_area=fov, m_exp=m_exp_run)
                                        area_deg_all.append(
                                            area_deg_arr[j_1][i])
                                        telescope_time1_all.append(
                                            telescope_time1)
                                        telescope_time2_all.append(
                                            telescope_time2)
                                        if (telescope_time1 > time_per_night_sec):
                                            cadence_matrix[m][i][k].append(0)
                                            prob1_all.append(0)
                                            prob2_all.append(0)
                                        elif (telescope_time2 > time_per_night_sec):
                                            cadence_matrix[m][i][k].append(0)
                                            prob1_all.append(0)
                                            prob2_all.append(0)

                                        else:
                                            if j_1 == j_2:
                                                if (telescope_time1+telescope_time2) > time_per_night_sec:
                                                    cadence_matrix[m][i][k].append(
                                                        0)
                                                    prob1_all.append(0)
                                                    prob2_all.append(0)
                                                else:
                                                    later_second_detection = get_detection_later(
                                                        fractions_all[m_2][k_2][i][j_2], fractions_all_later[m_2][k_2][i][j_2], time_delays[j_1], telescope_time1, kn_calc.delta_mjd_later)
                                                    cadence_matrix[m][i][k].append(
                                                        100*(float(fractions_all[m_1][k_1][i][j_1])/100.0)*(float(later_second_detection)/100.0))
                                                    prob1_all.append(
                                                        float(fractions_all[m_1][k_1][i][j_1]))
                                                    prob2_all.append(
                                                        later_second_detection)
                                            else:
                                                cadence_matrix[m][i][k].append(
                                                    100*(float(fractions_all[m_1][k_1][i][j_1])/100.0)*(float(fractions_all[m_2][k_2][i][j_2])/100.0))
                                                prob1_all.append(
                                                    float(fractions_all[m_1][k_1][i][j_1]))
                                                prob2_all.append(
                                                    float(fractions_all[m_2][k_2][i][j_2]))

                                    j = j+1
                            k = k+1
                m = m+1

        for k in range(0, len(filters_comb)):

            fig, ax = plt.subplots(figsize=(10, 15))

            ax = sns.heatmap(cadence_matrix[k][-1], annot=True, xticklabels=day_delays_comb,
                             yticklabels=exp_comblegend, linewidths=.5, vmin=0, vmax=100, ax=ax)

            # xticklabels
            plt.xlabel("Observing night (Days after merger)",
                       fontsize=14)  # days after merger

            plt.ylabel("Exposure times", fontsize=9)  # days after merger
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)

            plt.savefig(plot_name_+filters_comb[k]+"_cadence.png")
            plt.close()

        cadence_matrix = np.array(cadence_matrix)

        best_probs = []
        prob_labels = []
        best_filters = []
        best_exposure01 = []
        best_exposure02 = []
        best_time_delays01 = []
        best_time_delays02 = []
        best_area = []

        telescope_btim01 = []
        telescope_btim02 = []
        best_areadegree = []

        probs_all = []
        area_all = []
        area_all_deep = []
        filters_all = []
        exposure01_all = []
        exposure02_all = []
        exposure01_all_deep = []
        exposure02_all_deep = []
        time_delays01_all = []
        time_delays02_all = []

        cadence_full = cadence_matrix.copy()
        for m in range(0, len(cadence_matrix)):
            for i in range(0, len(cadence_matrix[m])):
                for k in range(0, len(cadence_matrix[m][i])):
                    for j in range(0, len(cadence_matrix[m][i][k])):
                        probs_all.append(cadence_matrix[m][i][k][j])
                        area_all.append(second_loop[i])
                        # area_all_deg.append(areas[i])
                        if m_exp_run == True:
                            if area_covered_deep[i] in warning_area_deep:
                                # area_all_deg_deep.append(0.0)
                                area_all_deep.append(area_covered_deep[i])
                                exposure01_all_deep.append(0.0)
                                exposure02_all_deep.append(0.0)

                            else:
                                # area_all_deg_deep.append(areas_deep[i])
                                area_all_deep.append(areas_cov_deep[i])
                                exposure01_all_deep.append(
                                    exp_comb_number_deep[k][0])
                                exposure02_all_deep.append(
                                    exp_comb_number_deep[k][1])

                        exposure01_all.append(exp_comb_number[k][0])
                        exposure02_all.append(exp_comb_number[k][1])

                        time_delays01_all.append(day_delays_comb_number[j][0])
                        time_delays02_all.append(day_delays_comb_number[j][1])
                        filters_all.append(filters_comb[m])

        f = open(plot_name_+"_allconfig.csv", 'w')
        f.close()
        f = open(plot_name_+"_allconfig.csv", 'a')

        f.write('# This is the '+kntype +
                ' component in a '+teff_kind+' night \n')
        if m_exp_run == False:
            pd.DataFrame({"Detection Probability": probs_all,
                          second_loop_legend: area_all,
                          "Filter_comb": filters_all,
                          "Exposure01": exposure01_all,
                          "Exposure02": exposure02_all,
                          "Observation01": time_delays01_all,
                          "Observation02": time_delays02_all,
                          "Telescope_time01": telescope_time1_all,
                          "Telescope_time02": telescope_time2_all,
                          "Region_coverage_deg": area_deg_all,
                          "Deprob1": prob1_all,
                          "Detprob2": prob2_all}).to_csv(f,  index=False)
        else:

            pd.DataFrame({"Detection Probability": probs_all,
                          second_loop_legend: area_all,
                          second_loop_legend+"_deep": area_all_deep,
                          "Filter_comb": filters_all,
                          "Exposure01": exposure01_all,
                          "Exposure01_deep": exposure01_all_deep,
                          "Exposure02": exposure02_all,
                          "Exposure02_deep": exposure02_all_deep,
                         "Observation01": time_delays01_all,
                         "Observation02": time_delays02_all,
                         "Telescope_time01": telescope_time1_all,
                         "Telescope_time02": telescope_time2_all,
                         "Region_coverage_deg": area_deg_all,
                         "Region_coverage_deg_deep": area_deg_all_deep,
                         "Deprob1": prob1_all,
                         "Detprob2": prob2_all}).to_csv(f,  index=False)

        f.close()

    if map_mode == 0:
        probability_coverage = 0.90
        bands_plot = ['g', 'r', 'i', 'z']
        fractions_all_tp = np.transpose(fractions_all, (0, 1, 3, 2))
        for m in range(0, len(bands_plot)):
            for k in range(0, len(exposure_times_calc)):
                ax1 = sns.heatmap(fractions_all_tp[m][k]*probability_coverage, annot=True,
                                  xticklabels=second_loop, yticklabels=day_delays, linewidths=.5, vmin=0, vmax=100)
                # days after merger
                plt.ylabel("Days after merger", fontsize=14)
                # days after merger
                plt.xlabel(second_loop_legend, fontsize=14)
                plt.xticks(fontsize=10)
                plt.yticks(fontsize=10)
                plt.title(teff_kind+" KN type: "+kntype+" " +
                          str(bands_plot[m])+" exp: "+str(exposure_times_calc[k]), fontsize=10)
                plt.savefig(
                    plot_name_+bands_plot[m]+"exp"+str(exposure_times_calc[k])+"_distvsdelay.png")
                plt.clf()
                plt.close()

        fractions_dist = np.transpose(fractions_all, (1, 3, 0, 2))
        for k in range(0, len(exposure_times_calc)):
            for j in range(0, len(time_delays)):
                plt.clf()
                fig = plt.figure()
                ax = fig.add_subplot(111)

                for m in range(0, len(bands_plot)):
                    ax.plot(second_loop, fractions_dist[k][j][m]*probability_coverage, color=tableau20[int(
                        2*m)], label=str(bands_plot[m]), ls='-.', linewidth=1, alpha=0.7)
                plt.xlabel(second_loop_legend, fontsize=14)
                plt.ylabel(" Detection Probability", fontsize=14)
                plt.legend(loc='lower left', fontsize=9)
                plt.title(teff_kind+" type: "+kntype+" Days:" +
                          str(day_delays[j])+" exp: "+str(exposure_times_calc[k]), fontsize=10)
                plt.savefig(
                    plot_name_+str(time_delays[j])+"exp"+str(exposure_times_calc[k])+"_dist.png")
