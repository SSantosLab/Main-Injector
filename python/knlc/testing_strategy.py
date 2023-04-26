from scipy.stats import norm, uniform
import numpy as np
import os
import time
import cProfile

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


def get_model_weights(kn_weight_type="uniform", kn_type='red', info_file=os.path.join(os.getenv("DESGW_DIR"), "knlc", "knsed_info.txt"), kn_weight_sigma=1.0):

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

if __name__ == '__main__':
    start = time.time()
    model, kn_ids = get_model_weights()
    end = time.time() - start
    print(f'{end:.2f} seconds!')
    print(model)