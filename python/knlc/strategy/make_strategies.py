import matplotlib
from matplotlib import colors
matplotlib.use('agg')
import fits_cat as fc
import os
import numpy as np
import seaborn as sns
from photo_lib import *
import pandas as pd
from scipy.interpolate import interpn
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import healpy as hp
from mpl_toolkits.mplot3d import Axes3D
import sys
from math import sqrt
from argparse import ArgumentParser
import glob
import pandas as pd
import logging as log
import traceback


def get_map_distance(map: str,
                     savedir: str =''):
    """Gets the map distance"""

    dtypes = [np.float64 for i in range(4)]
    pb, distmu, distsigma, distnorm = hp.read_map(map,
                                                  field=range(4),
                                                  dtype=dtypes)
    

    check_distmu = np.logical_not(np.isinf(distmu))
    pb_check = pb[check_distmu]
    distsigma_check = distsigma[check_distmu]
    distmu_check = distmu[check_distmu]

    check_distsigma = np.logical_not(np.isinf(distsigma_check))
    pb_check = pb_check[check_distsigma]
    distsigma_check = distsigma_check[check_distsigma]
    distmu_check = distsigma_check[check_distsigma]

    distmu_check_average= np.average(distmu_check,weights=pb_check)
    distsigma_check_average= np.average(distsigma_check,weights=pb_check)

    idx_sort_up = np.argsort(pb)[::-1]
    NSIDE = hp.npix2nside(pb.shape[0])
    resolution = hp.nside2pixarea(NSIDE, degrees=True)
    id_c = 0
    sum_full = 0
    id_full = 0
    area_max = np.sum(pb)

    while (sum_full < 0.9) and (id_full < len(idx_sort_up)):
        this_idx = idx_sort_up[id_full]
        sum_full = sum_full + pb[this_idx]
        id_full = id_full+1
        total_area = id_full*resolution

    if savedir!='':
        np.save(savedir,[distmu_check_average,distsigma_check_average])
    return distmu_check_average,distsigma_check_average





if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('--map-path',
                   help='Path to fits_flattened skymaps.')
    
    p.add_argument('--strategy-path',
                   help='Path to output csv files from strategy code.')
    
    p.add_argument('--model-label',
                   help='Kilonova model used for strategy. Default is bright and blue.',
                   default='Bright and Blue')

    p.add_argument('--constrain-time',
                   action='store_true',
                   help='make constrain time plots.',
                   default=True)
    p.add_argument('--no-plots',
                   action='store_true',
                   help='If set, won\'t make plots',
                   default=False)

    args = p.parse_args()
    map_path = args.map_path
    map_path_info = map_path
    strategy_dir = args.strategy_path
    constrain_time = args.constrain_time
    no_plots = args.no_plots
    global model_label
    model_label=args.model_label
    event_list_all=glob.glob(strategy_dir+'*_allconfig.csv')
    evaluate_strategies=True
    plot_individual_strategies=True
    m_exp=True

    nsns_limit=190
    nsbh_limit=330
    gw_o3_nsns_dl=[40,160] 
    nsns_names=['GW170817 BNS','GW190425 BNS']
    gw_o3_nsbh_dl=[240,290]
    nsbh_names=['GW190814 NSBH', 'GW200115 NSBH']
    nsns_markers=["P","v"]
    nsbh_markers=["s","o"]

    if m_exp == False:
        ares_depth_bins=[
            np.array([0.65,0.7,0.75,0.8,0.85,0.9])+0.025,
            [45.0,75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0]
        ]

        xdvsw=[0.7,0.75,0.8,0.85,0.9]
        ydvsw=[60,90,120,200,300,600,1200,2400,3600]

    if m_exp==True:
        dualexp_depth_bins=[
            [75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0,6825.0],
            [45.0,75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0]
        ]

        ares_depth_bins=[
            np.array([0.65,0.75,0.85,0.95]),
            [45.0,75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0]
        ]

        ares_depth_bins_deep=[
            np.array([0.25,0.35,0.55,0.75,0.85]),
            [75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0,6825.0]
        ]

        area_area_deep_bins=[
            np.array([0.65,0.75,0.85,0.95]),
            np.array([0.25,0.35,0.55,0.75,0.85])
        ]

        x_area_deep=[0.3,0.5,0.7,0.8]
        y_area=[0.7,0.8,0.9]
        y_depth_deep=[1.5,2.0,3.3,5.0,10.0,20.0,40.0,60.0,90.0]
        depth_outer=[1.0,1.5,2.0,3.3,5.0,10.0,20.0,40.0,60.0]

        xdvsw=[0.7,0.8,0.9]
        ydvsw=[60,90,120,200,300,600,1200,2400,3600]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    low_budget=True #mbnmo_timett10_bright  mbnmo_timett10_bright
    sufix='mbnmo_timett10_bright_final_reader'#'mrnmo_timett5_bright'#'mbnmo_narrow_timett5_bright'#'mbnmo_timett10_bright'#'mbnmo_timett10_bright_final_reader'#'mrnmo_timett5_bright'#'mbnmo_timett10_bright'#'mrnmo_timett5'#'mbnmo_time_20220919_v2_lowttbudget_consistency_testTTuniform'#'mbnmo_timett10_bright'##'mbnmo_time_bbtt5'#'mrnmo_time_20220629_v2_lowttbudget_tt15'#'mbnmo_time_20220629_v2_lowttbudget_tt5_narrow'#'mbnmo_time_20220503_v2'#'mbnmo_time_20220503_v2_narrow_vg'#'mbnmo_time_20220503_v2'#'mbnmo_time_20220503_v2_narrow_vg'#'mbnmo_time_20220503_v2'#20220323'#'srnmo_time_'#'mbnmo_time'#'sbnmo_time_'#'mrnmo_time_'#'srnmo_time_'#'_mexp'#'mbnmo_time_'


    _2d_map_dist=True
    use_ref_mult=True #false if you dont want the multiple exposure catalog to include the reference strategy

    #========= in case you are running multiple exposure times and want to plot the reference strategy together
    sufix_nomulti_strategy='_notmoony_blue__allconfig.csv'
    # path to not notmoony strategy
    path_nomulti_Strategy = strategy_dir

    #========================================================== bright plot

    plot_bright_night_strategy=False
    plot_bright_ref_st=False
    sufix_bright_strategy='_notmoony_blue__allconfig.csv'
    path_bright_Strategy=strategy_dir
    path_bright_ref_Strategy=''
    sufix_bright_ref_strategy='_notmoony_blue__allconfig.csv'
    #========================================

    strategy_csv_cat=strategy_dir+sufix+'.csv'

    old_strategy_arr=[[90.0,90.0,'gi',0.5,1.5,0.90],
                      [90.0,90.0,'zz',0.5,1.5,0.90],
                      [90.0,90.0,'zi',0.5,1.5,0.90],
                      [90.0,90.0,'iz',0.5,1.5,0.90],
                      [90.0,90.0,'ii',0.5,1.5,0.90]]

    ltt_config=[0.05,0.10,0.15] # ltt_config[1]changed for a test
    TT_max=[]
    TT_low=[]
    obs1=[]
    obs2=[]
    exp1=[]
    exp2=[]
    TT=[]
    TT1=[]
    TT2=[]
    probdet1=[]
    probdet2=[]
    strategy_type=[]
    bands=[]
    distance_all=[]
    distance_all_err=[]
    area_deg90_all=[]
    allstrategy_prob=[]
    prob_area=[]
    prob_area_deg=[]

    if m_exp==True:
        exp1_deep=[]
        exp2_deep=[]
        prob_area_deep=[]
        prob_area_deg_deep=[]  
        
    event_names=[]
    probs_low_tt=[]
    probs_top_=[]
    probs_old=[]

    if evaluate_strategies == False:
        event_list_all=[]

    for i in range(0,len(event_list_all[:2])):
        map_file_aux = os.path.basename(event_list_all[i])
        map_id = map_file_aux.split('_')[0]
        map_file = map_id+'.fits.gz'

        print(f"this is event {map_id}")
        try: 
            print(f"map info: {map_path_info+map_id}.npy")
            distance,distance_err=np.load(map_path_info+map_id+'.npy')
        except:
            print("Did not find npy with event info. "+\
                  "So, I am Calculating Map distance for "+map_file)
            distance,distance_err=get_map_distance(map_path+map_file,
                                                   savedir=map_path_info+map_id+'.npy')
            distance,distance_err=np.load(map_path_info+map_id+'.npy')

        #This is not working because area_deg_info
        #is splitted between outer and inner

        if _2d_map_dist==True:
            try:
                print('Loading informaation for 90\% and 50\% credible region')
                area_deg_info, area_deg_deep_info, resolution = np.load(map_path+f'/weights/{map_id}ac0.9ad0.5info.npy')
                use_map_info=map_path+"lowres/"+map_id+"_mapinfo.npy"
            except Exception as e:
                # continue
                print("Exception found!")
                print(e)
                print(traceback.format_exc())
        # area_deg_info = -1
        all_df=pd.read_csv(event_list_all[i],comment='#')

        print('Outer 0.9 in sq-deg:', area_deg_info)
        print('Inner 0.5 in sq-deg:', area_deg_deep_info)

        if constrain_time:
            print('Constraining time')
            observation01 = all_df.Observation01.values
            observation02 = all_df.Observation02.values
            observation02 - observation01
            strategy_time_delays = observation02 - observation01
            all_df=all_df[np.logical_or(strategy_time_delays > 0.6,
                                        strategy_time_delays < 0.4) ]

        if (m_exp==True) and (use_ref_mult==True):
            file_old=map_file_aux.split('_')[0]+sufix_nomulti_strategy
            file_old=path_nomulti_Strategy+file_old
            try:
                print("Opening file ",file_old)
                df_all_old=pd.read_csv(file_old,comment='#')
            except:
                print("Single exposure file not found: ",file_old)
                continue   
            if constrain_time:
                df_all_old_observation01 = df_all_old['Observation01'].values
                df_all_old_observation02 = df_all_old['Observation02'].values
                strategy_time_delays_nomulti = df_all_old_observation02 -\
                                                df_all_old_observation01
                
                df_all_old = df_all_old[
                    np.logical_or(strategy_time_delays_nomulti > 0.6,
                                  strategy_time_delays_nomulti < 0.4)
                ]

            telescope_time_01 = df_all_old["Telescope_time01"].values
            telescope_time_02 = df_all_old["Telescope_time02"].values
            total_telescope_time_aux_old = telescope_time_01 +\
                                            telescope_time_02
            df_all_old=df_all_old[total_telescope_time_aux_old > 0.0]
        else:
            df_all_old=all_df

        if plot_bright_night_strategy:
            file_bright=map_id+sufix_bright_strategy
            file_bright=path_bright_Strategy+file_bright
            print("Plot bright night strategy")
            print(file_bright)

            if plot_bright_ref_st:
                file_bright_ref=map_id+sufix_bright_ref_strategy
                file_bright_ref=path_bright_ref_Strategy+file_bright_ref
                try:
                    print("Opening file ",file_bright_ref)
                    df_all_bright_ref=pd.read_csv(file_bright_ref,comment='#')
                except:
                    print("Bright night strategy reference not file not found: ",file_bright_ref)
                    continue   

            try:
                df_all_bright=pd.read_csv(file_bright,comment='#')
            except:
                print("Bright night strategy not file not found: ",file_bright)
                continue

            if constrain_time:
                df_all_bright_obs01 = df_all_bright["Observation01"].values
                df_all_bright_obs02 = df_all_bright["Observation02"].values
                strategy_time_delays_bright = df_all_bright_obs02 -\
                                              df_all_bright_obs01
                df_all_bright=df_all_bright[
                    np.logical_or(strategy_time_delays_bright > 0.6,
                    strategy_time_delays_bright < 0.4)
                ]

                if plot_bright_ref_st:
                    df_all_bright_ref_obs01 = df_all_bright_ref["Observation01"].values
                    df_all_bright_ref_obs02 = df_all_bright_ref["Observation02"].values
                    strategy_time_delays_bright_ref = df_all_bright_ref_obs02 -\
                                                      df_all_bright_ref_obs01
                    df_all_bright_ref = df_all_bright_ref[
                        np.logical_or(strategy_time_delays_bright_ref > 0.6,
                        strategy_time_delays_bright_ref < 0.4)
                    ]
            
            if plot_bright_ref_st:
            # THIS SHOULDNT BE NECESSARY. CHECK IT.
                df_all_bright_ref_tt01 = df_all_bright_ref["Telescope_time01"].values
                df_all_bright_ref_tt02 = df_all_bright_ref["Telescope_time02"].values
                total_telescope_time_aux_bright_ref = df_all_bright_ref_tt01 +\
                                                      df_all_bright_ref_tt02
                df_all_bright_ref = df_all_bright_ref[
                    total_telescope_time_aux_bright_ref > 0.0
                ]
            
            df_all_bright_tt01 = df_all_bright["Telescope_time01"].values
            df_all_bright_tt02 = df_all_bright["Telescope_time02"].values
            total_telescope_time_aux_bright = df_all_bright_tt01 +\
                                              df_all_bright_tt02
            df_all_bright=df_all_bright[
                total_telescope_time_aux_bright > 0.0
            ]

        df_all_tt01 = all_df["Telescope_time01"].values
        df_all_tt02 = all_df["Telescope_time02"].values
        total_telescope_time= df_all_tt01 + df_all_tt02
        all_df = all_df[total_telescope_time > 0.0]
        prob_all = all_df["Detection Probability"].values
        prob_top_test=max(prob_all)-0.01

        all_df_top = (all_df[all_df["Detection Probability"]>prob_top_test]
                      .copy()
                      .reset_index(drop=True)
        )
        observation1_all=all_df_top["Observation01"].values
        observation2_all=all_df_top["Observation02"].values
        exposure1_all=all_df_top["Exposure01"].values
        exposure2_all=all_df_top["Exposure02"].values
        region_all=all_df_top["Region Coverage"].values
        region_all_deg=all_df_top["Region_coverage_deg"].values
        probdet1_all=all_df_top["Deprob1"].values#top["Detection Probability"].values
        probdet2_all=all_df_top["Detprob2"].values
        telescope_time1_all=all_df_top["Telescope_time01"].values
        telescope_time2_all=all_df_top["Telescope_time02"].values
        bands_all=all_df_top["Filter_comb"].values
        
        if m_exp:
            exp1_deep_all=all_df_top["Exposure01_deep"].values
            exp2_deep_all=all_df_top["Exposure02_deep"].values
            prob_area_deg_all=all_df_top["Region Coverage_deep"].values
            prob_area_deg_deep_all=all_df_top["Region_coverage_deg_deep"].values
        prob_top=max(all_df["Detection Probability"].values)
        if prob_top==0:
            print("undetected event for any strategy -- skipping")
            continue

        all_df_top_tt01 = all_df_top["Telescope_time01"].values
        all_df_top_tt02 = all_df_top["Telescope_time02"].values
        total_telescope_time_top = all_df_top_tt01 + all_df_top_tt02
        
        min_exposure_idx = np.argmin(exposure1_all)
        tt_maxprob=total_telescope_time_top[min_exposure_idx]
        obs1_maxprob=observation1_all[min_exposure_idx]
        obs2_maxprob=observation2_all[min_exposure_idx]
        exposure1_maxprob=exposure1_all[min_exposure_idx]
        exposure2_maxprob=exposure2_all[min_exposure_idx]
        area_maxprob=region_all[min_exposure_idx]
        area_maxprob_deg=region_all_deg[min_exposure_idx]
        bands_maxprob=bands_all[min_exposure_idx]
        probdet1_maxprob=probdet1_all[min_exposure_idx]
        probdet2_maxprob=probdet2_all[min_exposure_idx]
        tt1_maxprob=telescope_time1_all[min_exposure_idx]
        tt2_maxprob=telescope_time1_all[min_exposure_idx]
        
        if m_exp:
            exp1_deep_maxprob=exp1_deep_all[min_exposure_idx]
            exp2_deep_maxprob=exp2_deep_all[min_exposure_idx]
            prob_area_deg_maxprob=prob_area_deg_all[min_exposure_idx]
            prob_area_deg_deep_maxprob=prob_area_deg_deep_all[min_exposure_idx]

        allstrategy_prob.append(prob_top)
        TT.append(tt_maxprob)
        TT1.append(tt1_maxprob)
        TT2.append(tt2_maxprob)
        obs1.append(obs1_maxprob)
        obs2.append(obs2_maxprob)
        exp1.append(exposure1_maxprob)
        exp2.append(exposure2_maxprob)
        bands.append(bands_maxprob)
        prob_area.append(area_maxprob)
        prob_area_deg.append(area_maxprob_deg)
        distance_all.append(distance)
        distance_all_err.append(distance_err)
        probdet1.append(probdet1_maxprob)
        probdet2.append(probdet2_maxprob)
        area_deg90_all.append(area_deg_info)

        if m_exp==True:
            if exp2_deep_maxprob < exposure2_maxprob:
                if exp2_deep_maxprob> 0.0 :
                    print("Found exp2 deep < exp 2")
                    sys.exit()
            exp1_deep.append(exp1_deep_maxprob)
            exp2_deep.append(exp2_deep_maxprob)
            prob_area_deep.append(prob_area_deg_maxprob)
            prob_area_deg_deep.append(prob_area_deg_deep_maxprob)
        strategy_type.append("Top")
        event_names.append(map_file_aux) 

        print("===== top probability")
        print (prob_top)
        print(tt_maxprob)

        det_prob = all_df["Detection Probability"].values
        prob_top10 = prob_top - (ltt_config[1] * prob_top)
        df_ltt = all_df[det_prob > prob_top10]
        
        df_ltt_tt01 = df_ltt["Telescope_time01"].values    
        df_ltt_tt02 = df_ltt["Telescope_time02"].values    
        total_telescope_time_ltt = df_ltt_tt01 + df_ltt_tt02

        print(df_ltt[:10])
        print()
        if area_deg_info > 199.0 and distance > 275.0 and plot_individual_strategies:
            print('here')
            total_telescope_timeltt_hour = np.divide(total_telescope_time_ltt,
                                                     60 * 60)
            
            df_ltt_= df_ltt.sort_values(by = "Total_TT")
            df_ltt_.to_csv(map_id+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+"_sortTT.csv")
            df_ltt_= df_ltt_.sort_values(by = "Detection Probability")
            df_ltt_.to_csv(map_id+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+"_sortProb.csv")
            
       
            df_ltt_ = df_ltt.sort_values(by = "Total_TT", ascending=False)
            probdet1_ltt_aux=df_ltt_["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt_aux=df_ltt_["Detprob2"].values
            ltt_plot=min(total_telescope_time_ltt)/(60*60)
            total_telescope_timeltt_=df_ltt_["Total_TT"].values
            pdet1_ltt=probdet1_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
            pdet2_ltt=probdet2_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
            highlightpoints_=[[probdet1_maxprob,probdet2_maxprob],[pdet1_ltt,pdet2_ltt]]
            highlightpoints_label_=['Top','Low Telescope Time']

        #     plt.clf()
        #     df_ltt_prob=all_df[all_df["Detection Probability"].values > 2.0]
        #     df_ltt_prob = df_ltt_prob.sort_values(by = "Deprob1", ascending=True)
        #     total_telescope_timeltt_prob=np.add(df_ltt_prob["Telescope_time01"].values,df_ltt_prob["Telescope_time02"].values) 
        #     total_telescope_timeltt_hour_prob=np.divide(total_telescope_timeltt_prob,60*60)
        #     df_ltt_prob["Total_TT"]=total_telescope_timeltt_hour_prob
        #     probdet1_ltt_aux=df_ltt_prob["Deprob1"].values#top["Detection Probability"].values
        #     probdet2_ltt_aux=df_ltt_prob["Detprob2"].values
        #     index_x=np.arange(0,len(probdet1_ltt_aux))
        #     ltt_plot=min(total_telescope_timeltt)/(60*60)
        #     total_telescope_timeltt_prob=df_ltt_prob["Total_TT"].values
        #     tt_1pass=df_ltt_prob["Telescope_time01"].values
        #     tt_1pass=np.divide(tt_1pass,60*60)
            
        #     fig = plt.figure()
        #     ax1 = fig.add_subplot(111)
        #     cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
        #     normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
        #     out=ax1.scatter(index_x,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=20)
        #     timecbar=plt.colorbar(out,ax=ax1)
        #     timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
        #     ax1.set_xlabel("Observational Set index",fontsize=16 )
        #     #plt.xlabel('ra')
        #     ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
        #     ax1.set_title(model_label, fontsize=16)
        #     plt.savefig(map_id+'probdet_index_'+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+'.png')

        #     fig = plt.figure()
        #     ax1 = fig.add_subplot(111)
        #     cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
        #     normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
        #     out=ax1.scatter(tt_1pass,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=10)
        #     timecbar=plt.colorbar(out,ax=ax1)
        #     timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
        #     ax1.set_xlabel('Telescope Time of first pass (Hours)',fontsize=16 )
        #     ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
        #     ax1.set_title(model_label, fontsize=16)
        #     plt.savefig(map_id+'probdet_tt_'+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+'.png')
        #     plt.clf()

        # ## Parei AQUI
        # if area_deg_info >150.0 and area_deg_info <200.0 and  distance >124.0 and  distance <176.0 and plot_individual_strategies==True:
        #     print("Making individual Plots")
        #     tt_edges=np.arange(0,16.5,0.5)
        #     shift_bin=0.0#0.1
        #     total_telescope_timeltt_hour=np.divide(total_telescope_timeltt,60*60)
        #     df_ltt["Total_TT"]=total_telescope_timeltt_hour
        #     df_ltt_ = df_ltt.sort_values(by = "Total_TT")
        #     df_ltt_.to_csv(map_id+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+"_sortTT.csv")
        #     df_ltt_ = df_ltt.sort_values(by = "Detection Probability")
        #     df_ltt_.to_csv(map_id+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+"_sortProb.csv")

        #     plt.clf()

        #     fig = plt.figure()
        #     ax1 = fig.add_subplot(111)
        #     tt_hist, bin_edges_tt = np.histogram(total_telescope_timeltt_hour,bins=tt_edges,density=False)
        #     tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        #     ax1.fill_between(tt_center+(shift_bin),tt_hist, step="mid", alpha=0.4, facecolor="indigo")
        #     ax1.step(tt_center+(shift_bin),tt_hist,alpha=0.8,where='mid', c="indigo")#label=types_of_strategy_name[i]#types_of_strategy_name[i]) #hatch='\    \'#width=tt_edges[1]-
        #     plt.xlabel(r'Telescope Time for top 10 per cent observational configurations (hours)',fontsize=16)
        #     plt.ylabel(r'Sets of observational configurations for the event',fontsize=16)
        #     plt.legend()
        #     plt.savefig(map_id+'TT_hist_'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png')
        #     df_ltt_ = df_ltt.sort_values(by = "Total_TT", ascending=False)
        #     probdet1_ltt_aux=df_ltt_["Deprob1"].values#top["Detection Probability"].values
        #     probdet2_ltt_aux=df_ltt_["Detprob2"].values
        #     ltt_plot=min(total_telescope_timeltt)/(60*60)
        #     total_telescope_timeltt_=df_ltt_["Total_TT"].values
        #     pdet1_ltt=probdet1_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
        #     pdet2_ltt=probdet2_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
        #     highlightpoints_=[[probdet1_maxprob,probdet2_maxprob],[pdet1_ltt,pdet2_ltt]]
        #     highlightpoints_label_=['Top','Low Telescope Time']
        #     create_color_dist_scatter(xdata=probdet1_ltt_aux,ydata=probdet2_ltt_aux,zdata=total_telescope_timeltt_,bin_edges_x=np.arange(0,16.1,0.5),bin_edges_y=np.arange(0,101,20),xlims=[60,92],ylims=[60,92],outname=map_id+'probs1_probs2'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png', zlevels=[90,80,70,60,50,40], colorzlevels=['Indigo','Purple','DarkViolet','MediumOrchid','Plum','Thistle','Lavender'],markers_sc=["o","v","s","P","*","X","D"], color_hist='Indigo', plot_weights=False, log_scale=False,highlightpoints=highlightpoints_,usecolormap='RdYlBu', highlightpoints_label=highlightpoints_label_,x_label='Discovery Probability $1^{st}$',y_label='Discovery Probability $2^{nd}$', x_ticks=None, top_hist_only=True, hist_labelx='Total Telescope Time (Hours)',hist_labely='Obs sets')
        #     plt.clf()

        #     df_ltt_prob=all_df[all_df["Detection Probability"].values > 2.0]
        #     df_ltt_prob = df_ltt_prob.sort_values(by = "Deprob1", ascending=True)
        #     total_telescope_timeltt_prob=np.add(df_ltt_prob["Telescope_time01"].values,df_ltt_prob["Telescope_time02"].values) 
        #     total_telescope_timeltt_hour_prob=np.divide(total_telescope_timeltt_prob,60*60)
        #     df_ltt_prob["Total_TT"]=total_telescope_timeltt_hour_prob
        #     probdet1_ltt_aux=df_ltt_prob["Deprob1"].values#top["Detection Probability"].values
        #     probdet2_ltt_aux=df_ltt_prob["Detprob2"].values
        #     index_x=np.arange(0,len(probdet1_ltt_aux))
        #     ltt_plot=min(total_telescope_timeltt)/(60*60)
        #     total_telescope_timeltt_prob=df_ltt_prob["Total_TT"].values
        #     tt_1pass=df_ltt_prob["Telescope_time01"].values
        #     tt_1pass=np.divide(tt_1pass,60*60)
        #     fig = plt.figure()
        #     ax1 = fig.add_subplot(111)
        #     cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
        #     normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
        #     out=ax1.scatter(index_x,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=20)
        #     timecbar=plt.colorbar(out,ax=ax1)
        #     timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
        #     ax1.set_xlabel("Observational Set index",fontsize=16 )
        #     ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
        #     ax1.set_title(model_label, fontsize=16)
        #     plt.savefig(map_id+'probdet_index_'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png')
        #     plt.clf()

        #     fig = plt.figure()
        #     ax1 = fig.add_subplot(111)
        #     cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
        #     normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
        #     out=ax1.scatter(tt_1pass,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=10)
        #     timecbar=plt.colorbar(out,ax=ax1)
        #     timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
        #     ax1.set_xlabel('Telescope Time of first pass (Hours)',fontsize=16 )
        #     ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
        #     ax1.set_title(model_label, fontsize=16)
        #     plt.savefig(map_id+'probdet_tt_'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png')
        #     plt.clf()
