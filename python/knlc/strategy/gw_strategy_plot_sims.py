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
"""
Pipeline to ....
"""

def plt_style():
    plt.rcParams.update({
                        'lines.linewidth':1.0,
                        'lines.linestyle':'-',
                        'lines.color':'black',
                        'font.family':'serif',
                        'font.weight':'normal',
                        'font.size':13.0,
                        'text.color':'black',
                        'text.usetex':True,
                        'axes.edgecolor':'black',
                        'axes.linewidth':1.0,
                        'axes.grid':False,
                        'axes.titlesize':'x-large',
                        'axes.labelsize':'x-large',
                        'axes.labelweight':'normal',
                        'axes.labelcolor':'black',
                        'axes.formatter.limits':[-4,4],
                        'xtick.major.size':7,
                        'xtick.minor.size':4,
                        'xtick.major.pad':8,
                        'xtick.minor.pad':8,
                        'xtick.labelsize':'medium',
                        'xtick.minor.width':1.0,
                        'xtick.major.width':1.0,
                        'ytick.major.size':7,
                        'ytick.minor.size':4,
                        'ytick.major.pad':8,
                        'ytick.minor.pad':8,
                        'ytick.labelsize':'medium',
                        'ytick.minor.width':1.0,
                        'ytick.major.width':1.0,
                        'legend.numpoints':1,
                        #'legend.fontsize':'x-large',
                        'legend.shadow':False,
                        'legend.frameon':False})

    return 
def get_map_distance(map_,savedir=''):

    pb,distmu,distsigma,distnorm = hp.read_map(map_,field=range(4), dtype=[np.float64,np.float64,np.float64,np.float64])#clecio dtype='numpy.float64'
    NSIDE=hp.npix2nside(pb.shape[0])
    pb_check=pb[np.logical_not(np.isinf(distmu))]
    distsigma_check=distsigma[np.logical_not(np.isinf(distmu))]
    distmu_check=distmu[np.logical_not(np.isinf(distmu))]
            
    pb_check=pb_check[np.logical_not(np.isinf(distsigma_check))]
    distmu_check=distmu_check[np.logical_not(np.isinf(distsigma_check))]
    distsigma_check=distsigma_check[np.logical_not(np.isinf(distsigma_check))]


    
    distmu_check_average= np.average(distmu_check,weights=pb_check)
    distsigma_check_average= np.average(distsigma_check,weights=pb_check)

    idx_sort = np.argsort(pb)
    idx_sort_up = list(reversed(idx_sort))
    resolution=hp.nside2pixarea(NSIDE,degrees=True)
    #sum_ = 0
    id_c = 0
    sum_full=0
    id_full=0
    area_max=sum(pb)
            
    while (sum_full<0.9) and (id_full <len(idx_sort_up)) :
        this_idx = idx_sort_up[id_full]
        sum_full = sum_full+pb[this_idx]
        id_full = id_full+1
        total_area=id_full*resolution

    if savedir!='':
        np.save(savedir,[distmu_check_average,distsigma_check_average])
    return distmu_check_average,distsigma_check_average


def create_heatmap(
        x,
        y,
        nCut,
        outname,
        feature_name, 
        x_ticks,
        y_ticks,
        plot_num=True,
        show_norm=False,
        method='sum_norm',
        z_=None,
        font_size=20,ticks_size=16,
        label_size=0,
        precision=".2f",
        silent=True,
        remove_zeros=False,
        add_error=False,
        cbarlabel='',
        tit=''
):

    data = np.transpose([x,y])#(size = (nSamples, 2))
    data = pd.DataFrame(data)
    if type(z_)!=type(None):
        data['n']=z_
    else:
        data['n'] = np.ones(len(x))
    cuts = pd.DataFrame({str(feature_name[feature]): pd.cut(data[feature], nCut[feature]) for feature in [0, 1]})
    if method=='sum_norm' or method=='sum':
        sums = data.join(cuts).groupby( list(cuts) ).sum()#mean()#
    if method=='mean':
        sums = data.join(cuts).groupby( list(cuts) ).median()
    if method=='median':
        sums = data.join(cuts).groupby( list(cuts) ).median()
    if add_error==True:
        sums_q15=data.join(cuts).groupby( list(cuts) ).quantile(q=0.15)
        sums_q87=data.join(cuts).groupby( list(cuts) ).quantile(q=0.87)
        sums_q15 = sums_q15.unstack(level = 0)
        sums_q15 = sums_q15.iloc[::-1]
        sums_q87 = sums_q87.unstack(level = 0)
        sums_q87 = sums_q87.iloc[::-1] 
        data_q=(sums_q87['n'].values - sums_q15['n'].values)/2.
    if silent==False:
       print('========== heatmap not silent')
       raw_data=data.join(cuts).groupby( list(cuts) ) 
       print(type(raw_data))
       print(type(raw_data['n']))
       print(raw_data['n'])
       print('=============')
       print('============= Quantile')
       print('=============')
       sums_q=data.join(cuts).groupby( list(cuts) ).quantile(q=0.5)#print(raw_data['n'].head((2,1)))
       print(sums_q['n'])
       print('============= vis median')
       print('=============')
       sums_median=data.join(cuts).groupby( list(cuts) ).count()#print(raw_data['n'].head((2,1)))
       print(sums_median['n'])
       sums=sums_median
       print('========== END of heatmap not silent')

    sums = sums.unstack(level = 0) # Use level 0 to put 0Bin as columns.

    # Reverse the order of the rows as the heatmap will print from top to bottom.
    sums = sums.iloc[::-1]
    y_ticks=np.flip(y_ticks)
    sums['n'] = sums['n'].replace(np.nan, 0)
    if method=='sum_norm':
        sums['n']=sums['n']/ (sum(sum(sums['n'].values)))
        if add_error==True:
            data_q=data_q/(sum(sum(sums['n'].values)))
    plt.clf()
    
    if remove_zeros==True:
        data_array=sums['n'].values
        annotation_=[]
        for l in range(0,len(data_array)):
            annotation_line=[]
            for m in range(0,len(data_array[0])):
                if data_array[l][m] <0.01:
                    ann_aux=' - '
                elif add_error:
                    if round(data_q[l][m],int(precision))==0.0:
                       data_q[l][m]=round(sqrt(data_array[l][m]),int(precision))
                    if data_array[l][m]+data_q[l][m] > 100:
                        data_q[l][m]=100-data_array[l][m]
                    if data_array[l][m]-data_q[l][m] < 0.0:
                       data_q[l][m]=data_array[l][m]
 
                    if int(precision)!=0:
                        
                        ann_aux=str(round(data_array[l][m],int(precision)))+'±'+str(round(data_q[l][m],int(precision)))
                    else:
                        ann_aux=str(int(round(data_array[l][m],int(precision))))+'±'+str(int(round(data_q[l][m],int(precision)))) 
                elif int(precision)==0:
                    ann_aux=str(int(round(data_array[l][m],int(precision))))  
                else: 
                    ann_aux=str(round(data_array[l][m],int(precision)))
                annotation_line.append(ann_aux)
            annotation_.append(annotation_line)
        precision=":^"
    else:
        annotation_=plot_num
    if add_error==True:    
        plt.figure(figsize=(16, 10))
    g=sns.heatmap(sums['n'], linewidths=.5, xticklabels=x_ticks, yticklabels=y_ticks, annot=annotation_, fmt=precision, annot_kws={"size": font_size}) 
    g.set_xticklabels(g.get_xmajorticklabels(), fontsize = ticks_size)
    g.set_yticklabels(g.get_ymajorticklabels(), fontsize = ticks_size)
    if label_size==0:
        label_size=ticks_size
    g.set_xlabel(g.get_xlabel(),fontsize=label_size)
    g.set_ylabel(g.get_ylabel(),fontsize=label_size)
    if cbarlabel!='':
        cbar = g.collections[0].colorbar
        cbar.ax.tick_params(labelsize=label_size)
        g.collections[0].colorbar.set_label(cbarlabel)
        g.figure.axes[-1].yaxis.label.set_size(20)
    if tit!='':
        g.set_title(tit, fontsize=20)

    plt.tight_layout()
    plt.savefig(outname)
    if silent==False:
        sys.exit()
    return 0


def create_color_dist_scatter(xdata,ydata,zdata=None,bin_edges_x=np.arange(1,250,20),bin_edges_y=np.arange(1,200,20),xlims=[0,250],ylims=[0,180],outname="area_dist_color.png", zlevels=[90,80,70,60,50,40], colorzlevels=['Indigo','Purple','DarkViolet','MediumOrchid','Plum','Thistle','Lavender'],markers_sc=["o","v","s","P","*","X","D"], color_hist='Indigo', plot_weights=False, log_scale=False,highlightpoints=None,usecolormap=None, highlightpoints_label='',x_label='Luminosity Distance (Mpc)',y_label='Area (sq-deg)', x_ticks=[0,50,100,150,200,250,300,350],top_hist_only=False, hist_labelx='',hist_labely=''):

    fig_ = plt.figure(figsize=(8, 6))      #(8,6)   

    if top_hist_only==False:
        grid = plt.GridSpec(2, 2, hspace=0.0, wspace=0.0, 
                        left=0.1, right=0.95, bottom=0.1, top=0.95,
                        height_ratios=[1, 4], width_ratios=[4,1])
    else:
        grid = plt.GridSpec(2, 1, hspace=0.23, wspace=0.0, 
                        left=0.1, right=0.95, bottom=0.15, top=0.95,height_ratios=[1, 4])#height_ratios=[4], width_ratios=[1,1]# hspace =0.18 right/top 0.95 left/bottom =0.08
    
    if top_hist_only==False:
        ax_main = fig_.add_subplot(grid[1,0])
        ax_hist_x = fig_.add_subplot(grid[0, 0], sharex=ax_main)
    
        ax_hist_y = fig_.add_subplot(grid[1, 1], sharey=ax_main)
    else:
        ax_main = fig_.add_subplot(grid[1,0])
        ax_hist_x = fig_.add_subplot(grid[0, 0])

    if (plot_weights==True) and (top_hist_only==False): #density=False, weights=None
        ax_hist_x.hist(xdata, bin_edges_x,
                   density=True, weights=zdata,
                   color=color_hist, alpha=0.8, 
                   histtype='stepfilled')

        ax_hist_y.hist(ydata, bin_edges_y, 
                   color=color_hist, alpha=0.8,
                   density=True, weights=zdata, 
                   histtype='stepfilled',
                   orientation='horizontal')        
    elif (top_hist_only==False):
        ax_hist_x.hist(xdata, bin_edges_x, 
                   color=color_hist, alpha=0.8, 
                   histtype='stepfilled')
    


        ax_hist_y.hist(ydata, bin_edges_y, 
                   color=color_hist, alpha=0.8,
                   histtype='stepfilled',
                   orientation='horizontal')
    else:
        ax_hist_x.hist(zdata, bin_edges_x, 
                   color=color_hist, alpha=0.8, 
                   histtype='stepfilled')
        ax_hist_x.set_xlabel(hist_labelx,fontsize=16, labelpad=-2)
        ax_hist_x.set_ylabel(hist_labely, fontsize=18, labelpad=-1.5)
        ax_hist_x.tick_params(axis='both', which='major', labelsize=14)
        ax_hist_x.tick_params(axis='both', which='minor', labelsize=14)
       

    if type(usecolormap)==type(None):
        if len(zlevels)>0:
            for i in range(0,len(zlevels)+1):
                if i==0:
                    plot_x=xdata[zdata>zlevels[i]]#np.logical_and
                    plot_y=ydata[zdata>zlevels[i]]
                    zlabel='P>'+str(zlevels[i])+'% '
                elif i==len(zlevels):
                    plot_x=xdata[zdata<zlevels[i-1]]#np.logical_and
                    plot_y=ydata[zdata<zlevels[i-1]]
                    zlabel='P<'+str(zlevels[i-1])+'% '
                else:
                    plot_x=xdata[np.logical_and(zdata>zlevels[i],zdata<zlevels[i-1])]#np.logical_and
                    plot_y=ydata[np.logical_and(zdata>zlevels[i],zdata<zlevels[i-1])]
                    zlabel=str(zlevels[i])+'%<P<'+str(zlevels[i-1])+'%'

  
                ax_main.scatter(plot_x,plot_y,c=colorzlevels[i],marker=markers_sc[i],label=zlabel, s=40)
        else:
            ax_main.scatter(xdata,ydata,c=colorzlevels[0],marker='+',s=15)
    else:
        cmscatter =plt.cm.get_cmap(usecolormap) #'RdYlBu' ##plt.cm.rainbow
        normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
        out=ax_main.scatter(xdata,ydata,c=zdata,cmap=cmscatter, norm=normscatter, marker='o',s=60)
        timecbar=plt.colorbar(out,ax=ax_main, pad=0.01)
        timecbar.set_label(hist_labelx,size=16 , labelpad=-2)

    if top_hist_only==False:
        ax_main.set_xlabel(x_label,fontsize=18)
        #plt.xlabel('ra')
        ax_main.set_ylabel(y_label, fontsize=18)
    else:
        ax_main.set_xlabel(x_label,fontsize=18, labelpad=-2)
        #plt.xlabel('ra')
        ax_main.set_ylabel(y_label, fontsize=18,  labelpad=-1.5)

    xmin,xmax=xlims
    ymin,ymax=ylims  
    ax_main.set_xlim((xmin, xmax))
    ax_main.set_ylim((ymin, ymax))

    if type(highlightpoints)!=type(None):
        listy=['--','-.',':','-']
        alp=[0.5,0.5,0.3,0.15]
        colors_highlight=['Indigo',"dimgray","purple"]
        xmin_,xmax_=xlims
        ymin_,ymax_=ylims  
        for i in range(0,len(highlightpoints)):
            
            ax_main.hlines(y=highlightpoints[i][1],xmin=xmin_,xmax=highlightpoints[i][0],ls=listy[i], alpha=alp[i],colors=colors_highlight[i], label=highlightpoints_label[i],lw=4)
            ax_main.vlines(x=highlightpoints[i][0],ymin=ymin_,ymax=highlightpoints[i][1], ls=listy[i], alpha=alp[i],colors=colors_highlight[i],lw=4)

    if top_hist_only==False:
        plt.setp(ax_hist_x.get_xticklabels(), visible=False)
        plt.setp(ax_hist_x.get_yticklabels(), visible=False)
        plt.setp(ax_hist_y.get_yticklabels(), visible=False)
        plt.setp(ax_hist_y.get_xticklabels(), visible=False)
    ax_hist_x.set_title(model_label, fontsize=18)
    if type(x_ticks)!=type(None):
        ax_main.set_xticks(x_ticks)
    ax_main.legend(frameon=False, fontsize=16)
    
    ax_main.tick_params(axis='both', which='major', labelsize=14)
    ax_main.tick_params(axis='both', which='minor', labelsize=14)
    if log_scale==True:
        ax_main.set_yscale('log')
    plt.savefig(outname, dpi=400)
    return 0

def correct_expt_deep(exptime):

    exptime=float(exptime)
    if exptime==60.0:
        return 90.0

    if exptime==90.0:
        return 120.0

    if exptime==120.0:
        return 200.0
    if exptime==200.0:
        return 300.0
    if exptime==300.0:
        return 600.0
    if exptime==600.0:
        return 1200.0
    if exptime==1200.0:
        return 2400.0
    if exptime==2400.0:
        return 3600.0
    if exptime==0.0:
        return 0.0
    else:
        print('I am exiting. This is the exptime ')
        print(exptime)
        sys.exit()


tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]


markers20=["*","D",".","v","^","<",">","1","2","3","4","8","s","p","P",",","o"]


for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


#================ START 
#===================================================================

if __name__ == '__main__':

    FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
    log.basicConfig(format=FORMAT, level=log.INFO)

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
    
    args = p.parse_args()

    map_path = args.map_path
    map_path_info = map_path
    strategy_dir = args.strategy_path
    constrain_time = args.constrain_time
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
    gw_o3_nsbh_dl=[240,290]#,370]
    nsbh_names=['GW190814 NSBH', 'GW200115 NSBH']#,'GW190426_152155 NSBH'] #'GW200115_042309 NSBH'
    nsns_markers=["P","v"]
    nsbh_markers=["s","o"]

    if m_exp==False:
        ares_depth_bins=[np.array([0.65,0.7,0.75,0.8,0.85,0.9])+0.025,[45.0,75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0]]

    if m_exp==True:
        dualexp_depth_bins=[[75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0,6825.0],[45.0,75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0]]
        ares_depth_bins=[np.array([0.65,0.75,0.85,0.95]),[45.0,75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0]]
        ares_depth_bins_deep=[np.array([0.25,0.35,0.55,0.75,0.85]),[75.0,105.0,135.0,265.0,335.0,865.0,1535.0,3225.0,3975.0,6825.0]] #0.45
        area_area_deep_bins=[np.array([0.65,0.75,0.85,0.95]),np.array([0.25,0.35,0.55,0.75,0.85])]
        x_area_deep=[0.3,0.5,0.7,0.8]#0.4
        y_area=[0.7,0.8,0.9]
        y_depth_deep=[1.5,2.0,3.3,5.0,10.0,20.0,40.0,60.0,90.0]#[90,120,200,300,600,1200,2400,3600,5400]
        depth_outer=[1.0,1.5,2.0,3.3,5.0,10.0,20.0,40.0,60.0]#[60,90,120,200,300,600,1200,2400,3600]

    if m_exp==False:
        xdvsw=[0.7,0.75,0.8,0.85,0.9]
        ydvsw=[60,90,120,200,300,600,1200,2400,3600]
    else:
        xdvsw=[0.7,0.8,0.9]
        ydvsw=[60,90,120,200,300,600,1200,2400,3600]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    constrain_time=True#True
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

    #FIXME
    if m_exp==True:
        exp1_deep=[]
        exp2_deep=[]
        prob_area_deep=[]
        prob_area_deg_deep=[]  
        
    event_names=[]
    probs_low_tt=[]
    probs_top_=[]
    probs_old=[]

    if evaluate_strategies==False:
        event_list_all=[]

    for i in range(0,len(event_list_all)):
    # for i in range(3):
        map_file_aux = os.path.basename(event_list_all[i])
        map_id = map_file_aux.split('_')[0]
        map_file = map_id+'.fits.gz'

        print("this is event "+str(i+1))
        try: 
            print(f"map info: {map_path_info+map_id}.npy")
            distance,distance_err=np.load(map_path_info+map_id+'.npy')
        except:
            print("Did not find npy with event info. "+\
                  "So, I am Calculating Map distance for "+map_file)
            distance,distance_err=get_map_distance(map_path+map_file,
                                                   savedir=map_path_info+map_id+'.npy')

        if _2d_map_dist==True:
            use_map_weights_info=map_path+"weights/"+map_id+"ac"+str(0.9)+"info.npy"
            use_map_info=map_path+"lowres/"+map_id+"_mapinfo.npy"

            try:
                area_deg_info,resolution=np.load(use_map_weights_info)
            except:
                area_deg_info=-1
                log.warning("Did not find npy with are info for "+\
                      use_map_weights_info)     
        else:
            area_deg_info=-1
        #top=pd.read_csv(event_list_top[i],comment='#')
        all_df=pd.read_csv(event_list_all[i],comment='#')

        if constrain_time:
            observation01 = all_df.Observation01.values
            observation02 = all_df.Observation02.values
            strategy_time_delays = np.add(-1 * observation01, observation02)
            all_df=all_df[np.logical_or(strategy_time_delays > 0.6,
                                        strategy_time_delays < 0.4) ]

        if (m_exp==True) and (use_ref_mult==True):
            file_old=map_file_aux.split('_')[0]+sufix_nomulti_strategy
            file_old=path_nomulti_Strategy+file_old
            try:
            #print("Opening file ",file_old)
                df_all_old=pd.read_csv(file_old,comment='#')
            except:
                print("Single exposure file not found: ",file_old)
                continue   
            if constrain_time:
                strategy_time_delays_nomulti=np.add(-1*np.array(df_all_old["Observation01"].values),df_all_old["Observation02"].values) 
                df_all_old=df_all_old[np.logical_or(strategy_time_delays_nomulti > 0.6,strategy_time_delays_nomulti < 0.4) ]
            #continue
            total_telescope_time_aux_old=np.add(df_all_old["Telescope_time01"].values,df_all_old["Telescope_time02"].values) #
            df_all_old=df_all_old[total_telescope_time_aux_old>0.0 ] # THIS SHOULDNT BE NECESSARY. CHECK IT.
        else:
            df_all_old=all_df

        if plot_bright_night_strategy==True:
            file_bright=map_id+sufix_bright_strategy
            file_bright=path_bright_Strategy+file_bright
            #print(file_bright)
            if plot_bright_ref_st==True:
                file_bright_ref=map_id+sufix_bright_ref_strategy
                file_bright_ref=path_bright_ref_Strategy+file_bright_ref
                try:
                #print("Opening file ",file_old)
                    df_all_bright_ref=pd.read_csv(file_bright_ref,comment='#')
                except:
                    print("Bright night strategy reference not file not found: ",file_bright_ref)
                    continue   

            try:
                df_all_bright=pd.read_csv(file_bright,comment='#')
            except:
                print("Bright night strategy not file not found: ",file_bright)
                continue

            if constrain_time==True:
                strategy_time_delays_bright=np.add(-1*np.array(df_all_bright["Observation01"].values),df_all_bright["Observation02"].values) 
                df_all_bright=df_all_bright[np.logical_or(strategy_time_delays_bright > 0.6,strategy_time_delays_bright < 0.4) ]
                if plot_bright_ref_st==True:
                    strategy_time_delays_bright_ref=np.add(-1*np.array(df_all_bright_ref["Observation01"].values),df_all_bright_ref["Observation02"].values) 
                    df_all_bright_ref=df_all_bright_ref[np.logical_or(strategy_time_delays_bright_ref > 0.6,strategy_time_delays_bright_ref < 0.4) ]
            if plot_bright_ref_st==True:
            # THIS SHOULDNT BE NECESSARY. CHECK IT. ## THIS SHOULDNT BE NECESSARY. CHECK IT.  ## THIS SHOULDNT BE NECESSARY. CHECK IT.
                total_telescope_time_aux_bright_ref=np.add(df_all_bright_ref["Telescope_time01"].values,df_all_bright_ref["Telescope_time02"].values) # Telescope_time01,Telescope_time02
                df_all_bright_ref=df_all_bright_ref[total_telescope_time_aux_bright_ref>0.0 ] # THIS SHOULDNT BE NECESSARY. CHECK IT.
            
            total_telescope_time_aux_bright=np.add(df_all_bright["Telescope_time01"].values,df_all_bright["Telescope_time02"].values) #
            df_all_bright=df_all_bright[total_telescope_time_aux_bright>0.0 ] # THIS SHOULDNT BE NECESSARY. CHECK IT.
            #continue

        total_telescope_time_aux_ini=np.add(all_df["Telescope_time01"].values,all_df["Telescope_time02"].values) # Telescope_time01,Telescope_time02
        all_df=all_df[total_telescope_time_aux_ini>0.0 ] # THIS SHOULDNT BE NECESSARY. CHECK IT.
        total_telescope_time=np.add(all_df["Telescope_time01"].values,all_df["Telescope_time02"].values)
        prob_all=all_df["Detection Probability"].values#top["Detection Probability"].values
        prob_top_test=max(prob_all)-0.01#-1.5
        all_df_top=all_df[all_df["Detection Probability"]>prob_top_test].copy().reset_index(drop=True)
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
        if m_exp==True:
            exp1_deep_all=all_df_top["Exposure01_deep"].values
            exp2_deep_all=all_df_top["Exposure02_deep"].values
            prob_area_deg_all=all_df_top["Region Coverage_deep"].values
            prob_area_deg_deep_all=all_df_top["Region_coverage_deg_deep"].values
        prob_top=max(all_df["Detection Probability"].values)
        if prob_top==0:
            print("undetected event for any strategy -- skipping")
            continue

        total_telescope_time_top=np.add(all_df_top["Telescope_time01"].values,all_df_top["Telescope_time02"].values)
        tt_maxprob=total_telescope_time_top[np.where(exposure1_all==min(exposure1_all))][0]
        obs1_maxprob=observation1_all[np.where(exposure1_all==min(exposure1_all))][0]
        obs2_maxprob=observation2_all[np.where(exposure1_all==min(exposure1_all))][0]
        exposure1_maxprob=exposure1_all[np.where(exposure1_all==min(exposure1_all))][0]
        exposure2_maxprob=exposure2_all[np.where(exposure1_all==min(exposure1_all))][0]
        area_maxprob=region_all[np.where(exposure1_all==min(exposure1_all))][0]
        area_maxprob_deg=region_all_deg[np.where(exposure1_all==min(exposure1_all))][0]
        bands_maxprob=bands_all[np.where(exposure1_all==min(exposure1_all))][0]
        probdet1_maxprob=probdet1_all[np.where(exposure1_all==min(exposure1_all))][0]
        probdet2_maxprob=probdet2_all[np.where(exposure1_all==min(exposure1_all))][0]
        tt1_maxprob=telescope_time1_all[np.where(exposure1_all==min(exposure1_all))][0]
        tt2_maxprob=telescope_time1_all[np.where(exposure1_all==min(exposure1_all))][0]
        if m_exp==True:
            exp1_deep_maxprob=exp1_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
            exp2_deep_maxprob=exp2_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
            prob_area_deg_maxprob=prob_area_deg_all[np.where(exposure1_all==min(exposure1_all))][0]
            prob_area_deg_deep_maxprob=prob_area_deg_deep_all[np.where(exposure1_all==min(exposure1_all))][0]

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
        
        df_ltt=all_df[all_df["Detection Probability"].values > (prob_top-(ltt_config[1]*prob_top))]
        total_telescope_timeltt=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
        if area_deg_info >199.0  and  distance > 275.0 and plot_individual_strategies==True: #and distance < :

            print("Making individual Plots")
            
            tt_edges=np.arange(0,16.5,0.5)
            shift_bin=0.0#0.1
            total_telescope_timeltt_hour=np.divide(total_telescope_timeltt,60*60)
            df_ltt["Total_TT"]=total_telescope_timeltt_hour

            df_ltt_= df_ltt.sort_values(by = "Total_TT")
            df_ltt_.to_csv(map_id+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+"_sortTT.csv")

            df_ltt_= df_ltt_.sort_values(by = "Detection Probability")
            df_ltt_.to_csv(map_id+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+"_sortProb.csv")

            plt.clf()

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            tt_hist, bin_edges_tt = np.histogram(total_telescope_timeltt_hour,bins=tt_edges,density=False)
            tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

            ax1.fill_between(tt_center+(shift_bin),tt_hist, step="mid", alpha=0.4, facecolor="indigo")
            ax1.step(tt_center+(shift_bin),tt_hist,alpha=0.8,where='mid', c="indigo")#label=types_of_strategy_name[i]#types_of_strategy_name[i]) #hatch='\    \'#width=tt_edges[1]-

            plt.xlabel(r'Telescope Time for top 10 per cent observational configurations (hours)',fontsize=16)
            plt.ylabel(r'Sets of observational configurations for the event',fontsize=16)
            plt.legend()
            plt.savefig(map_id+'TT_hist_'+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+'.png')

            df_ltt_ = df_ltt.sort_values(by = "Total_TT", ascending=False)
            probdet1_ltt_aux=df_ltt_["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt_aux=df_ltt_["Detprob2"].values
            ltt_plot=min(total_telescope_timeltt)/(60*60)
            total_telescope_timeltt_=df_ltt_["Total_TT"].values
            pdet1_ltt=probdet1_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
            pdet2_ltt=probdet2_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
            highlightpoints_=[[probdet1_maxprob,probdet2_maxprob],[pdet1_ltt,pdet2_ltt]]
            highlightpoints_label_=['Top','Low Telescope Time']
            create_color_dist_scatter(xdata=probdet1_ltt_aux,ydata=probdet2_ltt_aux,zdata=total_telescope_timeltt_,bin_edges_x=np.arange(0,16.1,0.5),bin_edges_y=np.arange(0,101,20),xlims=[60,92],ylims=[60,92],outname=map_id+'probs1_probs2'+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+'.png', zlevels=[90,80,70,60,50,40], colorzlevels=['Indigo','Purple','DarkViolet','MediumOrchid','Plum','Thistle','Lavender'],markers_sc=["o","v","s","P","*","X","D"], color_hist='Indigo', plot_weights=False, log_scale=False,highlightpoints=highlightpoints_,usecolormap='RdYlBu', highlightpoints_label=highlightpoints_label_,x_label='Discovery Probability $1^{st}$',y_label='Discovery Probability $2^{nd}$', x_ticks=None, top_hist_only=True, hist_labelx='Total Telescope Time (Hours)',hist_labely='Obs sets')

            plt.clf()
            df_ltt_prob=all_df[all_df["Detection Probability"].values > 2.0]
            df_ltt_prob = df_ltt_prob.sort_values(by = "Deprob1", ascending=True)
            total_telescope_timeltt_prob=np.add(df_ltt_prob["Telescope_time01"].values,df_ltt_prob["Telescope_time02"].values) 
            total_telescope_timeltt_hour_prob=np.divide(total_telescope_timeltt_prob,60*60)
            df_ltt_prob["Total_TT"]=total_telescope_timeltt_hour_prob
            probdet1_ltt_aux=df_ltt_prob["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt_aux=df_ltt_prob["Detprob2"].values
            index_x=np.arange(0,len(probdet1_ltt_aux))
            ltt_plot=min(total_telescope_timeltt)/(60*60)
            total_telescope_timeltt_prob=df_ltt_prob["Total_TT"].values
            tt_1pass=df_ltt_prob["Telescope_time01"].values
            tt_1pass=np.divide(tt_1pass,60*60)
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
            normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
            out=ax1.scatter(index_x,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=20)
            timecbar=plt.colorbar(out,ax=ax1)
            timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
            ax1.set_xlabel("Observational Set index",fontsize=16 )
            #plt.xlabel('ra')
            ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
            ax1.set_title(model_label, fontsize=16)
            plt.savefig(map_id+'probdet_index_'+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+'.png')

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
            normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
            out=ax1.scatter(tt_1pass,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=10)
            timecbar=plt.colorbar(out,ax=ax1)
            timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
            ax1.set_xlabel('Telescope Time of first pass (Hours)',fontsize=16 )
            ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
            ax1.set_title(model_label, fontsize=16)
            plt.savefig(map_id+'probdet_tt_'+"highdist"+str(round(distance,0))+"_higharea"+str(round(area_deg_info,0))+'.png')
            plt.clf()

        if area_deg_info >150.0 and area_deg_info <200.0 and  distance >124.0 and  distance <176.0 and plot_individual_strategies==True:
            print("Making individual Plots")
            tt_edges=np.arange(0,16.5,0.5)
            shift_bin=0.0#0.1
            total_telescope_timeltt_hour=np.divide(total_telescope_timeltt,60*60)
            df_ltt["Total_TT"]=total_telescope_timeltt_hour
            df_ltt_ = df_ltt.sort_values(by = "Total_TT")
            df_ltt_.to_csv(map_id+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+"_sortTT.csv")
            df_ltt_ = df_ltt.sort_values(by = "Detection Probability")
            df_ltt_.to_csv(map_id+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+"_sortProb.csv")

            plt.clf()

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            tt_hist, bin_edges_tt = np.histogram(total_telescope_timeltt_hour,bins=tt_edges,density=False)
            tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

            ax1.fill_between(tt_center+(shift_bin),tt_hist, step="mid", alpha=0.4, facecolor="indigo")
            ax1.step(tt_center+(shift_bin),tt_hist,alpha=0.8,where='mid', c="indigo")#label=types_of_strategy_name[i]#types_of_strategy_name[i]) #hatch='\    \'#width=tt_edges[1]-
            plt.xlabel(r'Telescope Time for top 10 per cent observational configurations (hours)',fontsize=16)
            plt.ylabel(r'Sets of observational configurations for the event',fontsize=16)
            plt.legend()
            plt.savefig(map_id+'TT_hist_'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png')
            df_ltt_ = df_ltt.sort_values(by = "Total_TT", ascending=False)
            probdet1_ltt_aux=df_ltt_["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt_aux=df_ltt_["Detprob2"].values
            ltt_plot=min(total_telescope_timeltt)/(60*60)
            total_telescope_timeltt_=df_ltt_["Total_TT"].values
            pdet1_ltt=probdet1_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
            pdet2_ltt=probdet2_ltt_aux[np.where(total_telescope_timeltt_==min(total_telescope_timeltt_))][0]
            highlightpoints_=[[probdet1_maxprob,probdet2_maxprob],[pdet1_ltt,pdet2_ltt]]
            highlightpoints_label_=['Top','Low Telescope Time']
            create_color_dist_scatter(xdata=probdet1_ltt_aux,ydata=probdet2_ltt_aux,zdata=total_telescope_timeltt_,bin_edges_x=np.arange(0,16.1,0.5),bin_edges_y=np.arange(0,101,20),xlims=[60,92],ylims=[60,92],outname=map_id+'probs1_probs2'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png', zlevels=[90,80,70,60,50,40], colorzlevels=['Indigo','Purple','DarkViolet','MediumOrchid','Plum','Thistle','Lavender'],markers_sc=["o","v","s","P","*","X","D"], color_hist='Indigo', plot_weights=False, log_scale=False,highlightpoints=highlightpoints_,usecolormap='RdYlBu', highlightpoints_label=highlightpoints_label_,x_label='Discovery Probability $1^{st}$',y_label='Discovery Probability $2^{nd}$', x_ticks=None, top_hist_only=True, hist_labelx='Total Telescope Time (Hours)',hist_labely='Obs sets')
            plt.clf()

            df_ltt_prob=all_df[all_df["Detection Probability"].values > 2.0]
            df_ltt_prob = df_ltt_prob.sort_values(by = "Deprob1", ascending=True)
            total_telescope_timeltt_prob=np.add(df_ltt_prob["Telescope_time01"].values,df_ltt_prob["Telescope_time02"].values) 
            total_telescope_timeltt_hour_prob=np.divide(total_telescope_timeltt_prob,60*60)
            df_ltt_prob["Total_TT"]=total_telescope_timeltt_hour_prob
            probdet1_ltt_aux=df_ltt_prob["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt_aux=df_ltt_prob["Detprob2"].values
            index_x=np.arange(0,len(probdet1_ltt_aux))
            ltt_plot=min(total_telescope_timeltt)/(60*60)
            total_telescope_timeltt_prob=df_ltt_prob["Total_TT"].values
            tt_1pass=df_ltt_prob["Telescope_time01"].values
            tt_1pass=np.divide(tt_1pass,60*60)
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
            normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
            out=ax1.scatter(index_x,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=20)
            timecbar=plt.colorbar(out,ax=ax1)
            timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
            ax1.set_xlabel("Observational Set index",fontsize=16 )
            ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
            ax1.set_title(model_label, fontsize=16)
            plt.savefig(map_id+'probdet_index_'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png')
            plt.clf()

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            cmscatter =plt.cm.get_cmap('RdYlBu') #'RdYlBu' ##plt.cm.rainbow
            normscatter = colors.BoundaryNorm(np.arange(1, 16, 1), cmscatter.N)
            out=ax1.scatter(tt_1pass,probdet1_ltt_aux,c=total_telescope_timeltt_prob,cmap=cmscatter, norm=normscatter, marker='o',s=10)
            timecbar=plt.colorbar(out,ax=ax1)
            timecbar.set_label('Total Telescope Time (Hours)',size=12 , labelpad=-3)
            ax1.set_xlabel('Telescope Time of first pass (Hours)',fontsize=16 )
            ax1.set_ylabel("Discovery Probability $1^{st}$", fontsize=16)
            ax1.set_title(model_label, fontsize=16)
            plt.savefig(map_id+'probdet_tt_'+"lowdist"+str(round(distance,0))+"_lowarea"+str(round(area_deg_info,0))+'.png')
            plt.clf()


        prob_ltt=df_ltt["Detection Probability"].values
        prob_ltt_sel=prob_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))]
        ltt=min(total_telescope_timeltt)

        observation1_ltt=df_ltt["Observation01"].values
        observation2_ltt=df_ltt["Observation02"].values
        exposure1_ltt=df_ltt["Exposure01"].values
        exposure2_ltt=df_ltt["Exposure02"].values
        region_ltt=df_ltt["Region Coverage"].values
        region_ltt_deg=df_ltt["Region_coverage_deg"].values
        filters_ltt=df_ltt["Filter_comb"].values
        probdet1_ltt=df_ltt["Deprob1"].values#top["Detection Probability"].values
        probdet2_ltt=df_ltt["Detprob2"].values
        telescope_time1_ltt=df_ltt["Telescope_time01"].values
        telescope_time2_ltt=df_ltt["Telescope_time02"].values

        if m_exp==True:
            exp1_deep_ltt=df_ltt["Exposure01_deep"].values
            exp2_deep_ltt=df_ltt["Exposure02_deep"].values
            prob_area_deg_ltt=df_ltt["Region Coverage_deep"].values
            prob_area_deg_deep_ltt=df_ltt["Region_coverage_deg_deep"].values

        obs1_ltt=observation1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        obs2_ltt=observation2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        exp1_ltt=exposure1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        if len(exposure1_ltt)!=len(exposure2_ltt):
            print(len(exposure1_ltt)-len(exposure2_ltt) )
            print("inconsistency in catalog in telecope time. Exiting...")
            sys.exit()
        exp2_ltt=exposure2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        bands_ltt=filters_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        area_ltt=region_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        area_ltt_deg=region_ltt_deg[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

        pdet1_ltt=probdet1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        pdet2_ltt=probdet2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        tt1_ltt=telescope_time1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        tt2_ltt=telescope_time2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

        if m_exp==True:
            exp1_deep_lttprob=exp1_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            exp2_deep_lttprob=exp2_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            prob_area_deg_lttprob=prob_area_deg_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            prob_area_deg_deep_lttprob=prob_area_deg_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

        prob_top_ltt=prob_ltt_sel[0]
        TT.append(ltt)
        TT1.append(tt1_ltt)
        TT2.append(tt2_ltt)
        obs1.append(obs1_ltt)
        obs2.append(obs2_ltt)
        exp1.append(exp1_ltt)
        exp2.append(exp2_ltt)
        bands.append(bands_ltt)
        prob_area.append(area_ltt)
        prob_area_deg.append(area_ltt_deg)
        distance_all.append(distance)
        distance_all_err.append(distance_err)
        strategy_type.append("Telescope Time 10%")
        TT_low.append(ltt)
        allstrategy_prob.append(prob_ltt_sel[0])
        event_names.append(map_file_aux)
        probdet1.append(pdet1_ltt)
        probdet2.append(pdet2_ltt)
        #if 2d_map_dist==True:
        area_deg90_all.append(area_deg_info)

        if m_exp==True:
            exp1_deep.append(exp1_deep_lttprob)
            exp2_deep.append(exp2_deep_lttprob)
            prob_area_deep.append(prob_area_deg_lttprob)
            prob_area_deg_deep.append(prob_area_deg_deep_lttprob)

        # ====================== half nights
        half_night_time=60*60*4
        df_ltt=all_df[all_df["Detection Probability"].values > (prob_top_ltt-(ltt_config[0]*prob_top_ltt))]
        total_telescope_timeltt_hn_aux=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)

        df_ltt=df_ltt[total_telescope_timeltt_hn_aux < (ltt+(ltt_config[2]*ltt))]
        total_telescope_timeltt_hn=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
        strategy_time_delays_hn=np.add(-1*np.array(df_ltt["Observation01"].values),df_ltt["Observation02"].values) 
        df_ltt=df_ltt[np.logical_not(np.logical_and(total_telescope_timeltt_hn > half_night_time,strategy_time_delays_hn < 0.4)) ]
        df_ltt=df_ltt[df_ltt["Telescope_time01"].values < half_night_time]
        df_ltt=df_ltt[df_ltt["Telescope_time02"].values < half_night_time]

        #stopped here
        prob_ltt_hn=df_ltt["Detection Probability"].values
        if len(prob_ltt_hn)==0:
            TT.append(0.0)
            TT1.append(0.0)
            TT2.append(0.0)
            obs1.append(0.0)
            obs2.append(0.0)
            exp1.append(0.0)
            exp2.append(0.0)
            bands.append("gg")
            prob_area.append(0.0)
            prob_area_deg.append(-1.0)
            distance_all.append(distance)
            distance_all_err.append(distance_err)
            strategy_type.append("Half Nights")
            allstrategy_prob.append(-1.0)
            event_names.append(map_file_aux)
            probdet1.append(-1.0)
            probdet2.append(-1.0)
            #if 2d_map_dist==True:
            area_deg90_all.append(area_deg_info)
            if m_exp==True:
                exp1_deep.append(0.0)
                exp2_deep.append(0.0)
                prob_area_deep.append(0.0)
                prob_area_deg_deep.append(-1.0)
        else:

            prob_ltt_sel_hn=max(prob_ltt_hn)#prob_ltt[np.where(total_time_discovery==min(total_time_discovery))]
            total_telescope_timeltt_half=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
            #total_time_discovery=df_ltt["Observation02"].values
            observation1_ltt=df_ltt["Observation01"].values
            observation2_ltt=df_ltt["Observation02"].values
            exposure1_ltt=df_ltt["Exposure01"].values
            exposure2_ltt=df_ltt["Exposure02"].values
            region_ltt=df_ltt["Region Coverage"].values
            region_ltt_deg=df_ltt["Region_coverage_deg"].values
            filters_ltt=df_ltt["Filter_comb"].values
            probdet1_ltt=df_ltt["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt=df_ltt["Detprob2"].values
            telescope_time1_ltt=df_ltt["Telescope_time01"].values
            telescope_time2_ltt=df_ltt["Telescope_time02"].values


            if m_exp==True:
                exp1_deep_ltt=df_ltt["Exposure01_deep"].values
                exp2_deep_ltt=df_ltt["Exposure02_deep"].values
                prob_area_deg_ltt=df_ltt["Region Coverage_deep"].values
                prob_area_deg_deep_ltt=df_ltt["Region_coverage_deg_deep"].values


            ltt_half=total_telescope_timeltt_half[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            obs1_ltt=observation1_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            obs2_ltt=observation2_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            exp1_ltt=exposure1_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            if len(exposure1_ltt)!=len(exposure2_ltt):
                print(len(exposure1_ltt)-len(exposure2_ltt) )
                print("inconsistency in catalog in half nights. Exiting...")
                sys.exit()
            exp2_ltt=exposure2_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            bands_ltt=filters_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            area_ltt=region_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            area_ltt_deg=region_ltt_deg[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            #ltt_half=total_telescope_timeltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]

            pdet1_ltt=probdet1_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            pdet2_ltt=probdet2_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            tt1_ltt=telescope_time1_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
            tt2_ltt=telescope_time2_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]


            if m_exp==True:
                exp1_deep_lttprob=exp1_deep_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
                exp2_deep_lttprob=exp2_deep_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
                prob_area_deg_lttprob=prob_area_deg_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]
                prob_area_deg_deep_lttprob=prob_area_deg_deep_ltt[np.where(prob_ltt_hn==prob_ltt_sel_hn)][0]

            TT.append(ltt_half)
            TT1.append(tt1_ltt)
            TT2.append(tt2_ltt)
            obs1.append(obs1_ltt)
            obs2.append(obs2_ltt)
            exp1.append(exp1_ltt)
            exp2.append(exp2_ltt)
            bands.append(bands_ltt)
            prob_area.append(area_ltt)
            prob_area_deg.append(area_ltt_deg)
            distance_all.append(distance)
            distance_all_err.append(distance_err)
            strategy_type.append("Half Nights")
            allstrategy_prob.append(prob_ltt_sel_hn)
            event_names.append(map_file_aux)
            probdet1.append(pdet1_ltt)
            probdet2.append(pdet2_ltt)
            #if 2d_map_dist==True:
            area_deg90_all.append(area_deg_info)
            if m_exp==True:
                exp1_deep.append(exp1_deep_lttprob)
                exp2_deep.append(exp2_deep_lttprob)
                prob_area_deep.append(prob_area_deg_lttprob)
                prob_area_deg_deep.append(prob_area_deg_deep_lttprob)

        #========= Early discovery    ===========================

        if low_budget==True:
            df_ltt=all_df[all_df["Detection Probability"].values > (prob_top_ltt-(ltt_config[0]*prob_top_ltt))]
        else:
            df_ltt=all_df[all_df["Detection Probability"].values > (prob_top-(ltt_config[1]*prob_top))]
        
        print("prob_top=",prob_top," prob_ltt=",prob_top_ltt," probltt-0.05=",prob_top_ltt-(ltt_config[0]*prob_top_ltt)," ltt=",ltt)
        total_telescope_timeltt_aux=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
        print(len(total_telescope_timeltt_aux))
        if low_budget==True:
            df_ltt=df_ltt[total_telescope_timeltt_aux < (ltt+(ltt_config[2]*ltt))]
        print(len(total_telescope_timeltt))
        
        total_time_discovery=df_ltt["Observation02"].values#np.add(df_ltt["Observation01"].values,df_ltt["Observation02"].values)
        total_telescope_timeltt=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
        prob_ltt=df_ltt["Detection Probability"].values
        prob_ltt_sel_d=prob_ltt[np.where(total_time_discovery==min(total_time_discovery))]

        observation1_ltt=df_ltt["Observation01"].values
        observation2_ltt=df_ltt["Observation02"].values
        exposure1_ltt=df_ltt["Exposure01"].values
        exposure2_ltt=df_ltt["Exposure02"].values
        region_ltt=df_ltt["Region Coverage"].values
        region_ltt_deg=df_ltt["Region_coverage_deg"].values
        filters_ltt=df_ltt["Filter_comb"].values
        probdet1_ltt=df_ltt["Deprob1"].values#top["Detection Probability"].values
        probdet2_ltt=df_ltt["Detprob2"].values
        telescope_time1_ltt=df_ltt["Telescope_time01"].values
        telescope_time2_ltt=df_ltt["Telescope_time02"].values

        if m_exp==True:
            exp1_deep_ltt=df_ltt["Exposure01_deep"].values
            exp2_deep_ltt=df_ltt["Exposure02_deep"].values
            prob_area_deg_ltt=df_ltt["Region Coverage_deep"].values
            prob_area_deg_deep_ltt=df_ltt["Region_coverage_deg_deep"].values

        obs1_ltt=observation1_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        obs2_ltt=observation2_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        exp1_ltt=exposure1_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        if len(exposure1_ltt)!=len(exposure2_ltt):
            print(len(exposure1_ltt)-len(exposure2_ltt) )
            print("inconsistency in catalog in early discovery. Exiting...")
            sys.exit()
        exp2_ltt=exposure2_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        bands_ltt=filters_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        area_ltt=region_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        area_ltt_deg=region_ltt_deg[np.where(total_time_discovery==min(total_time_discovery))][0]
        tt_ed=total_telescope_timeltt[np.where(total_time_discovery==min(total_time_discovery))][0]

        pdet1_ltt=probdet1_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        pdet2_ltt=probdet2_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        tt1_ltt=telescope_time1_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
        tt2_ltt=telescope_time2_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]

        if m_exp==True:
            exp1_deep_lttprob=exp1_deep_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
            exp2_deep_lttprob=exp2_deep_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
            prob_area_deg_lttprob=prob_area_deg_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]
            prob_area_deg_deep_lttprob=prob_area_deg_deep_ltt[np.where(total_time_discovery==min(total_time_discovery))][0]

        TT.append(tt_ed)
        TT1.append(tt1_ltt)
        TT2.append(tt2_ltt)
        obs1.append(obs1_ltt)
        obs2.append(obs2_ltt)
        exp1.append(exp1_ltt)
        exp2.append(exp2_ltt)
        bands.append(bands_ltt)
        prob_area.append(area_ltt)
        prob_area_deg.append(area_ltt_deg)
        distance_all.append(distance)
        distance_all_err.append(distance_err)
        strategy_type.append("Early discovery")
        allstrategy_prob.append(prob_ltt_sel_d[0])
        event_names.append(map_file_aux)
        probdet1.append(pdet1_ltt)
        probdet2.append(pdet2_ltt)
        #if 2d_map_dist==True:
        area_deg90_all.append(area_deg_info)
        if m_exp==True:
            exp1_deep.append(exp1_deep_lttprob)
            exp2_deep.append(exp2_deep_lttprob)
            prob_area_deep.append(prob_area_deg_lttprob)
            prob_area_deg_deep.append(prob_area_deg_deep_lttprob)
        #=================== late discovery

        if low_budget==True:
            df_ltt=all_df[all_df["Detection Probability"].values > (prob_top_ltt-(ltt_config[0]*prob_top_ltt))]
        else:
            df_ltt=all_df[all_df["Detection Probability"].values > (prob_top-(ltt_config[1]*prob_top))]
        
        total_telescope_timeltt_aux=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
        print(len(total_telescope_timeltt_aux))
        if low_budget==True:
            df_ltt=df_ltt[total_telescope_timeltt_aux < (ltt+(ltt_config[2]*ltt))]
        print(len(total_telescope_timeltt))
        late_start_discovery=df_ltt["Observation01"].values#np.add(df_ltt["Observation01"].values,df_ltt["Observation02"].values)
        df_ltt=df_ltt[late_start_discovery > 0.6]
        late_start_discovery_2=df_ltt["Observation02"]
        df_ltt=df_ltt[late_start_discovery_2 > 0.6]
        total_time_discovery=df_ltt["Observation02"].values
        total_telescope_timeltt=np.add(df_ltt["Telescope_time01"].values,df_ltt["Telescope_time02"].values)
        prob_ltt=df_ltt["Detection Probability"].values
        observation1_ltt=df_ltt["Observation01"].values
        observation2_ltt=df_ltt["Observation02"].values
        exposure1_ltt=df_ltt["Exposure01"].values
        exposure2_ltt=df_ltt["Exposure02"].values
        region_ltt=df_ltt["Region Coverage"].values
        region_ltt_deg=df_ltt["Region_coverage_deg"].values
        filters_ltt=df_ltt["Filter_comb"].values
        probdet1_ltt=df_ltt["Deprob1"].values#top["Detection Probability"].values
        probdet2_ltt=df_ltt["Detprob2"].values
        telescope_time1_ltt=df_ltt["Telescope_time01"].values
        telescope_time2_ltt=df_ltt["Telescope_time02"].values

        if m_exp==True:
            exp1_deep_ltt=df_ltt["Exposure01_deep"].values
            exp2_deep_ltt=df_ltt["Exposure02_deep"].values
            prob_area_deg_ltt=df_ltt["Region Coverage_deep"].values
            prob_area_deg_deep_ltt=df_ltt["Region_coverage_deg_deep"].values

        if len(prob_ltt)==0:
            TT.append(0.0)
            TT1.append(0.0)
            TT2.append(0.0)
            obs1.append(0.0)
            obs2.append(0.0)
            exp1.append(0.0)
            exp2.append(0.0)
            bands.append("gg")
            prob_area.append(0.0)
            prob_area_deg.append(-1.0)
            distance_all.append(distance)
            distance_all_err.append(distance_err)
            strategy_type.append("Late discovery")
            allstrategy_prob.append(-1.0)
            event_names.append(map_file_aux)
            probdet1.append(-1.0)
            probdet2.append(-1.0)
            #if 2d_map_dist==True:
            area_deg90_all.append(area_deg_info)
            if m_exp==True:
                exp1_deep.append(0.0)
                exp2_deep.append(0.0)
                prob_area_deep.append(0.0)
                prob_area_deg_deep.append(-1.0)
        else:

            prob_ltt_sel_d=prob_ltt[np.where(total_time_discovery==max(total_time_discovery))]
            obs1_ltt=observation1_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            obs2_ltt=observation2_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            exp1_ltt=exposure1_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            exp2_ltt=exposure2_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            bands_ltt=filters_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            area_ltt=region_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            area_ltt_deg=region_ltt_deg[np.where(total_time_discovery==max(total_time_discovery))][0]
            ltt=total_telescope_timeltt[np.where(total_time_discovery==max(total_time_discovery))][0]


            probdet1_ltt=probdet1_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            probdet2_ltt=probdet2_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            tt1_ltt=telescope_time1_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
            tt2_ltt=telescope_time2_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]


            if m_exp==True:
                exp1_deep_lttprob=exp1_deep_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
                exp2_deep_lttprob=exp2_deep_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
                prob_area_deg_lttprob=prob_area_deg_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]
                prob_area_deg_deep_lttprob=prob_area_deg_deep_ltt[np.where(total_time_discovery==max(total_time_discovery))][0]

            TT.append(ltt)
            TT1.append(tt1_ltt)
            TT2.append(tt2_ltt)
            obs1.append(obs1_ltt)
            obs2.append(obs2_ltt)
            exp1.append(exp1_ltt)
            exp2.append(exp2_ltt)
            bands.append(bands_ltt)
            prob_area.append(area_ltt)
            prob_area_deg.append(area_ltt_deg)
            distance_all.append(distance)
            distance_all_err.append(distance_err)
            strategy_type.append("Late discovery")
            allstrategy_prob.append(prob_ltt_sel_d[0])
            event_names.append(map_file_aux)
            probdet1.append(probdet1_ltt)
            probdet2.append(probdet2_ltt)

            area_deg90_all.append(area_deg_info)
            if m_exp==True:
                exp1_deep.append(exp1_deep_lttprob)
                exp2_deep.append(exp2_deep_lttprob)
                prob_area_deep.append(prob_area_deg_lttprob)
                prob_area_deg_deep.append(prob_area_deg_deep_lttprob)

        #============================ bright night strategy =================================:
    
        if plot_bright_night_strategy==True:
            
            prob_all_bright=df_all_bright["Detection Probability"].values#top["Detection Probability"].values
            prob_top_test_bright=max(prob_all_bright)-0.01
            df_all_bright_top=df_all_bright[df_all_bright["Detection Probability"]>prob_top_test_bright].copy().reset_index(drop=True)

            observation1_all=df_all_bright_top["Observation01"].values
            observation2_all=df_all_bright_top["Observation02"].values
            exposure1_all=df_all_bright_top["Exposure01"].values
            exposure2_all=df_all_bright_top["Exposure02"].values
            region_all=df_all_bright_top["Region Coverage"].values
            region_all_deg=df_all_bright_top["Region_coverage_deg"].values
            probdet1_all=df_all_bright_top["Deprob1"].values#top["Detection Probability"].values
            probdet2_all=df_all_bright_top["Detprob2"].values
            telescope_time1_all=df_all_bright_top["Telescope_time01"].values
            telescope_time2_all=df_all_bright_top["Telescope_time02"].values
            bands_all=df_all_bright_top["Filter_comb"].values

            if m_exp==True:
                exp1_deep_all=df_all_bright_top["Exposure01_deep"].values
                exp2_deep_all=df_all_bright_top["Exposure02_deep"].values
                prob_area_deg_all=df_all_bright_top["Region Coverage_deep"].values
                prob_area_deg_deep_all=df_all_bright_top["Region_coverage_deg_deep"].values

            prob_top_bright=max(df_all_bright["Detection Probability"].values)
            total_telescope_time_top=np.add(df_all_bright_top["Telescope_time01"].values,df_all_bright_top["Telescope_time02"].values)
            tt_maxprob=total_telescope_time_top[np.where(exposure1_all==min(exposure1_all))][0]
            obs1_maxprob=observation1_all[np.where(exposure1_all==min(exposure1_all))][0]
            obs2_maxprob=observation2_all[np.where(exposure1_all==min(exposure1_all))][0]
            exposure1_maxprob=exposure1_all[np.where(exposure1_all==min(exposure1_all))][0]
            exposure2_maxprob=exposure2_all[np.where(exposure1_all==min(exposure1_all))][0]
            area_maxprob=region_all[np.where(exposure1_all==min(exposure1_all))][0]
            area_maxprob_deg=region_all_deg[np.where(exposure1_all==min(exposure1_all))][0]
            bands_maxprob=bands_all[np.where(exposure1_all==min(exposure1_all))][0]
            probdet1_maxprob=probdet1_all[np.where(exposure1_all==min(exposure1_all))][0]
            probdet2_maxprob=probdet2_all[np.where(exposure1_all==min(exposure1_all))][0]
            tt1_maxprob=telescope_time1_all[np.where(exposure1_all==min(exposure1_all))][0]
            tt2_maxprob=telescope_time1_all[np.where(exposure1_all==min(exposure1_all))][0]
            if m_exp==True:
                exp1_deep_maxprob=exp1_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
                exp2_deep_maxprob=exp2_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
                prob_area_deg_maxprob=prob_area_deg_all[np.where(exposure1_all==min(exposure1_all))][0]
                prob_area_deg_deep_maxprob=prob_area_deg_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
        
            allstrategy_prob.append(prob_top_bright)
        
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
                exp1_deep.append(exp1_deep_maxprob)
                exp2_deep.append(exp2_deep_maxprob)
                prob_area_deep.append(prob_area_deg_maxprob)
                prob_area_deg_deep.append(prob_area_deg_deep_maxprob)
            strategy_type.append("Top_bright")
            event_names.append(map_file_aux)

        # ============================== LTT for Bright night =========================

            df_ltt_bright=df_all_bright[df_all_bright["Detection Probability"].values > (prob_top_bright-(ltt_config[1]*prob_top_bright))]
            total_telescope_timeltt=np.add(df_ltt_bright["Telescope_time01"].values,df_ltt_bright["Telescope_time02"].values)
            prob_ltt=df_ltt_bright["Detection Probability"].values
            prob_ltt_sel=prob_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))]
            ltt=min(total_telescope_timeltt)

            observation1_ltt=df_ltt_bright["Observation01"].values
            observation2_ltt=df_ltt_bright["Observation02"].values
            exposure1_ltt=df_ltt_bright["Exposure01"].values
            exposure2_ltt=df_ltt_bright["Exposure02"].values
            region_ltt=df_ltt_bright["Region Coverage"].values
            region_ltt_deg=df_ltt_bright["Region_coverage_deg"].values
            filters_ltt=df_ltt_bright["Filter_comb"].values
            probdet1_ltt=df_ltt_bright["Deprob1"].values#top["Detection Probability"].values
            probdet2_ltt=df_ltt_bright["Detprob2"].values
            telescope_time1_ltt=df_ltt_bright["Telescope_time01"].values
            telescope_time2_ltt=df_ltt_bright["Telescope_time02"].values

            if m_exp==True:
                exp1_deep_ltt=df_ltt_bright["Exposure01_deep"].values
                exp2_deep_ltt=df_ltt_bright["Exposure02_deep"].values
                prob_area_deg_ltt=df_ltt_bright["Region Coverage_deep"].values
                prob_area_deg_deep_ltt=df_ltt_bright["Region_coverage_deg_deep"].values

            obs1_ltt=observation1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            obs2_ltt=observation2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            exp1_ltt=exposure1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            if len(exposure1_ltt)!=len(exposure2_ltt):
                print(len(exposure1_ltt)-len(exposure2_ltt) )
                print("inconsistency in catalog in telecope time. Exiting...")
                sys.exit()
            exp2_ltt=exposure2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            bands_ltt=filters_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            area_ltt=region_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            area_ltt_deg=region_ltt_deg[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

            pdet1_ltt=probdet1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            pdet2_ltt=probdet2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            tt1_ltt=telescope_time1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            tt2_ltt=telescope_time2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]


            if m_exp==True:
                exp1_deep_lttprob=exp1_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
                exp2_deep_lttprob=exp2_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
                prob_area_deg_lttprob=prob_area_deg_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
                prob_area_deg_deep_lttprob=prob_area_deg_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

            prob_top_ltt=prob_ltt_sel[0]
            TT.append(ltt)
            TT1.append(tt1_ltt)
            TT2.append(tt2_ltt)
            obs1.append(obs1_ltt)
            obs2.append(obs2_ltt)
            exp1.append(exp1_ltt)
            exp2.append(exp2_ltt)
            bands.append(bands_ltt)
            prob_area.append(area_ltt)
            prob_area_deg.append(area_ltt_deg)
            distance_all.append(distance)
            distance_all_err.append(distance_err)
            strategy_type.append("Telescope Time 10%_bright")
            TT_low.append(ltt)
            allstrategy_prob.append(prob_ltt_sel[0])
            event_names.append(map_file_aux)
            probdet1.append(pdet1_ltt)
            probdet2.append(pdet2_ltt)
            #if 2d_map_dist==True:
            area_deg90_all.append(area_deg_info)

            if m_exp==True:
                exp1_deep.append(exp1_deep_lttprob)
                exp2_deep.append(exp2_deep_lttprob)
                prob_area_deep.append(prob_area_deg_lttprob)
                prob_area_deg_deep.append(prob_area_deg_deep_lttprob)

            ######## reference strategy bright ========================================================================
            if plot_bright_ref_st==True:
                prob_ostrategy_bright=-1.0
                for j in range(0,len(old_strategy_arr)):
                    old_strategy=old_strategy_arr[j]
                    df_ostrategy_bright=df_all_bright_ref[(df_all_bright_ref["Exposure01"].values==old_strategy[0]) & (df_all_bright_ref["Exposure02"].values ==old_strategy[1]) & (df_all_bright_ref["Filter_comb"].values==old_strategy[2]) & (df_all_bright_ref["Observation01"].values ==old_strategy[3]) & (df_all_bright_ref["Observation02"].values==old_strategy[4]) & (df_all_bright_ref["Region Coverage"].values==old_strategy[5])].copy().reset_index(drop=True)
                    if df_ostrategy_bright["Detection Probability"].values[0] > prob_ostrategy_bright:
                        prob_ostrategy_bright=df_ostrategy_bright["Detection Probability"].values
                        TT_old_sel=np.add(df_ostrategy_bright["Telescope_time01"].values,df_ostrategy_bright["Telescope_time02"].values)[0]
                        TT1_old_sel=df_ostrategy_bright["Telescope_time01"].values[0]
                        TT2_old_sel=df_ostrategy_bright["Telescope_time02"].values[0]
                        obs1_old_sel=old_strategy[3]
                        obs2_old_sel=old_strategy[4]
                        exp1_old_sel=old_strategy[0]
                        exp2_old_sel=old_strategy[1]
                        bands_old_sel=old_strategy[2]
                        prob_area_old_sel=old_strategy[5]
                        prob_area_deg_old_sel=df_ostrategy_bright["Region_coverage_deg"].values[0]
                        probdet1_old_sel=df_ostrategy_bright["Deprob1"].values[0]
                        probdet2_old_sel=df_ostrategy_bright["Detprob2"].values[0]
                TT.append(TT_old_sel)
                TT1.append(TT1_old_sel)
                TT2.append(TT2_old_sel)
                obs1.append(obs1_old_sel)
                obs2.append(obs2_old_sel)
                exp1.append(exp1_old_sel)
                exp2.append(exp2_old_sel)
                bands.append(bands_old_sel)
                prob_area.append(prob_area_old_sel)
                prob_area_deg.append(prob_area_deg_old_sel)
                allstrategy_prob.append(prob_ostrategy_bright[0])
                distance_all.append(distance) 
                distance_all_err.append(distance_err)
                strategy_type.append("Reference_bright")
                event_names.append(map_file_aux)

                probdet1.append(probdet1_old_sel)
                probdet2.append(probdet2_old_sel)

            #if 2d_map_dist==True:
                area_deg90_all.append(area_deg_info)

                if m_exp==True:
                    exp1_deep.append(0.0)
                    exp2_deep.append(0.0)
                    prob_area_deep.append(0.0)
                    prob_area_deg_deep.append(0.0)

        #=========================================== Uniform exposures
        if (m_exp==True) and (use_ref_mult==False):
            continue
        
        prob_all_old=df_all_old["Detection Probability"].values#top["Detection Probability"].values
        prob_top_test_old=max(prob_all_old)-0.01
        df_all_old_top=df_all_old[df_all_old["Detection Probability"]>prob_top_test_old].copy().reset_index(drop=True)
        observation1_all=df_all_old_top["Observation01"].values
        observation2_all=df_all_old_top["Observation02"].values
        exposure1_all=df_all_old_top["Exposure01"].values
        exposure2_all=df_all_old_top["Exposure02"].values
        region_all=df_all_old_top["Region Coverage"].values
        region_all_deg=df_all_old_top["Region_coverage_deg"].values
        probdet1_all=df_all_old_top["Deprob1"].values#top["Detection Probability"].values
        probdet2_all=df_all_old_top["Detprob2"].values
        telescope_time1_all=df_all_old_top["Telescope_time01"].values
        telescope_time2_all=df_all_old_top["Telescope_time02"].values
        bands_all=df_all_old_top["Filter_comb"].values
        prob_top_single=max(df_all_old["Detection Probability"].values)
        total_telescope_time_top=np.add(df_all_old_top["Telescope_time01"].values,df_all_old_top["Telescope_time02"].values)
        tt_maxprob=total_telescope_time_top[np.where(exposure1_all==min(exposure1_all))][0]
        obs1_maxprob=observation1_all[np.where(exposure1_all==min(exposure1_all))][0]
        obs2_maxprob=observation2_all[np.where(exposure1_all==min(exposure1_all))][0]
        exposure1_maxprob=exposure1_all[np.where(exposure1_all==min(exposure1_all))][0]
        exposure2_maxprob=exposure2_all[np.where(exposure1_all==min(exposure1_all))][0]
        area_maxprob=region_all[np.where(exposure1_all==min(exposure1_all))][0]
        area_maxprob_deg=region_all_deg[np.where(exposure1_all==min(exposure1_all))][0]
        bands_maxprob=bands_all[np.where(exposure1_all==min(exposure1_all))][0]
        probdet1_maxprob=probdet1_all[np.where(exposure1_all==min(exposure1_all))][0]
        probdet2_maxprob=probdet2_all[np.where(exposure1_all==min(exposure1_all))][0]
        tt1_maxprob=telescope_time1_all[np.where(exposure1_all==min(exposure1_all))][0]
        tt2_maxprob=telescope_time1_all[np.where(exposure1_all==min(exposure1_all))][0]
        if m_exp==True:
            exp1_deep_maxprob=0.0 #exp1_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
            exp2_deep_maxprob=0.0 #exp2_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
            prob_area_deg_maxprob=0.0 #prob_area_deg_all[np.where(exposure1_all==min(exposure1_all))][0]
            prob_area_deg_deep_maxprob=0.0 #prob_area_deg_deep_all[np.where(exposure1_all==min(exposure1_all))][0]
        
        allstrategy_prob.append(prob_top_single)
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
            exp1_deep.append(exp1_deep_maxprob)
            exp2_deep.append(exp2_deep_maxprob)
            prob_area_deep.append(prob_area_deg_maxprob)
            prob_area_deg_deep.append(prob_area_deg_deep_maxprob)
        strategy_type.append("Top_single")
        event_names.append(map_file_aux) 

    #================================== LTT for uniform exp time
        if prob_top_single>0.01:
            df_ltt_s=df_all_old[df_all_old["Detection Probability"].values > (prob_top_single-(ltt_config[1]*prob_top_single))]
        else:
            df_ltt_s=df_all_old
        
        total_telescope_timeltt=np.add(df_ltt_s["Telescope_time01"].values,df_ltt_s["Telescope_time02"].values)
        prob_ltt=df_ltt_s["Detection Probability"].values
        print(file_old)
        
        prob_ltt_sel=prob_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))]
        ltt=min(total_telescope_timeltt)

        observation1_ltt=df_ltt_s["Observation01"].values
        observation2_ltt=df_ltt_s["Observation02"].values
        exposure1_ltt=df_ltt_s["Exposure01"].values
        exposure2_ltt=df_ltt_s["Exposure02"].values
        region_ltt=df_ltt_s["Region Coverage"].values
        region_ltt_deg=df_ltt_s["Region_coverage_deg"].values
        filters_ltt=df_ltt_s["Filter_comb"].values
        probdet1_ltt=df_ltt_s["Deprob1"].values#top["Detection Probability"].values
        probdet2_ltt=df_ltt_s["Detprob2"].values
        telescope_time1_ltt=df_ltt_s["Telescope_time01"].values
        telescope_time2_ltt=df_ltt_s["Telescope_time02"].values
        obs1_ltt=observation1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        obs2_ltt=observation2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        exp1_ltt=exposure1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        if len(exposure1_ltt)!=len(exposure2_ltt):
            print(len(exposure1_ltt)-len(exposure2_ltt) )
            print("inconsistency in catalog in telecope time. Exiting...")
            sys.exit()
        exp2_ltt=exposure2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        bands_ltt=filters_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        area_ltt=region_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        area_ltt_deg=region_ltt_deg[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

        pdet1_ltt=probdet1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        pdet2_ltt=probdet2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        tt1_ltt=telescope_time1_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
        tt2_ltt=telescope_time2_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

        if m_exp==True:
            exp1_deep_lttprob=0.0#exp1_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            exp2_deep_lttprob=0.0#exp2_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            prob_area_deg_lttprob=0.0#prob_area_deg_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]
            prob_area_deg_deep_lttprob=0.0#prob_area_deg_deep_ltt[np.where(total_telescope_timeltt==min(total_telescope_timeltt))][0]

        prob_top_ltt=prob_ltt_sel[0]
        TT.append(ltt)
        TT1.append(tt1_ltt)
        TT2.append(tt2_ltt)
        obs1.append(obs1_ltt)
        obs2.append(obs2_ltt)
        exp1.append(exp1_ltt)
        exp2.append(exp2_ltt)
        bands.append(bands_ltt)
        prob_area.append(area_ltt)
        prob_area_deg.append(area_ltt_deg)
        distance_all.append(distance)
        distance_all_err.append(distance_err)
        strategy_type.append("Telescope Time 10% single")
        TT_low.append(ltt)
        allstrategy_prob.append(prob_ltt_sel[0])
        event_names.append(map_file_aux)
        probdet1.append(pdet1_ltt)
        probdet2.append(pdet2_ltt)
        area_deg90_all.append(area_deg_info)

        if m_exp==True:
            exp1_deep.append(exp1_deep_lttprob)
            exp2_deep.append(exp2_deep_lttprob)
            prob_area_deep.append(prob_area_deg_lttprob)
            prob_area_deg_deep.append(prob_area_deg_deep_lttprob)

        prob_ostrategy=-1.0
        for j in range(0,len(old_strategy_arr)):
            old_strategy=old_strategy_arr[j]
            df_ostrategy=df_all_old[(df_all_old["Exposure01"].values==old_strategy[0]) & (df_all_old["Exposure02"].values ==old_strategy[1]) & (df_all_old["Filter_comb"].values==old_strategy[2]) & (df_all_old["Observation01"].values ==old_strategy[3]) & (df_all_old["Observation02"].values==old_strategy[4]) & (df_all_old["Region Coverage"].values==old_strategy[5])].copy().reset_index(drop=True)
            
            # print('prob ostrategy:', prob_ostrategy)
            # print(df_ostrategy.columns)
            # print(df_ostrategy[df_ostrategy.columns[0]])#.values[0])

            if (df_ostrategy["Detection Probability"].values[0] > prob_ostrategy).all():
                prob_ostrategy=df_ostrategy["Detection Probability"].values
                TT_old_sel=np.add(df_ostrategy["Telescope_time01"].values,df_ostrategy["Telescope_time02"].values)[0]
                TT1_old_sel=df_ostrategy["Telescope_time01"].values[0]
                TT2_old_sel=df_ostrategy["Telescope_time02"].values[0]
                obs1_old_sel=old_strategy[3]
                obs2_old_sel=old_strategy[4]
                exp1_old_sel=old_strategy[0]
                exp2_old_sel=old_strategy[1]
                bands_old_sel=old_strategy[2]
                prob_area_old_sel=old_strategy[5]
                prob_area_deg_old_sel=df_ostrategy["Region_coverage_deg"].values[0]
                probdet1_old_sel=df_ostrategy["Deprob1"].values[0]
                probdet2_old_sel=df_ostrategy["Detprob2"].values[0]
        TT.append(TT_old_sel)
        TT1.append(TT1_old_sel)
        TT2.append(TT2_old_sel)
        obs1.append(obs1_old_sel)
        obs2.append(obs2_old_sel)
        exp1.append(exp1_old_sel)
        exp2.append(exp2_old_sel)
        bands.append(bands_old_sel)
        prob_area.append(prob_area_old_sel)
        prob_area_deg.append(prob_area_deg_old_sel)
        allstrategy_prob.append(prob_ostrategy[0])
        distance_all.append(distance) 
        distance_all_err.append(distance_err)
        strategy_type.append("Reference")
        event_names.append(map_file_aux)
        probdet1.append(probdet1_old_sel)
        probdet2.append(probdet2_old_sel)
        area_deg90_all.append(area_deg_info)

        if m_exp==True:
            exp1_deep.append(0.0)
            exp2_deep.append(0.0)
            prob_area_deep.append(0.0)
            prob_area_deg_deep.append(0.0)
        color_plot_opt='indianred' 
        color_plot_os=tableau20[0]

    if evaluate_strategies==True:   
        TT=np.array(TT)#/(60.0*60.0)
        exp1=np.array(exp1)
        exp2=np.array(exp2)
        obs1=np.array(obs1)
        obs2=np.array(obs2)
        strategy_type=np.array(strategy_type)
        allstrategy_prob=np.array(allstrategy_prob)

        if (m_exp==True):

            strategy_dict={'Telescope Time': TT,'Observation1': np.array(obs1), 'Observation2': np.array(obs2), 'Strategy': np.array(strategy_type), 'Detection Prob': allstrategy_prob, 'Integrated Prob Area': np.array(prob_area), 'Filters': np.array(bands),'Exposure1': np.array(exp1), 'Exposure2': np.array(exp2), 'Distance':   np.array(distance_all), 'Coverage_deg': prob_area_deg,'TT/area': np.divide(TT,prob_area_deg), 'Prob01': probdet1  , 'Prob02': probdet2, 'TTime1': TT1, 'TTime2': TT2, 'Exposure1_deep':  np.array(exp1_deep), 'Exposure2_deep': np.array(exp2_deep), 'Integrated Prob Area deep': np.array(prob_area_deep), 'Coverage_deg_deep': np.array(prob_area_deg_deep), 'Area90_deg': np.array(area_deg90_all), 'Event_ID':event_names}

        else:
            strategy_dict={'Telescope Time': TT,'Observation1': np.array(obs1), 'Observation2': np.array(obs2), 'Strategy': np.array(strategy_type), 'Detection Prob': allstrategy_prob, 'Integrated Prob Area': np.array(prob_area), 'Filters': np.array(bands),'Exposure1': np.array(exp1), 'Exposure2': np.array(exp2), 'Distance':   np.array(distance_all), 'Coverage_deg': prob_area_deg,'TT/area': np.divide(TT,prob_area_deg), 'Prob01': probdet1  , 'Prob02': probdet2, 'TTime1': TT1, 'TTime2': TT2, 'Area90_deg': np.array(area_deg90_all)}

        strategy_df=pd.DataFrame.from_dict(strategy_dict)
        strategy_df.to_csv(strategy_csv_cat, index=False)  
    else:
        print('Loading strategy catalog file')
        strategy_df=pd.read_csv(strategy_csv_cat)

    print('Strategy after pair')
    print('Strategy after second pair')
    plt.clf()

    frac_label="Fraction of Simulated Events"

    if m_exp==False or use_ref_mult==True:

        if  plot_bright_night_strategy==True: #

            types_of_strategy_name=["Ref Strategy","Early Discovery","Low Telescope Time","Top", "Late discovery","Half Nights","Low Telescope Time - Bright"] #"Telescope Time 5%"
            types_of_strategy=["Reference","Early discovery","Telescope Time 10%", "Top", "Late discovery", "Half Nights", "Telescope Time 10%_bright"] #"Te
        else:
            types_of_strategy_name=["Ref Strategy","Early Discovery","Low Telescope Time","Top", "Late discovery","Half Nights"] #"Telescope Time 5%"
            types_of_strategy=["Reference","Early discovery","Telescope Time 10%", "Top", "Late discovery", "Half Nights"] #"Te
        color_strategy=["darkred","olive","indigo","dimgray","peru",'palegreen','indigo']
        ls_strategy=['dashdot','dotted','dashed','solid',(0, (3, 1, 1, 1, 1, 1)),'solid','dotted'] #  (0, (1, 10)) lossely dot and lossely dashdot
        hist_color=["lightgray"]

    else:
        types_of_strategy_name=["ED","TT", "Top", "LD"] #"Telescope Time 5%"
        types_of_strategy=["Early discovery","Telescope Time 10%", "Top", "Late discovery"] #"Telescope Time 5%"

    if  ((plot_bright_night_strategy==True) and (plot_bright_ref_st==True)): #

        types_of_strategy_bright=["Telescope Time 10%","Telescope Time 10%_bright", "Reference", "Reference_bright" ]
        types_of_strategy_name_bright=["Low Telescope Time","Low Telescope Time - Bright", "Ref Strategy", "Ref Strategy - Bright"]#["O4 proposed","O3 achieved"]
        tableau20_bright=["indigo",'indigo','darkred', 'darkred'] #'indigo'
        ls_bright=['solid','dotted','dashdot', 'dashed']
    elif plot_bright_night_strategy==True:
        types_of_strategy_bright=["Telescope Time 10%","Telescope Time 10%_bright", "Reference" ]
        types_of_strategy_name_bright=["Low Telescope Time","Low Telescope Time - Bright", "Ref Strategy"]#["O4 proposed","O3 achieved"]
        tableau20_bright=["indigo",'indigo','darkred', 'darkred'] #'indigo'
        ls_bright=['solid','dotted','dashdot', 'dashed']
    else:
        types_of_strategy_bright=[]
    if m_exp==True:
        #types_of_strategy_multi=["Top","O3 REF",]
        types_of_strategy_multi=["Telescope Time 10%","Reference","Telescope Time 10% single"]
        types_of_strategy_name_multi=["Low Telescope Time","Ref Strategy","Low Telescope Time - uniform"]#["O4 proposed","O3 achieved"]
        tableau20_multi=['indigo','darkred','indigo'] #'indigo'
        ls_multi=['solid','dashdot','dotted']


        #types_of_strategy_multi=["Top","O3 REF",]
        types_of_strategy_multi_two=["Top","Reference","Late discovery"]
        types_of_strategy_name_multi_two=["Top","Ref Strategy","Late discovery"]#["O4 proposed","O3 achieved"]
        tableau20_multi_two=['indigo','darkred','peru']
        ls_multi_two=['solid','dashdot',(0, (3, 1, 1, 1, 1, 1))]

        types_of_strategy_multi_two_tt=["Top","Reference","Telescope Time 10%"]
        types_of_strategy_name_multi_two_tt=["Top","Ref Strategy","Low Telescope Time"]#["O4 proposed","O3 achieved"]
        tableau20_multi_two_tt=["dimgray",'darkred','indigo']
        ls_multi_two_tt=['solid','dashdot','dashed']

    else:
        types_of_strategy_multi=[]
        types_of_strategy_name_multi=[]#["O4 proposed","O3 achieved"]
        tableau20_multi=[]
        ls_multi=[]

    types_of_strategy_restrict=["Telescope Time 10%","Reference"]#["Top","Reference"]#,"Late discovery"]
    types_of_strategy_name_restrict=["Low Telescope Time","Ref Strategy"]#,"Late discovery"]#["O4 proposed","O3 achieved"]
    tableau20_restrict=['indigo','darkred',"peru"]
    ls_restrict=['solid','dashdot',(0, (3, 1, 1, 1, 1, 1))]


    fig = plt.figure()
    ax = fig.add_subplot(111)

    #print('Strategy 01')
    markers_=['*','D','.']
    bin_edges_ref=np.arange(0,340,8)
    for i in range(0,len(types_of_strategy)):
        print('Strategy '+str(types_of_strategy[i]))
        print(strategy_df[:10])
        print(strategy_df.Strategy.values[:10])
        print(types_of_strategy)
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        print('Strategy after selection')
        total_prob_before_cut=len(strategy_plt_dict['Detection Prob'].values)
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<300.0)].copy().reset_index(drop=True)

        if m_exp==True:
            if types_of_strategy[i]!="Reference":
                strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg_deep'].values>0.0)].copy().reset_index(drop=True)
        objects_after=len(strategy_plt_dict['Distance'].values)

        if ((types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery")):
                total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
                strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
                total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
                print('=======+++++ '+types_of_strategy[i]+' lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))

        total_prob=strategy_plt_dict['Detection Prob'].values
        total_prob_after_cut=len(total_prob)
        print('+++++ events after coverage cut '+str(float(total_prob_after_cut)/float(total_prob_before_cut)))
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)')
    plt.axvline(x=nsns_limit, ls='dashed', c='r')
    plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')

    plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="darkred", marker="P", label="NSNS")
    plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="indigo", marker="s", label="NSBH")
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 380)
    plt.legend(loc='upper right')
    plt.ylabel(r'Confirmation Probability (2 X)')
    #plt.savefig('strategy_distance_all_'+sufix+'.png',dpi=400)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(111)

    #print('Strategy 01')
    markers_=['*','D','.']
    bin_edges_ref=np.arange(0,350,10)
    area_cut=300.0
    for i in range(0,len(types_of_strategy_restrict)):
        #print('Strategy '+str(types_of_strategy[i]))
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy_restrict[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)


        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)

        if (types_of_strategy_restrict[i])=="Half Nights" or (types_of_strategy_restrict[i]=="Late discovery") or (types_of_strategy_restrict[i]=="Late Discovery"):
                total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
                strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
                total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
                print('=======+++++ '+types_of_strategy_restrict[i]+' lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        

        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        


        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_restrict[i], label=types_of_strategy_name_restrict[i], color=tableau20_restrict[i])
        plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    ls_sty=['dotted',(0, (3, 5, 1, 5, 1, 5))]
    for j in range(0,len(gw_o3_nsns_dl)):
        #plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])
        plt.axvline(x=gw_o3_nsns_dl[j], label=nsns_names[j], ls=ls_sty[j] ,c="mediumblue")

    for j in range(0,len(gw_o3_nsbh_dl)):
        #plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])
        plt.axvline(x=gw_o3_nsbh_dl[j], label=nsbh_names[j], ls=ls_sty[j] ,c="orange")

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.title(model_label,  fontsize=16)
    plt.legend(loc='lower left', frameon=False)#'upper right' #'center left'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.savefig('strategy_distance_allfirstprob_'+sufix+'.png',dpi=400)
    plt.clf()



    for i in range(0,len(types_of_strategy_bright)):
        #print('Strategy '+str(types_of_strategy[i]))
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])

        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy_bright[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)

        if types_of_strategy_bright[i]=="Half Nights":
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))



        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_bright[i], label=types_of_strategy_name_bright[i], color=tableau20_bright[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_bright[i], label=types_of_strategy_name_bright[i]+' 2$^{nd}$ pass', color=tableau20_bright[i], alpha=0.3)
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    #nsns_markers=["P","v"]
    #for j in range(0,len(gw_o3_nsns_dl)):
    #    plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])

    #for j in range(0,len(gw_o3_nsbh_dl)):
    #    plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.legend(loc='lower left', frameon=False)#'upper right'#'center left'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.title(model_label, fontsize=16)
    plt.savefig('strategy_distance_allfirst_bright_prob_'+sufix+'.png',dpi=400)
    plt.clf()






    for i in range(0,len(types_of_strategy_multi)):
        #print('Strategy '+str(types_of_strategy[i]))
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])

        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy_multi[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)


        if types_of_strategy_multi[i]=="Half Nights":
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))

        if types_of_strategy_multi[i]=="Telescope Time 10% single":
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Uniform exp lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))


        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_multi[i], label=types_of_strategy_name_multi[i], color=tableau20_multi[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    #nsns_markers=["P","v"]
    #for j in range(0,len(gw_o3_nsns_dl)):
    #    plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])

    #for j in range(0,len(gw_o3_nsbh_dl)):
    #    plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.legend(loc='center left', frameon=False)#'upper right'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.title(model_label, fontsize=16)
    plt.savefig('strategy_distance_allfirst_multi_prob_'+sufix+'.png',dpi=400)
    plt.clf()



    for i in range(0,len(types_of_strategy_multi_two)):
        #print('Strategy '+str(types_of_strategy[i]))
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])

        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy_multi_two[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)

        if types_of_strategy_multi_two[i]=="Half Nights" or types_of_strategy_multi_two[i]=="Late discovery"  or (types_of_strategy_multi_two[i]=="Late Discovery"): 
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))

    #    if  types_of_strategy_multi_two[i]=="Late Discovery":
    #        total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
    #        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
    #        total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
    #        print('=======+++++ Late Discovery lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))

        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_multi_two[i], label=types_of_strategy_name_multi_two[i]+' 1$^{st}$ pass', color=tableau20_multi_two[i])


        bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_multi_two[i], label=types_of_strategy_name_multi_two[i]+' 2$^{nd}$ pass', color=tableau20_multi_two[i], alpha=0.5)
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    #nsns_markers=["P","v"]
    for j in range(0,len(gw_o3_nsns_dl)):
        plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])

    for j in range(0,len(gw_o3_nsbh_dl)):
        plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.legend(loc='center left', frameon=False)#'upper right'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.savefig('strategy_distance_allfirst_multi_prob_twopass'+sufix+'.png',dpi=400)
    plt.clf()





    for i in range(0,len(types_of_strategy_multi_two_tt)):
        #print('Strategy '+str(types_of_strategy[i]))
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])
        #print('THIS IS THE STRAGEGYYYYYYYYYYYYYYY:  '+types_of_strategy_multi[i])

        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy_multi_two_tt[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)
        if types_of_strategy_multi_two_tt[i]=="Half Nights":
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))

        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_multi_two_tt[i], label=types_of_strategy_name_multi_two_tt[i]+' 1$^{st}$ pass', color=tableau20_multi_two_tt[i])


        bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_multi_two_tt[i], label=types_of_strategy_name_multi_two_tt[i]+' 2$^{nd}$ pass', color=tableau20_multi_two_tt[i], alpha=0.5)
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    #nsns_markers=["P","v"]
    for j in range(0,len(gw_o3_nsns_dl)):
        plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])

    for j in range(0,len(gw_o3_nsbh_dl)):
        plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.legend(loc='center left', frameon=False)#'upper right'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.savefig('strategy_distance_allfirst_multi_prob_twopass_tt_'+sufix+'.png',dpi=400)
    plt.clf()






    for i in range(0,len(types_of_strategy_restrict)):
        #print('Strategy '+str(types_of_strategy[i]))
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy_restrict[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)
        if types_of_strategy_restrict[i]=="Half Nights":
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_restrict[i], label=types_of_strategy_name_restrict[i]+' 1$^{st}$ pass', color=tableau20_restrict[i])
        bin_center_2,bin_mean_2,bin_std_2=fc.make_bin(dist,second_prob,bin_edges_ref)
        plt.plot(bin_center_2, bin_mean_2, lw=2, ls=ls_restrict[i], label=types_of_strategy_name_restrict[i]+' 2$^{nd}$ pass', color=tableau20_restrict[i], alpha=0.5)


        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4 ", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    #nsns_markers=["P","v"]
    #for j in range(0,len(gw_o3_nsns_dl)):
    #    plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])

    #for j in range(0,len(gw_o3_nsbh_dl)):
    #    plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.legend(loc='lower left', frameon=False)#'upper right' #'lower left'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    plt.title(model_label, fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.savefig('strategy_distance_allfirst_secondprob_'+sufix+'.png',dpi=400)
    plt.clf()






    for i in range(0,len(types_of_strategy)):
        #print('Strategy '+str(types_of_strategy[i]))
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)
        if ((types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery")):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values<200.0)].copy().reset_index(drop=True)
            #strategy_plt_dict['Distance'].values
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
        total_prob=strategy_plt_dict['Detection Prob'].values
        first_prob=strategy_plt_dict['Prob01'].values
        second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls=ls_strategy[i], label=types_of_strategy_name[i], color=color_strategy[i])
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls=ls_restrict[i], label=types_of_strategy_name_restrict[i]+' second pass', color=tableau20_restrict[i], alpha=0.5)


        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20_restrict[i], alpha=0.3)


        #bin_center,bin_mean,bin_std=fc.make_bin(dist,first_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='dashed', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,second_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='-.', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.1)


        #plt.scatter(dist,total_prob, c=tableau20[i],marker=markers_[i])
        #plt.scatter([distance],[prob_ostrategy], c=tableau20[1], marker=markers_[i])
        #plt.scatter([distance],[prob_ltt_sel[0]], c=tableau20[2], marker='.')
    #plt.legend()
    plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)', fontsize=16)
    plt.axvline(x=nsns_limit, ls='dashed', c='r') 
    plt.axhline(y=90, ls='-.', alpha=0.3)
    plt.axhline(y=68, ls='-.', alpha=0.3)
    plt.axhline(y=50, ls='-.', alpha=0.3)
    plt.axhline(y=20, ls='-.', alpha=0.3)
    #plt.text(nsns_limit+5, 10, "NSNS O4", fontsize=10, rotation='vertical', color='r', alpha=0.7)
    plt.text(nsns_limit+5, 10, "BNS O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    plt.text(nsbh_limit+5, 10, "NSBH O4", fontsize=14, rotation='vertical', color='r', alpha=0.7)
    #plt.axvline(x=240, ls='-.', c='gray')
    #plt.text(245, 20, "GW190814", fontsize=10, rotation='vertica, color='gray', alpha=0.7)
    #plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')


    #nsns_markers=["P","v"]
    #for j in range(0,len(gw_o3_nsns_dl)):
    #    plt.scatter(x=[gw_o3_nsns_dl[j]], y=[0.0], c="black", marker=nsns_markers[j], label=nsns_names[j])

    #for j in range(0,len(gw_o3_nsbh_dl)):
    #    plt.scatter(x=[gw_o3_nsbh_dl[j]], y=[0.0], c="black", marker=nsbh_markers[j], label=nsbh_names[j])

    #plt.scatter(x=gw_o3_nsns_dl, y=[0.0 for i in gw_o3_nsns_dl], c="black", marker="P", label="NSNS")
    #plt.scatter(x=gw_o3_nsbh_dl, y=[0.0 for i in gw_o3_nsbh_dl], c="black", marker="s", label="NSBH")
    #plt.axhline(y=90, ls='-.', alpha=0.3)
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.xlim(10, 350)
    plt.legend(loc='lower left', frameon=False)#'upper right'
    plt.ylabel(r'Discovery Probability (%)', fontsize=16)
    plt.title(model_label, fontsize=16)
    #plt.axes().set_aspect('equal')
    plt.savefig('strategy_distance_allfirst_strategiesprob_'+sufix+'.png',dpi=400)
    plt.clf()







    #========== old plots








    for i in range(0,len(types_of_strategy)):
        #print('Strategy '+str(types_of_strategy[i]))
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        #print('Strategy after selection')
        objects_before=len(strategy_plt_dict['Distance'].values)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))

        ratio_ab=float(objects_after)/float(objects_before)
        total_prob_area=strategy_plt_dict['Integrated Prob Area'].values
        depth2=strategy_plt_dict['Exposure2'].values
        depth1=strategy_plt_dict['Exposure1'].values
        area90=strategy_plt_dict['Area90_deg'].values
        distance=strategy_plt_dict['Distance'].values

        #second_prob=strategy_plt_dict['Prob02'].values
        #dist=strategy_plt_dict['Distance'].values
        
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)

        #x='Exposure1', y='Integrated Prob Area'
        create_heatmap(x=total_prob_area,y=depth2,nCut=ares_depth_bins,outname=types_of_strategy[i]+sufix+'deepvswide.png',feature_name=['Integrated Prob Area','Exposure time (second pass)'], x_ticks=xdvsw,y_ticks=ydvsw, cbarlabel=frac_label, tit=model_label, font_size=22, ticks_size=20)
        create_heatmap(x=total_prob_area,y=depth1,nCut=ares_depth_bins,outname=types_of_strategy[i]+sufix+'deepvswideexp01.png',feature_name=['Integrated Prob Area','Exposure time (min)'], x_ticks=xdvsw,y_ticks=depth_outer, cbarlabel=frac_label, tit=model_label, font_size=22,ticks_size=15 ,label_size=20)

        filters=strategy_plt_dict['Filters'].values
        filters=[list(filts) for filts in filters]
        filter_1,filter_2=np.transpose(filters)
        #strategy_plt_dict['Filter01']=filter_1
        #strategy_plt_dict['Filter02']=filter_2
        for j in range(0,len(filter_1)):
            if filter_1[j]=='g':
                filter_1[j]='1'
            if filter_1[j]=='r':
                filter_1[j]='2'
            if filter_1[j]=='i':
                filter_1[j]='3'
            if filter_1[j]=='z':
                filter_1[j]='4'

            if filter_2[j]=='g':
                filter_2[j]='1'
            
            if filter_2[j]=='r':
                filter_2[j]='2'
            if filter_2[j]=='i':
                filter_2[j]='3'
            if filter_2[j]=='z':
                filter_2[j]='4'

        filter_1plot=filter_1.astype('float')
        filter_2plot=filter_2.astype('float')
        filter_bins=[np.array([1,2,3,4,5])-0.5,np.array([1,2,3,4,5])-0.5]
        filter_ticks=['g','r','i','z']    
        #h=sns.histplot(data=strategy_plt_dict, x='Filter01', y='Filter02', cbar=True, cbar_kws=dict(shrink=.75))
        create_heatmap(x=filter_1plot,y=filter_2plot,nCut=filter_bins,outname=types_of_strategy[i]+sufix+'filters.png',feature_name=['Filter 1','Filter 2'], x_ticks=filter_ticks,y_ticks=filter_ticks, cbarlabel=frac_label, tit=model_label, font_size=22,ticks_size=20)

        
        prob1det=strategy_plt_dict['Prob01'].values
        prob2det=strategy_plt_dict['Prob02'].values
        bins_prob=[[5,15,25,35,45,55,65,75,85,95,105],[5,15,25,35,45,55,65,75,85,95,105]]
        prob_ticks=[10,20,30,40,50,60,70,80,90,100]
        create_heatmap(x=prob1det,y=prob2det,nCut=bins_prob,outname=types_of_strategy[i]+sufix+'prob1_prob2.png',feature_name=['Prob 1 detection','Prob 2 detection'],    x_ticks=prob_ticks,y_ticks=prob_ticks, cbarlabel=frac_label, tit=model_label)
        plt.clf()


        if m_exp==True:
            exp1_deep=strategy_plt_dict['Exposure1_deep'].values
            exp2_deep=strategy_plt_dict['Exposure2_deep'].values
            total_prob_area_deep=strategy_plt_dict['Integrated Prob Area deep'].values  
            create_heatmap(x=exp1_deep,y=depth1,nCut=dualexp_depth_bins,outname=types_of_strategy[i]+sufix+'deepinvdeepout_1.png',feature_name=['Exposure 1 Inner','Exposure1'], x_ticks=y_depth_deep,y_ticks=depth_outer, cbarlabel=frac_label, tit=model_label, font_size=20,ticks_size=18)
            if i==2:
                #print(types_of_strategy[i])
                #print('Comparing the exp2_deep with exp')
                #print (np.transpose([exp2_deep[:10],depth2[:10]]))
                #print (len(exp2_deep))
                #print (len(depth2))   
                for k in range(0,len(exp2_deep)):
                    if (exp2_deep[k]< depth2[k]):
                        if exp2_deep[k] >0.0:
                            print('Found inconsintent pair')
                            print ([exp2_deep[k],depth2[k]])
                            sys.exit()
        
            create_heatmap(x=exp2_deep,y=depth2,nCut=dualexp_depth_bins,outname=types_of_strategy[i]+sufix+'deepinvdeepout_2.png',feature_name=['Exposure 2 Inner (min)','Exposure2 (min)'], x_ticks=y_depth_deep,y_ticks=depth_outer, cbarlabel=frac_label, tit=model_label, font_size=22,ticks_size=15 ,label_size=20) # font_size=20,ticks_size=18
            #FIXME
            create_heatmap(x=total_prob_area,y=total_prob_area_deep,nCut=area_area_deep_bins,outname=types_of_strategy[i]+sufix+'areainvsareaout_2.png',feature_name=['Integrated Prob Area','Integrated Prob Area Inner'], x_ticks=y_area,y_ticks=x_area_deep, cbarlabel=frac_label, tit=model_label, font_size=22,ticks_size=20)

            create_heatmap(x=total_prob_area_deep,y=exp1_deep,nCut=ares_depth_bins_deep,outname=types_of_strategy[i]+sufix+'deepinvareadeepout_1.png',feature_name=['Integrated Prob Area Inner','Exposure Time Inner (min)'], x_ticks=x_area_deep,y_ticks=y_depth_deep, show_norm=True, cbarlabel=frac_label, tit=model_label, font_size=22,ticks_size=15 ,label_size=20)

            create_heatmap(x=total_prob_area_deep,y=exp2_deep,nCut=ares_depth_bins_deep,outname=types_of_strategy[i]+sufix+'deepinvareadeepout_2.png',feature_name=['Integrated Prob Area Inner','Exposure Time 2 Inner (min)'], x_ticks=x_area_deep,y_ticks=y_depth_deep, cbarlabel=frac_label, tit=model_label, font_size=22,ticks_size=15 ,label_size=20)

            area90=strategy_plt_dict['Area90_deg'].values
            distance=strategy_plt_dict['Distance'].values
            #total_TT=strategy_plt_dict['Telescope Time'].values
            total_TT=np.divide(np.array(strategy_plt_dict['Telescope Time'].values).astype('float'),60*60)
            _bins_area_dist=[np.array([0.0,62.5,87.5,125,175.0,225.0,275.0,325.0,340.0]),[0,50,100,150,300.0]] #0.45 #200,250, #,37.5,
            x_dist=[50,75,100,150,200,250,300,350] #25
            x_area=[25,75,125,225]#200,250,275 

            print("================== ======================================")
            print("Last Heatmap")
            create_heatmap(x=distance,y=area90,nCut=_bins_area_dist,outname=types_of_strategy[i]+sufix+'dist_area_telescope_time.png',feature_name=['Luminosity Distance (Mpc)','90% area (deg)'], x_ticks=x_dist,y_ticks=x_area,z_=total_TT, method='median',font_size=22,ticks_size=20,precision="1", silent=True, remove_zeros=True ,add_error=False , cbarlabel="Average Telescope Time", tit=model_label)
            create_heatmap(x=distance,y=area90,nCut=_bins_area_dist,outname=types_of_strategy[i]+sufix+'dist_area_prob.png',feature_name=['Luminosity Distance (Mpc)','90% area (deg)'],  x_ticks=x_dist,y_ticks=x_area, z_=prob1det, method='median',font_size=24, ticks_size=20,precision="0", silent=True, remove_zeros=True,add_error=False, cbarlabel="Average Discovery Probability", tit=model_label)

            #create_heatmap(x,y,nCut,outname,feature_name, x_ticks,y_ticks, plot_num=True, show_norm=False,method='sum_norm',z_=None) 
            print("=================================")
        if _2d_map_dist==True:
            strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
            #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<250.0)].copy().reset_index(drop=True)
            #print(' ================= objects lost due to area cut was ' ,ratio_ab )
            distance_all=strategy_plt_dict['Distance'].values
            areas_90=strategy_plt_dict['Area90_deg'].values
            prob1det=strategy_plt_dict['Prob01'].values 
            print('Maximum probability:', prob1det.max())
            print('Minimum probability:', prob1det.min())
            #data_hist=pd.DataFrame({"90% area (deg)": areas_90,"Distance (Mpc)": distance_all, "Detection Probability" :prob1det})

            #h=sns.jointplot(data=data_hist, x="Distance (Mpc)",y="90% area (deg)", xlim=[1,200], ylim=[0,200] , hue="Detection Probability")#cbar=True

            #plt.savefig('GW_2D_scatter_areadist'+types_of_strategy[i]+'.png',dpi=400)    

            plt.clf()

            create_color_dist_scatter(xdata=distance_all,ydata=areas_90,zdata=prob1det,bin_edges_x=np.arange(1,362,20),bin_edges_y=np.arange(1,250,20),xlims=[-20,360],ylims=[0,250],outname='GW_PROBSIMS_2D_scatter_areadist'+types_of_strategy[i]+sufix+'.png', zlevels=[90,80,70,60],plot_weights=True)
            plt.clf()
            #create_color_dist_scatter(xdata=distance_all,ydata=areas_90,bin_edges_x=np.arange(1,400,20),bin_edges_y=np.arange(1,1000,20),xlims=[-20,400],ylims=[5,1000],outname='GW_SIMS_2D_scatter_areadist'+types_of_strategy[i]+sufix+'.png', zlevels=[],plot_weights=False, log_scale=True) #xscale('log')
            plt.clf()


        #'Exposure1_deep':  exp1_deep, 'Exposure2_deep': exp2_deep, 'Integrated Prob Area deep': prob_area_deep, 'Coverage_deg_deep': prob_area_deg_deep


    #for i in range(0,len(types_of_strategy)):
    #    strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
    #    total_prob=strategy_plt_dict['TT/area'].values
        #first_prob=strategy_plt_dict['TT/area'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        #dist=strategy_plt_dict['Distance'].values

        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #ax1.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #ax1.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)
    plt.clf()
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<area_cut)].copy().reset_index(drop=True)

        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))


        total_prob=strategy_plt_dict['TT/area'].values
        #first_prob=strategy_plt_dict['Prob01'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)


    plt.legend()
    #plt.yticks(np.arange(0, 100+1, 10))
    plt.xlabel(r'Luminosity Distance (MPC)')
    plt.axvline(x=nsns_limit, ls='dashed', c='r')
    plt.axvline(x=nsbh_limit, ls='-.', c='r')
    #plt.axvline(x=250, ls='-.', c='r')
    #plt.axhline(y=68, ls='-.', alpha=0.3)
    #plt.axhline(y=50, ls='-.', alpha=0.3)
    #plt.axhline(y=20, ls='-.', alpha=0.3)
    plt.ylabel(r'Telescope Time/Area')
    #plt.savefig('strategy_distance_tt_'+sufix+'.png',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,10,0.5)
    shift_bin=30
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values<120.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        total_prob=strategy_plt_dict['TT/area'].values
        total_prob=np.divide(1,total_prob)
        #first_prob=strategy_plt_dict['Prob01'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)

        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0]  fill=False

    plt.xlabel(r'Area/Telescope Time ($deg^2$/s)')
    plt.legend()
    #plt.savefig('tt_hist_lowD'+sufix+'.png')
    #plt.close('all')

    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,10,0.5)
    shift_bin=30
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values>120.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        total_prob=strategy_plt_dict['TT/area'].values
        total_prob=np.divide(1,total_prob)
        #first_prob=strategy_plt_dict['Prob01'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)

        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 

    plt.xlabel(r'Area/Telescope Time ($deg^2$/s)')
    plt.legend()
    #plt.savefig('tt_hist_highD'+sufix+'.png')




    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(1,55,3)
    shift_bin=0.5#30
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        #strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values>120.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        total_prob=np.divide(strategy_plt_dict['TT/area'].values,60*60)
        total_prob=np.divide(np.ones(len(total_prob)).astype('float'),total_prob.astype('float'))
        total_prob_nan=np.isnan(total_prob)
        print('greater than 0 in  tt/area is '+str(len(total_prob > 0.0)))
        print('NAns found in  area/tt is '+str(sum(total_prob_nan))+' of '+ str(len(total_prob)) )
        total_prob_plot=total_prob[np.logical_not(total_prob_nan)]
        #first_prob=strategy_plt_dict['Prob01'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)

        tt_hist, bin_edges_tt = np.histogram(total_prob_plot,bins=tt_edges,density=False)#bins=15
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 

    plt.xlabel(r'Area/Telescope Time ($deg^2$/hour)')
    plt.legend()
    #plt.savefig('tt_hist_allD'+sufix+'.png')





    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,16.5,0.5)
    shift_bin=0.1
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<area_cut)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values>190.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        total_prob=strategy_plt_dict['Telescope Time'].values
        total_prob=np.divide(strategy_plt_dict['Telescope Time'].values,60*60)
        #first_prob=strategy_plt_dict['Prob01'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)

        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 
    plt.xlabel(r'Telescope Time (hours)')
    plt.legend()
    #plt.savefig('tt_full_hist_highD'+sufix+'.png')



    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,16.5,0.5)
    shift_bin=0.1
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<area_cut)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values<190.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        total_prob=strategy_plt_dict['Telescope Time'].values
        total_prob=np.divide(strategy_plt_dict['Telescope Time'].values,60*60)
        #first_prob=strategy_plt_dict['Prob01'].values
        #second_prob=strategy_plt_dict['Prob02'].values
        dist=strategy_plt_dict['Distance'].values
        
        #bin_center,bin_mean,bin_std=fc.make_bin(dist,total_prob,bin_edges_ref)
        #plt.plot(bin_center, bin_mean, lw=2, ls='solid', label=types_of_strategy_name[i], color=tableau20[i])
        #plt.fill_between(bin_center, bin_mean+bin_std, bin_mean-bin_std, facecolor=tableau20[i], alpha=0.3)

        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 
    plt.xlabel(r'Telescope Time (hours)')
    plt.legend()
    #plt.savefig('tt_full_hist_lowD'+sufix+'.png')


    #===================================== novel histograms =====================================
    plt.clf()
    if m_exp==True:
        types_of_strategy=list(types_of_strategy)
        types_of_strategy.append('Top_single')
        types_of_strategy.append("Telescope Time 10% single")
        types_of_strategy_name=list(types_of_strategy_name)
        color_strategy=list(color_strategy)
        types_of_strategy_name.append("Top uniform exposure time")
        types_of_strategy_name.append("Low Telescope Time uniform")
        color_strategy.append("dimgray")
        color_strategy.append("indigo")


    tb_50tt=[]
    tb_90tt=[]
    tb_100tt=[]
    tb_n=[]
    n_events_50=[]
    n_events_90=[]

    tb_50ttarea=[]
    tb_90ttarea=[]
    tb_100ttarea=[]


    for i in range(0,len(types_of_strategy)):
        plt.clf()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        strategy_plt_top=strategy_df[(strategy_df['Strategy'].values=='Top')].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)
        strategy_plt_top=strategy_plt_top[(strategy_plt_top['Area90_deg'].values<area_cut)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late discovery") or (types_of_strategy[i]=="Late Discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        strategy_plt_dict_area=strategy_plt_dict.sort_values(by = 'Area90_deg', ascending=True) #        df_ltt_= df_ltt_.sort_values(by = "Detection Probability")
        total_tt_area=np.divide(strategy_plt_dict_area['Telescope Time'].values,60*60)
        n50_objects=int(len(total_tt_area)*0.5)
        n90_objects=int(len(total_tt_area)*0.9)
        n100_objects=len(total_tt_area)
        tb_50ttarea.append(sum(total_tt_area[0:n50_objects])/n50_objects)
        tb_90ttarea.append(sum(total_tt_area[0:n90_objects])/n90_objects)
        tb_100ttarea.append(sum(total_tt_area)/n100_objects) 
        total_prob=np.divide(strategy_plt_dict['Telescope Time'].values,60*60)
        total_prob_top=np.divide(strategy_plt_top['Telescope Time'].values,60*60)
        strategy_plt_dict_area=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
        dist=strategy_plt_dict['Distance'].values
        median_tt=np.percentile(total_prob,50)
        _90_tt=np.percentile(total_prob,90)
        _100_tt=np.percentile(total_prob,100)
        n_events_tt=len(total_prob)
        n_events_50.append(len(total_prob[total_prob<median_tt]))
        n_events_90.append(len(total_prob[total_prob<_90_tt]))
        print(types_of_strategy[i])
        if len(total_prob[total_prob<median_tt])==0:
            print(types_of_strategy[i])
            tb_50tt.append(0.0)
        else:
            tb_50tt.append(sum(total_prob[total_prob<median_tt])/len(total_prob[total_prob<median_tt]))
        tb_90tt.append(sum(total_prob[total_prob<_90_tt])/len(total_prob[total_prob<_90_tt]))
        tb_100tt.append(sum(total_prob[total_prob<_100_tt])/len(total_prob[total_prob<_100_tt]))
        tb_n.append(n_events_tt)
        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        tt_hist_top, bin_edges_tt_top = np.histogram(total_prob_top,bins=tt_edges,density=False)
        tt_center_top=bin_edges_tt_top[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)
        if types_of_strategy[i] != 'Top':
            ax1.fill_between(tt_center_top+(shift_bin),tt_hist_top, step="mid", alpha=0.4, facecolor=hist_color[0])
            ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.5,where='mid', c=hist_color[0], label='Top')#types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 

        ax1.fill_between(tt_center+(i*shift_bin),tt_hist, step="mid", alpha=0.4, facecolor=color_strategy[i])
        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=color_strategy[i], label=types_of_strategy_name[i])#types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-

        xmin_,xmax_=ax1.get_xlim()
        xmin_=0.0

        plt.axvline(x=median_tt, ls='dashed', c='r')
        if types_of_strategy[i] !="Early discovery":
            plt.text(median_tt+0.2, 25, "50%", fontsize=14, rotation='vertical', color='r', alpha=0.7)
        else:
            plt.text(median_tt+0.07, 25, "50%", fontsize=14, rotation='vertical', color='r', alpha=0.7) 
        plt.axvline(x=_90_tt, ls='dashed', c='r')
        plt.text(_90_tt+0.2, 25, "90%", fontsize=14, rotation='vertical', color='r', alpha=0.7)
        plt.xlabel(r'Telescope Time per event (hours)',fontsize=16)
        plt.ylabel(r'Number of simulated events',fontsize=16)
        plt.title(model_label, fontsize=16)
        plt.legend()
        ax1.set_xlim((0.1, xmax_))
        plt.savefig('tt_full_hist_lowD'+sufix+'_'+types_of_strategy[i]+'.png')






    tt_edges=np.arange(0,5.5,0.49)
    shift_bin=0.0#0.05
    for i in range(0,len(types_of_strategy)):
        plt.clf()

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        strategy_plt_top=strategy_df[(strategy_df['Strategy'].values=='Top')].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<area_cut)].copy().reset_index(drop=True)
        strategy_plt_top=strategy_plt_top[(strategy_plt_top['Coverage_deg'].values<area_cut)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
            print('=======+++++ Half Night lost '+str(float(total_hn_after_cut)/float(total_hn_before_cut))+' ',str(total_hn_after_cut),' of ', str(total_hn_before_cut))
        if (types_of_strategy[i]=="Late Discovery"):
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Late discovery"):
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
        total_prob=strategy_plt_dict['Observation2'].values
        total_prob_top=strategy_plt_top['Observation2'].values
        dist=strategy_plt_dict['Distance'].values
        median_tt=np.percentile(total_prob,50)
        _90_tt=np.percentile(total_prob,90)
        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)
        tt_hist_top, bin_edges_tt_top = np.histogram(total_prob_top,bins=tt_edges,density=False)
        tt_center_top=bin_edges_tt_top[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        print(types_of_strategy[i]+' bins and hist')
        print(tt_hist)
        print(bin_edges_tt)  
        if types_of_strategy[i] != 'Top':
            ax1.fill_between(tt_center_top,tt_hist_top, step="mid", alpha=0.5, facecolor=hist_color[0])
            ax1.step(tt_center+(shift_bin),tt_hist,alpha=0.3,where='mid', c=hist_color[0], label='Top')#types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 

        ax1.fill_between(tt_center+(shift_bin),tt_hist, step="mid", alpha=0.4, facecolor=color_strategy[i])
        ax1.step(tt_center_top+(shift_bin),tt_hist,alpha=0.8,where='mid', c=color_strategy[i], label=types_of_strategy_name[i])#types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-

        xmin_,xmax_=ax1.get_xlim()
        xmin_=0.0
        plt.xlabel(r'Observing day of Confirmation (second pass)',fontsize=16)
        plt.ylabel(r'Number of simulated events',fontsize=16)
        plt.legend()
        plt.title(model_label, fontsize=16)
        ax1.set_xlim((0.1, xmax_))
        plt.savefig('obs2_hist_'+sufix+'_'+types_of_strategy[i]+'.png')

    #=================================================================================
    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,16.5,0.5)
    shift_bin=0.1
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<area_cut)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late Discovery") or (types_of_strategy[i]=="Late discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
        total_prob=strategy_plt_dict['Telescope Time'].values
        total_prob=np.divide(strategy_plt_dict['Telescope Time'].values,60*60)
        dist=strategy_plt_dict['Distance'].values
        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 
    plt.xlabel(r'Telescope Time (hours)')
    plt.legend()

    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,5.5,0.5)
    shift_bin=0.1
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values<120.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late Discovery") or (types_of_strategy[i]=="Late discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
        total_prob=strategy_plt_dict['Observation2'].values
        dist=strategy_plt_dict['Distance'].values
        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 
    plt.xlabel(r'Telescope Time (hours)')
    plt.legend()

    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,4.5,0.5)
    shift_bin=0.1
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Distance'].values>120.0)].copy().reset_index(drop=True)
        if (types_of_strategy[i]=="Half Nights") or (types_of_strategy[i]=="Late Discovery") or (types_of_strategy[i]=="Late discovery"):
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
        total_prob=strategy_plt_dict['Observation2'].values
        dist=strategy_plt_dict['Distance'].values
        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)

        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 
    plt.xlabel(r'Days After Burst')
    plt.legend()
    plt.clf()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    tt_edges=np.arange(0,4.5,0.5)
    shift_bin=0.1
    for i in range(0,len(types_of_strategy)):
        strategy_plt_dict=strategy_df[(strategy_df['Strategy'].values==types_of_strategy[i])].copy().reset_index(drop=True)
        strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Coverage_deg'].values<100.0)].copy().reset_index(drop=True)
        if types_of_strategy[i]=="Half Nights":
            total_hn_before_cut=len(strategy_plt_dict['Detection Prob'].values)
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
            total_hn_after_cut=len(strategy_plt_dict['Detection Prob'].values)
        if types_of_strategy[i]=="Late Discovery":
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)
        if types_of_strategy[i]=="Late discovery":
            strategy_plt_dict=strategy_plt_dict[(strategy_plt_dict['Detection Prob'].values>0.0)].copy().reset_index(drop=True)

        total_prob=strategy_plt_dict['Observation2'].values
        dist=strategy_plt_dict['Distance'].values
        tt_hist, bin_edges_tt = np.histogram(total_prob,bins=tt_edges,density=False)
        tt_center=bin_edges_tt[:-1]+(float(tt_edges[1]-tt_edges[0])/2.0)
        ax1.step(tt_center+(i*shift_bin),tt_hist,alpha=0.8,where='mid', c=tableau20[i], label=types_of_strategy_name[i]) #hatch='\\'#width=tt_edges[1]-tt_edges[0] 

    plt.xlabel(r'Days After Burst')
    plt.legend()

    for i in range(0,len(types_of_strategy)):

        print(types_of_strategy[i])
        print('These are the telescope times percentiles 50,90,100 and the number of events it contains')
        print(round(tb_50ttarea[i],1))
        print(round(tb_90ttarea[i],1))
        print(round(tb_100ttarea[i],1))
        print(round(n_events_50[i],1))
        print(round(n_events_90[i],1))
        print(tb_n[i])
