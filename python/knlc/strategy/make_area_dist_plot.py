import matplotlib
import fits_cat as fc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

def create_color_dist_scatter(xdata,
                              ydata,
                              zdata=None,
                              bin_edges_x=np.arange(1,250,20),
                              bin_edges_y=np.arange(1,200,20),
                              xlims=[0,250],
                              ylims=[0,180],
                              outname="area_dist_color.png",
                              zlevels=[90,80,70,60,50,40],
                              colorzlevels=['Indigo','Purple','DarkViolet','MediumOrchid','Plum','Thistle','Lavender'],
                              markers_sc=["o","v","s","P","*","X","D"],
                              color_hist='Indigo',
                              plot_weights=False,
                              log_scale=False):
    
    fig_ = plt.figure(figsize=(8, 6))         

    grid = plt.GridSpec(2, 2, hspace=0.0, wspace=0.0, 
                        left=0.1, right=1.0, bottom=0.1, top=1.0,
                        height_ratios=[1, 4], width_ratios=[4,1])

    # main axis: predicted statistic vs. true statistic
    ax_main = fig_.add_subplot(grid[1,0])
    
    # hist axis [x]: histogram of data on x axis on top
    ax_hist_x = fig_.add_subplot(grid[0, 0], sharex=ax_main)
    
    # hist axis [y]: histogram of data on y axis on right
    ax_hist_y = fig_.add_subplot(grid[1, 1], sharey=ax_main)

    if plot_weights==True: #density=False, weights=None
        #print("====== Plotting distance weights")
        ax_hist_x.hist(xdata, bin_edges_x,
                   density=True, weights=zdata,
                   color=color_hist, alpha=0.8, 
                   histtype='stepfilled')

        ax_hist_y.hist(ydata, bin_edges_y, 
                   color=color_hist, alpha=0.8,
                   density=True, weights=zdata, 
                   histtype='stepfilled',
                   orientation='horizontal')        
    else:
        ax_hist_x.hist(xdata, bin_edges_x, 
                   color=color_hist, alpha=0.8, 
                   histtype='stepfilled')



        ax_hist_y.hist(ydata, bin_edges_y, 
                   color=color_hist, alpha=0.8,
                   histtype='stepfilled',
                   orientation='horizontal')
        
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
        #if log_scale==True:
        #    ydata=np.log10(ydata)
        ax_main.scatter(xdata,ydata,c=colorzlevels[0],marker='+',s=15)

    ax_main.set_xlabel('Luminosity Distance (Mpc)',fontsize=18)
    #plt.xlabel('ra')
    ax_main.set_ylabel('Area (sq-deg)', fontsize=18)
    xmin,xmax=xlims
    ymin,ymax=ylims  
    ax_main.set_xlim((xmin, xmax))
    #ax_main.set_ylim((0, 3))   
    ax_main.set_ylim((ymin, ymax))


    plt.setp(ax_hist_x.get_xticklabels(), visible=False)
    plt.setp(ax_hist_x.get_yticklabels(), visible=False)
    plt.setp(ax_hist_y.get_yticklabels(), visible=False)
    plt.setp(ax_hist_y.get_xticklabels(), visible=False)
    #ax_hist_y.set_xticks([50,100])
    #yticks(np.arange(0, 100+1, 10))
    ax_main.set_xticks([0,50,100,150,200,250,300,350])
    ax_main.legend(frameon=False, fontsize=16)
    ax_main.tick_params(axis='both', which='major', labelsize=14)
    ax_main.tick_params(axis='both', which='minor', labelsize=14)
    if log_scale==True:
        ax_main.set_yscale('log')
    plt.savefig(outname, dpi=400)
    plt.clf()


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='output tsv file from ligo-skymap-stats')
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='output image plot')
    parser.add_argument('--log-scale',
                        type=bool,
                        help='if True, output plot is in logscale.'
                        )
    
    args = parser.parse_args()
    event_file = args.input
    output_plot = args.output
    log_scale = args.log_scale

    matplotlib.use('agg')
    events = pd.read_csv(event_file,sep     ='\t',header=1)
    area = events['area(90)'].values
    dist = events['dist(50)'].values

    create_color_dist_scatter(xdata=dist,
                              ydata=area,
                              bin_edges_x=np.arange(1,400,20),
                              bin_edges_y=np.arange(1,1000,20),
                              xlims=[-20,200],
                              ylims=[5,300],
                              outname=output_plot,
                            #   zlevels=[90,80,70,60,50,40],
                            zlevels=[],
                              plot_weights=False,
                              log_scale=log_scale)

    # create_color_dist_scatter(xdata=dist,
    #                           ydata=area,
    #                           bin_edges_x=np.arange(1,400,20),
    #                           bin_edges_y=np.arange(1,1000,20),
    #                           xlims=[-20,400],
    #                           ylims=[5,1000],
    #                           outname='GW_SIMS_2D_scatter_areadist.png',
    #                           zlevels=[],
    #                           plot_weights=False,
    #                           log_scale=True)

