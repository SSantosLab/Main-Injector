import numpy as np; 
import matplotlib.pyplot as plt; 

license="""
   Copyright (C) 2014 James Annis

   This program is free software; you can redistribute it and/or modify it
   under the terms of version 3 of the GNU General Public License as
   published by the Free Software Foundation.

   More to the points- this code is science code: buggy, barely working,
   with little or no documentation. Science code in the the alpine fast 
   & light style.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""

# tru cmap = gray_r
def plot(x, y, vals, gridsize=160, save=False, file="", vmin="", vmax="", cmap="jet", sum=False) :
    if vmin > vals.max() : vmin = ""
    if vmax < vals.min() : vmax = ""
    if vmin == "" and vmax == "":
        plt.hexbin(x,y,vals,gridsize=gridsize, cmap=cmap);  
    elif vmax == "" :
        plt.hexbin(x,y,vals,gridsize=gridsize, vmin=vmin, cmap=cmap);  
    elif vmin == "" :
        plt.hexbin(x,y,vals,gridsize=gridsize, vmax=vmax, cmap=cmap);  
    else :
        if sum:
            print "\t sum!"
            plt.hexbin(x,y,vals,gridsize=gridsize, vmin=vmin, vmax=vmax, cmap=cmap, reduce_C_function=np.sum);  
        else :
            plt.hexbin(x,y,vals,gridsize=gridsize, vmin=vmin, vmax=vmax, cmap=cmap);  
    plt.axes().set_aspect('equal'); 
    plt.colorbar(shrink=0.5,pad=0.03); 
    plt.axes().set_frame_on(False); 
    plt.axes().set_xticks([]); 
    plt.axes().set_yticks([])
    if save :
        plt.savefig(file, bbox_inches='tight')

