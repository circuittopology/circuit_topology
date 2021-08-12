"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that creates a topological relations matrix plot for a single chain
"""
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import warnings

def matrix_plot(mat,protid):

    #create custom colormap
    newcolors = np.array([[218/255, 219/255, 228/255,1], #grey - 
                      [131/255, 139/255, 197/255,1],    #purple S
                      [172/255,200/255,247/255,1],       #blue P
                      [174/255,213/255,129/255,1],         #mint green P-1
                      [186/255, 155/255, 201/255,1],        #red purple - X
                      [172/255, 200/255, 247/255,1],        #blue P
                      [174/255, 213/255, 129/255,1],        #mint green
                      [131/255,139/255, 197/255,1]])        #purple S
    newcmp = ListedColormap(newcolors)

    fig, ax = plt.subplots()
    color = plt.get_cmap(newcmp, 8)

    #plot data
    pngmat = plt.imshow(mat,cmap=color,vmin = np.min(mat)-.5, vmax = np.max(mat)+.5)
    ax.set_title(protid)
    ax.tick_params(labelleft = False,labelbottom = False,bottom = False,left= False)
    cbar = fig.colorbar(pngmat)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cbar.ax.set_yticklabels(['-','S','P','P-1','X','CP','CP-1','CS'])
