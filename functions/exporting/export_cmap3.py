"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

For transforming Residue contact map indices to a contact map and exporting that to a csv file.
"""
import numpy as np 
import pandas as pd 

def export_cmap3(index,protid,numbering):
    cmap = np.zeros([len(numbering),len(numbering)],dtype='int')
    
    for row in index:
        x = row[0]
        y = row[1]
        cmap[x][y] = 1
        cmap[y][x] = 1

    df = pd.DataFrame(cmap)
    df.to_csv('results/circuit_diagram/'+protid+'_cmap3.csv',header = False,index = False) 
    print(f'Succesfully saved {protid}_cmap3.csv')