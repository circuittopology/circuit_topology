"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

For exporting a topological relations matrix to a csv. 
"""
import numpy as np
import pandas as pd

def export_mat(index,mat,protid):
    if mat.shape == (0,0):
        print(f'Error - mat is empty. Failed to export {protid}')

    if np.shape(index)[1] == 2:

        d = dict({0:'-',1:'S',2:'P',3:'P-1',4:'C',5:'CP',6:'CP-1',7:'CS'})
        c = np.vectorize(d.get)(mat.astype(int))

        names = []
        for i in range(len(index)):
            names.append(f'[{index[i,0]} - {index[i,1]}]')
            
        df = pd.DataFrame(c,columns=names,index=names)
        df.to_csv('results/matrix/'+protid+'_mat.csv',sep=';')
        print(f'Succesfully saved {protid}_mat.csv')

    elif np.shape(index)[1] == 4:

        d = dict({0:'-',1:'P',2:'S',3:'C',4:'I',5:'T',6:'L'})
        c = np.vectorize(d.get)(mat.astype(int))

        names = []
        for i in range(len(index)):
            names.append(f'[{index[i,0]}({index[i,2]}) - {index[i,1]}({index[i,3]})]')
        
        df = pd.DataFrame(c,columns=names,index=names)
        df.to_csv('results/matrix/'+protid+'_mat.csv',sep=';')
        print(f'Succesfully saved {protid}_mat.csv')