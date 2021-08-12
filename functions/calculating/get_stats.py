"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function for calculating the percentage of entangled contacts further along the diagonal
"""
import numpy as np  

def get_stats(mat):
    if mat.shape == (0,0):
        print('Error - mat empty')
        return np.array([0])
    #Calculates amount of overlapping(entangled), P and X, contacts across the diagonal
    entangled = np.zeros([len(mat),1])
    if mat.max() == 7:
        for i in range(0,len(mat)-1):
            diag = np.diag(mat,k=i)
            entangled[i] = 1 - (sum(diag == 1)+ sum(diag == 7))/len(diag)
    else:
        for i in range(0,len(mat)-1):
            diag = np.diag(mat,k=i)
            entangled[i] = 1 - (sum(diag == 2)/len(diag))

    return entangled