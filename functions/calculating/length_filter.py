"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function applying a length filter to an existing residue contact map index
"""
import numpy as np

def length_filter(index,distance,mode = '<'):
    if index.shape == (0,):
        print('Error - index empty')
        return index
    #checks the indices if they meet the requirement for the lenght filtering
    new_index = []
    if mode == '<':
        for i in index:
            dis = abs(i[0]-i[1])
            if dis <= distance:
                new_index.append(list(i))
    if mode == '>':
        for i in index:
            dis = abs(i[0]-i[1])
            if dis >= distance:
                new_index.append(list(i))
    return np.array(new_index)