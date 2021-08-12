"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that filters out residue contacts that are within a specified secondary structures. 
Uses Secondary structure produced in Stride_secondary_struc.py
"""
import numpy as np

def secondary_struc_filter(index,struc,filtered_structures = ['H','E'],ss_elements = ['H','E','B','b','G','b']):
    
    #STRIDE
    # H - Alpha-Helix
    # B - Isolated Beta-Bridge
    # b - Isolated Beta-Bridge
    # G - 3-10 Helix
    # I - Pi helix
    # T - Turn
    # C - Coil

    #DSSP
    # H - Alpha-Helix
    # B - Isolated Beta-Bridge
    # E - Strand
    # G - 3-10 Helix
    # I - Pi helix
    # T - Turn
    # S - Bend

    struc_length = len(struc)

    nstruc = 1
    struc_id = np.zeros([struc_length],dtype='int')
    #counts the number of structures within a protein
    for i in range(0,struc_length):
        if struc[i] in ss_elements:
            struc_id[i] = nstruc
            if i == struc_length-1:
                nstruc = nstruc + 1
            elif struc[i+1] != struc[i]:
                nstruc = nstruc + 1
    nstruc = nstruc - 1
    
    remove = []
    #checks whether a contact is between two residues in the same secondary structure
    for num,i in enumerate(index):
        res1 = i[0]
        res2 = i[1]
        if struc_id[res1] == struc_id[res2] and struc[res1] in filtered_structures:
            remove.append(num)

    index_filtered = np.delete(index,remove,0)

    return index_filtered,struc_id