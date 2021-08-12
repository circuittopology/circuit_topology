"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that imports a STRIDE secondary structure and creates a Segment-segment based contact map.
Note! does not produce a contact map but rather an index of the non-zero values in that contact map.
"""

import numpy as np
from Bio.PDB import PDBParser,Selection,NeighborSearch
from collections import Counter

def secondary_struc_cmap(
                        chain,
                        sequence,structure,
                        cutoff_distance = 4.5,
                        cutoff_numcontacts = 10,
                        exclude_neighbour= 3,
                        ss_elements = ['H','E','B','b','G']):

    #unpacks chain object into atom and residues
    atom_list = Selection.unfold_entities(chain,'A')
    res_list = Selection.unfold_entities(chain,'R')

    #extracts relevant protein information
    res_names, numbering = [], []
    for res in res_list:
        res_names.append(res.get_resname())
        numbering.append(res.get_id()[1])

    numbering =  np.array(numbering)
    res_range = np.array(range(len(numbering)))

    assert len(structure) == len(numbering), f'PDB file and Secondary structure map do not match!\n {chain.get_parent().get_parent().id} - PDB: {len(res_list)} Residues VS. STRIDE: {len(sequence)} Residues. '
    #Performs atom neighbour search within set distance
    ns = NeighborSearch(atom_list)
    all_neighbours = ns.search_all(cutoff_distance,'A')

    #Counts the amount of separate secondary structure present in the protein
    struc_length = len(structure)
    segment = np.zeros([struc_length],dtype='int')
    nseg = 1
    for i in range(struc_length):
        if structure[i] in ss_elements:
            segment[i] = nseg
            if i == struc_length:
                nseg += 1
            elif structure[i+1] != structure[i]:
                nseg += 1
    nseg -= 1

    #uses atom contacts to create a secondary structure contact map
    index_list = []
    for atompair in all_neighbours:
        res1 = res_range[numbering == atompair[0].get_parent().id[1]][0]
        res2 = res_range[numbering == atompair[1].get_parent().id[1]][0]

        if abs(res1-res2) > exclude_neighbour:
            if segment[res1] != 0 and segment[res2] != 0 and segment[res1] != segment[res2]:
                index_list.append((segment[res1]-1,segment[res2]-1))
                
    index_list.sort()
    count = Counter(index_list)
    index = [values for values in count if count[values] >= cutoff_numcontacts]

    return np.array(index),segment
        
    