"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function for calculating the absolute and relative contact order

input   = chain object, cutoff distance in angstrom
output  = chainlength, contact order, relative contact order
"""

import numpy as np
from Bio.PDB import Selection,NeighborSearch

def contact_order(chain,cutoff_distance):

    #Unpack chain object into atoms and residues
    atom_list = Selection.unfold_entities(chain,'A')
    res_list = Selection.unfold_entities(chain,'R')

    #obtain residue ID's 
    numbering = [res.get_id()[1] for res in res_list]
    numbering = np.array(numbering)
    segment = np.array(range(len(numbering)))

    #Look for neighbouring atoms within cutoff distance
    ns = NeighborSearch(atom_list)
    all_neighbours = ns.search_all(cutoff_distance,'A')

    dist = 0
    count = 0

    #Determine residue locations of atom contacts
    for atompair in all_neighbours:
        res1 = segment[numbering == atompair[0].get_parent().id[1]][0]
        res2 = segment[numbering == atompair[1].get_parent().id[1]][0]

        if abs(res1-res2) > 0:
            count += 1
            dist += abs(res1-res2)
    
    return len(numbering),dist/count,dist/count/len(numbering)