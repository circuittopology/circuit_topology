"""
Created on Mon March 28 17:00:09 2022

@author: DuaneM

Function for creating a Residue-Residue based contact map for either a single chain or a whole model.
Contacts are based on the Centre of Geometry of the residues.

Note! this does not produce a contact map but a matrix of the non-zero values in that contact map. 
"""

from Bio.PDB import MMCIFParser,Selection
from scipy.spatial.distance import pdist, squareform
import numpy as np


def get_cmap_cog(chain,cutoff_distance = 4.5,exclude_neighbour=3):

    atom_list = Selection.unfold_entities(chain,"A")

    #Make list of the atom information
    residue_number = np.zeros(len(atom_list),dtype='int')
    coords = np.zeros([len(atom_list),3])
    name = []


    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        name.append(atom.get_name())

    numbering = list(range(residue_number[0],residue_number[-1]+1))
    residue_number = residue_number - numbering[0] 

    duplicate = np.zeros(len(atom_list),dtype='int')
    for i in range(1,len(duplicate)):
        if residue_number[i] == residue_number[i-1] and name[i] == name[i-1]:
            duplicate[i] = 1


    residue_number = residue_number[np.where(duplicate != 1)]
    coords = coords[np.where(duplicate != 1)]

    cog = []
    for i in set(residue_number):
        cog.append(list(np.mean(coords[residue_number == i],axis=0)))

    cmap = squareform(pdist(cog))
    cmap = (cmap < cutoff_distance) * 1

    index1 = np.transpose(np.triu(cmap).nonzero())

    index = []
    for i in index1:
        if abs(i[0]-i[1]) > exclude_neighbour:
            index.append(list(i))

    return np.array(index),numbering
    

