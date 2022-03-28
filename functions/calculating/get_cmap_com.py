"""
Created on Mon March 28 17:00:09 2022

@author: DuaneM

Function for creating a Residue-Residue based contact map for either a single chain or a whole model.
Contacts are based on the Centre of Mass of the residues.
Note! this does not produce a contact map but a matrix of the non-zero values in that contact map. 
"""

from Bio.PDB import MMCIFParser,Selection
from scipy.spatial.distance import pdist, squareform
import numpy as np

def get_cmap_com(chain,cutoff_distance=7,exclude_neighbour=3):
    
    atom_list = Selection.unfold_entities(chain,"A")

    masses = {'H':1.008,'C':12.011,'O':15.999,'N':14.007,'S':32.066,'P':30.974}

    #Make list of the atom information
    residue_number = np.zeros(len(atom_list),dtype='int')
    coords = np.zeros([len(atom_list),3])
    name = []
    atom_element = []

    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        name.append(atom.get_name())
        atom_element.append(atom.element)

    numbering = list(range(residue_number[0],residue_number[-1]+1))
    residue_number = residue_number - numbering[0] 
    atom_element = np.array(atom_element)

    duplicate = np.zeros(len(atom_list),dtype='int')
    for i in range(1,len(duplicate)):
        if residue_number[i] == residue_number[i-1] and name[i] == name[i-1]:
            duplicate[i] = 1


    residue_number = residue_number[np.where(duplicate != 1)]
    coords = coords[np.where(duplicate != 1)]
    atom_element = atom_element[np.where(duplicate != 1)]

    coordsw = []
    m = []

    for num,i in enumerate(coords):
        coordsw.append(i*masses[atom_element[num]])
        m.append(masses[atom_element[num]])

    coordsw = np.array(coordsw)
    m = np.array(m)

    com = []
    for i in set(residue_number):
        coordw = coordsw[residue_number == i]
        mw = m[residue_number == i]
        com.append(sum(coordw)/sum(mw))

    cmap = squareform(pdist(com))
    cmap = (cmap < cutoff_distance) * 1

    index1 = np.transpose(np.triu(cmap).nonzero())

    index = []
    for i in index1:
        if abs(i[0]-i[1]) > exclude_neighbour:
            index.append(list(i))

    return np.array(index),numbering,com
