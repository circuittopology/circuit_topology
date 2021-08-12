"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function for applying an energy filter on an existing Residue contact map

input   =   indices of contact map, array of residue ID's, list of residue names,
            4 letter protein id, energy filter type
output  =   indices of filtered contact map, protein ID 
"""
import numpy as np

def energy_cmap(index,numbering,res_names,protid,potential_sign = '-'):
    
    #Checking if contacts are present
    if index.shape == (0,):
        print('Energy filtering failed - Empty index')
        return index,protid

    #amino acid keys
    n_amino = 20
    aamap = np.array(['CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR','ALA','GLY','THR','SER','GLN','ASN','GLU','ASP','HIS','ARG','LYS','PRO'])

    y = index[:,0]
    x = index[:,1]

    numbering = np.array([numbering]).T

    res_names = np.array(res_names)
    
    resnames1 = res_names[y]
    resnames2 = res_names[x]

    #transform txt file of potential matrix into array
    potential_matrix = np.loadtxt('input_files/matrix_potential.txt','object')
    nan_matrix = potential_matrix == 'n'
    nan_matrix = nan_matrix * 1
    potential_matrix[potential_matrix == 'n'] = 100
    potential_matrix = potential_matrix.astype('float')


    for i in range(0,n_amino):
        for j in range(0,n_amino):
            if potential_matrix[i][j] == 0:
                potential_matrix[i][j] = 10
                potential_matrix[j][i] = 10
            if nan_matrix[i][j] == 1:
                potential_matrix[i][j]=potential_matrix[j][i]

    #linking resnames to res keys and finding score in potential_matrix
    energy_cmap = np.zeros([len(numbering),len(numbering)],'float')

    for i in range(0,len(x)):
            score = potential_matrix[aamap == resnames1[i],aamap == resnames2[i]][0]
            energy_cmap[x[i]][y[i]] = score
            energy_cmap[y[i]][x[i]] = score

    #determine type of filtering
    if potential_sign == '+':
        energy_cmap = (energy_cmap > 0) * 1
    elif potential_sign == '-':
        energy_cmap = (energy_cmap < 0) * 1
    else:
        raise TypeError

    #transform into indices    
    energy_cmap = np.transpose(np.nonzero(np.triu(energy_cmap)))
    protid = f'{protid}_{potential_sign}_ef'

    return energy_cmap,protid

