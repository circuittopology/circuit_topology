"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function editing PDB files, this removes heteroatoms.
Note! Only run this in PyMol
"""

from glob import glob
import os 


for file in glob('/Users/duanemoes/Desktop/Desktop/RP1/Code/python/circuit_topology/input_files/cif/*.cif'):
    head,name = os.path.split(file)
    cmd.load(file)
    cmd.remove('hetatm')
    cmd.create(name,'all')
    cmd.save('/Users/duanemoes/Desktop/Desktop/RP1/Code/python/circuit_topology/input_files/pdb/'+name[:-4]+'.pdb',name)
    cmd.delete("all")
    
    

