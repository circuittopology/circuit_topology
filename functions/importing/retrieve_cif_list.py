"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that downloads all the cif files specified in the protlist.txt in the Input_files map
"""
from Bio.PDB import PDBList
import sys

def retrieve_cif_list():
    server = PDBList(server='ftp://ftp.wwpdb.org', pdb='input_files', obsolete_pdb=None ,verbose=True)
    pdb_list = open('input_files/protlist.txt','r')
    content = pdb_list.read().split()
    pdb_list.close()

    server.download_pdb_files(content,pdir="input_files/cif",file_format='mmCif', overwrite=True,obsolete= False)
