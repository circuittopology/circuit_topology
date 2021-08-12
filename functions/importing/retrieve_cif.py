"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that downloads single cif file specified in argument
"""
from Bio.PDB import PDBList
import sys

def retrieve_cif(prot_id):
    server = PDBList(server='ftp://ftp.wwpdb.org', pdb='input_files', obsolete_pdb=None ,verbose=True)
    server.retrieve_pdb_file(prot_id,pdir="input_files/cif",file_format='mmCif', overwrite=True,obsolete= False)
