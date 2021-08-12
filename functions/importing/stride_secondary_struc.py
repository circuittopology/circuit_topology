"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that reads STRIDE files, that are downloaded externally and generates a secondary structure.
"""

# H - Alpha-Helix
# B - Isolated Beta-Bridge
# b - Isolated Beta-Bridge
# G - 3-10 Helix
# I - Pi helix
# T - Turn
# C - Coil

def stride_secondary_struc(stride_file):

    stride_file = 'input_files/stride/' + stride_file
    
    f = open(stride_file,'r')

    chain = 1
    struc = []
    seq = []
    for line in f:
        if line[0:3] == 'CHN':
            if chain != 1:
                break
            else:
                chain += 1
        if line[0:3] == 'SEQ':
            seq.append(line[10:60])
        if line[0:3] == 'STR':
            struc.append(line[10:60])
        if line[0:3] == 'LOC':
            break
            
    seq = ''.join(seq).strip()
    struc = ''.join(struc)

    if len(struc) > len(seq):
            struc = struc[0:len(seq)]

    struc = struc.replace('b','B')

    return struc,seq
    