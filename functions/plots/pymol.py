import numpy as np
from string import Template

def pymol(
        protid,
        index,
        numbering,
        chain1):
    


    protid = f'{protid[:4]}_{chain1}'
    filterindex = index
    
    indextransformed = []
    
    for i in filterindex:
        res1 = numbering[int(i[0])]
        res2 = numbering[int(i[1])]
        indextransformed.append([res1,res2])
    
    res1 = '+'.join([str(sublist[0]) for sublist in indextransformed])
    res2 = '+'.join([str(sublist[1]) for sublist in indextransformed])

    d = {
    'chain1' : chain1,
    'res1' : res1,
    'res2' : res2}

    with open('functions/plots/template.txt', 'r') as f:
        src = Template(f.read())
        result = src.substitute(d)


    with open(f'{protid}.pml','w') as f:
        f.write(result)
