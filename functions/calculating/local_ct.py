import numpy as np

def local_ct(index,mat,numbering):
    localct = {}
    for i in range(len(numbering)):
        localct[i] = {"P":0,"IP":0,"X":0,"S":0}
    for i in range(len(numbering)):
        d = np.nonzero(index == i)
        localct[i]['IP'] = len(np.nonzero(np.sum((mat[d[0],:] == 5) | (mat[d[0],:] == 2),0))[0])
        localct[i]['P'] = len(np.nonzero(np.sum((mat[d[0],:] == 3) | (mat[d[0],:] == 6),0))[0])
        localct[i]['X'] = len(np.nonzero(np.sum(mat[d[0],:] == 4,0))[0])
        localct[i]['S'] = len(np.nonzero(np.sum((mat[d[0],:] == 7) | (mat[d[0],:] == 1),0))[0])
    
    return localct