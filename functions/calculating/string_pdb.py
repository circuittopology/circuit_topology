"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function that uses Anatoly's circuit theory to search for circuits within a Residue contact map
"""
import numpy as np
import copy

def string_pdb(index,numbering,threshold):

    if index.shape == (0,):
        print('Error - index is empty')
        segnums = 0
        segends = 0
        meanlength = 0
        return segnums ,meanlength,segends
        
    full_index = np.array(sorted(np.concatenate((index,np.flip(index,axis=1)),axis=0).tolist()))

    y,x = np.array(full_index)[:,0],np.array(full_index)[:,1]

    S = np.zeros([len(y),len(y)],dtype='int')
    cnt = np.ones([len(numbering)],dtype='int')

    for i in range(0,len(y)):
        a1 = np.nonzero(y==x[i])[0]
        n = cnt[x[i]]
        cnt[x[i]] = cnt[x[i]] + 1
        S[i][a1[n-1]] = 1
        
    d = copy.deepcopy(y)
    y,x = np.nonzero(S)
    a1 = np.minimum(x,y)
    s1 = copy.deepcopy(a1)
    s2 = copy.deepcopy(a1) 

    i = 0

    while s2.size != 0:
        s1[s1==min(s2)] = i
        
        s2 = np.delete(s2, s2 == min(s2))
        i = i + 1
        
    b1 = (a1 == x)

    parentheses = np.zeros(len(d),dtype='int')
    j = 0

    for i in range(max(d)+1):

        c1 = d[d==i]
        c2 = b1[d==i]
        c3 = s1[d==i]

        if len(c1) > 1:
            s1[d==i] = np.concatenate([np.flipud(c3[c2]),np.flipud(c3[~c2])])
            parentheses[d==i] = j
            j = j + 1
            
    #delete long-range contacts

    s2 = copy.deepcopy(s1)
    dell = np.zeros(len(s1),dtype='int')

    for i in range(0,len(s1)):
        if (index[s1[i]][1]-index[s1[i]][0]) > threshold:
            dell[i] = 1
    s2 = np.delete(s2,dell != 0)

    #look for circuits
    n = np.zeros([len(s2)],dtype='int')
    cnt = np.zeros([max(s2)+1],dtype='int')
    j = 1

    for i in range(len(s2)):
        n[i] = j
        cnt[s2[i]] = cnt[s2[i]] + 1
        if np.mean(cnt[cnt>0]) == 2:
            j = j + 1
            cnt = np.zeros([max(s2)+1],dtype='int')
        i = i + 1

        #combine consecutive contacts (e.g., 223344)
    num = [np.count_nonzero(n==i) for i in range(1,max(n)+1)]

    c = np.zeros([0,0],dtype='int')

    for i in range(0,len(num)):
        if num[i] > 2:
            c = np.concatenate((c,np.zeros([num[i]],dtype='int')),axis=None)

        else:
            c = np.concatenate((c,np.ones([num[i]],dtype='int')),axis=None)
            
    d = n * c

    for i in range(0,len(d)):
        if d[i] > 0:
            if i > 0 and d[i-1] > 0:
                d[i] = d[i-1]

    n[d>0] = d[d>0]  

    n2 = np.unique(n)

    n = np.array([np.where(i==n2)[0][0] for i in n])   

    segends = []
    for i in range(0,max(n)+1):
        try:
            segends.append(index[s2[np.max(np.where(n==i))]][1])
        except ValueError:
            continue

    segends = np.array(segends) + 1     
    segnums = max(n)+1
    length = np.array(segends[0])
    length = np.append(length,np.diff(segends))
    meanlength = np.around(np.mean(length),3)
    standard_dev = np.std(length,ddof=1)
    threshold_length = meanlength - standard_dev / 2

    #delete those circuits that are too small
    index_circuit = np.zeros(segnums,dtype='int')
    for i in range(len(length)):
        if segnums != 1 and length[i] <= threshold_length:
            index_circuit[i] = 1

    segnums = int(segnums - np.sum(index_circuit))
    segends = segends.tolist()
    
    return segnums,meanlength,segends