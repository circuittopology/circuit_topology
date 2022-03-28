"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function creating a topological relationship matrix for a Residue contact map, 
using either a single chain or a whole model.
"""
import numpy as np

def get_matrix(index,protid):
    
    if index.shape == (0,):
        print('Error - index empty')
        mat = np.zeros((len(index), len(index)),dtype = 'int')
        psc = [protid,0,0,0]
        return mat,psc
    #Determines whether index came from model or single chain
    if np.shape(index)[1] == 2:

        #create a numerical and character matrix based on the amount of nonzero values found in the previous function
        mat = np.zeros((len(index), len(index)),dtype = 'int')

        #Change the values based on the type of connection
        

        P = 0
        S = 0
        X = 0

        for x in range(0,len(index)):
            i = index[x,0]
            j = index[x,1]
            for y in range(x+1,len(index)):
                k = index[y,0]
                l = index[y,1]
                #series
                if (j < k):
                    S=S+1
                    mat[x, y]=1
                    mat[y, x]=1

                #parallel    
                elif (i>k and j<l):
                    P=P+1
                    mat[x, y]=2
                    mat[y, x]=3
                
                #5: CP
                #6: CP-1    
                elif (i==k and j<l):
                    mat[x, y]=5
                    mat[y, x]=6
                    P += 1
            
                elif (i==k and l<j):
                    mat[x, y]=6
                    mat[y, x]=5
                    P += 1
                
                elif (k>i and j==l):
                    mat[x,y]=6
                    mat[y,x]=5
                    P += 1
                    
                elif(i>k and l==j):
                    mat[x,y]=5
                    mat[y,x]=6
                    P += 1
                #inverse parallel
                elif (k>i and l<j):
                    P += 1
                    mat[x, y]=3
                    mat[y, x]=2
                #CS
                elif (j ==k):
                    mat[x,y]=7
                    mat[y,x]=7
                    S += 1
                #Cross
                if (k>i and k<j and j<l):
                    X += 1
                    mat[x, y]=4
                    mat[y, x]=4
                elif (i>k and i< l and j> l):
                    X += 1
                    mat[x, y]=4
                    mat[y, x]=4
        total = sum([P,S,X])
        psc = [protid,P,S,X,round(P/total,3),round(S/total,3),round(X/total,3)]

        return mat,psc

    elif np.shape(index)[1] == 4:
        
        chainids = np.unique(index[:,2:])
        chainstats = {}
        for i in chainids:
            chainstats[i] = {'p':0,'s':0,'x':0,'i2':0,'i3':0,'i4':0,'t2':0,'t3':0,'l':0}
            
        mat = np.zeros((len(index), len(index)),dtype = 'int')
    
        P = 0
        S = 0
        X = 0
        I2 = 0
        I3 = 0
        I4 = 0
        T2 = 0
        T3 = 0
        L = 0
        
        
        for x in range(len(index)):
            chain1 = False
            
            i = index[x][0]
            j = index[x][1]
            chaini = index[x][2]
            chainj = index[x][3]
            
            if chaini == chainj:
                chain1 = True
                
            for y in range(x+1,len(index)):
                chain2 = False
                
                k = index[y][0]
                l = index[y][1]
                chaink = index[y][2]
                chainl = index[y][3]
                
                set1 = set([chaini,chainj])
                set2 = set([chaink,chainl])

                if chaink == chainl:
                    chain2 = True
                    
                if chain1 and chain2:
                    if chaini == chaink:
                        #series
                        if j < k:
                            S += 1
                            mat[x,y] = 2
                            mat[y,x] = 2
                            chainstats[chaink]['s'] += 1
                            
                        #parallel
                        elif k < i and j < l:
                            P += 1
                            mat[x,y] = 1
                            mat[y,x] = 1
                            chainstats[chaink]['p'] += 1
                                
                        elif i < k and l < j:
                            P += 1
                            mat[x,y] = 1
                            mat[y,x] = 1
                            chainstats[chaink]['p'] += 1
                            
                        elif (i==k and j<l):
                            mat[x, y]=1
                            mat[y, x]=1
                            P += 1
                            chainstats[chaink]['p'] += 1
                            
                        elif (i==k and l<j):
                            mat[x, y]=1
                            mat[y, x]=1
                            P += 1
                            chainstats[chaink]['p'] += 1
                            
                        elif (k>i and j==l):
                            mat[x,y]=1
                            mat[y,x]=1
                            P += 1
                            chainstats[chaink]['p'] += 1
                            
                        elif(i>k and l==j):
                            mat[x,y]=1
                            mat[y,x]=1
                            P += 1
                            chainstats[chaink]['p'] += 1
                        #CS
                        elif j == k:
                            S += 1
                            mat[x,y] = 2
                            mat[y,x] = 2
                            chainstats[chaink]['s'] += 1
                            
                        #Cross
                        if (k>i and k<j and j<l):
                            X += 1
                            mat[x, y]=3
                            mat[y, x]=3
                            chainstats[chaink]['x'] += 1
                            
                        elif (i>k and i< l and j> l):
                            X += 1
                            mat[x, y]=3
                            mat[y, x]=3
                            chainstats[chaink]['x'] += 1
                    #Independent
                    else:
                        I2 += 1
                        mat[x,y] = 4
                        mat[y,x] = 4
                        chainstats[chaink]['i2'] += 1
                        chainstats[chainj]['i2'] += 1

                elif chain1 and not set1.intersection(set2):
                    I3 += 1
                    mat[x,y] = 4
                    mat[y,x] = 4
                    chainstats[chaini]['i3'] += 1
                    chainstats[chaink]['i3'] += 1
                    chainstats[chainl]['i3'] += 1
                
                elif chain2 and not set1.intersection(set2):
                    I3 += 1
                    mat[x,y] = 4
                    mat[y,x] = 4
                    chainstats[chaini]['i3'] += 1
                    chainstats[chaink]['i3'] += 1
                    chainstats[chainl]['i3'] += 1

                #I - multiple chains
                elif not set1.intersection(set2):
                    I4 += 1
                    mat[x,y] = 4
                    mat[y,x] = 4
                    chainstats[chaini]['i4'] += 1
                    chainstats[chainj]['i4'] += 1
                    chainstats[chaink]['i4'] += 1
                    chainstats[chainl]['i4'] += 1
                
                #T
                elif chain1 and set1.intersection(set2) :
                    T2 += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                    chainstats[chaini]['t2'] += 1
                    chainstats[list(set2-set1)[0]]['t2'] += 1
                    
                elif chain2 and set1.intersection(set2):
                    T2 += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                    chainstats[chaink]['t2'] += 1
                    chainstats[list(set1-set2)[0]]['t2'] += 1
                elif ~chain1 and ~chain2 and len(set1.intersection(set2)) == 1:
                    T3 += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                    chainstats[list(set1.intersection(set2))[0]]['t3'] += 1
                    chainstats[list(set2 - set1)[0]]['t3'] += 1
                    chainstats[list(set1 - set2)[0]]['t3'] += 1
                #L
                elif ~chain1 and ~chain2 and set1 == set2:
                    L += 1
                    mat[x,y] = 6
                    mat[y,x] = 6
                    chainstats[list(set1.intersection(set2))[0]]['l'] += 1
                    chainstats[list(set1.intersection(set2))[1]]['l'] += 1
                else:
                    print('error - ',i,chaini,j,chainj,k,chaink,l,chainl)
                
        stats = [protid,P,S,X,I2,I3,I4,T2,T3,L]
        return mat,stats,chainstats    