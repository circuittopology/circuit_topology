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

        psc = [protid,P,S,X]

        return mat,psc

    elif np.shape(index)[1] == 4:

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
                        #parallel
                        elif k < i and j < l:
                            P += 1
                            mat[x,y] = 1
                            mat[y,x] = 1
                            
                        elif i < k and l < j:
                            P += 1
                            mat[x,y] = 1
                            mat[y,x] = 1
                            
                        elif (i==k and j<l):
                            mat[x, y]=1
                            mat[y, x]=1
                            P += 1
                
                        elif (i==k and l<j):
                            mat[x, y]=1
                            mat[y, x]=1
                            P += 1

                        elif (k>i and j==l):
                            mat[x,y]=1
                            mat[y,x]=1
                            P += 1

                        elif(i>k and l==j):
                            mat[x,y]=1
                            mat[y,x]=1
                            P += 1
                        #CS
                        elif j == k:
                            S += 1
                            mat[x,y] = 2
                            mat[y,x] = 2
                        #Cross
                        if (k>i and k<j and j<l):
                            X += 1
                            mat[x, y]=3
                            mat[y, x]=3
                        elif (i>k and i< l and j> l):
                            X += 1
                            mat[x, y]=3
                            mat[y, x]=3
                        #Independent
                    else:
                        I2 += 1
                        mat[x,y] = 4
                        mat[y,x] = 4

                elif chain1 and not set1.intersection(set2):
                    I3 += 1
                    mat[x,y] = 4
                    mat[y,x] = 4

                #I - multiple chains
                elif not set1.intersection(set2):
                    I4 += 1
                    mat[x,y] = 4
                    mat[y,x] = 4
                
                #T
                elif chain1 and set1.intersection(set2) :
                    T2 += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                elif chain2 and set1.intersection(set2):
                    T2 += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                elif ~chain1 and ~chain2 and len(set1.intersection(set2)) == 1:
                    T3 += 1
                    mat[x,y] = 5
                    mat[y,x] = 5
                #L
                elif ~chain1 and ~chain2 and set1 == set2:
                    L += 1
                    mat[x,y] = 6
                    mat[y,x] = 6
                else:
                    print('error - ',i,chaini,j,chainj,k,chaink,l,chainl)
                
        stats = [protid,P,S,X,I2,I3,I4,T2,T3,L]
        return mat,stats    