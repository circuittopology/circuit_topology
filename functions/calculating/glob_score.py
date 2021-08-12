# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 17:00:09 2021

@author: scalvinib

Function for calculating advanced globularity score for a protein
"""
import numpy as np

def glob_score(mat):
    #Normalizes the matrix to its maximum possible value (2)
    mat_norm=np.ones((mat.shape[0], mat.shape[1]))*2
    
    #Transforms the matrix so that it only contains parallel-like and cross contacts
    nu_mat=np.zeros((mat.shape[0], mat.shape[1]))
    nu_mat[mat==2]=1
    nu_mat[mat==5]=1
    nu_mat[mat==3]=1
    nu_mat[mat==6]=1
    nu_mat[mat==4]=2

    score_tot=0
    norm_score=0

    #loops over diagonals and calculates the actual score and the maximum possible score 
    #If it existed solely out of cross contacts (2)
    for diagonal in range(1, len(nu_mat)):
        sum_entangled=0   
        all_twos=0
        len_diag=len(nu_mat)-diagonal

        for i in range(len_diag):
            sum_entangled += np.sum( nu_mat[i,i+diagonal])
            all_twos += np.sum(mat_norm[i,i+diagonal])

        sum_entangled = sum_entangled/len_diag  
        all_twos = all_twos/len_diag 
        score = sum_entangled*(diagonal/(len(nu_mat)-1))
        score_twos = all_twos*(diagonal/(len(nu_mat)-1))

        score_tot=score_tot+score
        norm_score=norm_score+score_twos
    if norm_score == False:
        return float("nan")
    #normalizes the actual score with the max possible score
    score_tot=score_tot/norm_score

    return score_tot