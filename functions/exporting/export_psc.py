"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

For exporting amount of PSC contacts to a csv file. 
also checks whether data came from model or single chain contact map
"""
import pandas as pd

def export_psc(psclist):
    if type(psclist[0]) == str:
        if len(psclist) == 4:
            df = pd.DataFrame(psclist,columns=['protid','P','S','C'])
        else:
            df = pd.DataFrame(psclist,columns=['protid','P','S','C','I','T','L'])
    elif len(psclist[0]) == 4:
        df = pd.DataFrame(psclist,columns=['protid','P','S','C'])
    else:
        df = pd.DataFrame(psclist,columns=['protid','P','S','C','I','T','L'])

    df.to_csv('results/statistics/pscresults.csv')
    print(f'Succesfully exported pscresults.csv')