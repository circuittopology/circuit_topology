import pandas as pd

"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

Function for exporting information on protein circuits to a csv
"""
def export_circuit(circlist):
    df = pd.DataFrame(circlist,columns = ['protid','segnums','meanlength','segends'])
    df.to_csv('results/circuit/circuitlist.csv',index=False)
    print('Succesfully exported circuitlist.csv')