import matplotlib.pyplot as plt
import numpy as np

def local_topology_plot(index,mat,numbering,protid,siteid,relation):
        plt.ion()
        nseg = len(numbering)
        segment = list(range(0,nseg))  
        ax = plt.subplots()[1]
        plt.plot()
        plt.xlim(0,len(numbering)+1)
        plt.hlines(0,0,len(numbering)+1,colors = 'k',linestyles='dashed')
        ax.tick_params(labelleft = False)
        ax.set_xlabel('Residues')
        #ax.set_title(protid)
        plt.yticks([])
        for row in index:
                X = np.array([row[0],row[1]])
                r = (X[1]-X[0])/2
                center = np.mean(X)
                xx = np.arange(X[0],X[1]+0.001,0.01)
                np.seterr(invalid = 'ignore')
                yy = np.sqrt(r**2 - (xx - center)**2)
                where_nans = np.isnan(yy)
                yy[where_nans] = 0
                if siteid in row:
                        plt.plot(xx,yy,
                        color ='r',linewidth = 5,zorder= 100)
                else:
                        plt.plot(xx,yy,
                        color ='k')
        
        plt.ion()
        nseg = len(numbering)
        segment = list(range(0,nseg))  
        ax = plt.subplots()[1]
        plt.plot()
        plt.xlim(0,len(numbering)+1)
        plt.hlines(0,0,len(numbering)+1,colors = 'k',linestyles='dashed')
        ax.tick_params(labelleft = False)
        ax.set_xlabel('Residues')
        #ax.set_title(protid)
        plt.yticks([])

        if relation == 'X':
                related_contacts  = np.where(mat[np.where(index == siteid)[0],:]==4)[1]
        if relation == 'S':
                related_contacts = np.where((mat[np.where(index == siteid)[0],:]==1) |(mat[np.where(index == siteid)[0],:]==7))[1]
        if relation == 'P':
                related_contacts = np.where((mat[np.where(index == siteid)[0],:]==3) |(mat[np.where(index == siteid)[0],:]==6))[1]
        if relation == 'IP':
                related_contacts = np.where((mat[np.where(index == siteid)[0],:]==5) |(mat[np.where(index == siteid)[0],:]==2))[1]

        for num,row in enumerate(index):
                X = np.array([row[0],row[1]])
                r = (X[1]-X[0])/2
                center = np.mean(X)
                xx = np.arange(X[0],X[1]+0.001,0.01)
                np.seterr(invalid = 'ignore')
                yy = np.sqrt(r**2 - (xx - center)**2)
                where_nans = np.isnan(yy)
                yy[where_nans] = 0
                if num in related_contacts:
                        plt.plot(xx,yy,
                        color ='b',linewidth = 2,zorder= 100)
                else:
                        plt.plot(xx,yy,
                        color ='k')
                        
        plt.ion()
        nseg = len(numbering)
        segment = list(range(0,nseg))  
        ax = plt.subplots()[1]
        plt.plot()
        plt.xlim(0,len(numbering)+1)
        plt.hlines(0,0,len(numbering)+1,colors = 'k',linestyles='dashed')
        ax.tick_params(labelleft = False)
        ax.set_xlabel('Residues')
        #ax.set_title(protid)
        plt.yticks([])
        
        for row in index:
                X = np.array([row[0],row[1]])
                r = (X[1]-X[0])/2
                center = np.mean(X)
                xx = np.arange(X[0],X[1]+0.001,0.01)
                np.seterr(invalid = 'ignore')
                yy = np.sqrt(r**2 - (xx - center)**2)
                where_nans = np.isnan(yy)
                yy[where_nans] = 0
                if siteid in row:
                        plt.plot(xx,yy,
                        color ='r',linewidth = 5,zorder= 100)
                else:
                        plt.plot(xx,yy,
                        color ='lightgray')
                        
        for num,row in enumerate(index):
                X = np.array([row[0],row[1]])
                r = (X[1]-X[0])/2
                center = np.mean(X)
                xx = np.arange(X[0],X[1]+0.001,0.01)
                np.seterr(invalid = 'ignore')
                yy = np.sqrt(r**2 - (xx - center)**2)
                where_nans = np.isnan(yy)
                yy[where_nans] = 0
                if num in related_contacts:
                        plt.plot(xx,yy,
                        color ='b',linewidth = 2,zorder= 99)
                else:
                        plt.plot(xx,yy,
                        color ='lightgray')
                        




        

            