# -*- coding: utf-8 -*-
"""
  For recreating the MDS of Brimacombe et al. (2023):
  On the nature of structure in open empirical food webs 

  This script takes both the pairwise GCD11 matrix 
  (gcd11.csv) between all networks and information about 
  each network (metaData_for_mds.csv) and calculates the  
  mean pairwise GCD-11 between networks
"""

import pandas as pd
import networkx as nx
import os
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from os import listdir
from os.path import isfile, join

#Path to directory that you would like to evaluate
mypath = 'D:\\food_web_seperation\\for_upload'
os.chdir(mypath)

#the distance matrix (DGCD-13 matrix)
distanceMatrix =  pd.read_csv('gcd11.csv',index_col=0)

metaData = pd.read_csv("metaData_for_mds.csv")

#storage for the type of interaction
arrayTypeOfInteraction  = np.array([])

#i is the name of the file
for i in distanceMatrix.columns:
    
        #the name of the file
        nameOfFile = i
        
        #the name of the file with the .csv removed
        breakTheName = nameOfFile.split(".")
        
        #obtaining the row
        rowOfInterest = metaData[metaData["name"].str.match(breakTheName[0])]
        
        typeOfInteraction = rowOfInterest.iloc[0]["author"]
        
        if typeOfInteraction == "Alcorlo_et_al_2001":
            typeOfInteraction = int(0)
            
        if typeOfInteraction == "Angelini_et_al_2006":
            typeOfInteraction = int(1)
                
        if typeOfInteraction == "Angelini_et_al_2013":
            typeOfInteraction = int(2)
            
        if typeOfInteraction == "Baeta_et_al_2011":
            typeOfInteraction = int(3)
            
        if typeOfInteraction == "Beaver_1985":
            typeOfInteraction = int(4)
            
        if typeOfInteraction == "Closs_and_Lake_1994":
            typeOfInteraction = int(5)

        if typeOfInteraction == "Cohen_et_al_2003":
            typeOfInteraction = int(6)
            
        if typeOfInteraction == "Fryer_1959":
            typeOfInteraction = int(7)
            
        if typeOfInteraction == "Menge_and_Sutherland_1976":
            typeOfInteraction = int(8)
                
        if typeOfInteraction == "Parker_and_Huryn_2006":
            typeOfInteraction = int(9)
            
        if typeOfInteraction == "Stewart_and_Sprules_2011":
            typeOfInteraction = int(10)
            
        if typeOfInteraction == "Tavares-Cromar_and_Williams_1996":
            typeOfInteraction = int(11)
            
        if typeOfInteraction == "Thompson_and_Townsend_2003":
            typeOfInteraction = int(12)

        if typeOfInteraction == "Thompson_and_Townsend_2004":
            typeOfInteraction = int(13)        
        
        if typeOfInteraction == "Valiela_1974":
            typeOfInteraction = int(14)

        if typeOfInteraction == "One_network_per_publication":
            typeOfInteraction = rowOfInterest.iloc[0]["Primary_type"]
            
            if typeOfInteraction == "Aquatic":
                typeOfInteraction = int(15)
                
            if typeOfInteraction == "Aquatic and terrestrial":
                typeOfInteraction = int(16)

            if typeOfInteraction == "Terrestrial":
                typeOfInteraction = int(17)
            
        arrayTypeOfInteraction = np.append(arrayTypeOfInteraction, typeOfInteraction)

arrayTypeOfInteraction.astype(int)


###############
##2d plot for entropy
###############

embedding=MDS(n_components=2,dissimilarity="euclidean")
here = embedding.fit_transform(distanceMatrix)
from mpl_toolkits import mplot3d 
import numpy as np 
import matplotlib.pyplot as plt 

colors = ['xkcd:pastel purple','xkcd:sky',"xkcd:caramel","xkcd:butter",
              "xkcd:pale aqua","xkcd:pale mauve","xkcd:pale olive","xkcd:azure","xkcd:dust",
              "xkcd:bubblegum","xkcd:watermelon","xkcd:light peach","xkcd:pastel orange","xkcd:light burgundy",
              "xkcd:light periwinkle","xkcd:jade","xkcd:jade","xkcd:jade"] 

plt.rcParams['figure.figsize'] = [7, 7]
plt.rc('font', size=12)
ax = plt.subplot(111, aspect='equal')
nstd = 1

for i in np.unique(arrayTypeOfInteraction).astype(int):
    subset = here[arrayTypeOfInteraction == i]
    x = [row[0] for row in subset]
    y = [row[1] for row in subset]
            
    if i==0:
                
        Alcorlo_et_al_2001 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:purple', linewidths=0.5)
            
    elif i==1:
                
        Angelini_et_al_2006 = ax.scatter(y, x, marker="*", c=colors[i] ,s=30, edgecolors='xkcd:cerulean', linewidths=0.5)
                
    elif i==2:
                
        Angelini_et_al_2013 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:milk chocolate', linewidths=0.5)
                
    elif i==3:
                
        Baeta_et_al_2011 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:dark gold', linewidths=0.5)
                
    elif i==4:
                
        Beaver_1985 = ax.scatter(y, x, c=colors[i], s=25, edgecolors='xkcd:turquoise blue', linewidths=0.5)
                
    elif i==5:
                
        Closs_and_Lake_1994 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:light blue grey', linewidths=0.5)
                
    elif i==6:
                
        Cohen_et_al_2003 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:dusty green', linewidths=0.5)  
        
    elif i==7:
                
        Fryer_1959 = ax.scatter(y, x, c=colors[i], marker="*", s=30, edgecolors='xkcd:lightish blue', linewidths=0.5)  
        
    elif i==8:
        
        Menge_and_Sutherland_1976 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:dark taupe', linewidths=0.5)  
                
    elif i==9:
        
        Parker_and_Huryn_2006 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:barbie pink', linewidths=0.5)  
        
    elif i==10:
        
        Stewart_and_Sprules_2011 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:faded red', linewidths=0.5)  
        
    elif i==11:
            
        Tavares_Cromar_and_Williams_1996 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:pale rose', linewidths=0.5)  
    
    elif i==12:
                
        Thompson_and_Townsend_2003 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:yellowish orange', linewidths=0.5)

    elif i==13:
                
        Thompson_and_Townsend_2004 = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:wine', linewidths=0.5)

    elif i==14:
                
        Valiela_1974 = ax.scatter(y, x ,marker="v",c=colors[i], s=25,edgecolors='xkcd:blue/grey',linewidths=0.5)
                
    elif i==15:
                
        One_network_per_publication = ax.scatter(y, x, marker="*", c=colors[i], s=30, edgecolors='xkcd:teal', linewidths=0.5)
        
    elif i==16:
                    
        One_network_per_publication = ax.scatter(y, x, c=colors[i], s=25, edgecolors='xkcd:teal', linewidths=0.5)
        
    elif i==17:
                    
        One_network_per_publication = ax.scatter(y, x, marker="v", c=colors[i], s=25, edgecolors='xkcd:teal', linewidths=0.5)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
ax.legend([Alcorlo_et_al_2001,Angelini_et_al_2006,Angelini_et_al_2013,Baeta_et_al_2011,Beaver_1985,
           Closs_and_Lake_1994,Cohen_et_al_2003,Fryer_1959,Menge_and_Sutherland_1976,
           Parker_and_Huryn_2006,Stewart_and_Sprules_2011,Tavares_Cromar_and_Williams_1996,Thompson_and_Townsend_2003,
           Thompson_and_Townsend_2004,Valiela_1974,One_network_per_publication],
          ["Alcorlo et al. (2001)","Angelini et al. (2006)","Angelini et al. (2013)","Baeta et al. (2011)","Beaver (1985)",
           "Closs and Lake (1994)","Cohen et al. (2003)","Fryer (1959)","Menge and Sutherland (1976)",
           "Parker and Huryn (2006)","Stewart and Sprules (2011)","Tavares-Cromar and Williams (1996)",
           "Thompson and Townsend (2003)","Thompson and Townsend (2004)","Valiela (1974)","One network per publication"],
          bbox_to_anchor=(1, 0.5),loc="center left",fontsize='xx-small',markerscale=1.8)
plt.savefig("2d_test.pdf",bbox_inches='tight')
plt.show()
plt.draw()