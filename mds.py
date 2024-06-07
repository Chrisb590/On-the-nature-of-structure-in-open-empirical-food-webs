# -*- coding: utf-8 -*-
"""
  For recreating the MDS of Brimacombe et al. (2024):
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
distanceMatrix =  pd.read_csv('gcd11_with_gateway.csv',index_col=0)

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
            
        if typeOfInteraction == "Cattin_Blandenier_2004":
            typeOfInteraction = int(5)
            
        if typeOfInteraction == "Closs_and_Lake_1994":
            typeOfInteraction = int(6)

        if typeOfInteraction == "Cohen_et_al_2003":
            typeOfInteraction = int(7)
            
        if typeOfInteraction == "Digel_et_al_2014":
            typeOfInteraction = int(8)
            
        if typeOfInteraction == "Fryer_1959":
            typeOfInteraction = int(9)

        if typeOfInteraction == "Havens_1992":
            typeOfInteraction = rowOfInterest.iloc[0]["Primary_type"]
            
            if typeOfInteraction == "Aquatic":
                typeOfInteraction = int(10)
                
            if typeOfInteraction == "Aquatic and terrestrial":
                typeOfInteraction = int(11)    

        if typeOfInteraction == "Layer_et_al_2010":
            typeOfInteraction = int(12)
            
        if typeOfInteraction == "Legagneux_et_al_2014":
            typeOfInteraction = int(13)
            
        if typeOfInteraction == "Menge_and_Sutherland_1976":
            typeOfInteraction = int(14)
            
        if typeOfInteraction == "O_Gorman_et_al_2019":
            typeOfInteraction = int(15)
                
        if typeOfInteraction == "Parker_and_Huryn_2006":
            typeOfInteraction = int(16)
            
        if typeOfInteraction == "Piechnik_et_al_2008":
            typeOfInteraction = int(17)

        if typeOfInteraction == "Stewart_and_Sprules_2011":
            typeOfInteraction = int(18)
            
        if typeOfInteraction == "Tavares-Cromar_and_Williams_1996":
            typeOfInteraction = int(19)
            
        if typeOfInteraction == "Thompson_and_Townsend_2003":
            typeOfInteraction = int(20)

        if typeOfInteraction == "Thompson_and_Townsend_2004":
            typeOfInteraction = int(21)        
        
        if typeOfInteraction == "Valiela_1974":
            typeOfInteraction = int(22)
            
        if typeOfInteraction == "One_network_per_publication":
            typeOfInteraction = rowOfInterest.iloc[0]["Primary_type"]
            
            if typeOfInteraction == "Aquatic":
                typeOfInteraction = int(23)
                
            if typeOfInteraction == "Aquatic and terrestrial":
                typeOfInteraction = int(24)

            if typeOfInteraction == "Terrestrial":
                typeOfInteraction = int(25)
            
        arrayTypeOfInteraction = np.append(arrayTypeOfInteraction, typeOfInteraction)

arrayTypeOfInteraction.astype(int)


###############
##2d plot for entropy
###############

embedding=MDS(n_components=2,dissimilarity="precomputed")
here = embedding.fit_transform(distanceMatrix)
from mpl_toolkits import mplot3d 
import numpy as np 
import matplotlib.pyplot as plt 

colors = ['xkcd:pastel purple','xkcd:sky',"xkcd:caramel","xkcd:butter",
              "xkcd:pale aqua","xkcd:pale mauve","xkcd:sand","xkcd:azure","xkcd:dust",
              "xkcd:bubblegum","xkcd:watermelon","xkcd:watermelon","xkcd:light peach","xkcd:pastel orange",
              "xkcd:light burgundy","xkcd:light periwinkle","xkcd:light royal blue","xkcd:light grey","xkcd:iris",
              "xkcd:golden yellow","xkcd:seafoam green","xkcd:sage","xkcd:pale blue","xkcd:bottle green",
              "xkcd:bottle green","xkcd:bottle green"] 

plt.rcParams['figure.figsize'] = [7, 7]
plt.rc('font', size=15)
ax = plt.subplot(111, aspect='equal')
nstd = 1

for i in np.unique(arrayTypeOfInteraction).astype(int):
    subset = here[arrayTypeOfInteraction == i]
    x = [row[0] for row in subset]
    y = [row[1] for row in subset]
            
    if i==0:
                
        Alcorlo_et_al_2001 = ax.scatter(y, x , marker = "*", c=colors[i], s=60, edgecolors='xkcd:purple', linewidths=1.1)
            
    elif i==1:
                
        Angelini_et_al_2006 = ax.scatter(y, x , marker = "*", c=colors[i] ,s=60, edgecolors='xkcd:cerulean', linewidths=1.1)
                
    elif i==2:
                
        Angelini_et_al_2013 = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:milk chocolate', linewidths=1.1)
                
    elif i==3:
                
        Baeta_et_al_2011 = ax.scatter(y, x , marker = "*", c=colors[i], s=60, edgecolors='xkcd:dark gold', linewidths=1.1)
                
    elif i==4:
                
        Beaver_1985 = ax.scatter(y, x, c=colors[i], s=60, edgecolors='xkcd:turquoise blue', linewidths=1.1)
                
    elif i==5:
                
        Cattin_Blandenier_2004 = ax.scatter(y, x, marker = "v", c=colors[i], s=60, edgecolors='xkcd:light blue grey', linewidths=1.1)
                
    elif i==6:
                
        Closs_and_Lake_1994 = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:sand brown', linewidths=1.1)  
        
    elif i==7:
                
        Cohen_et_al_2003 = ax.scatter(y, x, marker = "*", c=colors[i] , s=60, edgecolors='xkcd:lightish blue', linewidths=1.1)  
        
    elif i==8:
        
        Digel_et_al_2014 = ax.scatter(y, x, marker = "v", c=colors[i], s=60, edgecolors='xkcd:dark taupe', linewidths=1.1)  
                
    elif i==9:
        
        Fryer_1959 = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:barbie pink', linewidths=1.1)  
        
    elif i==10:
        
        Havens_1992_aquatic = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:faded red', linewidths=1.1)  
        
    elif i==11:
        
        Havens_1992_aquatic_terrestrial = ax.scatter(y, x, c=colors[i], s=60, edgecolors='xkcd:faded red', linewidths=1.1)  
        
    elif i==12:
            
        Layer_et_al_2010 = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:pale rose', linewidths=1.1)  
    
    elif i==13:
                
        Legagneux_et_al_2014 = ax.scatter(y, x, marker = "v", c=colors[i], s=60, edgecolors='xkcd:brick orange', linewidths=1.1)

    elif i==14:
                
        Menge_and_Sutherland_1976 = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:wine', linewidths=1.1)

    elif i==15:
                
        O_Gorman_et_al_2019 = ax.scatter(y, x, marker = "*", c=colors[i], s=60,edgecolors='xkcd:blue/grey',linewidths=1.1)
        
    elif i==16:
                
        Parker_and_Huryn_2006 = ax.scatter(y, x, marker = "*", c=colors[i], s=60,edgecolors='xkcd:light navy',linewidths=1.1)

    elif i==17:
                
        Piechnik_et_al_2008 = ax.scatter(y, x ,c=colors[i], s=60,edgecolors='xkcd:steel grey',linewidths=1.1)

    elif i==18:
                
        Stewart_and_Sprules_2011 = ax.scatter(y, x, marker = "*", c=colors[i], s=60,edgecolors='xkcd:indigo blue',linewidths=1.1)       

    elif i==19:
                
        Tavares_Cromar_and_Williams_1996 = ax.scatter(y, x , marker = "*", c=colors[i], s=60,edgecolors='xkcd:mango',linewidths=1.1)   

    elif i==20:
                
        Thompson_and_Townsend_2003 = ax.scatter(y, x , marker = "*", c=colors[i], s=60,edgecolors='xkcd:seaweed green',linewidths=1.1)     
        
    elif i==21:
                
        Thompson_and_Townsend_2004 = ax.scatter(y, x , marker = "*", c=colors[i], s=60,edgecolors='xkcd:slate green',linewidths=1.1)  
        
    elif i==22:
                    
        Valiela_1974 = ax.scatter(y, x , marker = "v", c=colors[i], s=60,edgecolors='xkcd:sky blue',linewidths=1.1)  
                
    elif i==23:
                
        One_network_per_publication = ax.scatter(y, x, marker = "*", c=colors[i], s=60, edgecolors='xkcd:dark green', linewidths=1.1)
        
    elif i==24:
                    
        One_network_per_publication = ax.scatter(y, x, c=colors[i], s=60, edgecolors='xkcd:dark green', linewidths=1.1)
        
    elif i==25:
                    
        One_network_per_publication = ax.scatter(y, x, marker="v", c=colors[i], s=60, edgecolors='xkcd:dark green', linewidths=1.1)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
ax.legend([Alcorlo_et_al_2001,Angelini_et_al_2006,Angelini_et_al_2013,Baeta_et_al_2011,Beaver_1985,Cattin_Blandenier_2004,
           Closs_and_Lake_1994,Cohen_et_al_2003,Digel_et_al_2014,Fryer_1959,Havens_1992_aquatic,
           Havens_1992_aquatic_terrestrial,Layer_et_al_2010,Legagneux_et_al_2014,
           Menge_and_Sutherland_1976,O_Gorman_et_al_2019,Parker_and_Huryn_2006,Piechnik_et_al_2008,Stewart_and_Sprules_2011,
           Tavares_Cromar_and_Williams_1996,Thompson_and_Townsend_2003,Thompson_and_Townsend_2004,Valiela_1974,
           One_network_per_publication],
          ["Alcorlo et al. (2001)","Angelini et al. (2006)","Angelini et al. (2013)","Baeta et al. (2011)","Beaver (1985)",
           "Cattin Blandenier (2004)","Closs and Lake (1994)","Cohen et al. (2003)","Digel et al. (2014)","Fryer (1959)",
           "Havens_aqutic (1992)","Havens_aqutic_terrestrial (1992)","Layer et al. (2010)","Legagneux et al. (2014)",
           "Menge and Sutherland (1976)","O'Gorman et al. (2019)","Parker and Huryn (2006)","Piechnik et al. (2008)",
           "Stewart and Sprules (2011)","Tavares Cromar and Williams (1996)","Thompson and Townsend (2003)",
           "Thompson and Townsend (2004)","Valiela (1974)","One network per publication"],
          bbox_to_anchor=(1, 0.5),loc="center left",fontsize='xx-small',markerscale=1.8)
plt.savefig("2d_test.pdf",bbox_inches='tight')
plt.show()
plt.draw()