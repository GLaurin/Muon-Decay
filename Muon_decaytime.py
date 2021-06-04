# -*- coding: utf-8 -*-
"""
Muon decay, sparks chamber

Sascha Zakaib-Bernier
"""

"""
Muon Decay
Spark chamber collaboration
"""

#%%
import os
import numpy as np
from matplotlib import pyplot as plt

#%%
#Lecture des documents
file_prefixe = "acq"    
file_ext     = ".txt"
folder       = "C:\\Users\\sasch\\Desktop\\Stage été 2021\\data_muons_0602\\"
for i in range(1000000):
    file_id = str(i+1)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + file_prefixe + file_id + file_ext
#Formation du graphique
        data = np.loadtxt(file_path)
        x = data[:,0]
        y = data[:,1]
        #Filtre 1 (moyenne de chunks de données)
        chunk = 15  #larguer d'intégration
        y_integre = np.zeros(y.size-chunk)  #y intégré
        for j in range(chunk, y_integre.size):
            y_integre[j] = y[j:j+chunk].mean() #moyenne de l'intervalle
        #Filtre 2 (soustraction de différents chunks de données)
        d_diff = 5    #distance des chunks
        y_final = np.zeros(y.size-chunk-d_diff)     #tableau des y finaux, avec le filtre de soustraction
        for k in range(y_final.size):
            y_final[k] = y_integre[k+d_diff]-y_integre[k] 
#Seuil
        for l in range(y_final.size):
            if y_final[l] <= 0,0012:    #Si la valeur de la fonction est plus grande que le seuil, 
                #je voudrais que lr programme enregistre le fichier dans un dossier appart
        break


plt.plot(y)
plt.plot(y_integre)
plt.plot(y_final)

