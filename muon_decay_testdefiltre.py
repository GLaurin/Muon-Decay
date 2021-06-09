# -*- coding: utf-8 -*-
"""
Muon Decay

Sparks chamber collaboration
"""

#%%

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
#%%

def IntegrateData_FindPeaks(x,y, l_min=4, l_max=12, c_box=6, seuil=0.0027, figshow=True, saveID="C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\Muon-Decay\\Desintegrations\\Test avec nouvelle echelle 0806\\Test ", file_id=1, file_prefixe = "aq"):
    """
    Integration des donnees par intervale variant entre l_min et l_max
    Integration des donnees par intervale fixe de grandeur c_box
    Identification des pics inferieurs au seuil
    Enregistrement de la figure si saveID est donne
    Affichage de la figure si figshow est True
    """
    #
    bkg_mean = np.mean(y[0:200])
    
    
    l_box       = l_max                     # Largeur d'integration initiale
    y_integre_1 = np.zeros(y.size)          # y integre
    y_integre_2 = np.zeros(y.size)          # y integre variable
    box_std_p   = 0                         # Initialisation de l'ecart-type d'intervalle
    box_arr     = np.zeros(y.size)          # Dimensions des intervales d'integration

    # Ecart-type seuil du bruit de fond
    bkg_std     = np.zeros(l_max-l_min+1)
    for j in range(bkg_std.size):
        l_box      = l_min+j
        bkg_std[j] = np.concatenate((y[0:200][y[0:200].argsort()[:200-(l_box-l_box//2):-1]], y[0:200][y[0:200].argsort()[:l_box//2]])).std()

    ### Integration par intervalle ajustable
    for j in range(y.size):

        ## Valeur initiale
        box_std      = y[j:j+l_box].std()
        bkg_std_j    = bkg_std[l_box-l_min]

        y_integre_1[j]  = y[j:j+c_box].mean()

        ## Ajustement de l'intervalle d'integration
        if (box_std > box_std_p) and (box_std > bkg_std_j) and (l_box > l_min):
            l_box -= 1
        elif (box_std < box_std_p) and (l_box < l_max):
            l_box += 1
        
        ## Valeur ajustee
        y_integre_2[j] = y[j:j+l_box].mean()
        box_std_p      = y[j:j+l_box].std()
        box_arr[j]     = l_box

    peaks,_ = find_peaks(-y_integre_1, height=seuil-bkg_mean)
    if peaks.size ==2:      #Ici, j'essaie de dire au programme, s'il y a deux maxs qui dépassent le seuil, enregistre le fichier dans un dossier à part, mais ça ne fonctionne pas
        plt.savefig("C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\Muon-Decay\\Desintegrations\\Acquisitions avec desintegrations\\"+file_prefixe+str(file_id))
    if peaks.size>2: 
        print("Plus de 2 minima locaux identifiees. \nVerifiez le seuil.")  #Aussi, ça me dit tout le temps ça, mais quand je regarde les seuils, il y a aucune figure où il y a 3 pics...

    ## Figure
    if saveID or figshow:
        fig,ax = plt.subplots(figsize=(10,8))
        ax.plot(y,'k.', label="0")
        ax.plot(y_integre_1, c='b', label=c_box)
        #ax.plot(y_integre_2, c='r', label="variable")
        ax.legend(title="intervale d'integration",loc=3,framealpha=1)
        ax.set_xlim(200,1000)
        ax.axhline(bkg_mean, c='r')
        ax.axhline(bkg_mean-seuil, c='r')
        for i in range(peaks.size):
            ax.axvline(peaks[i],c='k',ls=':',alpha=0.8)
        
        if peaks.size==10:
            ## Zoom manuel sur le pic
            axin = ax.inset_axes([0.35,0.08,0.6,0.5])
            axin.plot(y,'k.')
            axin.plot(y_integre_1, c='b')
            axin.plot(y_integre_2, c='r')
            axin.set_xlim(peaks[-1]-50,peaks[-1]+50)
            axin.set_ylim(-0.04,0.01)
            ax.indicate_inset_zoom(axin, edgecolor="black")
            for i in range(peaks.size):
                axin.axvline(peaks[i],c='k',ls=':',alpha=0.8)
            axin.axhline(-seuil,c='r',ls=':',alpha=0.8)
            axin.text(570,-seuil*1.1,"seuil",va="top",color='r')

        if saveID:
            plt.savefig(saveID+str(f+1))
            
        
        
        if figshow:
            plt.show()

    return y_integre_1, y_integre_2, peaks, box_arr
#%%
#Analyse
file_prefixe = "aq"    
file_ext     = ".txt"
folder       = "C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\Muon-Decay\\acq_nouv_echelle_0806\\"
for f in range (2,9):
    file_id = str(f)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + file_prefixe + file_id + file_ext
    data = np.loadtxt(file_path)
    
    #Paramètres
    x = data[:,0]
    y = data[:,1]
    l_min = 4
    l_max = 12
    c_box = 6
    seuil = 0.0025
    figshow = True
    saveID = "C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\Muon-Decay\\Desintegrations\\Test avec nouvelle echelle 0806\\Test "
    
    #Appel à la fonction
    y_var, y_con, ip, box = IntegrateData_FindPeaks(x,y, l_min, l_max, c_box, seuil, figshow, saveID, file_id)
    