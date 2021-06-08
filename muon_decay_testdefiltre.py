# -*- coding: utf-8 -*-
"""
Muon Decay

Spark chamber collaboration
"""

#%%

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
#%%

def IntegrateData_FindPeaks(x,y, l_min=4, l_max=12, c_box=6, seuil=0.005, figshow=False, saveID="C:\\Users\\sasch\\Desktop\\Stage été 2021\\data_muons_0602\\Deux pics\\Test de filtre "):
    """
    Integration des donnees par intervale variant entre l_min et l_max
    Integration des donnees par intervale fixe de grandeur c_box
    Identification des pics inferieurs au seuil
    Enregistrement de la figure si saveID est donne
    Affichage de la figure si figshow est True
    """

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

    peaks,_ = find_peaks(-y_integre_2, height=seuil)
    if peaks.size>2: 
        print("Plus de 2 minima locaux identifiees. \nVerifiez le seuil.")

    ## Figure
    if saveID or figshow:
        fig,ax = plt.subplots(figsize=(10,8))
        ax.plot(y,'k.', label="0")
        ax.plot(y_integre_1, c='b', label=c_box)
        ax.plot(y_integre_2, c='r', label="variable")
        ax.legend(title="intervale d'integration",loc=3,framealpha=1)
        ax.set_xlim(200,600)
        for i in range(peaks.size):
            ax.axvline(peaks[i],c='k',ls=':',alpha=0.8)
        
        if peaks.size==2:
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
            plt.close()
        
        if figshow:
            plt.show()

    return y_integre_1, y_integre_2, peaks, box_arr
#%%
#Analyse
for f in range (10):
    data = np.loadtxt("C:\\Users\\sasch\\Desktop\\Stage été 2021\\data_muons_0602\\Deux pics\\testdefiltre_"+str(f+1)+".txt")
    x = data[:,0]
    y = data[:,1]
    y_var, y_con, ip, box = IntegrateData_FindPeaks(x,y)