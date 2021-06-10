# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 12:51:21 2021

@author: sasch
"""

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

#%%

def IntegrateData_FindPeaks(y, c_box=6, seuil_abs=0.005, figshow=False, saveID=None, figsave=True):
    
    """
    Integration des donnees par intervale fixe de grandeur c_box
    Identification des pics inferieurs au seuil
    Enregistrement de la figure si saveID est donne
    Affichage de la figure si figshow est True
    """

    y_integre = np.zeros(y.size)    # Donnees integre par intervale
    bkg_mean  = np.mean(y[0:200])   # Hauteur zero des donnees
    seuil     = seuil_abs-bkg_mean  # Seuil reel positif

    ## Integration
    for j in range(y.size):
        y_integre[j]  = y[j:j+c_box].mean()

    ## Identification des pics
    peaks,_ = find_peaks(-y_integre, height=seuil)
    ip = 0  
    while ip<peaks.size-1:
        dp = peaks[ip+1]-peaks[ip]
        if dp<3:
            peaks = np.delete(peaks,ip+1)
        else: ip += 1        

    ## Figure
    if figsave or figshow:
        fig,ax = plt.subplots(figsize=(10,8))
        ax.plot(y,'k.', label="Donnees brutes")
        ax.plot(y_integre, c='b', label="Integrees sur "+str(c_box))
        #ax.set_xlim(200,600)
        ax.axhline(-seuil,c='r',label="Seuil")
        for i in range(peaks.size):
            ax.axvline(peaks[i],c='k',ls=':',alpha=0.8)
        ax.legend(framealpha=1, loc=3)

        ## Zoom sur le deuxieme pic
        if peaks.size==2:
            axin = ax.inset_axes([0.35,0.08,0.6,0.5])
            axin.plot(y,'k.')
            axin.plot(y_integre, c='b')
            axin.set_xlim(peaks[-1]-50,peaks[-1]+50)
            ax.indicate_inset_zoom(axin, edgecolor="black")
            for i in range(peaks.size):
                axin.axvline(peaks[i],c='k',ls=':',alpha=0.8)
            axin.axhline(-seuil,c='r')
            axin.text(peaks[-1]+40,-seuil*1.1,"seuil",va="top",color='r')

        if peaks.size>1 and figsave:
            plt.savefig(saveID)

        if figshow:
            plt.show()
        else:
            plt.close()

    return y_integre, peaks

#%% Analyse

#Analyse
file_prefixe = "aq"    
file_ext     = ".txt"
folder       = "C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\Muon-Decay\\acq_nouv_echelle_0806\\"
N_data       = len(os.listdir(folder))-1

for i in range(N_data):
    file_id = str(i+2)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + file_prefixe + file_id + file_ext
    data = np.loadtxt(file_path)

    #Parametres
    x = data[:,0]
    y = data[:,1]
    c_box   = 4
    seuil   = 0.00255
    fshow   = False
    sid     = folder+"figures\\fig"+file_id+"_box"+str(c_box)+"_seuil"+str(seuil)[2:]+"_"
    fsave   = True

    y_integre, ip = IntegrateData_FindPeaks(y, c_box, figsave=fsave, seuil_abs=seuil, figshow=fshow, saveID=sid)