# -*- coding: utf-8 -*-
"""
Filtrage grâce à la tansformations de Fourier

Sparks chamber collaboration
Sascha Zakaib-Bernier
"""
#%%%
import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.signal import find_peaks

#%%
def filtre_fourier(x,y,n,dt,seuil = 0.06, dp_min=3, figsave = True, figshow = False, saveID = None):
    #Première transformation de fourier
    fhat = np.fft.fft (y, n)        #Tableau de transformée de data
    PSD  = fhat*np.conj(fhat)/n     #Tableau du spectre de densité de puissance des différentes fréquences du signal
    freq = (1/(dt*n))*np.arange(n)  #Tableau de vecteurs de fréquences du signal
    #Deuxième transformation de fourier pour trouver la périodicité
    PSD_norm = PSD-PSD.min()        #Procédures pour mettre à niveau les valeurs de l'axe des y du graphique, pour fhathat
    PSD_norm /= PSD_norm.max()      #Lire en haut
    fhathat  = np.fft.fft(PSD_norm[1:n//2],n//2-1)
    PSD_hat  = fhathat*np.conj(fhathat)/(n//2-1)
    
    #Trouver les désintégrations
    peaks,_ = find_peaks(PSD_hat, height=seuil)
    '''
    ip = 0  
    while ip<peaks.size-1:
        dp = peaks[ip+1]-peaks[ip]
        if dp<dp_min:
            peaks = np.delete(peaks,ip+1)
        else: ip += 1
     '''   
    if figsave or figshow:
        fig,ax = plt.subplots(3,1)
        ax[0].plot(x,y, ".k") 
        ax[0].set_xlabel("Temps (s)")
        ax[0].set_ylabel("Signal brut\nAmplitude (V)")
        ax[1].plot(freq[1:n//2], PSD_norm[1:n//2], color = "r")
        ax[1].set_xlabel("Fréquences (Hz)")
        ax[1].set_ylabel("Après la première transformation de Fourier")
        ax[2].plot(PSD_hat[1:PSD_hat.size//2], color = "b")
        ax[2].set_xlabel("Fréquences (Hz)")
        ax[2].set_ylabel("Après la deuxième transformation de Fourier")
        if peaks.size>1 and figsave:
            plt.savefig(saveID)
            print(peaks)
        if figshow:
            plt.show()
        else: 
            plt.close()
        
        return PSD_hat, peaks
        
        
#%%
#Lecture des acquisitions
file_prefixe = "acq"
file_ext     = ".txt"
folder       = "C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\Test filtre fourier\\"
N_data       = len(os.listdir(folder))-1


for i in range(N_data):
    file_id = str(i+18900)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + file_prefixe + file_id + file_ext
    data = np.loadtxt(file_path)

    #Paramètres
    x      = data[:,0]   #La première colonne des acquisitions représente le temps
    y      = data[:,1]   #La deuxième colonne représente l'amplitude du signal
    n      = len (x)     #Le nombre de points dans un acquisition (toujours 2500)
    dt     = 2500/x[-1]  #Le nombre de points par seconde
    seuil  = 0.065      #Le seuil pour le PSD_hat, pour trouver les désintégrations
    dp_min = 20
    fsave  = True        #Voulons-nous sauvegarder les figures?
    fshow  = False       #Voulons-nous montrer les figures?
    sid    = folder+"figures\\fig"+file_id+"_fourier_seuil"+str(seuil)[2:]+"_"     #Nom sous lequel seront enregistrées les figures où il y a désintégration
    
    
    filtre_fourier(x, y, n, dt, seuil, dp_min, figsave = fsave, figshow = fshow, saveID = sid)  #Appel à la fonction
    
