"""
Analyse de la desintgration des muons
Été 2021
"""

#%%

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy import fft

#%% Fonctions

def MuonCount(t,N0,lambd):
    return N0*np.exp(-lambd*t)

def IntegrateData_FFT_FindPeaks(y, dt, c_box=6, seuil_abs=0.005, figshow=False, saveID=None, figsave=True):
    
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

    ## Transformee de Fourier
    fhat   = np.fft.fft(y_integre, y.size)        # Transformee de Fourier
    PSD_f  = fhat*np.conj(fhat)/y.size            # Spectre de densite de puissance
    freq_f = (1/(dt*y.size))*np.arange(y.size)    # Frequences

    ## Transformee de Fourier de PSD
    ghat   = np.fft.fft(PSD_f, PSD_f.size)        # Transformee de Fourier
    PSD_g  = ghat*np.conj(ghat)/PSD_f.size        # Spectre de densite de puissance
    PSD_gnorm = PSD_g-PSD_g.min()
    PSD_gnorm = PSD_gnorm/PSD_gnorm.max()

    ## Identification des pics
    peaks,_ = find_peaks(-y_integre, height=seuil)
    peaks   = np.delete(peaks,np.where(peaks<240)[0])
    ip = 0  
    while ip<peaks.size-1:
        dp = peaks[ip+1]-peaks[ip]
        if dp<3:
            peaks = np.delete(peaks,ip+1)
        else: ip += 1        

    ## Figure
    if figsave or figshow:
        fig,ax = plt.subplots(3,1,figsize=(10,10))

        ## Donnees et integration
        ax[0].plot(y,'k.', label="Donnees brutes")
        ax[0].plot(y_integre, c='b', label="Integrees sur "+str(c_box))
        ax[0].set_xlim(0,2500)
        ax[0].axhline(-seuil,c='r',label="Seuil")
        for i in range(peaks.size):
            ax[0].axvline(peaks[i],c='k',ls=':',alpha=0.8)
        ax[0].legend(framealpha=1, loc=3)

        ## Zoom sur le deuxieme pic
        if peaks.size==2:
            axin = ax[0].inset_axes([0.35,0.08,0.6,0.5])
            axin.plot(y,'k.')
            axin.plot(y_integre, c='b')
            axin.set_xlim(240,peaks[-1]+50)
            ax[0].indicate_inset_zoom(axin, edgecolor="black")
            for i in range(peaks.size):
                axin.axvline(peaks[i],c='k',ls=':',alpha=0.8)
            axin.axhline(-seuil,c='r')
        
        ## FFT
        ax[1].plot(freq_f[1:y.size//2],PSD_f[1:y.size//2],'b')
        ax[1].set_ylabel('FFT1')
        ax[2].plot(PSD_gnorm[1:y.size//2],'b')
        ax[2].set_ylabel('FFT2')

        if peaks.size>0 and figsave:
            plt.savefig(saveID)
        if figshow:
            plt.show()
        else:
            plt.close()

    return y_integre, peaks

#%% Analyse

## Parametres pour la lecture des fichiers
    # Le folder ne devrait contenir que les donnees et un dossier "figures"
file_prefixe = "aq"
file_ext     = ".txt"
folder       = "C:\\Users\\lauri\\Documents\\Sessions et Stages\\2021-2É\\Muon decay\\acq_nouv_echelle_0806\\"
N_data       = len(os.listdir(folder))-1

t_decay = np.zeros(0)   # Temps de desintegration des evenements identifies

for i in range(N_data):
    file_id = str(i+2)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + file_prefixe + file_id + file_ext
    data = np.loadtxt(file_path)

    #Parametres
    x  = data[:,0]
    y  = data[:,1]
    dt = x[-1]/x.size
    c_box   = 4
    seuil   = 0.0024
    fshow   = 0
    fsave   = 1
    sid     = folder+"figures\\fig_2fft_"+file_id+"_box"+str(c_box)+"_seuil"+str(seuil)[2:]

    y_integre, ip = IntegrateData_FFT_FindPeaks(y, dt, c_box, figsave=fsave, seuil_abs=seuil, figshow=fshow, saveID=sid)
    if ip.size==2:
        t_decay = np.append(t_decay,x[ip[1]]-x[ip[0]])

## Compte de desintegration selon le temps
N, bins = np.histogram(t_decay, bins='auto')
t       = bins[:-1]+ 0.5*(bins[1:] - bins[:-1])
#plt.plot(t, N, 'or', alpha=1)

#popt, pcov = curve_fit(MuonCount,N,t)