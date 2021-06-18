# -*- coding: utf-8 -*-
"""
Muon decay time
Sparks chamber collaboration
"""
#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
import os
from scipy.optimize import curve_fit
#%%
def MuonCount(t,N0,lambd):
    return N0*np.exp(-lambd*t)

def main(x, y, seuil, figshow = False, figsave = False, saveID = None):
    
    peaks,_ = find_peaks(-y, height = seuil)
    peaks   = np.delete(peaks,np.where(peaks<240)[0])
    ip = 0  
    while ip<peaks.size-1:
        dp = peaks[ip+1]-peaks[ip]
        if dp<30:
            peaks = np.delete(peaks,ip+1)
        else: ip += 1
        
    if figsave or figshow:
        fig,ax = plt.subplots(figsize=(10,8))
        ax.plot(y,'k.', label="Donnees brutes")
        ax.set_xlim(0,2500)
        ax.axhline(-seuil,c='r',label="Seuil")
        for i in range(peaks.size):
            ax.axvline(peaks[i],c='k',ls=':',alpha=0.8)
        ax.legend(framealpha=1, loc=3)

        if peaks.size>1 and figsave:
            if saveID == None:
                print("ATTENTION!!! Donner un saveID")
            else:
                plt.savefig(saveID)
        if figshow:
            plt.show()
        else:
            plt.close()
    
    return peaks
#%%
file_prefixe = "acq"
file_ext     = ".txt"
folder       = "C:\\Users\\sasch\\Desktop\\Stage été 2021\\codes\\muon_decay_FFA-50-50\\"
N_data       = len(os.listdir(folder))-1
t_decay = np.zeros(0)   # Temps de desintegration des evenements identifies

for i in range(N_data):
    file_id = str(i+1)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + file_prefixe + file_id + file_ext
    data = np.loadtxt(file_path)
    
    try:
        x       = data[:,0]   #La première colonne des acquisitions représente le temps
        y       = data[:,1]   #La deuxième colonne représente l'amplitude du signal
    except:
        print(file_id)
        continue
    
    seuil   = 0.02
    fshow   = 0
    fsave   = 1
    sid     = folder+"Decays\\figure"+file_id+"_seuil"+str(seuil)[2:]
    
    ip = main(x, y, seuil, figshow = fshow, figsave = fsave, saveID = sid) #Appel à a fonction
    
    if ip.size==2:
        t_decay = np.append(t_decay,x[ip[1]]-x[ip[0]])
N, bins = np.histogram(t_decay, bins='auto')
t       = bins[:-1]+ 0.5*(bins[1:] - bins[:-1])
plt.plot(t, N, 'or', alpha=1)
plt.savefig(folder+"Decays\\Histogramme_désintégrations")
plt.close()

popt, pcov = curve_fit(MuonCount,N,t)