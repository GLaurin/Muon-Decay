# -*- coding: utf-8 -*-
"""
Muon decay time
Sparks chamber collaboration
"""
#%%

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import os
import argparse

#%%

def MuonCount(t,N0,tau):
    return N0*np.exp(-t/tau)

def main(x, y, seuil, figshow = False, figsave = False, saveID = None, saveID_biz = None):
    
    peaks,_ = find_peaks(-y, height = seuil)
    peaks   = np.delete(peaks,np.where(peaks<240)[0])
    ip = 0  
    while ip<peaks.size-1:
        dp = peaks[ip+1]-peaks[ip]
        if dp<300:
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
        
        if saveID == None or saveID_biz == None:
            print("ATTENTION!!! Donner un saveID ou un saveID_biz")
        else:
            if peaks.size == 2 and figsave:
                plt.savefig(saveID)
            if peaks.size > 2  and figsave:     #Si une donnée a plus que deux pics, elle sera enregistrée ailleurs
                plt.savefig(saveID_biz)
        if figshow:
            plt.show()
        else:
            plt.close()
    
    return peaks

#%% Analyse des donnees

parser  = argparse.ArgumentParser(description="Calculer le temps de vie du muon")
parser.add_argument("--fshow",              type = int, choices = [0,1],    default = 0,        help = "Voulez-vous voir les figures? Non = 0, oui = 1")
parser.add_argument("-f_p","--file_prefixe",type = str,                     default = "acq",    help = "Entrez les lettres qui précèdent le numéro de chaque acquisition (ex: acq)")
parser.add_argument("-f_e","--file_ext",    type = str,                     default = ".txt",   help = "Entrez le nom du type de fichier à analyser (ex: .txt, .png, ...)")
parser.add_argument("-s", "--seuil",        type = float,                   default = 0.03,     help = "Entrez le seuil pour filtrer les désintération")
parser.add_argument("--fsave",              type = int,                     default = 1,        help = "Voulez-vous enregistrer les figures? Non = 0, oui = 1")
parser.add_argument("-f","--folder",        type = str,                                         help = "Entrez le chemin où se trouve les fichiers à analyser sur votre ordinateur")
args = parser.parse_args()

N_data       = len(os.listdir(args.folder))
t_decay      = np.zeros(0)   # Temps de desintegration des evenements identifies

for i in range(N_data):
    file_id = str(i+1)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = args.folder + "\\" + args.file_prefixe + file_id + args.file_ext
    
    try:
        data = np.loadtxt(file_path)
        x = data[:,0]   #La première colonne des acquisitions représente le temps
        y = data[:,1]   #La deuxième colonne représente l'amplitude du signal
    except:
        print(file_id)
        continue
    #Paramètres
    sid     = args.folder+"\\Decays\\figure"+file_id+"_seuil"+str(args.seuil)[2:]
    sid_biz = args.folder+"\\Decays\\Figures aberrantes\\figure"+file_id+"_seuil"+str(args.seuil)[2:]   #Dossier dans lequel les figures avec 3 pics ou plus seront enregistrées
    
    ip = main(x, y, args.seuil, figshow = args.fshow, figsave = args.fsave, saveID = sid, saveID_biz = sid_biz) #Appel à a fonction
    
    if ip.size==2:
        print(i)
        t_decay = np.append(t_decay,x[ip[1]]-x[ip[0]])

#np.savetxt("t_decay_1.txt", t_decay) je ferai cette partie lundi
    

#%% Temps de desintegration

#t_decay_1 = np.loadtxt("t_decay_1.txt")
#t_decay_2 = np.loadtxt("t_decay_2.txt")
#t_decay_tot = np.append(t_decay_1,t_decay_2)

N, bins = np.histogram(t_decay, bins='auto')        # Histogramme des desintegrations
t       = bins[:-1]+ 0.5*(bins[1:] - bins[:-1])     # Domaine discret des donnees
t_lin   = np.linspace(t[0]*0.8, t[-1]*1.2)          # Domaine lineaire

## Aujstements
# Curve_fit
popt, pcov = curve_fit(MuonCount, t, N, p0=[N.sum()*10,2.2e-6], sigma=N**0.5)
print(f"Curve_fit: $\ttau$ = {popt[1]:.3e} $\pm$ {np.sqrt(np.diag(pcov))[1]:.2e}")

## Figure
plt.figure(figsize = (9,9))
plt.errorbar(t, N, yerr=N**0.5, marker='o', color='r', ls='', capsize=3, label="Données")
plt.plot(t_lin, MuonCount(t_lin,*popt),'k', label=f"Curve_fit: $\ttau$ = {popt[1]:.3e} $\pm$ {np.sqrt(np.diag(pcov))[1]:.2e}")
plt.legend()
plt.savefig(args.folder+"\\Decays\\Histogramme_désintégrations")
plt.close()

