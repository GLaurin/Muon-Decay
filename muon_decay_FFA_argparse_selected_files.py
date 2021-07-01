# -*- coding: utf-8 -*-
"""
Muon decay time
Sparks chamber collaboration
"""

#%% imports

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import os
import argparse

#%% Parsing

parser  = argparse.ArgumentParser(description="Calculer le temps de vie du muon")
parser.add_argument("--fshow",                  type = bool,    default = 0,      help = "Voulez-vous voir les figures? Non = 0, oui = 1")
parser.add_argument("-f_p","--file_prefixe",    type = str,     default = "acq",  help = "Entrez les lettres qui précèdent le numéro de chaque acquisition (ex: acq)")
parser.add_argument("-f_e","--file_ext",        type = str,     default = ".txt", help = "Entrez le nom du type de fichier à analyser (ex: .txt, .png, ...)")
parser.add_argument("-s", "--seuil",            type = float,   default = 0.03,   help = "Entrez le seuil pour filtrer les désintérations")
parser.add_argument("--dp_min",                 type = int,     default = 400,    help = "Entrez largeur indicielle minimale d'un pic pour filtrer les désintérations")
parser.add_argument("--fsave",                  type = int,     default = 1,      help = "Voulez-vous enregistrer les figures? Non = 0, oui = 1")
parser.add_argument("-tdsave", "--save_times",  type = bool,    default = 0,  help = "Enregistrement des temps de desintegration dans un fichier .txt? Non = 0, oui = 1")
parser.add_argument("-tdID", "--times_file_ID", type = str,     default = "",     help = "Nom du fichier des temps de desintegrations.")
#parser.add_argument()
parser.add_argument("-f","--folder",            type = str,                       help = "Entrez le chemin où se trouve les fichiers à analyser sur votre ordinateur")
parser.add_argument("-d","--date",              type = str,                       help = "Entrez la date de l'analyse avec des tirets (ex: 30-06)")
args = parser.parse_args()

print("Analyse en cours...")

#%% fonctions

def MuonCount(t,N0,tau, k):
    """
    Nombre de muons en fonction du temps; fonction decroissante exponentille.

    Parameters
    ----------
    t : array-like
        Interval de temps sur lequel evaluer le compte de muons
    N0 : float
        Nombre initial de muons
    tau : float
        Temps de vie moyen des muons

    Returns
    -------
    numpy.ndarray

    """
    return N0*np.exp(-t/tau)+k

def FindMuonDecay(x, y, seuil, dp_min, figshow = False, figsave = False, saveID = None, saveID_biz = None):
    """
    Fonction principale trouvant les evenements superieurs au seuil; deux pics suggere une desintegration du muon

    Parameters
    ----------
    x : array-like
        Domaine discret des donnees
    y : array-like; y.size == x.size
        Amplitudes des donnees
    seuil : float
        Valeur absolue du seuil de discrimination des evenements
    dp_min : float
        Largeur min d'un pic
    figshow : bool, optional
        Presente la figure des donnees avec seuil et pics identifies. The default is False.
    figsave : bool, optional
        Enregistrement de la figure s'il y a plus d'un evenement. The default is False.
    saveID : string, optional
        Nom donne a la figure si elle est enregistree. The default is None.
    saveID_biz : string, optional
        Nom donne a la figure si elle contient plus de deux pics identifies. The default is None.

    Returns
    -------
    peaks : np.ndarray
        Indices des pics identifies.

    """
    peaks,_ = find_peaks(-y, height = seuil)                            # Identification premiere
    peaks   = np.delete(peaks,np.where(peaks<240)[0])                   # Suppression de pic avant l'evenement declencheur
    
    ## Application d'une largeur minimale aux pics
    ip = 0
    while ip<peaks.size-1:
        dp = peaks[ip+1]-peaks[ip]
        if dp<dp_min:
            peaks = np.delete(peaks,ip+1)
        else: ip += 1
        
    ## Figures
    if figsave or figshow:
        
        fig,ax = plt.subplots(figsize=(10,8))
        ax.plot(y,'k.', label="Donnees brutes")
        ax.set_xlim(0,2500)
        ax.axhline(-seuil,c='r',label="Seuil")
        ax.axvline(peaks[0]+dp_min)
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

N_data  = len(os.listdir(args.folder+"\\Decays"))  # Nombre de fichiers a lire
t_decay = np.zeros(0)                   # Initialisation : temps de desintegration des evenements identifies

for i in range(N_data):
    
    ## Construction du nom de fichier
    file_id = os.listdir(args.folder+"\\Decays")[i][6:12]
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
    file_path = args.folder + "\\" + args.file_prefixe + file_id + args.file_ext
    
    ## Lecture du fichier
    try:
        data = np.loadtxt(file_path)
        x = data[:,0]   #La première colonne des acquisitions représente le temps
        y = data[:,1]   #La deuxième colonne représente l'amplitude du signal
    except:
        print(file_id)
        continue
        
    ## Construction du nom de figure
    sid     = args.folder+"\\Decays\\figure"+file_id+"_seuil"+str(args.seuil)[2:]+"_dpmin"+str(args.dp_min)+"_date"+str(args.date)
    sid_biz = args.folder+"\\Decays\\Figures aberrantes\\figure"+file_id+"_seuil"+str(args.seuil)[2:]+"_dpmin"+str(args.dp_min)+"_date"+str(args.date)   #Dossier dans lequel les figures avec 3 pics ou plus seront enregistrées
    
    ## Identification des pics du fichier
    ip = FindMuonDecay(x, y, args.seuil, args.dp_min, figshow = args.fshow, figsave = args.fsave, saveID = sid, saveID_biz = sid_biz) #Appel à a fonction
    
    ## MaJ des temps de desintegration si applicable
    if ip.size==2:
        print(f"{i}/{N_data}")
        t_decay = np.append(t_decay,x[ip[1]]-x[ip[0]])
    elif i%50 == 0:
        print(f"{i}/{N_data}")

## Enregistrement des temps si applicable
if args.save_times == True:
    np.savetxt(args.folder + args.times_file_ID + ".txt", t_decay)

#%% Temps de desintegration


#t_decay_1 = np.loadtxt("t_decay_1.txt")
#t_decay_2 = np.loadtxt("t_decay_2.txt")
#t_decay_tot = np.append(t_decay_1,t_decay_2)

N, bins = np.histogram(t_decay, bins='auto')        # Histogramme des desintegrations
t       = bins[:-1]+ 0.5*(bins[1:] - bins[:-1])     # Domaine discret des donnees
t_lin   = np.linspace(t[0]*0.8, t[-1]*1.2)          # Domaine lineaire

## Aujstements
popt, pcov = curve_fit(MuonCount, t, N, p0=[N.sum()*10,2.2e-6, 1], sigma=N**0.5)
N_fit = MuonCount(t, *popt)
chi2 = sum((N_fit-N)**2/(np.sqrt(N)))
chi2_norm = chi2/(N.size-3)
print("chi2_norm = "+ str(chi2_norm))

## Figure
plt.figure(figsize = (9,9))
plt.errorbar(t, N, yerr=N**0.5, marker='o', color='r', ls='', capsize=3, label="Données\nNombre de désintégrations = {t_decay.size}")
plt.plot(t_lin, MuonCount(t_lin,*popt),'k', label="Curve_fit:\nDate: {args.date}\n" + r"$\tau$" + f"= {popt[1]:.3e} $\pm$ {np.sqrt(np.diag(pcov))[1]:.2e}\nN0 = {popt[0]}\nk = {popt[2]}\n"+r"$\chi^2$"+f"= {chi2_norm}")
plt.legend()
plt.savefig(args.folder+"\\Decays\\Histogramme_désintégrations_seuil"+str(args.seuil)[2:]+"_dpmin"+str(args.dp_min)+"_date"+str(args.date))
plt.close()

