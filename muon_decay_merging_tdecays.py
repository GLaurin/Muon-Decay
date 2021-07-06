# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 16:39:54 2021

@author: sasch
"""
#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import argparse

#%% Parsing

'''
Les valeurs de seuil et de dp_min optimales lorsqu'on utilise le scintillateur 2 sont seuil = 0.045 et dp_min = 300. Pour le scintillateur 1, une recherche poussée du seuil et du dp_min n'a pas été effectuée.
'''

parser  = argparse.ArgumentParser(description="Calculer le temps de vie du muon")
parser.add_argument("-f_e","--file_ext",        type = str,     default = ".txt", help = "Entrez le nom du type de fichier à analyser (ex: .txt, .png, ...)")
parser.add_argument("-tdsave", "--save_times",  type = bool,    default = 0,      help = "Enregistrement des temps de desintegration dans un fichier .txt? Non = 0, oui = 1")
parser.add_argument("-fIDtot", "--file_ID_tot", type = str,     default = "",     help = "Entrez le nom que vous voulez donner au nouveau fichier t_decay au besoin (ex: t_decay_06-07)")
parser.add_argument("-fID1", "--file_ID_1",     type = str,                       help = "Entrez le nom du premier fichier t_decay (ex: t_decay_05-07)")
parser.add_argument("-fID2", "--file_ID_2",     type = str,                       help = "Entrez le nom du deuxième fichier t_decay (ex: t_decay_08-06)")
parser.add_argument("-f","--folder",            type = str,                       help = "Entrez le chemin où se trouve les fichiers à analyser sur votre ordinateur")
parser.add_argument("-d","--date",              type = str,                       help = "Entrez la date de l'analyse avec des tirets (ex: 30-06)")
args = parser.parse_args()

        
#%%
def MuonCount(t,N0,tau):
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
    return N0*np.exp(-t/tau)
#%%
t_decay_1 = np.loadtxt(f"{args.folder}\\{args.file_ID_1}{args.file_ext}")
t_decay_2 = np.loadtxt(f"{args.folder}\\{args.file_ID_2}{args.file_ext}")
print("Analyse en cours...")
t_decay   = np.concatenate((t_decay_1[:,0],t_decay_2[:,0]))
for i in range (t_decay.size):
    if float(t_decay[i]) <= 1e-6:
        try:
            t_decay = np.delete(t_decay, i)
        except:
            print("L'acquisition {i} est problématique")
            continue
if args.save_times:
    np.savetxt(f"{args.folder}\\{args.file_ID_tot}.txt", t_decay)

N, bins = np.histogram(t_decay, bins='auto')        # Histogramme des desintegrations
t       = bins[:-1]+ 0.5*(bins[1:] - bins[:-1])     # Domaine discret des donnees
t_lin   = np.linspace(t[0]*0.8, t[-1]*1.2)          # Domaine lineaire

## Ajustement de courbe et chi2
popt, pcov  = curve_fit(MuonCount, t, N, p0=[N.sum()*10,2.2e-6], sigma=N**0.5)
N_fit       = MuonCount(t, *popt)
chi2        = sum((N_fit-N)**2/(N))
chi2_norm   = chi2/(N.size-3)
print(f"chi2_norm = {chi2_norm}")

## Figure de l'histogramme
plt.figure(figsize = (9,9))
plt.plot(t_lin, MuonCount(t_lin,*popt),'k', label=f"Curve_fit:\nDate: {args.date}\n" + r"$\tau$" + f"= {popt[1]:.3e} $\pm$ {np.sqrt(np.diag(pcov))[1]:.2e}\nN0 = {popt[0]:.3e} $\pm$ {np.sqrt(np.diag(pcov))[0]:.2e}\n"+r"$\chi^2_\nu$"+f"= {round(chi2_norm, 3)}")
plt.errorbar(t, N, yerr=N**0.5, marker='o', color='r', ls='', capsize=3, label=f"Données\nNombre de désintégrations = {(t_decay).size}")
plt.legend()
plt.title(f"Merging of t_decay_1: {args.file_ID_1} and t_decay_2: {args.file_ID_2}")
plt.savefig(f"{args.folder}\\Histogramme_désintégrations_date{str(args.date)}_fichiers_{args.file_ID_tot}")
plt.close()