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
from math import pi

#%% Parsing

'''
Les valeurs de seuil et de dp_min optimales lorsqu'on utilise le scintillateur 2 sont seuil = 0.045 et dp_min = 300. Pour le scintillateur 1, le seuil devrait être 0.03 et le dp_min, 500.
'''

parser  = argparse.ArgumentParser(description="Calculer le temps de vie du muon")
parser.add_argument("-ch", "--clean_h",         type = bool,    default = 0,      help = "Voulez-vous que l'histogramme soit propre/présentable?")
parser.add_argument("--scint",                  type = int,     choices = [1,2],  help = "Les signaux à analyser proviennent du scintillateur 1 ou 2?")
parser.add_argument("-sel_f","--selected_files",type = bool,    default = 0,      help = "Voulez-vous analyser des fichiers préselectionnés? Oui = 1, non = 0.")
parser.add_argument("--fshow",                  type = bool,    default = 0,      help = "Voulez-vous voir les figures? Non = 0, oui = 1.")
parser.add_argument("-f_p","--file_prefixe",    type = str,     default = "acq",  help = "Entrez les lettres qui précèdent le numéro de chaque acquisition (ex: acq)")
parser.add_argument("-f_e","--file_ext",        type = str,     default = ".txt", help = "Entrez le nom du type de fichier à analyser (ex: .txt, .png, ...).")
parser.add_argument("-s", "--seuil",            type = float,   default = 0,      help = "Entrez le seuil pour filtrer les désintérations")
parser.add_argument("--dp_min",                 type = int,     default = 0,      help = "Entrez largeur indicielle minimale d'un pic pour filtrer les désintérations.")
parser.add_argument("--fsave",                  type = int,     default = 1,      help = "Voulez-vous enregistrer les figures? Non = 0, oui = 1")
parser.add_argument("-tdsave", "--save_times",  type = bool,    default = 1,      help = "Enregistrement des temps de desintegration dans un fichier .txt? Non = 0, oui = 1.")
parser.add_argument("-tdID_a", "--tdID_analyse",type = str,     default = "",     help = "Nom du fichier des temps de desintegrations après l'analyse.")
parser.add_argument("-tdID_m", "--tdID_merge",  type = str,     default = "",     help = "Nom du fichier des temps de desintegrations après la fusion.")
parser.add_argument("-tds", "--t_decays",       type = str,     nargs = "+",      help = "Entrez les noms des fichiers que vous voulez fusionner un après les autres avec un espace entre chaque nom (ex: t_decay_05-07-a t_decay_29-06).")
parser.add_argument("-fm", "--folder_merge",    type = str,     default = "",     help = "Entrez le chemin du dossier où se trouve les fichiers à fusionner sur votre ordinateur. Si vous ne voulez pas fusionner de fichiers, n'entrez rien pour cet argument.")
parser.add_argument("-fa","--folder_analyse",   type = str,     default = "",     help = "Entrez le chemin du dossier où se trouve les données à analyser sur votre ordinateur. Si vous ne voulez pas analyser de fichiers, n'entrez rien pour cet argument.")
parser.add_argument("-d","--date",              type = str,                       help = "Entrez la date de l'analyse avec des tirets (ex: 30-06).")
args = parser.parse_args()

'''
Si on veut faire une analyse, il faut absolument préciser:
    -fa, -d, --scint, -sel_f (default 0), -tdID_a
Si on veut faire une fusion de fichiers, il faut absolument préciser:
    -fm, -d, -tdID_m,-tds
'''

#%% fonctions

def MuonCount(t,N0,tau):
    """
    Nombre de muons en fonction du temps; fonction decroissante exponentielle.

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
    if (figshow or figsave) and peaks.size > 1:
        
        fig,ax = plt.subplots(figsize=(10,8))
        ax.plot(y,'k.', label="Données brutes")
        ax.set_xlim(0,2500)
        ax.axhline(-seuil,c='r',label=f"Seuil: {str(seuil)} V")
        ax.axvline(peaks[0]+dp_min, label=f"Distance minimale: {str(dp_min)} canaux")
            
        for i in range(peaks.size):
            ax.axvline(peaks[i],c='k',ls=':',alpha=0.8)
            
        ax.legend(framealpha=1, loc=3)
        
        if peaks.size == 2 and figsave:
            plt.savefig(saveID)
        
        if peaks.size > 2  and figsave:     #Si une donnée a plus que deux pics, elle sera enregistrée ailleurs
            plt.savefig(saveID_biz)
            
        if figshow:
            plt.show()
        else:
            plt.close()
    
    return peaks

def MakeHistogram (t_decay, t_decays, date, folder, seuil=0.03, dp_min=500, scint = None, tdID_analyse = "", tdID_merge = "", ch = 0):
    
    N, bins = np.histogram(t_decay, bins='auto')        # Histogramme des desintegrations
    t       = bins[:-1]+ 0.5*(bins[1:] - bins[:-1])     # Domaine discret des donnees
    t_lin   = np.linspace(t[0]*0.8, t[-1]*1.2)          # Domaine lineaire
    
    ## Ajustement de courbe et chi2
    popt, pcov  = curve_fit(MuonCount, t, N, p0=[N.sum()*10,2.2e-6], sigma=N**0.5)
    tau         = popt[1]
    i_tau       = np.sqrt(np.diag(pcov))[1] #d_tau pour delta tau
    N_fit       = MuonCount(t, *popt)
    chi2        = sum((N_fit-N)**2/(N))
    chi2_norm   = chi2/(N.size-3)
    print(f"chi2_norm = {chi2_norm}")
    
    ## Figure de l'histogramme
    BIGGER_SIZE = 12
    plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
    plt.figure(figsize = (10,10))
    plt.plot(t_lin, MuonCount(t_lin,*popt),'k', label="Courbe:\n  " + r"$\tau$" + f"= {tau:.2e} $\pm$ {i_tau:.1e}\n  "+r"N_0"+f" = {popt[0]:.1e} $\pm$ {np.sqrt(np.diag(pcov))[0]:.0e}\n  "+r"$\chi^2_\nu$"+f"= {round(chi2_norm, 1)}")
    plt.errorbar(t, N, yerr=N**0.5, marker='o', color='r', ls='', capsize=3, label=f"Données:\n  Nombre total de désintégrations = {(t_decay).size}\n  Date: {date}")
    plt.xlabel("Temps (s)")
    plt.ylabel("Nombre de désintégrations")
    plt.legend()
    if ch:
        plt.title("Nombre de désintégrations des muons cosmiques en fonction du temps\npermettant de calculer le temps de vie du muon au repos")
        plt.savefig(f"{folder}\\Histogramme_muons_désintégrations_date{date}")
    elif tdID_merge != "":
        files = ""
        if tdID_analyse != "":
            files = tdID_analyse+", "
        for i in range (len(t_decays)) :
            files = files +f"{t_decays[i]}, "
        plt.title(f"Merging of {files[:-2]}")
        plt.savefig(f"{folder}\\Histogramme_désintégrations_date{date}_fichier_{tdID_merge}")
    else:
        plt.title(f"Scintillateur numéro {scint} - seuil {seuil} - dp_min {dp_min}")
        plt.savefig(f"{folder}\\Decays\\Histogramme_désintégrations_seuil{str(seuil)[2:]}_dpmin{str(dp_min)}_date{str(date)}")
    plt.close()
    return tau, i_tau

#%% Analyse des donnees
if args.folder_analyse != "":
    print("Analyse des acquisitions en cours...")
    if args.scint == None:
        print("Erreur: entrez un numéro de scintillateur (ex: --scint 1). ")
    
    elif args.scint == 2:
        if args.seuil == 0:
            args.seuil = 0.045
        if args.dp_min == 0:
            args.dp_min = 300
    
    elif args.scint == 1:
        if args.seuil == 0:
            args.seuil = 0.03
        if args.dp_min == 0:
            args.dp_min = 500
    
    N_data  = len(os.listdir(args.folder_analyse))  # Nombre de fichiers a lire
    if args.selected_files:
        N_data = len(os.listdir(args.folder_analyse+"\\Decays"))
    t_decay_1 = np.zeros((0,3))                   # Initialisation : temps de desintegration des evenements identifies
    
    for i in range(N_data):
        ## Construction du nom de fichier
        if args.selected_files:
            file_id = os.listdir(args.folder_analyse+"\\Decays")[i][6:12]
        else:
            file_id = str(i+1)
            if len(file_id)<6: 
                file_id = (6-len(file_id))*"0"+file_id
                
        file_path = f"{args.folder_analyse}\\{args.file_prefixe}{file_id}{args.file_ext}"
        
        ## Lecture du fichier
        try:
            data = np.loadtxt(file_path)
            x = data[:,0]   #La première colonne des acquisitions représente le temps
            y = data[:,1]   #La deuxième colonne représente l'amplitude du signal
        except:
            print(f"{file_id} Empty")
            continue
        
        #Recentrer les données
        bruit = y[:200].mean()
        y = y - bruit
        
        ## Construction du nom de figure
        sid     = f"{args.folder_analyse}\\Decays\\figure{file_id}_seuil{str(args.seuil)[2:]}_date{str(args.date)}"
        sid_biz = f"{args.folder_analyse}\\Decays\\Figures aberrantes\\figure{file_id}_seuil{str(args.seuil)[2:]}_date{str(args.date)}"   #Dossier dans lequel les figures avec 3 pics ou plus seront enregistrées
        
        ## Identification des pics du fichier
        ip = FindMuonDecay(x, y, args.seuil, args.dp_min, figshow = args.fshow, figsave = args.fsave, saveID = sid, saveID_biz = sid_biz) #Appel à a fonction
        
        ## MaJ des temps de desintegration si applicable
        if ip.size==2:
            print(f"{i}/{N_data}")
            t_decay_1 = np.concatenate((t_decay_1, np.array([[ x[ip[1]] - x[ip[0]], int(file_id), y[ip[1]]]])))
        elif i%100 == 0:
            print(f"{i}/{N_data}")
    
    ## Enregistrement des temps si applicable
    if args.save_times:
        np.savetxt(f"{args.folder_analyse}\\{args.tdID_analyse}.txt", t_decay_1)
    
    tau, i_tau = MakeHistogram(t_decay_1[:,0], None, args.date, args.folder_analyse, args.seuil, args.dp_min, args.scint, args.tdID_analyse, args.tdID_merge, args.clean_h) #Appel à la fonction

#%% Merging t_decays

if args.folder_merge != "":
    print("Fusion des fichiers des temps de désintégrations en cours...")
    t_decay     = np.zeros(0)
    t_decays    = np.zeros((len(args.t_decays)), dtype=object)    #tableau des noms des fichiers t_decays
    
    for i in range(len(args.t_decays)):
        t_decays[i] = np.loadtxt(f"{args.folder_merge}\\{args.t_decays[i]}")
        
        if len(t_decays[i].shape) > 1:
            t_decays[i] = t_decays[i][:,0]
            
        t_decay = np.concatenate((t_decay,t_decays[i]))
        
    if args.folder_analyse != "":
        t_decay = np.concatenate((t_decay, t_decay_1[:,0]))
     
    t_decay = np.delete(t_decay, t_decay<1e-6)
    
    if args.save_times:
        np.savetxt(f"{args.folder_merge}\\{args.tdID_merge}.txt", t_decay)
        
    tau, i_tau = MakeHistogram(t_decay, args.t_decays, args.date, args.folder_merge, args.seuil, args.dp_min, args.scint, args.tdID_analyse, args.tdID_merge, args.clean_h)

#%% Trouver la constante de couplage de Fermi

#Valeurs
h   = 6.58211915e-25  # Constante de planck sur 2pi en GeV s
c   = 1   # Vitesse de la lumière
mu  = 105.6583745e-3

GF      = (pi/(mu**2*c**5))*(192*pi*h/(mu*tau))**0.5 #En fait, c'est le GF/(hc)**3
i_GF    = i_tau*((pi/mu**2*c**5)*(192*pi*h/(mu))**0.5)*tau**-1.5/2
print(f"La constante de couplage de Fermi expérimentale est {GF:.3e} ± {i_GF:.1e} GeV^-2.")

