"""
Muon Decay

Spark chamber collaboration
"""

#%%

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

#%% Analyse

seuil = 0.005
for i in [0]:
    data = np.loadtxt("C:\\Users\\lauri\\Documents\\Sessions et Stages\\2021-2É\\Muon Decay\\acq_test\\te465.txt")
    x = data[:,0]
    y = data[:,1]

    l_min       = 4                         # Largeur minimale d'integration
    l_max       = 16                        # Largeur maximale d'integration
    l_box       = l_max                     # Largeur d'integration initiale
    c_box       = 7                         # Largeur d'integration constante
    y_integre_1 = np.zeros(y.size)          # y integre
    y_integre_2 = np.zeros(y.size)          # y integre variable
    box_std_p   = 0                         # Initialisation de l'ecart-type d'intervalle
    bkg_std     = np.zeros(l_max-l_min+1)   # Ecart-type seuil du bruit de fond
    box_arr     = np.zeros(y.size)          # Dimensions des intervales d'integration

    for j in range(bkg_std.size):
        l_box      = l_min+j
        bkg_std[j] = np.concatenate((y[0:200][y[0:200].argsort()[:200-(l_box-l_box//2):-1]], y[0:200][y[0:200].argsort()[:l_box//2]])).std()

    ### Integration par intervalle ajustable
    for j in range(y.size):

        ## Valeur initiale
        y_integre_j  = y[j:j+l_box].mean()
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

    ## Figure
    fig,ax = plt.subplots(figsize=(10,8))
    ax.plot(y,'k.', label="0")
    ax.plot(y_integre_1, c='b', label=c_box)
    ax.plot(y_integre_2, c='r', label="variable")
    ax.legend(title="intervale d'integration")
    ax.set_xlim(200,600)
    for i in range(peaks.size):
        ax.axvline(peaks[i],c='k',ls=':',alpha=0.8)
    
    ## Zoom manuel sur le pic
    axin = ax.inset_axes([0.31,0.27,0.6,0.5])
    axin.plot(y,'k.')
    axin.plot(y_integre_1, c='b')
    axin.plot(y_integre_2, c='r')
    axin.set_xlim(500,580)
    axin.set_ylim(-0.04,0.01)
    ax.indicate_inset_zoom(axin, edgecolor="black")
    for i in range(peaks.size):
        axin.axvline(peaks[i],c='k',ls=':',alpha=0.8)

    plt.show()

# %%
