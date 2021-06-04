"""
Muon Decay

Spark chamber collaboration
"""

#%%

import numpy as np
from matplotlib import pyplot as plt

#%%

for i in [0]:
    data = np.loadtxt("C:\\Users\\lauri\\Documents\\Sessions et Stages\\2021-2É\\Muon Decay\\acq_test\\te465.txt")
    x = data[:,0]
    y = data[:,1]

    l_min     = 4                 # Largeur minimale d'integration
    l_max     = 16                # Largeur maximale d'integration
    l_box     = l_max             # Largeur d'integration initiale
    y_integre_1 = np.zeros(y.size)  # y integre
    y_integre_2 = np.zeros(y.size)  # y integre
    box_std_p = 0                 # Initialisation de l'ecart-type d'intervalle
    bkg_std = np.concatenate((y[0:200][y[0:200].argsort()[:-l_min//2:-1]], y[0:200][y[0:200].argsort()[:l_min//2]])).std()
    # Ecart-type seuil du bruit de fond

    ### Integration par intervalle ajustable
    for j in range(y.size):

        ## Valeur initiale
        y_integre_j  = y[j:j+l_box].mean()
        box_std      = y[j:j+l_box].std()

        y_integre_1[j]  = y[j:j+6].mean()

        ## Ajustement de l'intervalle d'integration
        if (box_std > box_std_p) and (box_std > bkg_std) and (l_box > l_min):
            l_box -= 1
        elif (box_std < box_std_p) and (l_box < l_max):
            l_box += 1
        
        ## Valeur ajustee
        y_integre_2[j] = y[j:j+l_box].mean()
        box_std_p    = y[j:j+l_box].std()

    ## Figure
    plt.figure(figsize=(10,10)) 
    plt.plot(y, label="donnees")
    plt.plot(y_integre_1 +0.25,label="constant avec 6")
    plt.plot(y_integre_2 +0.5, label="variable")
    plt.xlim(200,600)
    plt.legend()

    plt.show()

