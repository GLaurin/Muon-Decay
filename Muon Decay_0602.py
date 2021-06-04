"""
Muon Decay

Spark chamber collaboration
"""

#%%

import numpy as np
from matplotlib import pyplot as plt

#%%
N_data = 7
for i in [0]:
    data = np.loadtxt("C:\\Users\\lauri\\Documents\\Sessions et Stages\\2021-2É\\Muon Decay_rep\\acq_test\\te465.txt")
    x = data[:,0]
    y = data[:,1]

    l_integre = 15          # Largeur d'intégration
    y_integre = np.zeros(y.size-l_integre)  # Y intégré

    for j in range(y_integre.size):
        y_integre[j] = y[j:j+l_integre].mean()
    
    d_diff    = 15           # distance des blocs 
    y_final   = np.zeros(y.size-l_integre-d_diff)
    for k in range(y_final.size):
        y_final[k] = y_integre[k+d_diff]-y_integre[k]
    break

plt.plot(y)
#plt.plot(y_integre)
plt.plot(y_final)
plt.xlim(200,600)
