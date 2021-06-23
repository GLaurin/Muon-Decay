# -*- coding: utf-8 -*-
"""
Height of peaks produced by muons
Sparks chamber collaboration
"""
#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
import os
#%%
file_prefixe = "acq"
file_ext     = ".txt"
folder       = "C:\\Users\\ho57982\\Desktop\\Donées muons\\Data 10-06, 5mV et 1 microseconde"
N_data       = len(os.listdir(folder))
hauteur      = np.zeros(0)
for i in range(N_data):
    file_id = str(i+1)
    if len(file_id)<6: 
        file_id = (6-len(file_id))*"0"+file_id
        file_path = folder + "\\" + file_prefixe + file_id + file_ext
    try:               
        data = np.loadtxt(file_path)
        y = data[:,1]
    except:
        print(file_id)
        continue
    peaks,_ = find_peaks(-y, height = 0.03)
    if peaks.size==1:
        hauteur = np.append(hauteur, y[peaks[0]])
hauteur_moy = hauteur.mean()
print(hauteur_moy)
