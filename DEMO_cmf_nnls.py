#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 14:21:40 2020

@author: gilleschardon
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.io import loadmat
import damas


mat = loadmat("damasdata2D.mat")

# dictionary of sources
def dictionary(PX, PS, k):
    dx = PX[:, 0:1] - PS[:, 0:1].T
    dy = PX[:, 1:2] - PS[:, 1:2].T
    dz = PX[:, 2:3] - PS[:, 2:3].T

    d = np.sqrt(dx**2 + dy**2 + dz**2);
    
    D = np.exp( -1j * k * d) / d
    
    return D

Data = mat['Data']

# we select the 64 inner microphones (no fundamental reason, just to have
# an array of reasonable size)
N = 64;
Pmic = mat['Pmic']
Norms = Pmic[:, 0]**2 + Pmic[:, 1]**2
idx = np.argsort(Norms)[:N]

Pmic = Pmic[idx, :]

Data = Data[np.ix_(idx, idx)]

# source grid
Lx = 180
Ly = 60
xx = np.linspace(-2, 1, Lx+1)
yy = np.linspace(-1, 0, Ly+1)

xx = xx[:-1]
yy = yy[:-1]

# plotting grid
xxp = np.linspace(-2, 1, Lx+1)
yyp = np.linspace(-1, 0, Ly+1)

Xg, Yg = np.meshgrid(xx, yy);
Xp, Yp = np.meshgrid(xxp, yyp);

Z = 4.4;

# dictionary
D = dictionary(Pmic,
               np.vstack([Xg.ravel(), Yg.ravel(), np.ones(Lx*Ly)*Z]).T,
               mat['k'])


# beamforming
bf_map = np.real(  np.sum((D.conj().T @ Data) * D.T, 1))

# beamforming
DataDR = Data - np.diag(np.diag(Data))
bf_dr_map = np.real(  np.sum((D.conj().T @ DataDR) * D.T, 1))


# CMF-NNLS
cmf_map, _ = damas.cmf_nnls_lh(D, Data)
# CMF-NNLS with diagonal removal
cmf_map_dr, _ = damas.cmf_nnls_dr_lh(D, Data)
# CMF-NNLS with noise estimation
Dnoise = np.hstack([D, np.eye(D.shape[0])])
cmf_map_noise, _ = damas.cmf_nnls_lh(Dnoise, Data)

# DAMAS-NNLS
damas_nnls_map, _ = damas.damas_nnls_lh(D, Data)

damas_nnls_dr_map, _ = damas.damas_nnls_dr_lh(D, Data)

cmfDB = 10*np.log10(cmf_map)

def plotmap(pmap, name, dynrange):
    mapDB = 10*np.log10(pmap)
    plt.figure()
    m = np.max(mapDB)
    plt.pcolor(Xp, Yp, np.reshape(mapDB, [Ly, Lx]), cmap='hot', vmax=m, vmin=m-dynrange)
    current_cmap = plt.cm.get_cmap()
    current_cmap.set_bad(color='black')
    plt.axis('image')
    plt.title(name)
    ax = plt.gca()
    ax.set_facecolor('black')



#%%


plotmap(bf_map, 'Beamforming', 40)
plotmap(cmf_map, 'CMF-NNLS', 40)
plotmap(cmf_map_dr, 'CMF-NNLS DR', 40)
plotmap(damas_nnls_map, 'DAMAS-NNLS', 40)


plt.figure()
plt.stem(cmf_map_noise[Xg.size:])
plt.title('CMF-NNLS noise estimation - noise')


