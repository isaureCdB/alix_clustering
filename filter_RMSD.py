#!/usr/bin/env python3

import numpy as np

#filter_RMSD.py AAAr.npy AAA_trigo_mask.npy AAA_RMSD_mask.npy

coor = np.load(sys.argv[1])
mask1 = np.load(sys.argv[2])
cutoff = float(sys.argv[3])

def fit_multi_npy(a, ref):
    rotation, translation, RMSD = multifit(a, ref)
    rot = np.transpose(rotation, axes=(0,2,1))
    COM = a.sum(axis=1)/a.shape[1]
    centered = a - COM[:,None,:]
    rotated = np.einsum('...ij,...jk->...ik',centered,rot)
    fitted = rotated + COM[:,None,:]
    translated = fitted - translation[:,None,:]
    return translated, RMSD

n = coor.shape()[0]

for i in range(n):
    ref = coor[i]
    tofit = coor[mask1[i]]
    indices = np.argwhere(mask1[i])
    fitted, rmsd = fit_multi_npy(tofit, ref)
    mask1[indices] = rmsd < cutoff

np.save(sys.argv[4], mask1)