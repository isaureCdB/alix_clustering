#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 2 10:28:14 2021

@author: alixdelannoy
"""

import numpy as np
import argparse
import time

#Communication via le terminal de commande
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy', help="np array")
args = parser.parse_args()

#Import
structures = np.load(args.npy)
XXX=args.npy.split("r.npy")[0]


#Script récupéré de fit.py
def fit_xyz(atoms1, atoms2):
    assert len(atoms1.shape) == 2
    assert atoms1.shape == atoms2.shape

    # Must always center the two proteins to avoid affine transformations.
    com1 = atoms1.mean(axis=0)
    com2 = atoms2.mean(axis=0)
    atoms1 = atoms1 - com1 # not -=, because that modifies the coordinates in-line!
    atoms2 = atoms2 - com2

    # Initial residual, see Kabsch.
    e0 = np.einsum("ij,ij->", atoms1, atoms1) + np.einsum("ij,ij->", atoms2, atoms2)
    # equivalent to: (atoms1*atoms1).sum() + (atoms2*atoms2).sum()

    # This beautiful step provides the answer.  V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error.
    v, s, wt = np.linalg.svd(atoms1.T.dot(atoms2))

    # We just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.0 (and thus products are +/- 1.0 ).
    reflect = np.linalg.det(v) * np.linalg.det(wt)
    if reflect < 0:
        s[-1] = -s[-1]
        v[:,-1] = -v[:,-1]

    rotmat = v.dot(wt).T
    pivot = com1 - com2
    sd = e0 - 2.0 * s.sum()
    if sd < 0:
        rmsd = 0
    else:
        rmsd = np.sqrt(sd / len(atoms1))
    return rotmat, com1, pivot, rmsd

start_time = time.time()

rmsdarray=np.zeros((structures.shape[0],structures.shape[0]))
#la matrice des RMSD est symétrique, donc on ne calcule qu'un coté
for i in range(rmsdarray.shape[0]):
    for j in range(i):
        rmsdarray[j,i]=fit_xyz(structures[i],structures[j])[3]

filename = "rmsd"+args.npy.split("r.npy")[0]+".npy"
np.save(filename,rmsdarray)
print("La matrice des RMSD est enregistrée sous "+filename)
print("Temps d'éxécution")
print("--- %s seconds ---" % (time.time() - start_time))