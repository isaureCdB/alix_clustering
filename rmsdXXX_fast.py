#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 16:09:43 2021

@author: alixdelannoy
"""

import numpy as np
import argparse
import time

#Communication via le terminal de commande
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('rnpy',help="np array - pseudo-atom representation of XXX")
parser.add_argument('trigonpy', help="np array - intern coordinates representation of XXX")
parser.add_argument('seuilsnpy', help="np array - dictionnary of intern coordinates thresholds")
args = parser.parse_args()

structures=np.load(args.rnpy)
trigo=np.load(args.trigonpy)
seuils=np.load(args.seuilsnpy,allow_pickle='TRUE').item()

start_time = time.time()

#Calcul des matrices de similarités
similarity_d12=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_d23=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_d13=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_bb1=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_bb2=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_bb3=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_mu12=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_mu23=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_mu13=np.zeros((trigo.shape[0],trigo.shape[0]))

id_all=[]
for i in range(trigo.shape[0]):
    for j in range(i):
        #pour les distances on mesure la similarité par une différence
        similarity_d12[j,i]=abs(trigo[i,0]-trigo[j,0]) 
        similarity_d23[j,i]=abs(trigo[i,1]-trigo[j,1]) 
        similarity_d13[j,i]=abs(trigo[i,2]-trigo[j,2])
        #pour les angles on mesure la similarité par une différence des valeurs absolues
        similarity_bb1[j,i]=abs(abs(trigo[i,3])-abs(trigo[j,3]))
        similarity_bb2[j,i]=abs(abs(trigo[i,4])-abs(trigo[j,4]))
        similarity_bb3[j,i]=abs(abs(trigo[i,5])-abs(trigo[j,5]))
        similarity_mu12[j,i]=abs(abs(trigo[i,6])-abs(trigo[j,6]))
        similarity_mu23[j,i]=abs(abs(trigo[i,7])-abs(trigo[j,7]))
        similarity_mu13[j,i]=abs(abs(trigo[i,8])-abs(trigo[j,8]))
        if 0<similarity_d12[j,i]<=seuils['d12'] and 0<similarity_d23[j,i]<=seuils['d23'] and 0<similarity_d23[j,i]<=seuils['d13'] and 0<similarity_bb1[j,i]<=seuils['bb1'] and 0<similarity_bb2[j,i]<=seuils['bb2'] and 0<similarity_bb3[j,i]<=seuils['bb3'] and 0<similarity_mu12[j,i]<=seuils['mu12'] and 0<similarity_mu23[j,i]<=seuils['mu23'] and 0<similarity_mu13[j,i]<=seuils['mu13']:
            id_all.append((j,i))


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

rmsdarray=np.ones((structures.shape[0],structures.shape[0]))*2
for (i,j) in id_all:
    rmsdarray[i,j]=fit_xyz(structures[i],structures[j])[3]

#filename = "rmsd"+args.rnpy.split("r.npy")[0]+".npy"
np.save('rmsdAAA_fast.npy',rmsdarray)

print("Temps d'éxécution")
print("--- %s seconds ---" % (time.time() - start_time))