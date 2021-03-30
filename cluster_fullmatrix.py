#!/usr/bin/env python3

import sys
import numpy as np
from scipy.spatial.distance import pdist, squareform
from math import sqrt

comp = sys.argv[1] # boolean of compatible fragments
output_clust = npy_file.split(".npy")[0]+"-clust"+str(cutoff)
output_npy = npy_file.split(".npy")[0]+"-clust"+str(cutoff) + ".npy"

def write_clustfile(clust, clustfile):
  cf = open(clustfile, "w")
  for cnr, c in enumerate(clust):
    print("Cluster %d ->" % (cnr+1), end=' ', file=cf)
    for cc in c: print(cc, end=' ', file=cf)
    print("", file=cf)

clust = []
n = comp.shape()[0]
struc = range(n)
centers = []

#comp = np.array( [[1,0,0,0,0],[0,1,1,0,0],[0,1,1,1,0],[0,0,1,1,0],[0,0,0,0,1] ] , dtype=bool)

for steps in range(n*n):
    maxconnected, maxconnections = 0, 0
    Nconnections = comp.sum(axis=0)

    if Nconnections > maxconnections:
        maxconnected, maxconnections = i, Nconnections
    if maxconnections == 1 or len(struc) == 0:
        print("Converged after "+str(steps)+" steps")
        break
    centers.append(struc[maxconnected])
    clusters.append(np.argwhere(comp[maxconnected]))
    struc = [struc[a] for a in range(len(struc)) if rmsd[a][maxconnected] > cutoff]
    rmsd = [[ rmsdtot[a][b] for a in struc] for b in struc ]

centers = centers + struc
clusters = clusters + [ [i] for i in struc]



d = squareform(pdist(clust_struc[:nstruc].reshape(nstruc, nfloat), 'sqeuclidean'))
d2 = d<lim
clustered = 0
while clustered < nstruc:
  neigh = d2.sum(axis=0)
  heart = neigh.argmax()
  leaf = np.where(d2[heart])[0]
  for cs in leaf:
    d2[cs,:] = False
    d2[:, cs] = False
  leaf = [heart+1] + [v+1 for v in leaf if v != heart]
  clust.append(leaf)
  clustered += len(leaf)

write_clustfile(clust, output_clust)
clust_npy = coors[[c[0]-1 for c in clust]]
np.save(output_npy, clust_npy)
