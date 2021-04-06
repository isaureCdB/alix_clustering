#!/usr/bin/env python3

import sys
import numpy as np
from random import randrange

def write_clustfile(clust, clustfile):
  cf = open(clustfile, "w")
  for cnr, c in enumerate(clust):
    print("Cluster %d ->"%(cnr+1), file=cf, end=" ") 
    for cc in c: 
        print(cc, end = " ",file=cf)
    print("", file=cf)

#tests:
#comp = np.array( [[1,0,0,0,0],[0,1,1,0,0],[0,1,1,1,0],[0,0,1,1,0],[0,0,0,0,1] ] , dtype=bool)

comp = np.load(sys.argv[1]) # symetric mask of compatible pairs
np.fill_diagonal(comp, 1)
output_clust = sys.argv[2]

clusters = []
centers = []

n = comp.shape[0]
left = set(range(n))

#inv_comp = np.flip(comp, axis=1)
#comp = np.append(comp,inv_comp), axis=0
n_max = n
while n_max > 1:
    N_comp = comp.sum(axis=0)
    print(N_comp)
    n_max = np.max(N_comp)
    #select one of the fragments with the max number of neighbors
    x_max = np.where(N_comp == n_max)[0]
    i_rand = randrange(len(x_max))
    i_max = x_max[i_rand]
    if len(x_max) > 1:
        print("      %i frag with %i neighbors"%(len(x_max), n_max), file=sys.stderr)
    #
    clust = np.where(comp[i_max])[0]
    centers.append(i_max)
    clusters.append(clust)
    comp[clust] = False
    comp[:,clust] = False
    #
    left = left - set(clust)
    a, b = len(left), len(centers)
    print("%i clusters done, %i structures left"%(b,a), file=sys.stderr)

for i in left:
    clusters.append(np.array([i]))

write_clustfile(clusters, output_clust)

for c in centers:
    print(c)