#!/usr/bin/env python3

import sys
import numpy as np

def write_clustfile(clust, clustfile):
  cf = open(clustfile, "w")
  for cnr, c in enumerate(clust):
    print("Cluster %d ->"%(cnr+1), file=cf, end=" ") 
    for cc in c: 
        print(cc, end = " ",file=cf)
    print("", file=cf)

def get_max(comp):
    a = np.uniq(comp[:,0])

#tests:
#comp = np.array( [[1,0,0,0,0],[0,1,1,0,0],[0,1,1,1,0],[0,0,1,1,0],[0,0,0,0,1] ] , dtype=bool)

comp = np.load(sys.argv[1]) # symetric mask of compatible pairs
output_clust = sys.argv[2]

clusters = []
centers = []

n = comp.shape[0]
left = set(range(n))

#inv_comp = np.flip(comp, axis=1)
#comp = np.append(comp,inv_comp), axis=0

for steps in range(n*n):
    N_comp = comp.sum(axis=0)
    print(N_comp)
    n_max = np.max(N_comp)
    if n_max == 0:
        print("Converged after "+str(steps)+" steps", file=sys.stderr)
        break
    i_max = np.where(N_comp == n_max)[0][0]
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