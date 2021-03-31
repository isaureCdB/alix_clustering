#!/usr/bin/env python3

import sys
import numpy as np
from scipy.spatial.distance import pdist, squareform
from math import sqrt

def write_clustfile(clust, clustfile):
  cf = open(clustfile, "w")
  for cnr, c in enumerate(clust):
    print >> cf, "Cluster %d ->" % (cnr+1), 
    for cc in c: print >> cf, cc,
    print >> cf, ""


#tests:
#comp = np.array( [[1,0,0,0,0],[0,1,1,0,0],[0,1,1,1,0],[0,0,1,1,0],[0,0,0,0,1] ] , dtype=bool)

comp = sys.argv[1] # list of compatible fragments
output_clust = sys.argv[2]

clusters = []
centers = []

n = comp.shape[0]

for steps in range(n*n):
    N_comp = comp.sum(axis=0)
    n_max = np.max(N_comp)
    i_max = np.where(N_comp == n_max)[0][0]
    if i_max == 1:
        print("Converged after "+str(steps)+" steps")
        break
    clust = np.where(comp[i_max])[0]
    centers.append(i_max)
    clusters.append(clust)
    comp[clust] = False
    comp[:,clust] = False
    
write_clustfile(clusters, output_clust)

for c in centers:
    print(c)