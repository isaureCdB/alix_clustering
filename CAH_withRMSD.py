#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 11:56:32 2021

@author: alixdelannoy
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
import argparse

#Communication via le terminal de commande
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy', help="np array")
args = parser.parse_args()

#Import
data_matrix = np.load('rmsdAAA_fast.npy')
#Classification ascendante hiérachique avec seuil de 1Å
model = AgglomerativeClustering(distance_threshold=1,n_clusters=None,affinity='precomputed', linkage='complete').fit(data_matrix)

#fonction pour afficher un dendogramme
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

#Affichage du dendogramme
plt.figure()
plt.title('Hierarchical Clustering Dendrogram')
plot_dendrogram(model)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()

#Affichage du nombre de clusters obtenus
labels=model.labels_
n=model.n_clusters_
print("Nombre de clusters")
print(n)

        