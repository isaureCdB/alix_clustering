#!/usr/bin/env python3

import sys
import numpy as np

def get_threshold(rmsd_matrix, intern_coordinate, ndist, cutoff, percent):
    nstruc, ncoor = intern_coordinate.shape

    assert rmsd_matrix.shape == (nstruc, nstruc), "The RMSD matrix has not the size of the number of structures"

    intcoor_diff = np.zeros((nstruc, nstruc, ncoor))

    for m in range(ndist):
        intcoor_diff[:,:,m] = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m]) 
 
    for m in range(ndist,ncoor):
        dd = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m])
        intcoor_diff[:,:,m] = np.minimum(dd, 360-dd) 
 
    intcoor_diff_inf = intcoor_diff[(0 < rmsd_matrix) & (rmsd_matrix < cutoff)]
    
    thresh = []
    p = int(percent*len(intcoor_diff_inf))
    for c in range(ncoor):
        order_diff = np.sort(intcoor_diff_inf[:,c])
        thresh.append(order_diff[p])
    thresholds = np.array(thresh)

    return thresholds


def main():
    rmsd = np.load(sys.argv[1]) #rmsdAAA_test.npy
    intern_coordinate = np.load(sys.argv[2]) #AAA_test_trigo.npy
    ndist = int(sys.argv[3])
    cutoff = float(sys.argv[4]) # 1A RMSD cutoff
    percent = float(sys.argv[5]) # 0.999
    
    thresholds = get_threshold(rmsd, intern_coordinate, ndist, cutoff, percent)

    np.save(sys.argv[6], thresholds)

if __name__ == "__main__":
    main()
