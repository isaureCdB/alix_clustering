#!/usr/bin/env python3

import sys
import numpy as np

def get_threshold(rmsd_matrix, intern_coordinate, cutoff):
    nstruc, ncoor = intern_coordinate.shape

    assert rmsd_matrix.shape == (nstruc, nstruc), "The RMSD matrix has not the size of the number of structures"

    intcoor_diff = np.zeros((nstruc, nstruc, ncoor))

    for m in range(3):
        intcoor_diff[:,:,m] = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m]) 
 
    for m in range(3,ncoor):
        dd = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m])
        intcoor_diff[:,:,m] = np.minimum(dd, 360-dd) 
 
    intcoor_diff_inf = intcoor_diff[(0 < rmsd_matrix) & (rmsd_matrix < cutoff)]

    thresholds = np.max(intcoor_diff_inf, axis=0)
    return thresholds

def main():
    rmsd = np.load(sys.argv[1]) #rmsdAAA_test.npy
    intern_coordinate = np.load(sys.argv[2]) #AAA_test_trigo.npy
    cutoff = float(sys.argv[3]) # 1A RMSD cutoff

    thresholds = get_threshold(rmsd, intern_coordinate, cutoff)

    np.save(sys.argv[4], thresholds)

if __name__ == "__main__":
    main()
