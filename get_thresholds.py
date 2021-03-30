#!/usr/bin/env python3

import sys
import numpy as np

def get_threshold(rmsd_matrix, intern_coordinate, cutoff):

    nstruc = intern_coordinate.shape[0]
    try:
        assert rmsd_matrix.shape == (nstruc, nstruc)
        intcoor_diff = np.abs(intern_coordinate[:,None, :]-intern_coordinate[None, :, :])
        intcoor_diff_inf = intcoor_diff[rmsd_matrix < cutoff]
        thresholds = np.max(intcoor_diff_inf, axis=0)
        return thresholds
    except:
        print("The RMSD matrix has not the size of the number of structures.")

def main():
    rmsd = np.load(sys.argv[1]) #rmsdAAA_test.npy
    intern_coordinate = np.load(sys.argv[2]) #AAA_test_trigo.npy
    cutoff = float(sys.argv[3]) # 1A RMSD cutoff

    nstruc, ncoor = intern_coordinate.shape
    assert rmsd.shape == (nstruc, nstruc)

    intcoor_diff = np.abs(intern_coordinate[:,None, :]-intern_coordinate[None, :, :])

    intcoor_diff_inf = intcoor_diff[rmsd < cutoff]

    thresholds = np.max(intcoor_diff_inf, axis=0)

    np.save(sys.argv[4], thresholds)

if __name__ == "__main__":
    main()
