#!/usr/bin/env python3

import sys
import numpy as np

def get_threshold(rmsd_matrix, intern_coordinate, cutoff):
    nstruc, ncoor = intern_coordinate.shape()
    try:
        assert rmsd_matrix.shape() == (nstruc, nstruc)

        intcoor_diff = np.abs(intern_coordinate[:,None, :]-intern_coordinate[None, :, :])

        intcoor_diff_inf = intcoor_diff[rmsd_matrix < cutoff]

        thresholds = np.min(intcoor_diff_inf)

        return thresholds
    except:
        print("The RMSD matrix has not the size of the number of structures.")


def main():
    rmsd = np.load(sys.argv[1]) #rmsdAAA_test.npy
    intcoor = np.load(sys.argv[2]) #AAA_test_trigo.npy
    cutoff = float(sys.argv[3]) # 1A RMSD cutoff

    nstruc, ncoor = intcoor.shape()
    assert rmsd.shape() == (nstruc, nstruc)

    intcoor_diff = np.abs(intcoor[:,None, :]-intcoor[None, :, :])

    intcoor_diff_inf = intcoor_diff[rmsd < cutoff]

    thresholds = np.min(intcoor_diff_inf)

    np.save(sys.argv[4], thresholds)

if __name__ == "__main__":
    main()
