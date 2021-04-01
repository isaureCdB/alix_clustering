#!/usr/bin/env python3

import sys
import numpy as np

def get_threshold(mask_rmsd, intern_coordinate, ndist, cutoff, percent):
    nstruc, ncoor = intern_coordinate.shape

    assert mask_rmsd.shape == (nstruc, nstruc), "The RMSD matrix has not the size of the number of structures"

    intcoor_diff = np.zeros((nstruc, nstruc, ncoor))

    for m in range(ndist):
        intcoor_diff[:,:,m] = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m]) 
 
    for m in range(ndist,ncoor):
        dd = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m])
        intcoor_diff[:,:,m] = np.minimum(dd, 360-dd) 
    
    n = mask_rmsd.sum()
    diff_inf = np.zeros((n, ncoor+2))
    diff_inf[:,:ncoor] = intcoor_diff[mask_rmsd]

    diff_inf[:,ncoor] = diff_inf[:,:ndist].sum(axis=1)
    diff_inf[:,ncoor+1] = diff_inf[:,ndist:].sum(axis=1)

    print(ncoor, file=sys.stderr)

    if percent == 1:
        thresh = np.max(diff_inf, axis=0)
    else:
        thresh = []
        p = int(percent*len(diff_inf)) - 1
        for c in range(ncoor+2):
            order_diff = np.sort(diff_inf[:,c])
            thresh.append(order_diff[p])
    
    return thresh


def main():
    mask_rmsd = np.load(sys.argv[1]) #mask_rmsd.npy
    intern_coordinate = np.load(sys.argv[2]) #AAA_test_trigo.npy
    ndist = int(sys.argv[3])
    cutoff = float(sys.argv[4]) # 1A RMSD cutoff
    percent = float(sys.argv[5]) # 0.999
    
    thresholds = get_threshold(mask_rmsd, intern_coordinate, ndist, cutoff, percent)
    for t in thresholds:
        print("%.2f"%t)

if __name__ == "__main__":
    main()
