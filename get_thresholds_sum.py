#!/usr/bin/env python3

import sys
import numpy as np

def get_threshold_sum(mask_rmsd, intern_coordinate, ndist, cutoff, percent):
    nstruc, ncoor = intern_coordinate.shape

    assert mask_rmsd.shape == (nstruc, nstruc), "The RMSD matrix has not the size of the number of structures"

    intcoor_diff = np.zeros((nstruc, nstruc, 2))

    for m in range(ndist):
        delta = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m])
        intcoor_diff[:,:,0] = np.add(intcoor_diff[:,:,0], delta)
    
    for m in range(ndist,ncoor):
        delta = np.abs(intern_coordinate[:,None, m]-intern_coordinate[None, :, m])
        intcoor_diff[:,:,1] = np.add(intcoor_diff[:,:,1], np.minimum(delta, 360-delta))

    intcoor_diff_inf = intcoor_diff[mask_rmsd]
    
    if percent == 1:
        thresh = np.max(intcoor_diff_inf, axis=0)   
    else:
        thresh = []
        p = int(percent*len(intcoor_diff_inf)) - 1
        for c in range(2):
            order_diff = np.sort(intcoor_diff_inf[:,c])
            thresh.append(order_diff[p])
        
    return thresh


def main():
    mask_rmsd = np.load(sys.argv[1]) #mask_rmsd.npy
    intern_coordinate = np.load(sys.argv[2]) #AAA_test_trigo.npy
    ndist = int(sys.argv[3])
    cutoff = float(sys.argv[4]) # 1A RMSD cutoff
    percent = float(sys.argv[5]) # 0.999
    
    thresholds = get_threshold_sum(mask_rmsd, intern_coordinate, ndist, cutoff, percent)
    for t in thresholds:
        print("%.2f"%t)

if __name__ == "__main__":
    main()
