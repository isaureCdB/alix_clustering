#!/usr/bin/env python3

import numpy as np
import sys

def filter_intcoor(intcoor, thresholds, ndist):
    # x first coor are delta of distances
    # next coor are delta of angles

    nstruc, ncoor = intcoor.shape
    keep = np.ones((nstruc, nstruc), dtype=bool)

    for m in range(ndist):
        delta = np.abs(intcoor[:,None, m]-intcoor[None, :, m]) 
        d = delta < thresholds[m]
        keep = keep & d

    for m in range(ndist,ncoor):
        dd = np.abs(intcoor[:,None, m]-intcoor[None, :, m])
        delta = np.minimum(dd, 360-dd) 
        d = delta < thresholds[m]
        keep = keep & d
    
    return keep

def main():
    intcoor = np.load(sys.argv[1])
    thresholds = np.load(sys.argv[2])
    ndist = int(sys.argv[3])

    keep = filter_intcoor(intcoor, thresholds,ndist)
    np.save(sys.argv[4], keep)

if __name__ == "__main__":
    main()
