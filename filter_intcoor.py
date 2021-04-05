#!/usr/bin/env python3

import numpy as np
import sys

def filter_intcoor(intcoor, thresholds, ndist):
    # The ndist first coor are delta of distances,
    # next coor are delta of angles

    nstruc, ncoor = intcoor.shape
    assert nstruc < 2**16
    assert len(thresholds) == ncoor + 2

    keep = np.ones((nstruc, nstruc), dtype=bool)
    np.fill_diagonal(keep, 0)
    keep = np.triu(keep)

    dist_sum = np.zeros((nstruc,nstruc))
    for m in range(ndist):
        print("coor %i"%(m+1))
        diff = np.abs(intcoor[:,None, m]-intcoor[None, :, m])
        dist_sum += diff
        d = diff < thresholds[m]
        diff = None # to free memory
        keep = keep & d
        d = None
    
    d = dist_sum < thresholds[-2]
    dist_sum = None
    keep = keep & d
    d=None  
    
    ang_sum = np.zeros((nstruc,nstruc))
    for m in range(ndist,ncoor):
        print("coor %i"%(m+1))
        x = np.abs(intcoor[:,None, m]-intcoor[None, :, m])
        diff = np.minimum(x, 360-x)
        x = None
        ang_sum += diff
        d = diff  < thresholds[m]
        diff = None
        keep = keep & d
        d = None

    d = ang_sum < thresholds[-1]
    keep = keep & d
    d=None
    ang_sum = None

    print( 2*(keep.sum()) / (len(intcoor)**2-len(intcoor)) )
    print("all coor done")

    keep_ind = np.argwhere(keep)
    keep = None
    # uint16_t = Unsigned integer (0 to 65535)
    # max number of conformers is 65535
    # if neeeded, change to uint32
    keep_list = np.array(keep_ind, dtype=np.uint16)
    keep_ind = None
    return keep_list

def main():
    intcoor = np.load(sys.argv[1])
    thresholds = [float(l.strip()) for l in open(sys.argv[2]).readlines()]
    ndist = int(sys.argv[3])

    keep_list = filter_intcoor(intcoor, thresholds,ndist)
    np.save(sys.argv[4], keep_list)

if __name__ == "__main__":
    main()