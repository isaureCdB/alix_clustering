#!/usr/bin/env python3

import numpy as np
import sys

def filter_intcoor(intcoor, thresholds, ndist):
    # The ndist first coor are delta of distances,
    # next coor are delta of angles

    nstruc, ncoor = intcoor.shape

    keep = np.ones((nstruc, nstruc), dtype=bool)
    np.fill_diagonal(keep, 0)

    for m in range(ndist):
        print("coor %i"%(m+1))
        delta = np.abs(intcoor[:,None, m]-intcoor[None, :, m]) < thresholds[m]
        d = delta < thresholds[m]
        keep = keep & d
        d = []

    for m in range(ndist,ncoor):
        print("coor %i"%(m+1))
        dd = np.abs(intcoor[:,None, m]-intcoor[None, :, m])
        d = np.minimum(dd, 360-dd)  < thresholds[m]
        dd = []
        keep = keep & d
        d = []

    print("all coor done")
    # max number of expected pairs with rmsd potentially bellow cutoff
    # maxpair = 0.001 * ncoor * ncoor
    
    # uint16_t = Unsigned integer (0 to 65535)
    # max number of conformers is 65535
    # if neeeded, change to uint32_t
    keep_ind = np.argwhere(keep)
    keep = keep_ind[keep_ind[:,0]<keep_ind[:,1]]
    keep_list = np.array(keep, dtype=np.uint16)
    keep_ind = []
    return keep_list

def main():
    intcoor = np.load(sys.argv[1])
    thresholds = [float(l.strip()) for l in open(sys.argv[2]).readlines()]
    ndist = int(sys.argv[3])

    keep_list = filter_intcoor(intcoor, thresholds,ndist)
    np.save(sys.argv[4], keep_list)

if __name__ == "__main__":
    main()
