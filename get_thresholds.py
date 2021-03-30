#!/usr/bin/env python3

import numpy as np

rmsd = np.load(sys.argv[1]) #rmsdAAA_test.npy
intcoor = np.load(sys.argv[2]) #AAA_test_trigo.npy
cutoff = float(sys.argv[3]) # 1A RMSD cutoff

nstruc, ncoor = intcoor.shape()
assert rmsd.shape() == (nstruc, nstruc)

intcoor_diff = np.abs(intcoor[:,None, :]-intcoor[None, :, :]) 

intcoor_diff_inf = intcoor_diff[rmsd < cutoff]

thresholds = np.min(intcoor_diff_inf1, axis=0)

np.save(sys.argv[4], thresholds)