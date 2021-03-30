#!/usr/bin/env python3

import numpy as np
import sys

# 3 first coor are delta of  distances
# 6 next coor are delta of angles
intcoor = np.load(sys.argv[1])
thresholds = np.load(sys.argv[2])

nstruc, ncoor = intcoor.shape
keep = np.ones((nstruc, nstruc), dtype=bool)

for m in range(3):
    delta = np.abs(intcoor[:,None, m]-intcoor[None, :, m]) 
    d = delta < thresholds[m]
    keep = keep | d

for m in range(4,ncoor):
    dd = np.abs(intcoor[:,None, m]-intcoor[None, :, m])
    delta = np.mod(dd,360) 
    d = delta < thresholds[m]
    keep = keep | d


np.save(sys.argv[3], keep)