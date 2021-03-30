#!/usr/bin/env python3

import numpy as np
import sys

intcoor = np.load(sys.argv[1])
thresholds = np.load(sys.argv[2])

nstruc, ncoor = intcoor.shape
keep = np.ones((nstruc, nstruc), dtype=bool)

for m in range(ncoor):
    delta = np.abs(intcoor[:,None, m]-intcoor[None, :, m]) 
    d = delta < thresholds[m]
    keep = keep | d

np.save(sys.argv[3], keep)