#!/usr/bin/env python3

import numpy as np
import sys, random

coor = np.load(sys.argv[1])
k = int(sys.argv[2])

ncoor = coor.shape[0]
x = random.sample(range(ncoor), k)
sample = coor[x]

np.save(sys.argv[3], sample)