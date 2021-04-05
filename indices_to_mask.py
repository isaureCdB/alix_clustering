#!/usr/bin/env python3

import sys
import numpy as np

a = np.load(sys.argv[1])
a = np.array(a)                                                                                      

n = int(sys.argv[2])
#m = int(sys.argv[3])

b = np.zeros((n,n), dtype=bool)

c0=a[:,0]                                                                                          
c1=a[:,1]                                                                                          
b[c0,c1] = 1                                                                                       
b[c1,c0] = 1     
np.save(sys.argv[3], b)