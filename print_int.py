#!/usr/bin/env python3

import numpy as np
import sys

c = np.load(sys.argv[1])
a = int(sys.argv[2])-1

n, m = c.shape
print(n, m)

'''
for a in range(ndist):
    print("coor %i"%(a+1), file=sys.stderr )
    o = open("int%i.txt"%(a+1), "w")
    for i in range(n):
        for j in range(i+1,n):
            x = ((c[i,a]-c[j,a])**2)**0.5
            print("%.3f"%x, file=o)
    o.close()
'''
#for a in range(ndist, m):
print("coor %i"%(a+1), file=sys.stderr )
o = open("int%i.txt"%(a+1), "w")
for i in range(n):
    for j in range(i+1,n):
        x = ((c[i,a]-c[j,a])**2)**0.5
        xx = min(x, 360-x)
        print("%.3f"%xx, file=o)
o.close()
