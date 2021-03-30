#!/usr/bin/env python3

import numpy as np
import math, sys

def dist_at(coor, ind1, ind2):
    a = coor[:,ind1]
    b = coor[:,ind2]
    d = ((b-a)**2).sum(axis=1)
    dist = np.sqrt(d)
    return dist

def dihedral_angle(coord, a, b, c, d):
    #https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    b1 = coord[:,b]-coord[:,a]
    b2 = coord[:,c]-coord[:,b]
    b3 = coord[:,d]-coord[:,c]

    norm_b2 = np.sqrt( (b2**2).sum(axis=1) )
    unit_b2 = np.divide(b2, norm_b2)

    n1 = np.cross(b1, b2) # normal vector to plan containing a,b,c
    n2 = np.cross(b2, b3) # normal vector to plan containin b,c,d

    m1 = np.cross(n1, unit_b2)

    x = (n1 * n2).sum(axis=1) # dot product
    y = (m1 * n2).sum(axis=1) # dot product

    rdh = np.arctan2(y,x)
    ddh = rdh*180/np.pi

    return ddh

#coordinates of all conformers 
coor = np.load(sys.argv[1])

#indices of pseudo-atoms GP1_1, GX3/4_1, GX3/4_2, GX3/4_3, GS1_3
#AAA: 1 7 14 21 16
ind = [int(i)-1 for i in sys.argv[3:] ]

nstruc = coor.shape[0]

intcoor = np.zeros((nstruc, 4))

# dist base1-base2
intcoor[:,0] = dist_at(coor, ind[1], ind[2])
# dist base2-base3
intcoor[:,1] = dist_at(coor, ind[2], ind[3])
# dist base1-base3
intcoor[:,2] = dist_at(coor, ind[1], ind[3])
#dist ph1-sug3
intcoor[:,3] = dist_at(coor, ind[0], ind[4])

np.save(sys.argv[2], intcoor)