#!/usr/bin/env python3

import numpy as np
import math, sys

def distance(coor, ind1, ind2):
    a = coor[:,ind1]
    b = coor[:,ind2]
    d = ((b-a)**2).sum(axis=1)
    dist = np.sqrt(d)
    return dist

def dihedral(coord, a, b, c, d):
    #https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    b1 = coord[:,b]-coord[:,a]
    b2 = coord[:,c]-coord[:,b]
    b3 = coord[:,d]-coord[:,c]

    norm_b2 = np.sqrt( (b2**2).sum(axis=1) )
    unit_b2 = b2 / norm_b2[:, None]

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

#Nb of beads in each nucleotide
#AAA: 7 7 7
nbeads = [int(i) for i in sys.argv[3:] ]

nstruc = coor.shape[0]

intcoor = np.zeros((nstruc, 13))

gp1 = 0
gp2 = gp1 + nbeads[0]
gp3 = gp2 + nbeads[1]
gs11 = 1
gs12 = gs11 + nbeads[0]
gs13 = gs12 + nbeads[1]
gs21 = 2
gs22 = gs21 + nbeads[0]
gs23 = gs22 + nbeads[1]
gx11 = 3
gx12 = gx11 + nbeads[0]
gx13 = gx12 + nbeads[1]
gx21 = 4
gx22 = gx21 + nbeads[0]
gx23 = gx22 + nbeads[1]
gx1 = nbeads[0] -1
gx2 = gx1 + nbeads[0]
gx3 = gx2 + nbeads[1]

# dist base-base
intcoor[:,0] = distance(coor, gx1, gx2)
intcoor[:,1] = distance(coor, gx2, gx3)
intcoor[:,2] = distance(coor, gx1, gx3)
#dist ph1-sug3
intcoor[:,3] = distance(coor, gx1, gs13)
#bb angles
intcoor[:,4] = dihedral(coor, gp1, gs11, gp2, gs12)
intcoor[:,5] = dihedral(coor, gs11, gp2, gs12, gp3)
intcoor[:,6] = dihedral(coor, gp2, gs12, gp3, gs13)
#mu angles (base-base)
intcoor[:,7] = dihedral(coor, gx11, gs21, gs22, gx12)
intcoor[:,8] = dihedral(coor, gx12, gs22, gs23, gx13)
intcoor[:,9] = dihedral(coor, gx11, gs21, gs23, gx13)
#chi angles
intcoor[:,10] = dihedral(coor, gs11, gs21, gx11, gx21)
intcoor[:,11] = dihedral(coor, gs12, gs22, gx12, gx22)
intcoor[:,12] = dihedral(coor, gs13, gs23, gx13, gx23)

np.save(sys.argv[2], intcoor)