#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 23:03:08 2021

@author: alixdelannoy
"""

import numpy as np
import math
import argparse


############################
#Calcul de distance/d'angles
############################

def distance(x, y):
    """
    Distance between two points x and y of respective coordinates (x0,x1,x2) and (y0,y1,y2)
    """
    return math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)


def vectorproduct(v, w):
    """
    Vector product between two vectors v and w
    """
    return (v[1]*w[2]-v[2]*w[1], v[2]*w[0]-v[0]*w[2], v[0]*w[1]-v[1]*w[0])


def scalarproduct(v, w):
    """
    Scalar product between two vectors v and w
    """
    return (v[0]*w[0]+v[1]*w[1]+v[2]*w[2])


def norm(v):
    """
    Norm of a given vector v
    """
    return distance(v, (0, 0, 0))


def angle(a, b, c):
    """
    Angle formed by three points a, b and c
    """
    ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]]
    bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]
    abVec = scalarproduct(ab, ab)
    bcVec = scalarproduct(bc, bc)
    abNorm = [ab[0] / abVec, ab[1] / abVec, ab[2] / abVec]
    bcNorm = [bc[0] / bcVec, bc[1] / bcVec, bc[2] / bcVec]
    res = scalarproduct(abNorm, bcNorm)
    return math.acos(res)*180/math.pi


def dihedral_angle(a, b, c, d):
    """
    Returns the dihedral angle formed by four points a,b,c and d
    """
    ab = (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    bc = (b[0]-c[0], b[1]-c[1], b[2]-c[2])
    cd = (c[0]-d[0], c[1]-d[1], c[2]-d[2])
    normalVec1 = vectorproduct(ab, bc) #vecteur normal au plan formé par a,b,c
    normalVec2 = vectorproduct(bc, cd) #vecteur normal au plan formé par b,c,d

    if scalarproduct(bc, vectorproduct(normalVec1, normalVec2)) < 0:
        sc=scalarproduct(normalVec1, normalVec2)
        n=norm(normalVec1)*norm(normalVec2)
        r=sc/n
        #important pour éviter les problèmes d'arrondi: parfois sc et n sont très proches et arrondis on peut avoir un rapport qui n'est pas entre -1 et 1
        if r<-1:
            r=-1
        if r>1:
            r=1
        return math.acos(r)*180/math.pi

    else:
        sc=scalarproduct(normalVec1, normalVec2)
        n=norm(normalVec1)*norm(normalVec2)
        r=sc/n
        #important pour éviter les problèmes d'arrondi: parfois sc et n sont très proches et arrondis on peut avoir un rapport qui n'est pas entre -1 et 1
        if r<-1:
            r=-1
        if r>1:
            r=1
        return -math.acos(r)*180/math.pi


##############################
#Calcul de coordonées internes
##############################

def nucl(X):
    """
    Returns the number of pseudo-atoms in the pseudo-atom representation of a nucleotide X
    """
    if X=='A' or X=='G':
        return 7
    else:
        return 6

def basedistance12(conf, sequence):
    """
    Returns the distance between the "tip" of base 1 and 2
    """
    return distance(conf[nucl(sequence[0])-1], conf[nucl(sequence[0])+nucl(sequence[1])-1])

def basedistance23(conf, sequence):
    """
    Returns the distance between the "tip" of base 2 and 3
    """
    return distance(conf[nucl(sequence[0])+nucl(sequence[1])-1],conf[nucl(sequence[0])+nucl(sequence[1])+nucl(sequence[2])-1])

def basedistance13(conf, sequence):
    """
    Returns the distance between the "tip" of base 1 and 3
    """
    return distance(conf[nucl(sequence[0])-1],conf[nucl(sequence[0])+nucl(sequence[1])+nucl(sequence[2])-1])

def backbone(conf, sequence):
    """
    Returns the coordinates of the pseudo-atoms that are part of the backbone of a conformation (conf)
    """
    bbi = [0, 1, nucl(sequence[0]), nucl(sequence[0])+1, nucl(sequence[0])+nucl(sequence[1]), nucl(sequence[0])+nucl(sequence[1])+1]
    bb=[]
    for i in bbi:
        bb.append(conf[i])
    return bb

def backboneangles(conf, sequence):
    """
    Returns the three dihedral angles formed by pseudo-atoms along the backbone of a conformation (conf)
    """
    bb = backbone(conf, sequence)
    bba = []
    for i in range(3):
        bba.append(dihedral_angle(bb[i], bb[i+1], bb[i+2], bb[i+3]))
    return bba

def mu(conf, sequence):
    """
    Returns the torsion angles between bases of nucleotides 1 and 2, 2 and 3, and 1 and 3 of a pseudo-atom conformation (conf)
    """
    return [dihedral_angle(conf[3], conf[2], conf[nucl(sequence[0])+2], conf[nucl(sequence[0])+3]),
            dihedral_angle(conf[nucl(sequence[0])+3], conf[nucl(sequence[0])+2], conf[nucl(sequence[0])+nucl(sequence[1])+2], conf[nucl(sequence[0])+nucl(sequence[1])+3]),
            dihedral_angle(conf[3], conf[2],conf[nucl(sequence[0])+nucl(sequence[1])+2], conf[nucl(sequence[0])+nucl(sequence[1])+3])]

def return_internal_coordinates(atoms, sequence):
    d12 = basedistance12(atoms, sequence)
    d23 = basedistance23(atoms, sequence)
    d13 = basedistance13(atoms, sequence)
    bb1, bb2, bb3 = backboneangles(atoms, sequence)
    mu12, mu23, mu13 = mu(atoms, sequence)
    trigo = np.array([d12, d23, d13, bb1, bb2, bb3, mu12, mu23, mu13])

    return trigo

########################################
#Application sur le trinucléotide choisi
########################################

def main():

    #Communication via le terminal de commande
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('npy', help="np array")
    args = parser.parse_args()

    #Import
    structures = np.load(args.npy)
    sequence=args.npy.split("r.npy")[0]

    trigo = []
    for i in range(structures.shape[0]):
        trigo.append(return_internal_coordinates(structures[i], sequence))

    filename = args.npy.split("r.npy")[0]+"_trigo.npy"
    np.save(filename,np.array(trigo))

if __name__ == "__main__":
    main()
