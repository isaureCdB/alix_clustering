#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 15:21:36 2021

@author: alixdelannoy
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import argparse

#Communication via le terminal de commande
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('npy', help="np array of pseudo-atoms conformations")
args = parser.parse_args()

#Import 
XXXr = np.load(args.npy)
XXX=args.npy.split("r.npy")[0]


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
    if X=='C' or X=='U':
        return 6
    
def basedistance12(conf):
    """
    Returns the distance between the "tip" of base 1 and 2
    """
    return distance(conf[nucl(XXX[0])-1], conf[nucl(XXX[0])+nucl(XXX[1])-1])

def basedistance23(conf):
    """
    Returns the distance between the "tip" of base 2 and 3
    """
    return distance(conf[nucl(XXX[0])+nucl(XXX[1])-1],conf[nucl(XXX[0])+nucl(XXX[1])+nucl(XXX[2])-1])

def basedistance13(conf):
    """
    Returns the distance between the "tip" of base 1 and 3
    """
    return distance(conf[nucl(XXX[0])-1],conf[nucl(XXX[0])+nucl(XXX[1])+nucl(XXX[2])-1])


def backbone(conf):
    """
    Returns the coordinates of the pseudo-atoms that are part of the backbone of a conformation (conf)
    """
    bbi = [0, 1, nucl(XXX[0]), nucl(XXX[0])+1, nucl(XXX[0])+nucl(XXX[1]), nucl(XXX[0])+nucl(XXX[1])+1]
    bb=[]
    for i in bbi:
        bb.append(conf[i])
    return bb

def backboneangles(conf):
    """
    Returns the three dihedral angles formed by pseudo-atoms along the backbone of a conformation (conf)
    """
    bb = backbone(conf)
    bba = []
    for i in range(3):
        bba.append(dihedral_angle(bb[i], bb[i+1], bb[i+2], bb[i+3]))
    return bba


def chi_1(conf):
    """
    Returns the coordinates of the pseudo-atoms that form the angle between sugar and base of each nucleotide of a conformation (conf)
    Method 1 : GS1-GS2-GX1-GX2
    """
    chii = [1, 1+nucl(XXX[0]), 1+nucl(XXX[0])+nucl(XXX[1])]
    chi = []
    for i in chii:
        chi.append(conf[i])
        chi.append(conf[i+1])
        chi.append(conf[i+2])
        chi.append(conf[i+3])
    return(chi)

def chi_2(conf):
    """
    Returns the coordinates of the pseudo-atoms that form the angle between sugar and base of each nucleotide of a conformation (conf)
    Method 2 : GS1-GS2-GX1-GX4
    """
    chii=[1,2,3,nucl(XXX[0])-1,
          nucl(XXX[0])+1,nucl(XXX[0])+2,nucl(XXX[0])+3,nucl(XXX[0])+nucl(XXX[1])-1,
          nucl(XXX[0])+nucl(XXX[1])+1,nucl(XXX[0])+nucl(XXX[1])+2,nucl(XXX[0])+nucl(XXX[1])+3,nucl(XXX[0])+nucl(XXX[1])+nucl(XXX[2])-1]
    chi=[]
    for i in chii:
        chi.append(conf[i])
    return(chi)

def chi_3(conf):
    """
    Returns the coordinates of the pseudo-atoms that form the angle between sugar and base of each nucleotide of a conformation (conf)
    Method 3 : GS1-GS2-GX1-GX3
    """
    chii=[1,2,3,nucl(XXX[0])-2,
          nucl(XXX[0])+1,nucl(XXX[0])+2,nucl(XXX[0])+3,nucl(XXX[0])+nucl(XXX[1])-2,
          nucl(XXX[0])+nucl(XXX[1])+1,nucl(XXX[0])+nucl(XXX[1])+2,nucl(XXX[0])+nucl(XXX[1])+3,nucl(XXX[0])+nucl(XXX[1])+nucl(XXX[2])-2]
    chi=[]
    for i in chii:
        chi.append(conf[i])
    return(chi)

def chiangles(conf):
    """
    Returns the angle between sugar and base for each nucleotide of a pseudo-atom conformation (conf)
    """
    c = chi_1(conf)
    chiangles = []
    for i in [0, 4, 8]:
        chiangles.append(dihedral_angle(c[i], c[i+1], c[i+2], c[i+3]))
    return chiangles

def mu(conf):
    """
    Returns the torsion angles between bases of nucleotides 1 and 2, 2 and 3, and 1 and 3 of a pseudo-atom conformation (conf)
    """
    return [dihedral_angle(conf[3], conf[2], conf[nucl(XXX[0])+2], conf[nucl(XXX[0])+3]),
            dihedral_angle(conf[nucl(XXX[0])+3], conf[nucl(XXX[0])+2], conf[nucl(XXX[0])+nucl(XXX[1])+2], conf[nucl(XXX[0])+nucl(XXX[1])+3]),
            dihedral_angle(conf[3], conf[2],conf[nucl(XXX[0])+nucl(XXX[1])+2], conf[nucl(XXX[0])+nucl(XXX[1])+3])]




########################################
#Application sur le trinucléotide choisi
########################################
d12=[]
d23=[]
d13=[]
bb1=[]
bb2 = []
bb3 = []
chi1=[]
chi2 = []
chi3 = []
mu12=[]
mu23=[]
mu13=[]


for i in range(XXXr.shape[0]):
    bb1.append(backboneangles(XXXr[i])[0])
    bb2.append(backboneangles(XXXr[i])[1])
    bb3.append(backboneangles(XXXr[i])[2])
    chi1.append(chiangles(XXXr[i])[0])
    chi2.append(chiangles(XXXr[i])[1])
    chi3.append(chiangles(XXXr[i])[2])
    d12.append(basedistance12(XXXr[i]))
    d23.append(basedistance23(XXXr[i]))
    d13.append(basedistance23(XXXr[i]))
    mu12.append(mu(XXXr[i])[0])
    mu23.append(mu(XXXr[i])[1])
    mu13.append(mu(XXXr[i])[2])
    

#Tracé et enregistrement des distributions pour ce nucléotide
plt.figure(figsize=(150,200))

plt.subplot(431)
plt.hist(d12, color = 'blue',edgecolor='black',bins=100)
plt.xlabel('Distance interbase',fontsize=70)
plt.ylabel('Nombre de trinucléotides',fontsize=70)
plt.title('Distribution - Distance nucléotide 1/nucléotide 2',fontsize=90)

plt.subplot(432)
plt.hist(d23, color = 'blue',edgecolor='black',bins=100)
plt.xlabel('Distance interbase',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Distance nucléotide 2/nucléotide 3',fontsize=90)

plt.subplot(433)
plt.hist(d13, color = 'blue',edgecolor='black',bins=100)
plt.xlabel('Distance interbase',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Distance nucléotide 1/nucléotide 3',fontsize=90)

plt.subplot(434)
plt.hist(bb1, color = 'orangered',edgecolor='black',bins=100)
plt.xlabel('Backbone angle 1',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Backbone angle 1',fontsize=90)

plt.subplot(435)
plt.hist(bb2, color = 'orangered',edgecolor='black',bins=100)
plt.xlabel('Backbone angle 2',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Backbone angle 2',fontsize=90)

plt.subplot(436)
plt.hist(bb3, color = 'orangered',edgecolor='black',bins=100)
plt.xlabel('Backbone angle 3',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Backbone angle 3',fontsize=90)

plt.subplot(437)
plt.hist(chi1, color = 'green',edgecolor='black',bins=100)
plt.xlabel('Chi angle of first nucleotide',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Chi angle of first nucleotide',fontsize=90)

plt.subplot(438)
plt.hist(chi2, color = 'green',edgecolor='black',bins=100)
plt.xlabel('Chi angle of second nucleotide',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Chi angle of second  nucleotide',fontsize=90)

plt.subplot(439)
plt.hist(chi3, color = 'green',edgecolor='black',bins=100)
plt.xlabel('Chi angle of third nucleotide',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Chi angle of third nucleotide',fontsize=90)

plt.subplot(4,3,10)
plt.hist(mu12, color = 'yellow',edgecolor='black',bins=100)
plt.xlabel('Mu between first and second nucleotide',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Mu between first and second nucleotide',fontsize=90)

plt.subplot(4,3,11)
plt.hist(mu23, color = 'yellow',edgecolor='black',bins=100)
plt.xlabel('Mu between second and third nucleotide',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Mu between second and third nucleotide',fontsize=90)

plt.subplot(4,3,12)
plt.hist(mu13, color = 'yellow',edgecolor='black',bins=100)
plt.xlabel('Mu between first and third nucleotide',fontsize=70)
plt.ylabel('Nombre de trinucléotide',fontsize=70)
plt.title('Distribution - Mu between first and third nucleotide',fontsize=90)

plt.savefig(XXX+"_distriCI",edgecolor="black",transparent=True)
print("Plots have been succesfully saved.")


