#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 17:53:56 2021

@author: alixdelannoy
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse


#Change the path according to where you store your files if needed
rmsd=np.load('rmsdAAA.npy')
trigo=np.load('AAA_trigo.npy')

#Calcul des matrices de similarités
similarity_d12=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_d23=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_d13=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_bb1=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_bb2=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_bb3=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_mu12=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_mu23=np.zeros((trigo.shape[0],trigo.shape[0]))
similarity_mu13=np.zeros((trigo.shape[0],trigo.shape[0]))

for i in range(trigo.shape[0]):
    for j in range(i):
        #pour les distances on mesure la similarité par une différence
        similarity_d12[j,i]=abs(trigo[i,0]-trigo[j,0]) 
        similarity_d23[j,i]=abs(trigo[i,1]-trigo[j,1]) 
        similarity_d13[j,i]=abs(trigo[i,2]-trigo[j,2])
        #pour les angles on mesure la similarité par une différence des valeurs absolues
        similarity_bb1[j,i]=abs(abs(trigo[i,3])-abs(trigo[j,3]))
        similarity_bb2[j,i]=abs(abs(trigo[i,4])-abs(trigo[j,4]))
        similarity_bb3[j,i]=abs(abs(trigo[i,5])-abs(trigo[j,5]))
        similarity_mu12[j,i]=abs(abs(trigo[i,6])-abs(trigo[j,6]))
        similarity_mu23[j,i]=abs(abs(trigo[i,7])-abs(trigo[j,7]))
        similarity_mu13[j,i]=abs(abs(trigo[i,8])-abs(trigo[j,8]))


#Figures : RMSD en fonction de la similarité de coordonées internes
# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité de distance nucléotide 1/nucléotide 2',fontsize=70);
# plt.ylabel('RMSD',fontsize=70);
# plt.title('RMSD en fonction de la similarité de distance interbase',fontsize=80);
# plt.scatter(similarity_d12,rmsd);
# #plt.savefig('RMSD_Diff_d12',edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité de distance nucléotide 2/nucléotide 3',fontsize=70);
# plt.ylabel('RMSD',fontsize=70);
# plt.title('RMSD en fonction de la similarité de distance interbase',fontsize=80);
# plt.scatter(similarity_d23,rmsd);
# #plt.savefig('RMSD_Diff_d23',edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité de distance nucléotide 1/nucléotide 3',fontsize=70);
# plt.ylabel('RMSD',fontsize=70);
# plt.title('RMSD en fonction de la similarité de distance interbase',fontsize=80);
# plt.scatter(similarity_d13,rmsd);
# #plt.savefig('RMSD_Diff_d13',edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité entre les angles 1 du backbone',fontsize=70)
# plt.ylabel('RMSD',fontsize=70)
# plt.title("RMSD en fonction de la similarité BB1",fontsize=80)
# plt.scatter(similarity_bb1,rmsd);
# #plt.savefig("RMSD_AbsDiff_bb1",edgecolor="black",transparent=True)
# plt.show();
    
# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité entre les angles 2 du backbone',fontsize=70)
# plt.ylabel('RMSD',fontsize=70)
# plt.title("RMSD en fonction de la similarité BB2",fontsize=80)
# plt.scatter(similarity_bb2,rmsd);
# #plt.savefig("RMSD_AbsDiff_bb2",edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité entre les angles 3 du backbone',fontsize=70)
# plt.ylabel('RMSD',fontsize=70)
# plt.title("RMSD en fonction de la similarité BB3",fontsize=80)
# plt.scatter(similarity_bb3,rmsd);
# #plt.savefig("RMSD_AbsDiff_bb3",edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité entre les angles Mu12',fontsize=70)
# plt.ylabel('RMSD',fontsize=70)
# plt.title("RMSD en fonction de la similarité Mu12",fontsize=80)
# plt.scatter(similarity_mu12,rmsd);
# #plt.savefig("RMSD_AbsDiff_mu1",edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité entre les angles Mu23',fontsize=70)
# plt.ylabel('RMSD',fontsize=70)
# plt.title("RMSD en fonction de la similarité Mu23",fontsize=80)
# plt.scatter(similarity_mu23,rmsd);
# #plt.savefig("RMSD_AbsDiff_mu2",edgecolor="black",transparent=True)
# plt.show();

# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité entre les angles Mu13',fontsize=70)
# plt.ylabel('RMSD',fontsize=70)
# plt.title("RMSD en fonction de la similarité Mu13",fontsize=80)
# plt.scatter(similarity_mu13,rmsd);
# #plt.savefig("RMSD_AbsDiff_mu3",edgecolor="black",transparent=True)
# plt.show();

#Calculs des seuils "maximum"
idx=[]
for i in range(rmsd.shape[0]):
    for j in range(rmsd.shape[0]):
        #on prend tous les rmsd =! 0 (on ne veut pas le "bas" de la matrice symétrique)
        if 0<rmsd[i,j]<=1:
            idx.append((i,j))

dist={}
dist['d12']=[]
dist['d23']=[]
dist['d13']=[]
dist['bb1']=[]
dist['bb2']=[]
dist['bb3']=[]
dist['mu12']=[]
dist['mu23']=[]
dist['mu13']=[]
for (i,j) in idx:
    dist['d12'].append(similarity_d12[i,j])
    dist['d23'].append(similarity_d23[i,j])
    dist['d13'].append(similarity_d13[i,j])
    dist['bb1'].append(similarity_bb1[i,j])
    dist['bb2'].append(similarity_bb2[i,j])
    dist['bb3'].append(similarity_bb3[i,j])
    dist['mu12'].append(similarity_mu12[i,j])
    dist['mu23'].append(similarity_mu23[i,j])
    dist['mu13'].append(similarity_mu13[i,j])


seuils={}
seuils['d12']=max(dist['d12'])
seuils['d23']=max(dist['d23'])
seuils['d13']=max(dist['d13'])
seuils['bb1']=max(dist['bb1'])
seuils['bb2']=max(dist['bb2'])
seuils['bb3']=max(dist['bb3'])
seuils['mu12']=max(dist['mu12'])
seuils['mu23']=max(dist['mu23'])
seuils['mu13']=max(dist['mu13'])

np.save('thresholds_CI.npy', seuils) 

#Figures avec seuil 
# plt.figure(figsize=(40,40));
# plt.xlabel('Similarité de distance nucléotide 1/nucléotide 3',fontsize=70);
# plt.ylabel('RMSD',fontsize=70);
# plt.title('RMSD en fonction de la similarité de distance interbase',fontsize=80);
# plt.scatter(similarity_d13,rmsd);
# plt.axhline(y=1,color='red')
# plt.axvline(x=seuils['d13'],color='green')
# plt.show();


#Proportion de points avant le seuil 
# id_d12=[]
# id_d23=[]
# id_d13=[]
# id_bb1=[]
# id_bb2=[]
# id_bb3=[]
# id_mu12=[]
# id_mu23=[]
# id_mu13=[]

# for i in range(rmsd.shape[0]):
#     for j in range(rmsd.shape[0]):
#         if 0<similarity_d12[i,j]<=seuils['d12']:
#             id_d12.append((i,j))
#         if 0<similarity_d23[i,j]<=seuils['d23']:
#             id_d23.append((i,j))
#         if 0<similarity_d13[i,j]<=seuils['d13']:
#             id_d13.append((i,j))
#         if 0<similarity_bb1[i,j]<=seuils['bb1']:
#             id_bb1.append((i,j))
#         if 0<similarity_bb2[i,j]<=seuils['bb2']:
#             id_bb2.append((i,j))
#         if 0<similarity_bb3[i,j]<=seuils['bb3']:
#             id_bb3.append((i,j))
#         if 0<similarity_mu12[i,j]<=seuils['mu12']:
#             id_mu12.append((i,j))
#         if 0<similarity_mu23[i,j]<=seuils['mu23']:
#             id_mu23.append((i,j))
#         if 0<similarity_mu13[i,j]<=seuils['mu13']:
#             id_mu13.append((i,j))
 
# count=0
# for i in range(rmsd.shape[0]):
#     for j in range(rmsd.shape[0]):
#         if rmsd[i,j]!=0:
#             count+=1
    
# proportion={}
# proportion['d12']=len(id_d12)/(count)
# proportion['d23']=len(id_d23)/(count)
# proportion['d13']=len(id_d13)/(count)
# proportion['bb1']=len(id_bb1)/(count)
# proportion['bb2']=len(id_bb2)/(count)
# proportion['bb3']=len(id_bb3)/(count)
# proportion['mu12']=len(id_mu12)/(count)
# proportion['mu23']=len(id_mu23)/(count)
# proportion['mu13']=len(id_mu13)/(count)


