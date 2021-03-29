#!/usr/bin/env python3

import numpy as np

#Communication via le terminal de commande
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('rnpy',help="np array - pseudo-atom representation of XXX")
parser.add_argument('intcoornpy', help="np array - 9 intern coordinates representation of XXX")
parser.add_argument('seuilsnpy', help="np array - dictionnary of intern coordinates thresholds")
args = parser.parse_args()

structures = np.load(args.rnpy)
intcoor = np.load(args.intcoornpy)
seuils = np.load(args.seuilsnpy,allow_pickle='TRUE').item()

#Calcul des matrices de similarités
similarity_d12 = abs(intcoor[]-intcoor) 
similarity_d23[j,i]=abs(intcoor[i,1]-intcoor[j,1]) 
similarity_d13[j,i]=abs(intcoor[i,2]-intcoor[j,2])
#pour les angles on mesure la similarité par une différence des valeurs absolues
similarity_bb1[j,i]=abs(abs(intcoor[i,3])-abs(intcoor[j,3]))
similarity_bb2[j,i]=abs(abs(intcoor[i,4])-abs(intcoor[j,4]))
similarity_bb3[j,i]=abs(abs(intcoor[i,5])-abs(intcoor[j,5]))
similarity_mu12[j,i]=abs(abs(intcoor[i,6])-abs(intcoor[j,6]))
similarity_mu23[j,i]=abs(abs(intcoor[i,7])-abs(intcoor[j,7]))
similarity_mu13[j,i]=abs(abs(intcoor[i,8])-abs(intcoor[j,8]))
        if 0<similarity_d12[j,i]<=seuils['d12'] and 0<similarity_d23[j,i]<=seuils['d23'] and 0<similarity_d23[j,i]<=seuils['d13'] and 0<similarity_bb1[j,i]<=seuils['bb1'] and 0<similarity_bb2[j,i]<=seuils['bb2'] and 0<similarity_bb3[j,i]<=seuils['bb3'] and 0<similarity_mu12[j,i]<=seuils['mu12'] and 0<similarity_mu23[j,i]<=seuils['mu23'] and 0<similarity_mu13[j,i]<=seuils['mu13']:
            id_all.append((j,i))

