#!/usr/bin/env python3

import numpy as np

#Communication via le terminal de commande
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('rnpy',help="np array - pseudo-atom representation of XXX")
parser.add_argument('intcoornpy', help="np array - 9 intern coordinates representation of XXX")
parser.add_argument('thresholdsnpy', help="np array - dictionnary of intern coordinates thresholds")
parser.add_argument('filtered', help="np array - dictionnary of intern coordinates thresholds")

args = parser.parse_args()

structures = np.load(args.rnpy)
intcoor = np.load(args.intcoornpy)
thresholds = np.load(args.seuilsnpy,allow_pickle='TRUE').item()

ncoor = intcoor.shape[1]
keep = np.ones((ncoor, ncoor), dtype=bool)

for m in range(ncoor):
    delta = np.abs(intcoor[:,None, m]-intcoor[None, :, m]) 
    d = delta < thresholds[m]
    keep = keep | d

np.dump(keep, args.filtered)