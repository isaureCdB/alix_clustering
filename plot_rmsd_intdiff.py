#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import numpy as np

coor = np.load(sys.argv[1])
n, m = coor.shape
a = int(sys.argv[2])

npairs = (n**2-n)/2
diff = np.zeros(npairs)

[TO FINISH]

for i in range(n):
    for j in range(i+1,n):
        x = ((c[i,a]-c[j,a])**2)**0.5
        xx = 0.1*int(10*x)


x= [6, 6, 6, 6, 6, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 13, 13, 13, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 3, 3, 12, 12, 12, 16, 15, 15, 14, 11]
y= [2, 5, 6, 4, 1, 3, 10, 7, 9, 8, 8, 7, 1, 9, 6, 3, 7, 4, 1, 6, 3, 5, 2, 2, 5, 4, 1, 3, 6, 13, 11, 2, 4, 1, 3, 8, 6, 4, 5, 2, 7, 1, 3, 1, 11, 7, 8, 11, 15, 9, 10, 11]
z= [6, 6, 12, 6, 4, 4, 12, 2, 4, 2, 2, 4, 12, 6, 2, 2, 20, 20, 8, 6, 8, 16, 6, 10, 20, 16, 8, 4, 2, 2, 2, 18, 26, 16, 26, 6, 6, 2, 2, 2, 4, 2, 14, 4, 4, 2, 2, 2, 2, 2, 2, 2]
z = [pow(1.25,x) for x in z]

fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(x,y,s=z,c="black")
plt.plot([1,17],[1,17],'k-',lw=1)
plt.xlabel("Number of nucleotide in the loop.", color="black")
plt.ylabel("Number of nucleotide in contact with the protein in the loop.", color="black")

plt.savefig("hairpin_fig_contact.png")

