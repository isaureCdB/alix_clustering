#!/usr/bin/env python3

import numpy as np
import sys

coor = np.load(sys.argv[1])
keep_list = np.load(sys.argv[2])
cutoff = float(sys.argv[3])

def multifit(array_atoms1, atoms2):
  """
  Fits an array of atom sets (array_atoms1) onto an atom set (atoms2)
  """
  import numpy
  assert isinstance(array_atoms1, numpy.ndarray)
  assert isinstance(atoms2, numpy.ndarray)

  assert len(array_atoms1.shape) == 3
  assert len(atoms2.shape) == 2

  assert len(atoms2) > 0
  assert array_atoms1.shape[2] == atoms2.shape[1] == 3
  assert len(atoms2) == array_atoms1.shape[1], (atoms2.shape, array_atoms1.shape)
  L = len(atoms2)

  # must alway center the two proteins to avoid
  # affine transformations.  Center the two proteins
  # to their selections.
  COM1 = numpy.sum(array_atoms1,axis=1) / float(L)
  COM2 = numpy.sum(atoms2,axis=0) / float(L)

  array_atoms1 = array_atoms1 - COM1[:,numpy.newaxis, :]
  atoms2 = atoms2 - COM2

  # Initial residual, see Kabsch.
  E0 = numpy.sum( numpy.sum(array_atoms1 * array_atoms1,axis=1),axis=1) + numpy.sum( numpy.sum(atoms2 * atoms2,axis=0),axis=0)

  #
  # This beautiful step provides the answer.  V and Wt are the orthonormal
  # bases that when multiplied by each other give us the rotation matrix, U.
  # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
  d = numpy.einsum("ijk,jl->ikl", array_atoms1, atoms2)
  V, S, Wt = numpy.linalg.svd( d )

  # We already have our solution, in the results from SVD.
  # we just need to check for reflections and then produce
  # the rotation.  V and Wt are orthonormal, so their det's
  # are +/-1.0 (and thus products are +/- 1.0 ).
  reflect = numpy.linalg.det(V) * numpy.linalg.det(Wt)

  S[:,-1] *= reflect
  V[:,:,-1] *= reflect[:, numpy.newaxis]
  U = numpy.einsum('...ij,...jk->...ki', V, Wt)
  RMSD = E0 - (2.0 * S.sum(axis=1))
  RMSD = numpy.sqrt(abs(RMSD / L))
  return U, COM1-COM2, RMSD

def fit_multi_npy(a, ref):
    rotation, translation, RMSD = multifit(a, ref)
    rot = np.transpose(rotation, axes=(0,2,1))
    COM = a.sum(axis=1)/a.shape[1]
    centered = a - COM[:,None,:]
    rotated = np.einsum('...ij,...jk->...ik',centered,rot)
    fitted = rotated + COM[:,None,:]
    translated = fitted - translation[:,None,:]
    return translated, RMSD

n = coor.shape[0]
new_keep = []

for i in range(n):
    print(i)
    ref = coor[i]
    keep = keep_list[ keep_list[:,0] == i][:,1]
    tofit = coor[keep]
    fitted, rmsd = fit_multi_npy(tofit, ref)
    new = [[i, j] for j in keep[rmsd < cutoff]]
    new_keep.extend(new)

new_keep_list = np.array(new_keep)
np.save(sys.argv[4], new_keep_list)