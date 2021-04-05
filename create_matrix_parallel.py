#!/usr/bin/env python3

import sys

import numpy as np
from math import *
import rmsdlib
from multiprocessing import Pool, Queue

def read_pdb(pdb):
  atoms = []
  lines0 = open(pdb).readlines()
  lines = []
  extralines = []
  for l in lines0:
    if not l.startswith("ATOM"):
      extralines.append((len(lines), l))
      continue
    x = float(l[30:38])
    y = float(l[38:46])
    z = float(l[46:54])
    atoms.append((x,y,z))
    lines.append(l)
  return lines, np.array(atoms), extralines

def run(runarg):
  ref = runarg
  print(ref)
  lines_ref, atoms_ref, extralines_ref = read_pdb(mobiles[ref])
  results = np.zeros(len(mobiles))
  for ii in range(len(mobiles)):
      lines, atoms, extralines = read_pdb(mobiles[ii])
      rotmap, offset, rmsd = rmsdlib.fit(atoms_ref, atoms)
      results[ii] = round(rmsd, 3)

  return ref, results

def run2(runarg):
  ref = runarg
  print(ref)
  atoms_ref = test[ref]
  results = np.zeros(len(mobiles))
  for ii in range(len(mobiles)):
      if ii == ref:
          rmsd = 0
      else:
          atoms = test[ii]
          rotmap, offset, rmsd = rmsdlib.fit(atoms_ref, atoms)
      results[ii] = round(rmsd, 3)

  return ref, results

import sys
import argparse
a = argparse.ArgumentParser(prog="fit-multi.py")
a.add_argument("mobilelist")
a.add_argument("output")
a.add_argument("--np",type=int)
a.add_argument("--matrix", default=False)
args = a.parse_args()

mobiles = [l.strip().strip("\n") for l in open(args.mobilelist) if len(l.strip().strip("\n"))]

if args.matrix == False:
    runargs = []
    for i in range(len(mobiles)):
        runargs.append(i)

    pool = Pool(args.np)

    result = pool.map(run, runargs)
    # print(result)
    matrix_out = np.zeros((len(mobiles), len(mobiles)))

    for jj in range(len(result)):
        line = result[jj][0]
        matrix_out[line, :] = result[jj][1]

    np.save(args.output, matrix_out)
else:
    test= np.load(args.matrix)
    runargs = []
    for i in range(test.shape[0]):
        runargs.append(i)

    pool = Pool(args.np)

    result = pool.map(run2, runargs)
    # print(result)
    matrix_out = np.zeros((len(mobiles), len(mobiles)))

    for jj in range(len(result)):
        line = result[jj][0]
        matrix_out[line, :] = result[jj][1]

    np.save(args.output, matrix_out)
