
get_internal_coordinate.py AAAr.npy AAA_test_dist.npy 1 7 14 21 16 

############################################
# inputs
############################################
# AAA_test_trigo.npy = internal coordinates for subset of AAA fragments
# AAA_trigo.npy = internal coordinates for all AAA fragments
# rmsd_AAAr_test.npy = exact pairwise RMSDs of subset of AAA fragments

############################################
# get thresholds of delta in internal coordinate above which a pair has RMSD > 1A.
# output is a vector of length = number of internal coordinates.
# !!!! Use on a subset of fragments (< 10.000) !!! 
#
# usage: get_thresholds.py rmsd internal_coor Ndistances RMSD_cutoff percent_to_keep outp
get_thresholds.py rmsd_AAAr_test.npy AAA_4dist6angle.npy 4 1 0.999 threshold_AAA.npy

############################################
# get boolean mask of which pair has all delta values below thresholds.
# Use on the full set of fragments.
filter_intcoor.py coor_AAA_trigo.npy threshold_AAA.npy mask_AAA.npy

############################################
# filter by pairwise RMSDs the compatible fragments.
# get boolean mask of which pair has real RMSD < 1A
filter_RMSD.py AAAr.npy mask_AAA.npy 1 mask_AAA_rmsd.npy

############################################
#cluster with full boolean RMSD matrix
cluster_fullmatrix.py mask_AAA_rmsd.npy AAA.clust > AAA.centers