ndist=4         # number of ditances in the internal coordinates
percent=0.999   # percent accepted false negative
cutoff=1        # RMSD cutoff in A
coor=AAA-dr0.2r.npy
nsample=4771

############################################
# get threshold on a sample of conformers
############################################
# get thresholds of delta in internal coordinate above which a pair has RMSD > 1A.
# output is a vector of length = number of internal coordinates.
# !!!! Use on a subset of fragments (< 10.000) !!! 
sample.py $coor $nsample sample.npy
pairwise_rmsd.py sample.npy mask_rmsd_sample.npy

get_internal_coordinate.py sample.npy intcoor_sample.npy 7 7 7 
get_thresholds.py mask_rmsd_sample.npy intcoor_sample.npy $ndist $cutoff $percent > thresholds-${cutoff}A-$percent.txt

############################################
# apply on full set of conformers
############################################
# get boolean mask of which pair has all delta values below thresholds.
# Use on the full set of fragments.

get_internal_coordinate.py $coor intcoor.npy 7 7 7 
filter_intcoor.py intcoor.npy thresholds-${cutoff}A-$percent.txt $ndist mask_intcoor.npy

############################################
# filter by pairwise RMSDs the compatible fragments.
# get boolean mask of which pair has real RMSD < 1A
filter_RMSD.py $coor mask_intcoor.npy $cutoff mask_rmsd.npy

############################################
#cluster with full boolean RMSD matrix
cluster_fullmatrix.py mask_rmsd.npy clusters > centers