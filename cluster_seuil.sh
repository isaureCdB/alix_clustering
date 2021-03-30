############################################
# inputs
############################################
# AAA_test_trigo.npy = internal coordinates for subset of AAA fragments
# AAA_trigo.npy = internal coordinates for all AAA fragments
# rmsd_AAAr_test.npy = exact pairwise RMSD of subset of AAA fragments


# get thresholds of delta in internal coordinate above which a pair has RMSD > 1A.
# output is a vector of length = number of internal coordinates.
# !!!! Use on a subset of fragments (< 10.000) !!! 
get_thresholds.py rmsd_AAAr_test.npy AAA_test_trigo.npy 1 threshold_AAA_test_1A.npy


# get boolean mask of which pair has all delta values below thresholds.
# Use on the full set of fragments.
filter_intcoor.py AAAr.npy AAA_trigo.npy threshold_AAA_test_1A.npy AAA_trigo_mask.npy

# filter by pairwise RMSDs the compatible fragments
filter_RMSD.py AAAr.npy AAA_trigo_mask.npy 1 AAA_RMSD_mask.npy
