# Let's inspect the dataset origin. The barcodes have prefixes like L1_, L2_, L3_, L5_, E2L1_, E2L2_, ...
# This is a classic benchmarking dataset (Ding et al. 2020 Nature Biotechnology - systematic comparison of single-cell and single-nucleus RNA-sequencing methods, or similar PBMC dataset, e.g. 10x Chromium v2 / v3 etc.)
# Let's see how many samples/libraries there are:
prefixes = [x.split('_')[0] for x in adata.obs_names]
from collections import Counter
print(Counter(prefixes))