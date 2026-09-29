# Let's inspect the dataset metadata if any or see if this is Hao et al. 2021 (CITE-seq human PBMC)
# Note the barcodes: L3_, E2L7_, E2L2_, E2L4_, L5_, etc. - classic Hao et al. 2021 (PBMC 10k or multi-lane PBMC)!
import scanpy as sc
adata = sc.read_h5ad('raw_counts.h5ad')
print(adata.obs_names.map(lambda x: x.split('_')[0]).value_counts())