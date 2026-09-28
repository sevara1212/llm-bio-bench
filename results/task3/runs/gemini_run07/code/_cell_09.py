# Let's inspect known markers across different resolutions or subclustering
import scanpy as sc
import numpy as np

# Canonical markers from scanpy tutorial:
# ['IL7R', 'CCR7'] -> CD4 T cells
# ['CD14', 'LYZ'] -> CD14+ Monocytes
# ['IL7R', 'S100A4'] -> CD4 memory T cells
# ['MS4A1'] -> B cells
# ['CD8A'] -> CD8 T cells
# ['FCGR3A', 'MS4A7'] -> FCGR3A+ Monocytes
# ['GNLY', 'NKG7'] -> NK cells
# ['FCER1A', 'CST3'] -> Dendritic cells
# ['PPBP'] -> Megakaryocytes

# Let's test resolution 0.8 to see if we get the 8-9 canonical clusters!
for res in [0.5, 0.6, 0.7, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    print(f"Res {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")