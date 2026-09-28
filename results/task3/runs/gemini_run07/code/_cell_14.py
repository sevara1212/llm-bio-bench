# Let's inspect cluster 2 in resolution 0.6:
# Is it NK cells or CD8 T cells, or both?
# Let's check expression of CD8A, CD8B, GNLY, NKG7, NCAM1 in Cluster 2.
import numpy as np

for g in ['CD8A', 'CD8B', 'GNLY', 'NKG7', 'CD3D', 'CD14', 'FCGR3A', 'MS4A1', 'FCER1A', 'PPBP']:
    if g in adata.raw.var_names:
        expr = np.asarray(adata.raw[:, g].X.todense()).flatten()
        print(f"{g}:")
        for cl in sorted(adata.obs['leiden'].unique()):
            mask = (adata.obs['leiden'] == cl).values
            print(f"  cl {cl}: mean={expr[mask].mean():.3f}, pct={(expr[mask] > 0).mean():.2f}")