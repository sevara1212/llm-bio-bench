import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')

# Let's inspect raw unscaled expression for these clusters to see absolute expression levels
raw_df = pd.DataFrame(
    adata.raw.X[:, [adata.raw.var_names.get_loc(g) for g in ['CD3D', 'CD4', 'CD8A', 'CD8B', 'NCAM1', 'FCGR3A', 'KLRC1', 'KLRD1', 'TRDC', 'TRGC1', 'TRAC', 'FOXP3', 'IL2RA'] if g in adata.raw.var_names]].toarray(),
    columns=[g for g in ['CD3D', 'CD4', 'CD8A', 'CD8B', 'NCAM1', 'FCGR3A', 'KLRC1', 'KLRD1', 'TRDC', 'TRGC1', 'TRAC', 'FOXP3', 'IL2RA'] if g in adata.raw.var_names],
    index=adata.obs_names
)
raw_df['leiden_0.5'] = adata.obs['leiden_0.5']

for cl in [0, 2, 3, 4, 6, 8, 14, 15, 17]:
    sub = raw_df[raw_df['leiden_0.5'] == str(cl)].drop(columns='leiden_0.5')
    print(f"Cluster {cl} mean raw normalized expr:\n{(sub > 0).mean().round(2)}")