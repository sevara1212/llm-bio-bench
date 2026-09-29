# Let's inspect CD4, CD8A, CD8B, CD14, FCGR3A, NCAM1, MS4A1, CD19, CD3D, FOXP3, etc. across clusters
genes = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'CCR7', 'FOXP3', 'NKG7', 'GNLY', 'NCAM1', 
         'MS4A1', 'CD19', 'CD14', 'FCGR3A', 'HLA-DRA', 'CLEC9A', 'CLEC10A', 'IL3RA', 'MZB1', 
         'PPBP', 'HBB', 'MKI67', 'CD34', 'GZMK']

# Check which genes are present
genes = [g for g in genes if g in adata.var_names]

# Get mean expression by cluster
import scipy.sparse as sp
mean_df = pd.DataFrame(index=sorted(adata.obs['leiden'].unique(), key=int))
for g in genes:
    idx = adata.var_names.get_loc(g)
    vals = adata.X[:, idx]
    if sp.issparse(vals):
        vals = vals.toarray().flatten()
    mean_df[g] = [vals[adata.obs['leiden'] == c].mean() for c in mean_df.index]

print(mean_df[['CD3D', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'NKG7', 'GNLY', 'MS4A1', 'CD14', 'FCGR3A']])
print(mean_df[['CLEC9A', 'IL3RA', 'MZB1', 'PPBP', 'HBB', 'MKI67', 'GZMK']])