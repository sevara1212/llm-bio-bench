# Check all clusters from 0 to 20
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')
adata.raw = None  # drop the unlogged raw so rank_genes_groups uses log-normalized adata.X!
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')

result = adata.uns['rank_genes_groups']
groups = [str(i) for i in range(21)]
top_markers = pd.DataFrame(
    {group: result['names'][group][:10] for group in groups}
)
for col in top_markers.columns:
    print(f"Cluster {col}: {list(top_markers[col].values)}")