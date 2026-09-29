# Let's inspect Leiden 0.5 marker genes
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('processed.h5ad')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_markers = pd.DataFrame(
    {group: result['names'][group][:10] for group in groups}
)
print("Top 10 markers for Leiden 0.5:")
print(top_markers)