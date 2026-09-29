# Let's inspect cluster 0 and ribosomal genes.
# In single-cell PBMC analysis, ribosomal genes (RPS/RPL) and mitochondrial genes often dominate top marker lists if not excluded.
# Let's check mean expression of canonical markers across all clusters.
import pandas as pd
mean_exp = pd.DataFrame(index=adata_proc.obs['leiden_0.5'].cat.categories)
for gene in present_keys:
    if gene in adata_proc.raw.var_names:
        # compute mean log1p normalized expression in raw
        mean_exp[gene] = [adata_proc.raw[adata_proc.obs['leiden_0.5'] == cl, gene].X.mean() for cl in mean_exp.index]

print(mean_exp.round(2))