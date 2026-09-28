import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('adata_filtered.h5ad')

print("--- Res 0.8 ---")
sc.tl.rank_genes_groups(adata, groupby='leiden_0.8', method='wilcoxon')
result = adata.uns['rank_genes_groups']
for g in result['names'].dtype.names:
    genes = [result['names'][g][i] for i in range(10)]
    print(f"Cluster {g} (N={sum(adata.obs['leiden_0.8'] == g)}): {genes}")

print("\n--- Res 1.0 ---")
sc.tl.rank_genes_groups(adata, groupby='leiden_1.0', method='wilcoxon')
result1 = adata.uns['rank_genes_groups']
for g in result1['names'].dtype.names:
    genes = [result1['names'][g][i] for i in range(10)]
    print(f"Cluster {g} (N={sum(adata.obs['leiden_1.0'] == g)}): {genes}")

# Check key marker expression across clusters for res 0.8 and 1.0
key_genes = ['IL7R', 'CCR7', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 'CD8B', 'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3', 'PPBP', 'CD3D']
for g in key_genes:
    if g not in adata.raw.var_names:
        print(f"Missing {g}")