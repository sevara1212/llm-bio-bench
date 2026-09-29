import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('adata_final.h5ad')
adata_raw = adata.raw.to_adata()
adata_raw.obs['leiden'] = adata.obs['leiden'].values

genes_check = ['CD4','CD8A','CD8B','IL7R','CCR7','SELL','LEF1','TCF7','CD27','CCL5','CST7','GZMK','GZMA']
df = pd.DataFrame(adata_raw[:, genes_check].X.toarray(), columns=genes_check, index=adata_raw.obs_names)
df['leiden'] = adata_raw.obs['leiden'].values
means = df.groupby('leiden')[genes_check].mean()
print(means.loc[['4','6']].round(2))

# check cluster sizes and top genes for 4,6,9,13,14,16,17
result = adata.uns['rank_genes_groups']
for cl in ['4','6','9','13','14','16','17']:
    print(f"\nCluster {cl} top genes:", [result['names'][cl][i] for i in range(10)])
    print("n cells:", (adata.obs['leiden']==cl).sum())