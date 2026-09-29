# Check batch effect across samples
import matplotlib.pyplot as plt
adata.obs['sample'] = adata.obs_names.map(lambda x: x.split('_')[0])
# Let's check markers
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = pd.DataFrame({group: [result['names'][group][i] for i in range(10)] for group in groups})
print("Top genes per leiden 0.5 cluster:")
print(top_genes.iloc[:, :10])
print(top_genes.iloc[:, 10:])