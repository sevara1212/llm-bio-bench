# Check marker genes
markers = {
    'T_cells': ['CD3D', 'CD3E', 'CD3G'],
    'CD4_T': ['CD4', 'IL7R', 'CCR7'],
    'CD8_T': ['CD8A', 'CD8B'],
    'NK_cells': ['GNLY', 'NKG7', 'NCAM1'],
    'B_cells': ['MS4A1', 'CD19', 'CD79A'],
    'CD14_Mono': ['CD14', 'LYZ', 'S100A9', 'S100A8'],
    'FCGR3A_Mono': ['FCGR3A', 'MS4A7'],
    'DC': ['FCER1A', 'CST3', 'HLA-DQA1'],
    'Platelet': ['PPBP', 'PF4', 'GP9']
}

sc.tl.rank_genes_groups(adata, 'leiden_0.8', method='t-test')
import pandas as pd
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = pd.DataFrame({group: result['names'][group][:10] for group in groups})
print("Top 10 marker genes per cluster (res 0.8):")
print(top_genes)