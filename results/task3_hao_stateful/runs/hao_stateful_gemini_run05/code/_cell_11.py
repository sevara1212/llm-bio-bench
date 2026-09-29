# Let's inspect mean expression of canonical markers across clusters in leiden_0.3
df_expr = pd.DataFrame(index=adata_proc.obs['leiden_0.3'].cat.categories)
for gene in present_markers:
    df_expr[gene] = [adata_proc[adata_proc.obs['leiden_0.3'] == cl, gene].X.mean() for cl in df_expr.index]

pd.set_option('display.max_columns', 35)
pd.set_option('display.width', 1000)
print(df_expr)