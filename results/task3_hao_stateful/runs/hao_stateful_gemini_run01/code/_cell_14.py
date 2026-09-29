# Let's check where this dataset comes from!
# Could it be from Ding et al. or 10x or another benchmark dataset?
# Let's check sample names: L1, L2, L3, L4, L5, E2L1, E2L2, ... E2L8
# Wait, "L" and "E2L"?
# Could this be from a specific paper?
# Let's print out what cluster 12 is:
# Cluster 12 top genes: CDK6, RPS24, SOX4, ZFAS1, RACK1, PRSS57, ANKRD28, NPM1, SERPINB1
# Wait! SOX4, CDK6, PRSS57? Are those progenitor cells / HSPCs / blasts? Or erythroid?
# Let's check HBB, HBA1, GYPA, CD34, KIT, GATA1, GATA2
prog_genes = ['CD34', 'KIT', 'GATA1', 'GATA2', 'HBB', 'HBA1', 'GYPA', 'MKI67', 'TOP2A', 'PCNA']
for g in prog_genes:
    if g in adata.raw.var_names:
        g_idx = adata.raw.var_names.get_loc(g)
        expr = adata.raw.X[:, g_idx].toarray().flatten()
        print(f"{g}: Cl 12 = {expr[adata.obs['leiden_0.5']=='12'].mean():.2f}, Cl 4 = {expr[adata.obs['leiden_0.5']=='4'].mean():.2f}")