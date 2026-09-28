# Let's inspect markers at resolution 0.6, 0.8, and 1.0
for res in [0.6, 0.8, 1.0]:
    sc.tl.rank_genes_groups(adata, f'leiden_{res}', method='t-test')
    result = adata.uns['rank_genes_groups']
    print(f"=== Resolution {res} ===")
    for g in result['names'].dtype.names:
        print(f"Cluster {g} (n={(adata.obs[f'leiden_{res}']==g).sum()}):", list(result['names'][g][:6]))