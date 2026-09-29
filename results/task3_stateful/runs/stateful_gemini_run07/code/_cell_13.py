for g in marker_genes:
    if g in adata.raw.var_names:
        idx = adata.raw.var_names.get_loc(g)
        # adata.raw.X is csr_matrix
        vals = adata.raw.X[:, idx].toarray().flatten()
        print(f"{g:8s}:", [f"C{c}: {vals[adata.obs['leiden'].values == str(c)].mean():.2f}" for c in range(6)])