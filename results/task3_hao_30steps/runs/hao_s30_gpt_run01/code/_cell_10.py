import scanpy as sc, pandas as pd
raw=sc.read_h5ad('raw_counts.h5ad'); raw.var['mt']=raw.var_names.str.startswith('MT-');sc.pp.calculate_qc_metrics(raw,qc_vars=['mt'],inplace=True,log1p=False)
out=pd.read_csv('labels.csv').set_index('barcode')
x=raw.obs.join(out)
print(x[x.cell_type=='Low-quality cell'][['total_counts','n_genes_by_counts','pct_counts_mt']].describe())
print(x[x.cell_type=='Low-quality cell'].sort_values('total_counts',ascending=False)[['total_counts','n_genes_by_counts','pct_counts_mt']].head(20))
print(pd.read_csv('labels.csv').head());print('lines',sum(1 for _ in open('labels.csv')))