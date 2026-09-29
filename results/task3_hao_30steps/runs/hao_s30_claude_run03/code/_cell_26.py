print("""
==========================================
PBMC scRNA-seq Analysis Summary
==========================================

WORKFLOW:
1. QUALITY CONTROL
   - Loaded raw counts: 7,841 cells x 20,264 genes
   - Filtered: min_genes=200 per cell, min_cells=3 per gene
   - Removed cells with >15% mitochondrial reads (likely dying/stressed cells)
   - Removed likely doublets (>6000 genes/cell)
   - Result: 7,838 cells x 19,064 genes passed QC (3 cells removed)

2. NORMALIZATION
   - Library-size normalization (target_sum=1e4)
   - Log1p transformation
   - Selected 2,000 highly variable genes (Seurat flavor)

3. DIMENSIONALITY REDUCTION & CLUSTERING
   - Scaled data (max_value=10 clipping)
   - PCA: 50 components on HVGs
   - Neighbors graph: 15 neighbors, 30 PCs
   - Leiden clustering at resolution=0.5 -> 20 clusters
   - UMAP embedding for visualization

4. MARKER GENE IDENTIFICATION
   - Wilcoxon rank-sum test (rank_genes_groups) per cluster
   - Cross-referenced canonical PBMC markers:
     T cells: CD3D/E, CD4, CD8A/B, IL7R, GZMK, FOXP3
     NK cells: NKG7, GNLY, KLRD1, KLRF1, FCGR3A
     B cells: MS4A1, CD79A/B, IGHD, TCL1A (naive) vs BANK1 (memory)
     Plasma cells: MZB1, JCHAIN, TNFRSF17
     Monocytes: CD14, LYZ, S100A8/9 (classical) vs MS4A7, FCGR3A (non-classical)
     DCs: CST3, FCER1A, CD1C (cDC2) vs CLEC9A, BATF3 (cDC1); LILRA4, IL3RA (pDC)
     Progenitors: CD34, PRSS57
     Platelets: PPBP, PF4, NRGN
     Erythrocytes: HBB, HBA1

5. CELL TYPE ASSIGNMENT (20 clusters -> 20 labels)
   Cluster sizes range from 85 (Erythrocyte) to 1,213 (CD4+ Naive T)
   Full immune cell diversity captured: T cell subsets (naive/memory/
   effector/Treg/proliferating), NK subsets (CD56bright/dim), B cells
   (naive/memory/plasma), myeloid (CD14+/CD16+ mono, cDC1/cDC2/pDC),
   HSPCs, platelets, and rare erythrocyte contamination.

OUTPUT: labels.csv
   - Columns: barcode, cell_type
   - 7,841 rows (all original barcodes, matching raw data 1:1)
   - 3 cells that failed QC labeled "Low quality (filtered out)"
   - No missing values, no duplicate barcodes
==========================================
""")