# Single-Cell RNA-seq Analysis: PBMC 3k

Standard single-cell RNA-seq analysis workflow on the 10x Genomics PBMC 3k dataset using Scanpy (Python) and Seurat (R).

## Dataset

3,000 PBMCs from a healthy donor (10x Genomics). Pre-processed count matrix — no alignment step required.

## Analysis

- QC and cell filtering (mitochondrial content, gene counts)
- Normalization and highly variable gene selection
- PCA, neighbor graph construction, UMAP embedding
- Leiden clustering
- Marker gene identification and cell type annotation
- Differential expression between clusters

## How to run

```bash
pip install -r requirements.txt
jupyter notebook notebooks/pbmc3k_analysis.ipynb
```

For the Seurat analysis:
```R
Rscript scripts/seurat_pbmc3k.R
```

## Results

Figures and outputs saved to `results/`.
