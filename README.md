# Single-Cell RNA-seq Analysis

Single-cell RNA-seq workflows using Scanpy (Python) and Seurat (R) on 10x Genomics PBMC data.

## Analyses

### 1. Standard Workflow (`notebooks/pbmc3k_analysis.ipynb`)

Single-sample analysis of 2,700 PBMCs: QC, normalization, HVG selection, PCA, UMAP, Leiden clustering, marker gene identification, cell type annotation, and differential expression.

### 2. Multi-Sample Integration (`notebooks/multi_sample_integration.ipynb`)

Integration of PBMC 3k and PBMC 4k datasets (two different donors, ~7k cells total) using Harmony batch correction. Demonstrates handling batch effects, joint clustering across samples, cell type proportion comparison, and cross-donor differential expression within cell types.

### 3. Seurat Analysis (`scripts/seurat_pbmc3k.R`)

Same standard workflow in R/Seurat for comparison.

## Setup

```bash
conda env create -f environment.yml
conda activate scrnaseq
```

Or with pip:
```bash
pip install -r requirements.txt
```

## Run

```bash
jupyter notebook notebooks/pbmc3k_analysis.ipynb
jupyter notebook notebooks/multi_sample_integration.ipynb
Rscript scripts/seurat_pbmc3k.R
```

## Results

Figures saved to `results/`. Key outputs:

- QC violin plots and scatter plots
- UMAP embeddings (before/after batch correction)
- Cluster marker dotplots
- Cell type proportion comparison across donors
- Differential expression results
