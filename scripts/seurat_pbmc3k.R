library(Seurat)
library(ggplot2)

# Download PBMC 3k from 10x
data_url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
data_dir <- "../data"

if (!file.exists(file.path(data_dir, "filtered_gene_bc_matrices"))) {
  download.file(data_url, destfile = file.path(data_dir, "pbmc3k.tar.gz"))
  untar(file.path(data_dir, "pbmc3k.tar.gz"), exdir = data_dir)
}

pbmc_data <- Read10X(data.dir = file.path(data_dir, "filtered_gene_bc_matrices/hg19"))
pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)

# QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf("../results/seurat_qc_violin.pdf", width = 10, height = 4)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale and PCA
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pdf("../results/seurat_elbow.pdf", width = 6, height = 4)
ElbowPlot(pbmc, ndims = 50)
dev.off()

# Cluster
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

pdf("../results/seurat_umap_clusters.pdf", width = 7, height = 5)
DimPlot(pbmc, reduction = "umap", label = TRUE)
dev.off()

# Marker genes
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)

pdf("../results/seurat_heatmap.pdf", width = 12, height = 8)
DoHeatmap(pbmc, features = top5$gene) + NoLegend()
dev.off()

# Cell type annotation based on canonical markers
new_ids <- c("CD4 T", "CD14 Mono", "B", "CD8 T", "FCGR3A Mono", "NK", "DC", "Platelets")
names(new_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_ids)

pdf("../results/seurat_umap_celltypes.pdf", width = 7, height = 5)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

# DE between monocyte subtypes
mono_de <- FindMarkers(pbmc, ident.1 = "CD14 Mono", ident.2 = "FCGR3A Mono")
write.csv(mono_de, "../results/seurat_de_monocytes.csv")

saveRDS(pbmc, file = "../results/pbmc3k_seurat.rds")
cat("Done.\n")
