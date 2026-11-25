# Load packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(celldex)
#load built in data
pbmc_data <- Read10X(data.dir = "~/Bioinformatics/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
# Create Seurat Object
seu_obj <- CreateSeuratObject(
    counts = pbmc_data,
    project = "pbmc3k",
    min.cells = 3,
    min.features = 200
)

# What does data in a count matrix look like?
seu_obj

# Exploring Seurat object with "@"
seu_obj@assays
metadata <- seu_obj@meta.data

# Lets examine a few genes in the first thirty cells
pbmc_data[c("CD3D","TCL1A","MS4A1"), 1:30]

# Dense size
dense_size <- object.size(as.matrix(pbmc_data))
dense_size

# Sparse size
sparse_size <- object.size(pbmc_data)
sparse_size

dense_size / sparse_size

# Explore metadata
# nCount_RNA: Total UMI counts per cell (expect: 1,000-50,000)
# nFeature_RNA: Unique genes detected per cell (expect: 500-5,000)
# percent.mt: Mitochondrial gene percentage (expect: <5% for healthy cells)
# Doublet rate: Expected ~0.5-2% for our cell count
# Ambient RNA contamination: Expected <5% for well-prepared PBMCs
head(seu_obj@meta.data)

# Standard pre-processing workflow
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Mitochondrial
seu_obj[["percent_mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
head(seu_obj@meta.data)

# Ribosomal
seu_obj[["percent_ribo"]] <- PercentageFeatureSet(seu_obj, pattern = "^RP[SL]")
head(seu_obj@meta.data)

# Hemoglobin
seu_obj[["percent_hb"]] <- PercentageFeatureSet(seu_obj, pattern = "^HB[^(P|E|S)]")
head(seu_obj@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(seu_obj, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_ribo", "percent_hb"), ncol = 5)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtering
seu_obj <- subset(seu_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)

# Normalizing the data
seu_obj <- NormalizeData(seu_obj)
seu_obj

# Identification of highly variable features (feature selection)
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu_obj), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu_obj)
plot2 <- LabelPoints(plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all_genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj)

# Perform linear dimensional reduction
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
VizDimLoadings(seu_obj, dims = 1:2, reduction = "pca")
DimPlot(seu_obj, reduction = "pca") + NoLegend()
DimHeatmap(seu_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seu_obj, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(seu_obj)
DimHeatmap(seu_obj, dims = 7:10, cells = 500, balanced = TRUE)

# Cluster the cells
seu_obj <- FindNeighbors(seu_obj, dims = 1:10)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)

# non-linear dimensional reduction
seu_obj <- RunUMAP(seu_obj, dims = 1:10)
DimPlot(seu_obj, reduction = "umap")

