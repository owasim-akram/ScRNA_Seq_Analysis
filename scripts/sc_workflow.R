# Step 1. Download the data
# Data Source: https://www.nature.com/articles/s41598-020-64929-x
#options(timeout = 3600)

#download.file(
  #  url = "https://single-cell-transcriptomics.s3.eu-central-1.amazonaws.com/course_data.tar.gz",
 #   destfile = "data.tar.gz"
#)
# untar
#untar("data.tar.gz")

# remove compressed files
#file.remove("data.tar.gz")

# Step 2: Package Installation
install.packages("pak")
pak::pkg_install(c(
    "tidyverse",
    "Seurat",
    "patchwork",
    "Matrix",
    "celldex",
    "SingleR",
    "edgeR",
    "limma",
    "scuttle",
    "presto",
    "clusterProfiler",
    "ggiraph",
    "enrichplot",
    "clustree",
    "dittoSeq"
))
# Step 3: Load required packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(Matrix)
library(clustree)
library(celldex)
library(SingleR)
library(edgeR)
library(limma)
library(scuttle)
library(presto)
library(clusterProfiler)
library(ggiraph)
library(enrichplot)
library(dittoSeq)
# Step 4: Load the data
sample_info <- read.csv("data/sample_info.csv")

data_dirs <- file.path(
    "data",
    "count_matrices",
    sample_info$SampleName,
    "outs",
    "filtered_feature_bc_matrix"
)

names(data_dirs) <- gsub("_", "-", sample_info$SampleName)
data_dirs <- data_dirs[1:3]
# Sparse matrix
sparse_matrix <- Read10X(data.dir = data_dirs)
sparse_matrix[c("PECAM1", "CD8A", "TSPAN1"), 1:30]

# Step 5: Create seaurat object
seu_obj <- CreateSeuratObject(
    counts  = sparse_matrix,
    project = "pbmc",
    min.cells = 3,
    min.features = 100
)
seu_obj
dim(sparse_matrix)
dim(seu_obj)

# cell Metadata exploration
# seurat_object <- Seurat::function(seurat_object)
head(seu_obj@meta.data)

# Step 6. Perform QC
# Visualizing QC per cell and gene
FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(seu_obj, features = c("nCount_RNA",
                          "nFeature_RNA"))
# mitochondrial genes
seu_obj <- PercentageFeatureSet(seu_obj,
                            pattern = "^MT-",
                            col.name = "percent_mito")
# ribosomal genes
seu_obj <- PercentageFeatureSet(seu_obj,
                            pattern = "^RP[SL]",
                            col.name = "percent_ribo")
# hemoglobin genes (but not HBP)
seu_obj <- PercentageFeatureSet(seu_obj,
                            pattern = "^HB[^(P)]",
                            col.name = "percent_globin")

# Visualize
VlnPlot(seu_obj, features = c("percent_mito",
                          "percent_ribo",
                          "percent_globin"))
FeatureScatter(seu_obj,
               feature1 = "percent_globin",
               feature2 = "percent_ribo")
# Gene expression
most_expressed_boxplot <- function(object, ngenes = 20){

    # matrix of raw counts
    cts <- Seurat::GetAssayData(object, assay = "RNA", layer = "counts")

    # get percentage/cell
    cts <- t(cts)/colSums(cts)*100
    medians <- SparseArray::colMedians(cts)

    # get top n genes
    most_expressed <- order(medians, decreasing = T)[ngenes:1]
    most_exp_matrix <- as.matrix((cts[,most_expressed]))

    # prepare for plotting
    most_exp_df <- stack(as.data.frame(most_exp_matrix))
    colnames(most_exp_df) <- c("perc_total", "gene")

    # boxplot with ggplot2
    boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
        geom_boxplot() +
        coord_flip()
    return(boxplot)
}

most_expressed_boxplot(seu_obj, 20)

# Step 7. Filtering
seu_obj <- subset(seu_obj, subset = nFeature_RNA > 200 &
                  nFeature_RNA < 5000 &
                  percent_mito < 8)
# Visualize after filtering
VlnPlot(seu_obj, features = c("nFeature_RNA",
                          "percent_mito"))
# Save the filtered data
saveRDS(seu_obj, "data/filtered_data.rds")

# Step 8. Normalization and scaling
seu_obj <- readRDS("data/filtered_data.rds")

seu_obj <- NormalizeData(seu_obj)
# Inspect normalize data
GetAssayData(seu_obj)
GetAssayData(seu_obj, layer = "counts")
GetAssayData(seu_obj, layer = "data")

# Step 9: HVG
seu_obj <- FindVariableFeatures(seu_obj,
                            selection.method = "vst",
                            nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu_obj), 10)
top10

# Plot top10
vf_plot <- VariableFeaturePlot(seu_obj)

# Build the labeled plot
p <- LabelPoints(
    plot = vf_plot,
    points = top10,
    repel = TRUE,
    xnudge = 0,
    ynudge = 0
)

# Save it
ggplot2::ggsave(
    filename = "figures/variable_features.png",
    plot = p,
    width = 7,
    height = 5,
    dpi = 300,
    bg = "white"
)

# Scaling the data
seu_obj <- ScaleData(seu_obj)

# Save normalized data
saveRDS(seu_obj, "data/normalized_data.rds")

# Step 10: Dimensionality reduction
seu_obj <- readRDS("data/normalized_data.rds")

# Dimension reduction
seu_obj <- RunPCA(seu_obj)
DimPlot(seu_obj, reduction = "pca")
FeaturePlot(seu_obj, reduction = "pca", features = "percent_globin")
DimHeatmap(seu_obj, dims = 1:12, cells = 500, balanced = TRUE)

# Determine PCs
ElbowPlot(seu_obj, ndims = 40) +
    geom_vline(xintercept = 25,
               linetype = "dashed",
               color = "red",
               linewidth = 1)
# UMAP
seu_obj <- RunUMAP(seu_obj, dims = 1:25)
DimPlot(seu_obj, reduction = "umap")

# Save the data
saveRDS(seu_obj, "data/dimred_data.rds")

# Step 11. Integration
seu_obj <- readRDS("data/dimred_data.rds")

# Data integration
seu_obj[["RNA"]] <- split(seu_obj[["RNA"]], f = seu_obj$orig.ident)
seu_obj <- IntegrateLayers(
    object = seu_obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
)
# re-join layers after integration
seu_obj[["RNA"]] <- JoinLayers(seu_obj[["RNA"]])
seu_obj <- RunUMAP(seu_obj, dims = 1:25,
               n.components = 3, ,
               reduction = "integrated.cca",
               reduction.name = "umap")
DimPlot(seu_obj, reduction = "umap")

# Cleaned-up code with proper checks
seu_obj <- FindNeighbors(seu_obj, reduction = "integrated.cca", dims = 1:25)
seu_obj <- FindClusters(seu_obj, resolution = 0.5)
seu_obj <- RunUMAP(seu_obj, reduction = "integrated.cca", dims = 1:25,
               reduction.name = "umap.cca")
# 1. Mixing: shuffle so no sample hides another, then view each sample on its own
DimPlot(seu_obj, reduction = "umap.cca", group.by = "orig.ident", shuffle = TRUE)
DimPlot(seu_obj, reduction = "umap.cca", split.by = "orig.ident")
# 2. Cell types still separate? Use known markers
FeaturePlot(seu_obj, reduction = "umap.cca",
            features = c("CD3E", "MS4A1", "CD14", "HBB", "CD34", "NKG7"))

# 3. Numbers: what fraction of each cluster comes from each sample?
prop.table(table(seu_obj$seurat_clusters, seu_obj$orig.ident), margin = 1)

# save the data
saveRDS(seu_obj, "data/integrated_data.rds")

# Step 12. Clustering
seu_obj <- readRDS("data/integrated_data.rds")

# Find FindNeighbors
seu_obj <- FindNeighbors(seu_obj, dims = 1:25, reduction = "integrated.cca")

# FindClusters
seu_obj <- FindClusters(seu_obj, resolution = seq(0.1, 0.8, by=0.1))

head(seu_obj@meta.data)

# Cluster tree
clustree(seu_obj@meta.data[, grep("RNA_snn_res", colnames(seu_obj@meta.data))],
         prefix = "RNA_snn_res.")

# UMAP coloring each cell according to a cluster id
DimPlot(seu_obj, group.by = "RNA_snn_res.0.1")

# save the data
saveRDS(seu_obj, "data/clustered_data.rds")

# Step 13. Cell annotation
seu_obj <- readRDS("data/clustered_data.rds")

# SetIdent
SetIdent(seu_obj, value = seu_obj$RNA_snn_res.0.3)
DefaultAssay(seu_obj) <- "RNA"
FeaturePlot(seu_obj, "HBA1")

# Markers
tcell_genes <- c("IL7R", "LTB", "TRAC", "CD3D")
monocyte_genes <- c("CD14", "CST3", "CD68", "CTSS")

FeaturePlot(seu_obj, tcell_genes, ncol = 2)

VlnPlot(seu_obj,
        features = tcell_genes,
        ncol = 2)
# AddModuleScore
seu_obj <- AddModuleScore(seu_obj,
                      features = list(tcell_genes),
                      name = "tcell_genes")

head(seu_obj@meta.data)

FeaturePlot(seu_obj, "tcell_genes1")
VlnPlot(seu_obj, "tcell_genes1")

# Annotating cells according to cycling phase
cc.genes
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

seu_obj <- CellCycleScoring(
    seu_obj,
    s.features = s_genes,
    g2m.features = g2m_genes
)

head(seu_obj@meta.data)

DimPlot(seu_obj, group.by = "Phase")

# Cell type annotation using SingleR
ref <- celldex::NovershternHematopoieticData()
class(ref)
table(ref$label.main)

seu_SingleR <- SingleR(
    test = GetAssayData(seu_obj),
    ref = ref,
    labels = ref$label.main
)

head(seu_SingleR)

# Visualize singleR score quality scores
plotScoreHeatmap(seu_SingleR)

# plotDeltaDistribution
plotDeltaDistribution(seu_SingleR)

# annotate
singleR_labels <- seu_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

head(seu_obj@meta.data)
# UMAP
# Add SingleR annotations to Seurat metadata
seu_obj$SingleR_annot <- seu_SingleR$labels
# Plot UMAP
dittoDimPlot(seu_obj, "SingleR_annot", size = 0.7)
