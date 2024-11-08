# Load required packages 
library(monocle3)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(tidyverse)
library(patchwork)
library(Cairo)

# Set working path 
setwd("/mnt/sdb16t/pancancer_NKT/final/DNT_CD8_评分/DNT+CD8拟时序数据/拟时序_new")

# Read data
merged_CD8_DNT <- readRDS("/mnt/sdb16t/pancancer_NKT/final/DNT_CD8_评分/DNT+CD8拟时序数据/拟时序_new/merged_CD8_DNT.rds")
merged_CD8_DNT <- FindVariableFeatures(merged_CD8_DNT, selection.method = "vst", nfeatures = 2000)

# Normalize data 
merged_CD8_DNT <- NormalizeData(merged_CD8_DNT, normalization.method = "LogNormalize", scale.factor = 10000)

# Standardize all genes
all.genes <- rownames(merged_CD8_DNT)
merged_CD8_DNT <- ScaleData(merged_CD8_DNT, features = all.genes)

# Run PCA analysis
merged_CD8_DNT <- RunPCA(merged_CD8_DNT, features = VariableFeatures(object = merged_CD8_DNT))

# Perform batch correction using Harmony
merged <- RunHarmony(merged_CD8_DNT, "batch", plot_convergence = FALSE)

# Perform UMAP dimensionality reduction using Harmony results
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:15)

# Save the processed data (takes a long time, temporarily not saving)
# saveRDS(merged, file = "/mnt/sdb16t/pancancer_NKT/final/DNT_CD8_评分/DNT+CD8拟时序数据/拟时序_new/merged.rds")

# Obtain expression matrix
data.sample = merged

# Display the grouping of the leiden_renamed column
table(merged@meta.data$leiden_renamed)

# Convert leiden, celltype, and leiden_renamed columns to factors
merged$leiden <- as.factor(merged$leiden)
merged$celltype <- as.factor(merged$celltype)
merged$leiden_renamed <- as.factor(merged$leiden_renamed)

# Select the 2000 most variable features
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# Get the names of all genes
all.genes <- rownames(merged)

# Normalize all genes
merged <- ScaleData(merged, features = all.genes)

# Perform Principal Component Analysis (PCA)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

# Perform UMAP dimensionality reduction using the first 15 principal components
merged <- RunUMAP(merged, dims = 1:15)

# Define a custom set of 36 colors
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#F7F398', '#8C549C', '#585658', '#23452F', '#BD956A',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

# Plot the UMAP, colored by leiden_renamed
DimPlot(merged, reduction = "umap", group.by = "leiden_renamed", raster=FALSE, cols = my36colors)

# Save the UMAP image
ggsave("UMAP_plot.png")

# Use monocle3 for further analysis
data <- GetAssayData(data.sample, assay = "RNA", slot = "counts")

# Filter genes with zero expression, keeping genes expressed in at least 3 cells
data <- data[rowSums(data > 0) >= 3,]
dim(data)

# Set cell metadata and gene annotation
cell_metadata <- data.sample@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Construct a single-cell dataset
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Data preprocessing, selecting 100 principal components
cds <- preprocess_cds(cds, num_dim = 100)

# Perform dimensionality reduction using PCA
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# Plot UMAP
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "leiden_renamed")
p1
ggsave(filename = "UMAP1.pdf", height = 10, width = 10)

# Integrate UMAP coordinates
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data.sample, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

# Plot the integrated UMAP
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "leiden_renamed")
p2
ggsave(filename = "UMAP_integrated.pdf", height = 10, width = 10)

# UMAP colored by cell type
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "celltype")
ggsave(filename = "UMAP_celltype_integrated.pdf", height = 10, width = 10)

# Perform cell clustering analysis
cds <- cluster_cells(cds, resolution = 0.01, k = 20, random_seed = 18, verbose = T)

# Manually set clustering labels
cds@clusters$UMAP$clusters <- data.sample$leiden_renamed

# Plot clustering results
plot_cells(cds, color_cells_by = "partition")
ggsave(filename = "partition_plot.pdf", height = 10, width = 10)

p1 <- plot_cells(cds, group_cells_by = 'cluster')
p1
ggsave(filename = "cluster_plot.pdf", height = 10, width = 10)

# Perform pseudotime analysis
cds <- learn_graph(cds, verbose = T, use_partition = T, close_loop = F)
save(cds, file = "cds.Rdata")

# Plot pseudotime graph
p <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = T,
                label_leaves = T, label_branch_points = T, cell_size = 0.5, group_label_size = 4)
p
ggsave(filename = "pseudotime_plot.pdf", height = 10, width = 10)

# Plot UMAP colored by pseudotime
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = F)
p1
p2 <- plot_cells(cds, color_cells_by = "leiden_renamed", label_branch_points = FALSE, label_leaves = F) +
  scale_color_manual(values = my36colors)
p2 | p1

ggsave("monocle3_FigS3G2.pdf", width = 20, height = 8)
ggsave("monocle3_FigS3G2.png", width = 20, height = 8)

# Plot UMAP without legend, colored by pseudotime
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE) +
  NoLegend() + 
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot UMAP without legend, colored by cell type
p2 <- plot_cells(cds, color_cells_by = "leiden_renamed", show_trajectory_graph = FALSE, 
                 label_branch_points = FALSE, label_leaves = FALSE) +
  NoLegend() +
  theme_minimal() +
  scale_color_manual(values = my36colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Extract the legend
legend <- get_legend(
  plot_cells(cds, color_cells_by = "leiden_renamed", show_trajectory_graph = T) +
    scale_color_manual(values = my36colors)
)

# Combine UMAP plots and legend horizontally, with the legend on the right
final_plot <- plot_grid(p1, p2, legend, ncol = 3, rel_widths = c(1, 1, 0.2))

# Display the combined plot
print(final_plot)

# Save the final image
ggsave("monocle3_2.pdf", final_plot, width = 20, height = 8)
ggsave("monocle3_2.png", final_plot, width = 20, height = 8)

# Plot trajectory changes for the specified gene
Track_genes_sig2 <- c("CCR7",
                      "SELL",
                      "TCF7",
                      "LEF1")

# Plot gene expression trends
pdf("genes_in_pseudotime.pdf", width = 8, height = 3)
plot_genes_in_pseudotime(cds[Track_genes_sig2,], color_cells_by = "leiden_renamed", 
                         min_expr = 0.5, ncol = 3)
dev.off()

# Save the gene expression trends plot as a PNG for higher resolution display
CairoPNG("genes_in_pseudotime.png", width = 1600, height = 1200, res = 150)
# Plot the same gene expression trends
plot_genes_in_pseudotime(cds[Track_genes_sig2,], 
                         color_cells_by = "leiden_renamed", 
                         min_expr = 0.5, 
                         ncol = 3)

# Close the PNG device to complete image saving
dev.off()
save(cds, file = "cds_final.Rdata")