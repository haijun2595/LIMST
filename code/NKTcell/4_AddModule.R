# Load packages #
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(tidyverse)
library(stringr)
library(harmony)
library(ggplot2)
library(ggpubr)
library(AUCell) 
library(clusterProfiler)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsignif)
DNT <- readRDS("path/DNT.rds")
CD8 <- readRDS("path/CD8_cells_selected_samplegroup_clusters.rds")

Tcell <- merge(CD8, DNT)
Tcell@active.ident <- as.factor(Tcell@meta.data$leiden_renamed)
table(Tcell@active.ident)
Tcell <- NormalizeData(Tcell)
rm(list = ls()[ls() != "Tcell"])
saveRDS(Tcell, file = "/mnt/sdb16t/pancancer_NKT/final/DNT_CD8_评分/Tcell.rds")

# Scoring TEFF #

# Assign value #
scRNA <- Tcell
# TEFF scoring #
gene_list <- read.table("path/TEFF.txt", header = T)
table(scRNA@meta.data$leiden_renamed)
Idents(scRNA) <- "leiden_renamed"
# Convert genes to character format #
genes_vector <- as.character(gene_list$gene)
# AddModuleScore scoring #

DefaultAssay(scRNA) <- "RNA"

scRNA <- AddModuleScore(scRNA,
                        features = gene_list, 
                        ctrl = 100, 
                        name = "AddModuleScore")
## Rename metadata columns ####
meta = scRNA@meta.data
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == "AddModuleScore1"] <- "TEFF_Score"
table(scRNA@meta.data$leiden_renamed)
table(scRNA@meta.data$celltype)
# Extract required data #
data <- FetchData(scRNA, vars = c("leiden_renamed", "TEFF_Score"))
# Check data #

## Create violin plot and add boxplot ##
p <- ggplot(data, aes(x = leiden_renamed, y = TEFF_Score, fill = leiden_renamed)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  #facet_wrap(~leiden_renamed, scales = "free_y") +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "white") +
  labs(title = "TEFF_Score") +
  xlab("leiden_renamed") +
  ylab("TEFF_Score") +
  scale_fill_manual(values = c(
    'CD8_Tn' = '#E59CC4', 'CD8_Tm_CXCR4' = '#585658', 'CD8_Tm' = '#AB3282', 
    'CD8_Tm_CD44' = '#5F3D69', 'CD8_Tm_IL7R' = '#58A4C3', 'CD8_Teff' = '#E63863', 
    'CD8_Teff_DUSP2' = '#8C549C', 'CD8_Tex' = '#BD956A', 
    'CD8_Tex_PDCD1' = '#E0D4CA', 'CD8_Tex_LAG3' = '#C5DEBA', 
    'CD8_Tex_MKI67' = '#9FA3A8', 'DNT' = '#F7F398')) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),  
        axis.line = element_line(color = "black"),  
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))  
p <- p + geom_signif(comparisons = list(c("DNT", "CD8_Teff"), 
                                        c("DNT", "CD8_Teff_DUSP2")),
                     map_signif_level = TRUE, 
                     test = "wilcox.test", 
                     textsize = 3, vjust = 0.5, y_position = c(1.3, 1.2))

ggsave(filename = "TEFF_Score.pdf", plot = p, width = 6, height = 4, dpi = 300)
dev.off()

# Scoring TEX #
# rm(list=ls())

scRNA <- Tcell
# TEX scoring #
gene_list <- read.table("path/TEX.txt", header = T)
table(scRNA@meta.data$leiden_renamed)
Idents(scRNA) <- "leiden_renamed"
# Convert genes to character format #
genes_vector <- as.character(gene_list$gene)

DefaultAssay(scRNA) <- "RNA"

scRNA <- AddModuleScore(scRNA,
                        features = gene_list, 
                        ctrl = 100, 
                        name = "AddModuleScore")
## Rename metadata columns ####
meta = scRNA@meta.data
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == "AddModuleScore1"] <- "TEX_Score"

table(scRNA@meta.data$leiden_renamed)
table(scRNA@meta.data$celltype)
# Extract required data #
data <- FetchData(scRNA, vars = c("leiden_renamed", "TEX_Score"))
# Check data #
head(data)
library(ggplot2)
library(ggsignif)
## Scoring
## Create violin plot and add boxplot ##
p <- ggplot(data, aes(x = leiden_renamed, y = TEX_Score, fill = leiden_renamed)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  #facet_wrap(~leiden_renamed, scales = "free_y") +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "white") +
  labs(title = "TEX_Score") +
  xlab("leiden_renamed") +
  ylab("TEX_Score") +
  scale_fill_manual(values = c(
    'CD8_Tn' = '#E59CC4', 'CD8_Tm_CXCR4' = '#585658', 'CD8_Tm' = '#AB3282', 
    'CD8_Tm_CD44' = '#5F3D69', 'CD8_Tm_IL7R' = '#58A4C3', 'CD8_Teff' = '#E63863', 
    'CD8_Teff_DUSP2' = '#8C549C', 'CD8_Tex' = '#BD956A', 
    'CD8_Tex_PDCD1' = '#E0D4CA', 'CD8_Tex_LAG3' = '#C5DEBA', 
    'CD8_Tex_MKI67' = '#9FA3A8', 'DNT' = '#F7F398')) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),  
        axis.line = element_line(color = "black"),  
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))  
p <- p + geom_signif(comparisons = list(c("DNT", "CD8_Tex"), 
                                        c("DNT", "CD8_Tex_LAG3"), 
                                        c("DNT", "CD8_Tex_PDCD1"), c("DNT", "CD8_Tex_MKI67")
                                        ),
                     map_signif_level = TRUE, 
                     test = "wilcox.test", 
                     textsize = 3, vjust = 0.5, y_position = c(1.45, 1.3, 1.20, 1.10))

# Display plot
print(p)
ggsave(filename = "P_Tex_Score.pdf", plot = p, width = 6, height = 4, dpi = 300)
dev.off()
