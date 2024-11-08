library(Seurat)
library(SeuratData)
write.csv(t(as.matrix(DNT_samplegroup@assays$RNA@counts)),file = "DNT_samplegroupfor.scenic.data.csv")

####################################################
###################################
##################################
conda activate pyscenic
python
import os,sys
import loompy as lp
import numpy as np
import scanpy as sc
dir= '/path/'
x=sc.read_csv("/path/DNT_samplegroupfor.scenic.data.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create(dir+"sample.loom",x.X.transpose(),row_attrs,col_attrs);
import loompy
ds = loompy.connect(dir+"sample.loom")

#####################
cd /path
dir='/path/' 
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl 
input_loom=./sample.loom

#2.1 grn
pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$tfs 

#2.2 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20 \
--mask_dropouts

#2.3 AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 20


rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(devtools)
library(S4Arrays)
loom <- open_loom('out_SCENIC.loom') 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons") ##########33
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
head(regulonAUC)[1:3,1:3]
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)  
embeddings
close_loom(loom)
rownames(regulonAUC)
names(regulons)
library(SeuratData) 
DNT_samplegroup <- readRDS("/DNT_samplegroup.rds")
seurat.data = DNT_samplegroup
seurat_obj=seurat.data
table(seurat.data$sample_type)
seurat.data$group=seurat.data$sample_type
table(seurat.data$group)
##################################
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
identical(colnames(sub_regulonAUC), colnames(seurat.data))
cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           group = as.character(seurat.data$group))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        group = seurat.data$group)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = 'for_rss_and_visual.Rdata')

table(seurat.data$group)
dim(sub_regulonAUC)
sub_regulonAUC <- sub_regulonAUC[, match(colnames(seurat.data), colnames(sub_regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(seurat.data)) 
group_info <- seurat.data$group
if (is(sub_regulonAUC, "SummarizedExperiment")) {
  sub_regulonAUC_matrix <- assay(sub_regulonAUC)
} else {
  stop("sub_regulonAUC is not a SummarizedExperiment object")
}

avg_auc <- sub_regulonAUC_matrix %>% 
  as.data.frame() %>% 
  t() %>%  
  as.data.frame() %>%  
  mutate(group = group_info) %>%  
  group_by(group) %>%  
  summarize_all(mean)
avg_auc_matrix <- as.matrix(avg_auc[,-1])
rownames(avg_auc_matrix) <- avg_auc$group

new_order <- c("Met","PT","NT")  
new_order <- new_order[new_order %in% avg_auc$group]  
avg_auc_matrix <- avg_auc_matrix[new_order, ]
write.csv(avg_auc_matrix, file = "avg_auc_matrix.ori.csv", row.names = TRUE)

##############################
selectedResolution<-"group" 
cellsPerGroup <- split(rownames(cellTypes),cellTypes[,selectedResolution])
sub_regulonAUC<-sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
dim(sub_regulonAUC)
rss<-calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),selectedResolution])

rss=na.omit(rss)


pdf("PT_RSS_plot.pdf", width = 6, height = 6)
plotRSS_oneSet(rss,setName="PT")#clusterID
dev.off()


pdf("Met_RSS_plot.pdf", width = 6, height = 6)
plotRSS_oneSet(rss,setName="Met")#clusterID
dev.off()
