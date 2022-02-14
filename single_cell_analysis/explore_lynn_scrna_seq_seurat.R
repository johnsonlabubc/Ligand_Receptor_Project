# explore Rdata file seurat object containing human islet 
# scRNA-seq from Lynn lab


library(Seurat)
library(SeuratObject)
library(tidyverse)
library(cowplot)


# open the 4.6 GB RDS file containing the already processed data
filename <- file.choose()
seurat_obj <- readRDS(filename)


# view plot of UMAP showing clustering of different cell types
Seurat::DimPlot(seurat_obj)

ggsave("single_cell_analysis/figures/UMAP_human_islet_scRNAseq.png",
       width = 1600, height = 1200, units = "px")



# Get cell and feature names, and total numbers
head(colnames(seurat_obj))
head(Cells(seurat_obj))
# colnames are same as cell names


rownames(seurat_obj)
# row names appear to be HGNC gene symbols

seurat_obj$cell_type

ncol(seurat_obj)
nrow(seurat_obj)
# 68650 columns (which I think is the number of cells?)
# 3000 rows (which I think is the number of genes?)

# view different cell types
head(Idents(seurat_obj))

# increase memory limit & free up RAM used by R
# cause the DoHeatmap funciton was giving error: 
# "Error: cannot allocate vector of size 276 Kb"
gc()
memory.limit(size=1000000)

DoHeatmap(seurat_obj)
