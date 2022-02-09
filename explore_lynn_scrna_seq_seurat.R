# explore Rdata file seurat object containing human islet 
# scRNA-seq from Lynn lab



library(Seurat)
library(SeuratObject)




# open the 4.6 GB RDS file containing the already processed ata
filename <- file.choose()
seurat_obj <- readRDS(filename)


# view plot of UMAP showing clustering of different cell types
Seurat::DimPlot(seurat_obj)
