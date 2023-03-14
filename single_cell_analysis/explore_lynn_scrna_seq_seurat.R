# explore Rdata file seurat object containing human islet 
# scRNA-seq from Lynn lab


library(Seurat)
library(SeuratObject)
library(tidyverse)
library(cowplot)


# open the 4.6 GB RDS file containing the already processed data
filename <- file.choose()
seurat_obj <- readRDS(filename)

seurat_obj@meta.data

seurat_obj@assays$RNA

# view plot of UMAP showing clustering of different cell types
Seurat::DimPlot(seurat_obj, group.by = "cell_type",
                reduction = "umap") +
  labs(x = "UMAP 1",
       y = "UMAP 2") +
  # NoLegend() +
  #ggtitle("Unsorted SC-islets") +
  scale_color_viridis(discrete=TRUE) +
  ggtitle("human islets")

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/umap_h_islet.png",
       scale = 1.5)



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



############ dot plot of top ligrec ####################

# reopen ave counts of sorted 
seurat_ave_expr_df <- read.csv("single_cell_analysis/Islet_sc_rnaseq_francis/ave_counts_cell_types.tsv",
                                      sep = "\t")

# get final lig rec lists
ligands <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_hgnc_symbols.tsv",
                    sep = "\t") %>% 
  filter(type == "ligand")

receptors <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_hgnc_symbols.tsv",
                      sep = "\t") %>% 
  filter(type == "receptor")


seurat_ave_expr_lig <- seurat_ave_expr_df %>% 
  filter(hgnc_symbol %in% ligands$hgnc_symbol) %>% 
  arrange(desc(Beta)) %>% 
  head(25)

#plot ligands

DotPlot(seurat_obj,
        features = seurat_ave_expr_lig$hgnc_symbol,
        group.by = "cell_type",
        assay = "SCT",
        scale = TRUE) +
  scale_color_distiller(palette = "BuPu",
                        direction = 1) +
  RotatedAxis() +
  coord_flip() +
  scale_x_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/h_islet_dotplot_lig.png", 
       width = 5, height = 8,
       scale = 1)









