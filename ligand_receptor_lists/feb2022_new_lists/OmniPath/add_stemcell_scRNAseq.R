# add the stem cell dervied cellsscRNA-seq data from Franics' lab to our main spreadsheet
# data was provided as seurat objects by Meltem


# load libraries
library(tidyverse)
library(Seurat)
library(cowplot)
library(viridis)

############ load main datatables to append new data to ########################

# load genes list with all prior expression data and annotations
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_annot_expr.tsv", 
                     sep = "\t") %>% 
  # drop columns with outdated info that we are replacing
  select(-keep_in_list,
         -comments)

# get latest manual annotations from the google spreadsheet
manual_annot <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_annot_manual_May2022.tsv", 
                         sep = "\t") %>% 
  # keep just the columns we want from here
  dplyr::select(hgnc_symbol,
         #keep_in_list = keep_in_list..CORRECTED.DATA.,
         keep_in_list,
         comments)

gene_list_annot <- gene_list %>%
  left_join(manual_annot,
            by = "hgnc_symbol") %>% 
  # reorder columns to put keep in list and comments back in same spot
  relocate(keep_in_list,
           comments,
           .after = description)


########### explore seurat object of stem cell data ###############

######### Unsorted data set

# open the  RDS file containing the already processed data
# start with the unsorted dataset
filename <- file.choose()
seurat_obj <- readRDS(filename)

# view different cell types
head(Idents(seurat_obj))


# view plot of UMAP showing clustering of different cell types
# in human islet seurat object the grouping was done by cell_type
# but in this dataset there is no "cell_type" variable, and the cell identities
# are just a number from 1 through 11

levels(seurat_obj@meta.data$clusters)
# in metadata we see they have been grouped into clusters into:
# "Immature Beta cells"    "INS+/eGFP+"   "INS+/eGFP+/GCG+/SST+"  
# "Pancreatic Progenitors" "UNK" 

# UMAP plot
 Seurat::DimPlot(seurat_obj, group.by = "clusters",
                reduction = "umap") +
  labs(x = "UMAP 1",
       y = "UMAP 2") +
  # NoLegend() +
  #ggtitle("Unsorted SC-islets") +
  scale_color_viridis(discrete=TRUE) +
  ggtitle("S7 unsorted cells")

# plot with clusters labelled on plot
# 
# LabelClusters(plot = unsorted_umap,
#               id = "clusters",
#               size = 5)

# save plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/umap_unsorted_sc_2.png",
       scale = 1.5)

ncol(seurat_obj)
nrow(seurat_obj)
# 6875 columns (which I think is the number of cells?)
# 16328 rows (which I think is the number of genes?)


# view the details of the 3 assays
seurat_obj@assays
# they are identical in the number of cells and features (genes)
# confirms there are 6875 cells



# calculate ave expression
seurat_ave_expr <- AverageExpression(seurat_obj,
                                     # group by the predetermined cell types
                                     group.by = "clusters",
                                     # only do RNA, not SCT & integrated
                                     assays = "RNA", 
                                     # data is not log normalized
                                     # so don't do exponentiation
                                     slot = "counts")

View(seurat_ave_expr$RNA)

# convert to dataframe 
seurat_ave_expr_unsorted_df <- as.data.frame(seurat_ave_expr$RNA) %>% 
  # make rownames a column
  rownames_to_column(var = "hgnc_symbol") %>% 
  # rename columns
  rename(unsorted_immature_beta = "Immature Beta cells") %>% 
  rename(unsorted_ins_gfp = "INS+/eGFP+") %>%
  rename(unsorted_ins_gfp_gcg_sst = "INS+/eGFP+/GCG+/SST+") %>%
  rename(unsorted_pancreatic_proj = "Pancreatic Progenitors")


######### Repeat for sorted data set


# open the  RDS file containing the already processed data
# start with the sorted dataset
filename <- file.choose()
seurat_obj_sorted <- readRDS(filename)

# view different cell types
head(Idents(seurat_obj_sorted))


# view plot of UMAP showing clustering of different cell types
# in human islet seurat object the grouping was done by cell_type
# but in this dataset there is no "cell_type" variable, and the cell identities
# are just a number from 1 through 9

levels(seurat_obj_sorted@meta.data$clusters)
# in metadata we see they have been grouped into clusters into:
# "eBCs"                   "Immature Beta cells"    "INS+/eGFP+/GCG+/SST+"  
# "INS+/eGFP+/SST+"        "Pancreatic Progenitors"


# # UMAP plot
# Seurat::DimPlot(seurat_obj_sorted, group.by = "clusters") +
#   labs(x = "UMAP 1",
#        y = "UMAP 2") +
#   ggtitle("Sorted & aggregated SC-islets") +
#   scale_color_viridis(discrete=TRUE,
#                       # change eBCs to SC beta
#                       labels = c("SCβ-cells", 
#                                  "Immature Beta cells",
#                                  "INS+/eGFP+/GCG+/SST+",
#                                  "INS+/eGFP+/SST+",
#                                  "Pancreatic Progenitors"))



# UMAP plot
Seurat::DimPlot(seurat_obj_sorted, group.by = "clusters",
                                 reduction = "umap") +
  labs(x = "UMAP 1",
       y = "UMAP 2") +
 # NoLegend() +
  #ggtitle("Unsorted SC-islets") +
  scale_color_viridis(discrete=TRUE,
                      labels = c(expression(paste("SC", beta, "-cells")), 
                                 "Immature Beta cells",
                                 "INS+/eGFP+/GCG+/SST+",
                                 "INS+/eGFP+/SST+",
                                 "Pancreatic Progenitors")) +
  ggtitle("S7 sorted cells") +
  # re left align 
  theme(legend.text.align = 0)

# # plot with clusters labelled on plot
# 
# LabelClusters(plot = sorted_umap,
#               id = "clusters",
#               size = 5)

# save plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/umap_sorted_2.png",
       scale = 1.5)



ncol(seurat_obj_sorted)
nrow(seurat_obj_sorted)
# 5355 columns (which I think is the number of cells?)
# 15795 rows (which I think is the number of genes?)

# determine number of eBCs
sum(seurat_obj_sorted@meta.data$clusters == "eBCs")
# 891 eBCs


# calculate ave expression
seurat_ave_expr_sorted <- AverageExpression(seurat_obj_sorted,
                                     # group by the predetermined cell types
                                     group.by = "clusters",
                                     # only do RNA, not SCT & integrated
                                     assays = "RNA", 
                                     # data is not log normalized
                                     # so don't do exponentiation
                                     slot = "counts")

View(seurat_ave_expr_sorted$RNA)


# convert to dataframe 
seurat_ave_expr_sorted_df <- as.data.frame(seurat_ave_expr_sorted$RNA) %>% 
  # make rownames a column
  rownames_to_column(var = "hgnc_symbol") %>% 
  # rename columns
  rename(sorted_eBCs = "eBCs") %>% 
  rename(sorted_immature_beta = "Immature Beta cells") %>%
  rename(sorted_ins_gfp_gcg_sst = "INS+/eGFP+/GCG+/SST+") %>%
  rename(sorted_ins_gfp_sst = "INS+/eGFP+/SST+") %>% 
  rename(sorted_pancreatic_proj = "Pancreatic Progenitors")


## save average counts by cell type in scRNA-seq data

# unsorted
write_tsv(seurat_ave_expr_unsorted_df, 
          "single_cell_analysis/stem_cell_scrnaseq_francis/unsorted_ave_counts_clusters.tsv")

# sorted
write_tsv(seurat_ave_expr_sorted_df, 
          "single_cell_analysis/stem_cell_scrnaseq_francis/sorted_ave_counts_clusters.tsv")




############# dot plot of top expressed lig and rec ############################

# reopen ave counts of sorted 
seurat_ave_expr_sorted_df <- read.csv("single_cell_analysis/stem_cell_scrnaseq_francis/sorted_ave_counts_clusters.tsv",
                                      sep = "\t")

# get final lig rec lists
ligands <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_hgnc_symbols.tsv",
                                      sep = "\t") %>% 
  filter(type == "ligand")

receptors <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_hgnc_symbols.tsv",
                    sep = "\t") %>% 
  filter(type == "receptor")


seurat_ave_expr_sorted_lig <- seurat_ave_expr_sorted_df %>% 
  filter(hgnc_symbol %in% ligands$hgnc_symbol) %>% 
  arrange(desc(sorted_eBCs)) %>% 
  head(25)

 #plot receptors

DotPlot(seurat_obj_sorted,
        features = seurat_ave_expr_sorted_lig$hgnc_symbol,
        group.by = "clusters",
        assay = "SCT",
        scale = TRUE) +
   scale_color_distiller(palette = "BuPu",
                         direction = 1) +
  facet_wrap(facets = ligands) +
  RotatedAxis() +
  coord_flip() +
  scale_y_discrete(labels = c(expression(paste("SC", beta, "-cells")), 
                              "Immature Beta cells",
                              "INS+/eGFP+/GCG+/SST+",
                              "INS+/eGFP+/SST+",
                              "Pancreatic Progenitors")) +
  scale_x_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/sorted_sc_dotplot_lig.png", 
       width = 5, height = 8,
       scale = 1)




seurat_ave_expr_sorted_rec <- seurat_ave_expr_sorted_df %>% 
  filter(hgnc_symbol %in% receptors$hgnc_symbol) %>% 
  arrange(desc(sorted_eBCs)) %>% 
  head(25)

#plot receptors

DotPlot(seurat_obj_sorted,
        features = seurat_ave_expr_sorted_rec$hgnc_symbol,
        group.by = "clusters",
        assay = "SCT",
        scale = TRUE) +
  scale_color_distiller(palette = "BuGn",
                        direction = 1) +
  RotatedAxis() +
  coord_flip() +
  scale_y_discrete(labels = c(expression(paste("SC", beta, "-cells")), 
                              "Immature Beta cells",
                              "INS+/eGFP+/GCG+/SST+",
                              "INS+/eGFP+/SST+",
                              "Pancreatic Progenitors")) +
  scale_x_discrete(limits=rev) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/sorted_sc_dotplot_rec.png", 
       width = 5, height = 8,
       scale = 1)



###### append stem cell scRNA-seq sorted & unsorted ave counts to genes list ########
gene_list_annot_scdata <- left_join(gene_list_annot,
                              seurat_ave_expr_unsorted_df,
                              by = "hgnc_symbol") %>% 
  left_join(seurat_ave_expr_sorted_df,
            by = "hgnc_symbol") 




# replace all NA's in stemcell scRNA-seq data with zeros 
gene_list_annot_scdata <- mutate_at(gene_list_annot_scdata, 
                                    c("unsorted_immature_beta", 
                                      "unsorted_ins_gfp",
                                      "unsorted_ins_gfp_gcg_sst",
                                      "unsorted_pancreatic_proj",
                                      "sorted_eBCs",
                                      "sorted_immature_beta",
                                      "sorted_ins_gfp_gcg_sst",
                                      "sorted_ins_gfp_sst",
                                      "sorted_pancreatic_proj"), 
                                    ~replace(., is.na(.), 0.000))


View(gene_list_annot_scdata %>% select(hgnc_symbol, unsorted_immature_beta))
# check type of value of 0 vs 0.000 in 
# cause the replaced NAs do not show decimals while the original 0's do
typeof(gene_list_annot_scdata$unsorted_immature_beta)
typeof(gene_list_annot_scdata2$unsorted_immature_beta)




# will continue with analysis & ranks in another R file
# because first need to correct the human RNA-seq data


# save full spreadsheet with added stem cell scRNA-seq
write_tsv(gene_list_annot_scdata, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_stemcell.tsv")



