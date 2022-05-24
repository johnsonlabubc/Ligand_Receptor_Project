# add several analysis columns & create figures to explore the combined datasets


# install.packages("ggrepel")
# install.packages("pheatmap")
# install.packages("RColorBrewer")


# load libraries
library(tidyverse)
library(cowplot)
# for adding non-overlapping text labels to points on ggplot scatterplot
library(ggrepel) 
library(pheatmap)

# load receptors list with all prior expression data and annotations
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                      sep = "\t")
  

#### add ranks & analysis columns
gene_list <- gene_list %>% 
  mutate(sorted_eBCs_rank = dense_rank(desc(sorted_eBCs)),
         sorted_eBCs_specificity
         # had to add very small amounts to prevent infinity values
         # the values are so small that adding 0.01 screws up the data
         eBCs_human_beta_ratio = ((sorted_eBCs + 0.0000001) / (Beta + 0.0000001)),
         human_beta_eBCs_ratio = ((Beta + 0.0000001) / (sorted_eBCs + 0.0000001)))


## explore interesting filters and top genes

gene_list %>% 
  filter(keep_in_list == "Yes") %>% 
  select(hgnc_symbol,
         description,
         sorted_eBCs,
         Beta,
         human_beta_eBCs_ratio) %>% 
  arrange(desc(human_beta_eBCs_ratio)) %>% 
  View()

# plot human_beta_eBCs_ratio against huamn beta abundance
gene_list %>% 
  filter(keep_in_list == "Yes") %>% 
  ggplot(aes(x = log2(Beta), y = log2(human_beta_eBCs_ratio), label = hgnc_symbol)) +
  geom_point(colour = "dodgerblue2", size = 2) +
  labs(x = "Log2 beta cell RNA abundance",
       y = "Log2 beta cell / eBCs ratio") +
  # this uses the ggrepel to position
  geom_text_repel(data=subset(gene_list %>% filter(keep_in_list == "Yes"), 
                              log2(human_beta_eBCs_ratio) > 7),
                  size = 2.25,
                  max.overlaps = 15,
                  force = 0.5,
                  max.time = 2,
                  force_pull = 5) +
  theme_cowplot()




########## Heatmap of all data #############################################



genes_list_rowname <- gene_list
genes_list_rowname <- column_to_rownames(as.data.frame(genes_list_rowname), 
                                         var = "hgnc_symbol")

# heatmap of almost everything
genes_list_rowname %>% 
  filter(keep_in_list == "Yes") %>% 
  select(islet_tpm,
         islet_specificity,
         Immune,
         Duct,
         Endothelial,
         Pericytes,
         Acinar,
         INS.SST,
         INS.GCG,
         PPY,
         Delta,
         Beta,
         Alpha,
         beta_specificity,
         unsorted_immature_beta,
         unsorted_ins_gfp,
         unsorted_ins_gfp_gcg_sst,
         unsorted_pancreatic_proj,
         sorted_eBCs,
         sorted_immature_beta,
         sorted_ins_gfp_gcg_sst,
         sorted_ins_gfp_sst,
         sorted_pancreatic_proj, 
         sorted_sc_proteomics,
         unsorted_sc_proteomics,
         ND_islet_proteomics,
         T2D_islet_proteomics) %>%
  pheatmap(cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           clustering_method = "average",
           clustering_distance_cols = "euclidean",
           show_colnames = TRUE,
           show_rownames = FALSE,
           treeheight_col = 20, # changes the height of the col clustering lines
           fontsize = 10)


# heatmap of just beta cell stuff
genes_list_rowname %>% 
  filter(keep_in_list == "Yes") %>% 
  select(islet_tpm,
         islet_specificity,
         Beta,
         beta_specificity,
         sorted_eBCs) %>% 
#         sorted_sc_proteomics,
#         unsorted_sc_proteomics,
#         ND_islet_proteomics,
#         T2D_islet_proteomics) %>%
  pheatmap(cluster_rows = TRUE,
           cluster_cols = TRUE,
           scale = "row",
           clustering_method = "average",
           clustering_distance_cols = "euclidean",
           show_colnames = TRUE,
           show_rownames = FALSE,
           treeheight_col = 20, # changes the height of the col clustering lines
           fontsize = 10)

# heatmap of just ranks
genes_list_rowname %>% 
  filter(keep_in_list == "Yes") %>% 
  select(islet_tpm_rank,
         islet_specificity_rank,
         beta_rank,
         beta_specificity_rank,
         sorted_eBCs_rank, 
         sorted_sc_proteomics_rank,
         unsorted_sc_proteomics_rank,
         ND_islet_proteomics_rank,
         T2D_islet_proteomics_rank) %>%
  pheatmap(cluster_rows = TRUE,
           cluster_cols = FALSE,
           scale = "row",
           clustering_method = "average",
           clustering_distance_cols = "euclidean",
           show_colnames = TRUE,
           show_rownames = FALSE,
           treeheight_row = 0, # changes the height of the col clustering lines
           fontsize_col = 6,
           fontsize_row = 1,
           angle_col = 45)

