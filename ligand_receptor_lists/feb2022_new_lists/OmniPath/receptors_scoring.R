# scoring the 349 receptors based on aggregate rankings of all data of interest
# that we have compiled. Receptor scores will then be combined with ligand
# scores based on known interactions from omnipath data. The combined final 
# ligand scores will be used to decide which ligands to order of in vitro 
# screening.

library(tidyverse)
library(pheatmap)


############# load data ###############################


receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_bulkRNA.tsv", 
                          sep = "\t"))




############### prep proteomics scoring ############################



# replace NAs with 0 in proteomics so score function doesn't break
# we will use pooled sc proteomics to capture all proteins of interest
receptors_df$all_sc_proteomics_0 <- 
  lapply(receptors_df$all_sc_proteomics,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))



############# prep proteomics differential abundance data ###########################


# this returns fold change values only for proteins upregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
receptors_df <- receptors_df %>% 
# Fold change > 1 indicates higher abundance in stem cell  
mutate(proteomics_up_in_sc = ifelse(stemcell_NDislet_proteomics_FC > 1
                                    , ifelse(stemcell_NDislet_proteomics_adj_p_value < 0.05,
                                             # return fold change for scoring
                                             stemcell_NDislet_proteomics_FC,
                                             NA), 
                                    NA))
# replace NAs with 1 in proteomics_up_in_sc score function doesn't break
receptors_df$proteomics_up_in_sc_0 <- 
  lapply(receptors_df$proteomics_up_in_sc,
         replace_na, replace = 1)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))



################ prep scRNAseq differential expression data #########################
       

# this returns fold change values only for proteins upregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
# with minimum fold change of log2FC = 0.5 (41% change)
receptors_df <- receptors_df %>% 
  # Fold change > 1 indicates higher abundance in stem cell  
  mutate(scRNA_up_in_eBC = ifelse(eBC_beta_scRNA_log2FC > 0.5
                                      , ifelse(eBC_beta_scRNA_adj_p_value < 0.05,
                                               # return fold change for scoring
                                               eBC_beta_scRNA_log2FC,
                                               NA), 
                                      NA))
# replace NAs with 0 in scRNA_up_in_eBC so score function doesn't break
receptors_df$scRNA_up_in_eBC <- 
  lapply(receptors_df$scRNA_up_in_eBC,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))




################ prep bulk SC-islet RNA-seq melton data #########################

# replace NAs with 0 in sc_islet_bulk_rna so score function doesn't break
receptors_df$sc_islet_bulk_rna <- 
  lapply(receptors_df$sc_islet_bulk_rna,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))



# prep SCislet_Hislet_bulkRNA_Ave_log2FC 

# this returns fold change values only for genes upregulated in SC-islets
# Fold change > 1 indicates higher expression in SC-islet
# and provides a good number of hits
receptors_df <- receptors_df %>% 
  mutate(bulkRNA_up_in_SCislet = ifelse(SCislet_Hislet_bulkRNA_Ave_log2FC > 1,
                                           # return fold change for scoring
                                           SCislet_Hislet_bulkRNA_Ave_log2FC,
                                           NA))
# replace NAs with 0 so score function doesn't break
receptors_df$bulkRNA_up_in_SCislet <- 
  lapply(receptors_df$bulkRNA_up_in_SCislet,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))




############################ scoring function ####################################### 

# percent_rank() assigns a highest score of 1 and lowest score of 0. 
# The benefit of using percent_rank() over min_rank() is that min rank
# would have a different highest score for our receptor & ligands list, 
# because we have 422 ligands and 349 receptors

receptors_scores <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  
  # create scores for sorted sc's cell populations
  mutate(sorted_eBCs_score = percent_rank(sorted_eBCs),
         sorted_immature_beta_score = percent_rank(sorted_immature_beta),
         sorted_ins_gfp_gcg_sst_score = percent_rank(sorted_ins_gfp_gcg_sst),
         sorted_ins_gfp_sst_score = percent_rank(sorted_ins_gfp_sst),
         sorted_pancreatic_proj_score = percent_rank(sorted_pancreatic_proj)) %>% 
  
  # create scores for all unsorted sc's cell populations
  mutate(unsorted_immature_beta_score = percent_rank(unsorted_immature_beta),
         unsorted_ins_gfp_score = percent_rank(unsorted_ins_gfp),
         unsorted_ins_gfp_gcg_sst_score = percent_rank(unsorted_ins_gfp_gcg_sst),
         unsorted_pancreatic_proj_score = percent_rank(unsorted_pancreatic_proj)) %>% 
  
  # create scores for stem cell proteomics
  mutate(all_sc_proteomics_score = percent_rank(all_sc_proteomics_0)) %>% 
  
  ## create scores for proteomics differential abundance sc pooled vs ND islets
  mutate(proteomics_up_in_sc_score = percent_rank(proteomics_up_in_sc_0)) %>% 
  
  # create scores for scRNA-seq differential expression eBCs vs human beta cells
  mutate(scRNA_up_in_eBC_score = percent_rank(scRNA_up_in_eBC)) %>% 
  
  # create scores for SC-islet bulk RNA-seq melton data
  mutate(sc_islet_bulk_rna_score = percent_rank(sc_islet_bulk_rna)) %>% 
  
  # create scores for SC-islet bulk RNAseq log2FC vs H-islet ave of melton & lund
  mutate(bulkRNA_up_in_SCislet_score = percent_rank(bulkRNA_up_in_SCislet)) %>% 

  # create aggregate sorted scores
  rowwise() %>% 
  mutate(sorted_sc_mean_score = mean(c(sorted_eBCs_score,
                                       sorted_immature_beta_score,
                                       sorted_ins_gfp_gcg_sst_score,
                                       sorted_ins_gfp_sst_score,
                                       sorted_pancreatic_proj_score)),
         sorted_sc_max_score = max(c(sorted_eBCs_score,
                                     sorted_immature_beta_score,
                                     sorted_ins_gfp_gcg_sst_score,
                                     sorted_ins_gfp_sst_score,
                                     sorted_pancreatic_proj_score)),
         # create aggregate unsorted scores
         unsorted_sc_mean_score = mean(c(unsorted_immature_beta_score,
                                         unsorted_ins_gfp_score,
                                         unsorted_ins_gfp_gcg_sst_score,
                                         unsorted_pancreatic_proj_score)),
         unsorted_sc_max_score = max(c(unsorted_immature_beta_score,
                                       unsorted_ins_gfp_score,
                                       unsorted_ins_gfp_gcg_sst_score,
                                       unsorted_pancreatic_proj_score)),
         # create aggregate of sorted & unsorted scores
         combined_sc_max_score = max(c(sorted_sc_max_score,
                                       unsorted_sc_max_score)),
         
         # also create final aggregate score based on all data of interest
         aggregate_score = mean(c(sorted_eBCs_score,
                               combined_sc_max_score,
                               scRNA_up_in_eBC_score,
                               sc_islet_bulk_rna_score,
                               bulkRNA_up_in_SCislet_score,
                               all_sc_proteomics_score,
                               proteomics_up_in_sc_score))) %>% 

  ungroup() %>% 
  
  # rerank combined_sc_max_score & aggregate_score to rescale the data
  mutate(combined_sc_max_score = percent_rank(combined_sc_max_score),
         aggregate_score= percent_rank(aggregate_score)) %>% 
  
  dplyr::select(hgnc_symbol,
         sorted_eBCs_score,
         sorted_sc_mean_score,
         unsorted_sc_mean_score,
         sorted_sc_max_score,
         unsorted_sc_max_score,
         combined_sc_max_score,
         scRNA_up_in_eBC_score,
         sc_islet_bulk_rna_score,
         bulkRNA_up_in_SCislet_score,
         all_sc_proteomics_score,
         proteomics_up_in_sc_score,
         aggregate_score)
  




# compare max scores of sorted vs unsorted stem cell scRNA-seq
receptors_scores %>% 
  ggplot(aes(x = unsorted_sc_max_score, 
             y = sorted_sc_max_score)) +
  geom_point()
# comparison of unsorted_sc_max_score & sorted_sc_max_score shows they are 
# very strognly correlated, with a few outliers. Therefor it makes sense to 
# just combine them into 1 score for the purpose of the aggregate ranking. 
# In other words, take the max score across all 9 cell types between the 
# sorted and unsorted scRNA-seq



############## histogram of aggregate scores ########################


# just aggregate score
receptors_scores %>% 
  ggplot(aes(x = aggregate_score)) +
  geom_histogram(bins = 20,
                 fill = "#1F9274") +
  labs(x = "Aggregate Score",
       y = "Count") +
  ggtitle("Receptors") +
  theme_bw()

# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptor_scores_histo.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)




## histogram of all scores facet wrapped
receptors_scores_longer <- receptors_scores %>%
  dplyr::select(hgnc_symbol,
         sorted_eBCs_score,
         combined_sc_max_score,
         scRNA_up_in_eBC_score,
         sc_islet_bulk_rna_score,
         bulkRNA_up_in_SCislet_score,
         all_sc_proteomics_score,
         proteomics_up_in_sc_score,
         aggregate_score) %>% 
  pivot_longer(cols = c(sorted_eBCs_score,
               combined_sc_max_score,
               scRNA_up_in_eBC_score,
               sc_islet_bulk_rna_score,
               bulkRNA_up_in_SCislet_score,
               all_sc_proteomics_score,
               proteomics_up_in_sc_score,
               aggregate_score)) 
    
# reorder the name factor to change order of plots in the facet wrap
# and add labels for the facet labels on the plot
receptors_scores_longer$name<- factor(receptors_scores_longer$name, 
                                      levels = c("sorted_eBCs_score",
                                                 "combined_sc_max_score",
                                                 "scRNA_up_in_eBC_score",
                                                 "sc_islet_bulk_rna_score",
                                                 "bulkRNA_up_in_SCislet_score",
                                                 "all_sc_proteomics_score",
                                                 "proteomics_up_in_sc_score",
                                                 "aggregate_score"),
                                      labels = c("SCβ-cell scRNA-seq",
                                                 "Max SC-islet cell scRNA-seq",
                                                 "Up in SCβ-cell scRNA-seq",
                                                 "SC-islet bulk RNA-seq",
                                                 "Up in SC-islet bulk RNA-seq",
                                                 "SC-islet proteomics",
                                                 "Up in SC-islet proteomics",
                                                 "Overall Rank"))

receptors_scores_longer %>% 
  ggplot(aes(x = value)) +
  facet_wrap(~ name,
             scales = "free_y") +
  geom_histogram(bins = 20,
                 fill = "#1F9274") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  labs(x = "Rank",
       y = "Count") +
  theme_bw()


# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptor_scores/receptor_scores_histo_facet_4.JPG",
       device = "jpg",
       width = 2400,
       height = 1500,
       units = "px",
       scale = 0.8)




#################### heatmap using pheatmap #####################################


receptors_scores_rowname <- receptors_scores
receptors_scores_rowname <- column_to_rownames(as.data.frame(receptors_scores), 
                                         var = "hgnc_symbol")

# column legend table
# col_annotations <- data.frame(col_name = c("sorted_eBCs_score",
#                                            "combined_sc_max_score",
#                                            "scRNA_up_in_eBC_score",
#                                            "all_sc_proteomics_score",
#                                            "proteomics_up_in_sc_score",
#                                            "aggregate_score"),
#                               display_name = c("SCβ-cell scRNA",
#                                                "Max SC-islet scRNA",
#                                                "Up in SCβ-cell scRNA",
#                                                "SC-islet protein",
#                                                "Up in SC-islet protein",
#                                                "Final aggregate"))
# col_annotations_rowname <- col_annotations
# col_annotations_rowname <- column_to_rownames(as.data.frame(col_annotations), 
#                                                var = "col_name")


# heatmap of just ranks
receptors_scores_rowname %>% 
  dplyr::select(sorted_eBCs_score,
         combined_sc_max_score,
         scRNA_up_in_eBC_score,
         sc_islet_bulk_rna_score,
         bulkRNA_up_in_SCislet_score,
         all_sc_proteomics_score,
         proteomics_up_in_sc_score,
         aggregate_score) %>%
  arrange(desc(aggregate_score)) %>% 
  head(50) %>% 
  pheatmap(labels_col = c("SCβ-cell scRNA",
                          "Max SC-islet scRNA",
                          "Up in SCβ-cell scRNA",
                          "SC-islet bulk RNA",
                          "Up in SC-islet bulk RNA",
                          "SC-islet protein",
                          "Up in SC-islet protein",
                          "Final aggregate"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           scale = "column",
           show_colnames = TRUE,
           fontsize_col = 10,
           fontsize_row = 6,
           angle_col = 90)




################# pheatmap with flipped x and y axes ####################333

receptors_scores_flip_rowname <- column_to_rownames(as.data.frame(receptors_scores), 
                                                    var = "hgnc_symbol") %>% 
  dplyr::select(sorted_eBCs_score,
         combined_sc_max_score,
         scRNA_up_in_eBC_score,
         sc_islet_bulk_rna_score,
         bulkRNA_up_in_SCislet_score,
         all_sc_proteomics_score,
         proteomics_up_in_sc_score,
         aggregate_score) %>%
  arrange(desc(aggregate_score)) %>% 
  head(50) %>% 
  t()



# row labels table
# row_labels <- data.frame(col_name = c("sorted_eBCs_score",
#                                            "combined_sc_max_score",
#                                            "scRNA_up_in_eBC_score",
#                                            "all_sc_proteomics_score",
#                                            "proteomics_up_in_sc_score",
#                                            "aggregate_score"),
#                               display_name = c("SCβ-cell scRNA",
#                                                "Max SC-islet scRNA",
#                                                "Up in SCβ-cell scRNA",
#                                                "SC-islet protein",
#                                                "Up in SC-islet protein",
#                                                "Final aggregate"))
# row_labels_rowname <- row_labels
# row_labels_rowname <- column_to_rownames(as.data.frame(row_labels), 
#                                               var = "col_name")


# heatmap of scores
receptors_scores_flip_rowname %>% 
  pheatmap(labels_row = c("SCβ-cell scRNA-seq",
                         "Max SC-islet cell scRNA-seq",
                         "Up in SCβ-cell scRNA-seq",
                         "SC-islet bulk RNA-seq",
                         "Up in SC-islet bulk RNA-seq",
                         "SC-islet proteomics",
                         "Up in SC-islet proteomics",
                         "Overall rank"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_col = 9,
           fontsize_row = 10,
           angle_col = 90)


########### save table ######################################

# save new datatables
receptors_scores %>% 
  dplyr::select(hgnc_symbol,
                sorted_eBCs_score,
                combined_sc_max_score,
                scRNA_up_in_eBC_score,
                sc_islet_bulk_rna_score,
                bulkRNA_up_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_up_in_sc_score,
                aggregate_score) %>%
  arrange(desc(aggregate_score)) %>% 
write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_scores_2.tsv")
