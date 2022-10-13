# scoring the 422 ligands based on aggregate rankings of all data of interest
# for now just use the ligands data, and later will add scoring 
# based on ligand-receptor interactions

library(tidyverse)
library(pheatmap)


############# load data ###############################


ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_bulkRNA.tsv", 
                          sep = "\t"))


############### prep data for scoring ############################



# replace NAs with 0 in proteomics so score function doesn't break
# we will use pooled sc proteomics to capture all proteins of interest
ligands_df$all_sc_proteomics_0 <- 
  lapply(ligands_df$all_sc_proteomics,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_df <- as.data.frame(lapply(ligands_df, unlist))


############# prep proteomics differential abundance data ###########################


# this returns fold change values only for proteins downregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
ligands_df <- ligands_df %>% 
  # Fold change < 1 indicates higher abundance in human islet  
  mutate(proteomics_down_in_sc = ifelse(stemcell_NDislet_proteomics_FC < 1
                                      , ifelse(stemcell_NDislet_proteomics_adj_p_value < 0.05,
                                               # return fold change for scoring
                                               stemcell_NDislet_proteomics_FC,
                                               NA), 
                                      NA))
# replace NAs with 1 in proteomics_up_in_sc score function doesn't break
ligands_df$proteomics_down_in_sc_0 <- 
  lapply(ligands_df$proteomics_down_in_sc,
         replace_na, replace = 1)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_df <- as.data.frame(lapply(ligands_df, unlist))





################ prep scRNAseq differential expression data #########################


# this returns fold change values only for proteins downregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
# with minimum fold change of log2FC = 0.5 (41% change)
ligands_df <- ligands_df %>% 
  # Fold change < 0.5 indicates higher abundance in human beta cell  
  mutate(scRNA_down_in_eBC = ifelse(eBC_beta_scRNA_log2FC < -0.5
                                  , ifelse(eBC_beta_scRNA_adj_p_value < 0.05,
                                           # return fold change for scoring
                                           eBC_beta_scRNA_log2FC,
                                           NA), 
                                  NA))
# replace NAs with 0 in scRNA_down_in_eBC so score function doesn't break
ligands_df$scRNA_down_in_eBC <- 
  lapply(ligands_df$scRNA_down_in_eBC,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_df <- as.data.frame(lapply(ligands_df, unlist))




################ prep bulk SC-islet RNA-seq melton data #########################

# replace NAs with 0 in sc_islet_bulk_rna so score function doesn't break
ligands_df$sc_islet_bulk_rna <- 
  lapply(ligands_df$sc_islet_bulk_rna,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_df <- as.data.frame(lapply(ligands_df, unlist))




# prep SCislet_Hislet_bulkRNA_Ave_log2FC 

# this returns fold change values only for genes downregulated in SC-islets
# Fold change < 1 indicates lower expression in SC-islet
# and provides a good number of hits
ligands_df <- ligands_df %>% 
  mutate(bulkRNA_down_in_SCislet = ifelse(SCislet_Hislet_bulkRNA_Ave_log2FC < -1,
                                        # return fold change for scoring
                                        SCislet_Hislet_bulkRNA_Ave_log2FC,
                                        NA))
# replace NAs with 0 so score function doesn't break
ligands_df$bulkRNA_down_in_SCislet <- 
  lapply(ligands_df$bulkRNA_down_in_SCislet,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_df <- as.data.frame(lapply(ligands_df, unlist))





############################ scoring function ####################################### 

# percent_rank() assigns a highest score of 1 and lowest score of 0. 
# The benefit of using percent_rank() over min_rank() is that min rank
# would have a different highest score for our receptor & ligands list, 
# because we have 422 ligands and 349 receptors

ligand_scores <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  
  # create scores for sorted sc's cell populations
  mutate(sorted_eBCs_score = percent_rank(desc(sorted_eBCs)),
         sorted_immature_beta_score = percent_rank(desc(sorted_immature_beta)),
         sorted_ins_gfp_gcg_sst_score = percent_rank(desc(sorted_ins_gfp_gcg_sst)),
         sorted_ins_gfp_sst_score = percent_rank(desc(sorted_ins_gfp_sst)),
         sorted_pancreatic_proj_score = percent_rank(desc(sorted_pancreatic_proj))) %>% 
  
  # create scores for all unsorted sc's cell populations
  mutate(unsorted_immature_beta_score = percent_rank(desc(unsorted_immature_beta)),
         unsorted_ins_gfp_score = percent_rank(desc(unsorted_ins_gfp)),
         unsorted_ins_gfp_gcg_sst_score = percent_rank(desc(unsorted_ins_gfp_gcg_sst)),
         unsorted_pancreatic_proj_score = percent_rank(desc(unsorted_pancreatic_proj))) %>% 
  
  # create scores for stem cell proteomics
  mutate(all_sc_proteomics_score = percent_rank(desc(all_sc_proteomics_0))) %>% 
  
  ## create scores for proteomics differential abundance sc pooled vs ND islets
  mutate(proteomics_down_in_sc_score = percent_rank(desc(proteomics_down_in_sc_0))) %>% 
  
  # create scores for scRNA-seq differential expression eBCs vs human beta cells
  mutate(scRNA_down_in_eBC_score = percent_rank(desc(scRNA_down_in_eBC))) %>% 
  
  # create scores for SC-islet bulk RNA-seq melton data
  mutate(sc_islet_bulk_rna_score = percent_rank(desc(sc_islet_bulk_rna))) %>% 
  
  # create scores for SC-islet bulk RNAseq log2FC vs H-islet ave of melton & lund
  mutate(bulkRNA_down_in_SCislet_score = percent_rank(desc(bulkRNA_down_in_SCislet))) %>% 
  
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
                                  scRNA_down_in_eBC_score,
                                  sc_islet_bulk_rna_score,
                                  bulkRNA_down_in_SCislet_score,
                                  all_sc_proteomics_score,
                                  proteomics_down_in_sc_score))) %>% 
  
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
                scRNA_down_in_eBC_score,
                sc_islet_bulk_rna_score,
                bulkRNA_down_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_down_in_sc_score,
                aggregate_score)




############## histogram of aggregate scores ########################


# just aggregate score
ligand_scores %>% 
  ggplot(aes(x = aggregate_score)) +
  geom_histogram(bins = 20,
                 fill = "#36226B") +
  labs(x = "Aggregate Score",
       y = "Count") +
  ggtitle("Ligands") +
  theme_bw()

# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligand_scores/ligand_scores_histo.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)




## histogram of all scores facet wrapped
ligand_scores_longer <- ligand_scores %>%
  dplyr::select(hgnc_symbol,
                sorted_eBCs_score,
                combined_sc_max_score,
                scRNA_down_in_eBC_score,
                sc_islet_bulk_rna_score,
                bulkRNA_down_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_down_in_sc_score,
                aggregate_score) %>% 
  pivot_longer(cols = c(sorted_eBCs_score,
                        combined_sc_max_score,
                        scRNA_down_in_eBC_score,
                        sc_islet_bulk_rna_score,
                        bulkRNA_down_in_SCislet_score,
                        all_sc_proteomics_score,
                        proteomics_down_in_sc_score,
                        aggregate_score)) 

# reorder the name factor to change order of plots in the facet wrap
# and add labels for the facet labels on the plot
ligand_scores_longer$name<- factor(ligand_scores_longer$name, 
                                      levels = c("sorted_eBCs_score",
                                                 "combined_sc_max_score",
                                                 "scRNA_down_in_eBC_score",
                                                 "sc_islet_bulk_rna_score",
                                                 "bulkRNA_down_in_SCislet_score",
                                                 "all_sc_proteomics_score",
                                                 "proteomics_down_in_sc_score",
                                                 "aggregate_score"),
                                      labels = c("SCβ-cell scRNA-seq",
                                                 "Max SC-islet cell scRNA-seq",
                                                 "Down in SCβ-cell scRNA-seq",
                                                 "SC-islet bulk RNA-seq",
                                                 "Down in SC-islet bulk RNA-seq",
                                                 "SC-islet proteomics",
                                                 "Down in SC-islet proteomics",
                                                 "Overall Rank"))

ligand_scores_longer %>% 
  ggplot(aes(x = value)) +
  facet_wrap(~ name,
             scales = "free_y") +
  geom_histogram(bins = 20,
                 fill = "#36226B") +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  labs(x = "Rank",
       y = "Count") +
  theme_bw()


# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligand_scores/ligand_scores_histo_facet_2.JPG",
       device = "jpg",
       width = 2400,
       height = 1500,
       units = "px",
       scale = 0.8)




################# pheatmap with flipped x and y axes ####################333

ligand_scores_flip_rowname <- column_to_rownames(as.data.frame(ligand_scores), 
                                                    var = "hgnc_symbol") %>% 
  dplyr::select(sorted_eBCs_score,
                combined_sc_max_score,
                scRNA_down_in_eBC_score,
                sc_islet_bulk_rna_score,
                bulkRNA_down_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_down_in_sc_score,
                aggregate_score) %>%
  arrange(desc(aggregate_score)) %>% 
  head(50) %>% 
  t()


# heatmap of scores
ligand_scores_flip_rowname %>% 
  pheatmap(labels_row = c("SCβ-cell scRNA-seq",
                          "Max SC-islet cell scRNA-seq",
                          "Down in SCβ-cell scRNA-seq",
                          "SC-islet bulk RNA-seq",
                          "Down in SC-islet bulk RNA-seq",
                          "SC-islet proteomics",
                          "Down in SC-islet proteomics",
                          "Overall Rank"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_col = 9,
           fontsize_row = 10,
           angle_col = 90)


########### save table ######################################

# save new datatables
ligand_scores %>% 
  dplyr::select(hgnc_symbol,
                sorted_eBCs_score,
                combined_sc_max_score,
                scRNA_down_in_eBC_score,
                sc_islet_bulk_rna_score,
                bulkRNA_down_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_down_in_sc_score,
                aggregate_score) %>%
  arrange(desc(aggregate_score)) %>% 
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligand_scores_2.tsv")












