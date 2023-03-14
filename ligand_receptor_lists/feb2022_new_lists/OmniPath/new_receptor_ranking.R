# new receptor ranking 
# exploring 2 sided ranks for DE columns
# and ranking in a way that works for all 3 separate final rankings



library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)


############# load data ###############################


receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_meltonbalboa.tsv", 
                          sep = "\t")) %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) 




############### prep proteomics scoring ############################

# normalize to AA length
# load AA length uniprot data
# use raw list b/c it has more proteins than howards filtered list
uniprot_size <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/uniprot_size.txt",
                           sep = "\t", header = T) %>% 
  select(hgnc_symbol = Genes,
         AA_length = Length)


# compute abundance normalized to length
receptors_AA <- receptors_df %>% 
  left_join(uniprot_size) %>% 
  mutate(all_sc_proteomics_AAnorm = all_sc_proteomics/AA_length) %>% 
  # some genes have multiple AA size, so lets take the average
  select(hgnc_symbol,
         all_sc_proteomics_AAnorm) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(all_sc_proteomics_AAnorm = mean(all_sc_proteomics_AAnorm))


# add back to main table
receptors_df <- receptors_df %>% 
  left_join(receptors_AA,
            by = "hgnc_symbol")



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



# receptors_df %>% 
#   select(hgnc_symbol,
#          stemcell_NDislet_proteomics_FC,
#          stemcell_NDislet_proteomics_adj_p_value,
#          proteomics_up_in_sc_0) %>% 
#   mutate(rank = percent_rank(proteomics_up_in_sc_0)) %>% 
#   View()



# this returns fold change values only for proteins downregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
receptors_df <- receptors_df %>% 
  # Fold change < 1 indicates higher abundance in human islet  
  mutate(proteomics_down_in_sc = ifelse(stemcell_NDislet_proteomics_FC < 1
                                        , ifelse(stemcell_NDislet_proteomics_adj_p_value < 0.05,
                                                 # return fold change for scoring
                                                 stemcell_NDislet_proteomics_FC,
                                                 NA), 
                                        NA))
# replace NAs with 1 in proteomics_up_in_sc score function doesn't break
receptors_df$proteomics_down_in_sc_0 <- 
  lapply(receptors_df$proteomics_down_in_sc,
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



# this returns fold change values only for proteins downregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
# with minimum fold change of log2FC = 0.5 (41% change)
receptors_df <- receptors_df %>% 
  # Fold change < 0.5 indicates higher abundance in human beta cell  
  mutate(scRNA_down_in_eBC = ifelse(eBC_beta_scRNA_log2FC < -0.5
                                    , ifelse(eBC_beta_scRNA_adj_p_value < 0.05,
                                             # return fold change for scoring
                                             eBC_beta_scRNA_log2FC,
                                             NA), 
                                    NA))
# replace NAs with 0 in scRNA_down_in_eBC so score function doesn't break
receptors_df$scRNA_down_in_eBC <- 
  lapply(receptors_df$scRNA_down_in_eBC,
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



# prep SCislet_Hislet_bulkRNA_Melton_log2FC 

# this returns fold change values only for genes upregulated in SC-islets
# Fold change > 1 indicates higher expression in SC-islet
# and provides a good number of hits
# use comparison to just melton's own h-islet, not to lund
# b/c the comparison to lund produced discrepant data
receptors_df <- receptors_df %>% 
  mutate(bulkRNA_up_in_SCislet = ifelse(SCislet_Hislet_bulkRNA_Melton_log2FC > 1,
                                        # return fold change for scoring
                                        SCislet_Hislet_bulkRNA_Melton_log2FC,
                                        NA))
# replace NAs with 0 so score function doesn't break
receptors_df$bulkRNA_up_in_SCislet <- 
  lapply(receptors_df$bulkRNA_up_in_SCislet,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))


# this returns fold change values only for genes downregulated in SC-islets
# Fold change < 1 indicates lower expression in SC-islet
# and provides a good number of hits
# use comparison to just melton's own h-islet, not to lund
# b/c the comparison to lund produced discrepant data
receptors_df <- receptors_df %>% 
  mutate(bulkRNA_down_in_SCislet = ifelse(SCislet_Hislet_bulkRNA_Melton_log2FC < -1,
                                          # return fold change for scoring
                                          SCislet_Hislet_bulkRNA_Melton_log2FC,
                                          NA))
# replace NAs with 0 so score function doesn't break
receptors_df$bulkRNA_down_in_SCislet <- 
  lapply(receptors_df$bulkRNA_down_in_SCislet,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))



########## prep melton 2019 scRNA expression & DE data ################

## prep expression data
# replace NAs with 0  so score function doesn't break
receptors_df$sc_beta_scRNA_Melton2019 <- 
  lapply(receptors_df$sc_beta_scRNA_Melton2019,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))


## prep DE data
# this returns fold change values only for genes upregulated in stemcell
# so will score based on fold change
receptors_df <- receptors_df %>% 
  # Fold change > 1 indicates higher abundance in stem cell  
  mutate(scRNA_up_in_SCbeta_melton2019 = ifelse(SCbeta_Hbeta_scRNA_Melton2019_log2FC > 1,
                                                # return fold change for scoring
                                                SCbeta_Hbeta_scRNA_Melton2019_log2FC,
                                                NA))
# replace NAs with 0 in scRNA_up_in_eBC so score function doesn't break
receptors_df$scRNA_up_in_SCbeta_melton2019 <- 
  lapply(receptors_df$scRNA_up_in_SCbeta_melton2019,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))



## prep DE data
# this returns fold change values only for genes dowregulated in stemcell
# so will score based on fold change
receptors_df <- receptors_df %>% 
  # Fold change < 1 indicates higher abundance in human islet  
  mutate(scRNA_down_in_SCbeta_melton2019 = ifelse(SCbeta_Hbeta_scRNA_Melton2019_log2FC < -1,
                                                  # return fold change for scoring
                                                  SCbeta_Hbeta_scRNA_Melton2019_log2FC,
                                                  NA))
# replace NAs with 0 in scRNA_up_in_eBC so score function doesn't break
receptors_df$scRNA_down_in_SCbeta_melton2019 <- 
  lapply(receptors_df$scRNA_down_in_SCbeta_melton2019,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))




################### prep balboa 2022 scRNA DE data ##########################



# this returns fold change values only for proteins upregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
# with minimum fold change of log2FC = 0.5 (41% change)
receptors_df <- receptors_df %>% 
  # Fold change > 1 indicates higher abundance in stem cell  
  mutate(scRNA_up_in_SCbeta_balboa2022 = ifelse(SCbeta_Hbeta_scRNA_Balboa2022_log2FC > 0.5
                                                , ifelse(SCbeta_Hbeta_scRNA_Balboa2022_adj_p_value < 0.05,
                                                         # return fold change for scoring
                                                         SCbeta_Hbeta_scRNA_Balboa2022_log2FC,
                                                         NA), 
                                                NA))
# replace NAs with 0 in scRNA_up_in_eBC so score function doesn't break
receptors_df$scRNA_up_in_SCbeta_balboa2022 <- 
  lapply(receptors_df$scRNA_up_in_SCbeta_balboa2022,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))




# this returns fold change values only for proteins downregulated in stemcell
# with a p value of less than 0.05
# so will score based on fold change
# with minimum fold change of log2FC = 0.5 (41% change)
receptors_df <- receptors_df %>% 
  # Fold change < 0.5 indicates higher abundance in human beta cell  
  mutate(scRNA_down_in_SCbeta_balboa2022 = ifelse(SCbeta_Hbeta_scRNA_Balboa2022_log2FC < -0.5
                                                  , ifelse(SCbeta_Hbeta_scRNA_Balboa2022_adj_p_value < 0.05,
                                                           # return fold change for scoring
                                                           SCbeta_Hbeta_scRNA_Balboa2022_log2FC,
                                                           NA), 
                                                  NA))
# replace NAs with 0 in scRNA_up_in_eBC so score function doesn't break
receptors_df$scRNA_down_in_SCbeta_balboa2022 <- 
  lapply(receptors_df$scRNA_down_in_SCbeta_balboa2022,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))




# receptors_ranks <- receptors_df %>%
#   transmute(hgnc_symbol = hgnc_symbol,
#             all_sc_proteomics_0,
#             all_sc_proteomics_score = percent_rank(all_sc_proteomics_0),
#             all_sc_proteomics_score_rev = percent_rank(desc(all_sc_proteomics_0)),
#             all_sc_proteomics_score_1_minus = (1 - percent_rank(all_sc_proteomics_0)),
#             all_sc_proteomics_score_1_minus_rerank = percent_rank(1 - percent_rank(all_sc_proteomics_0)))
# 




############################ abundant rank function ####################################### 

# percent_rank() assigns a highest score of 1 and lowest score of 0. 
# The benefit of using percent_rank() over min_rank() is that min rank
# would have a different highest score for our receptor & ligands list, 
# because we have 422 ligands and 349 receptors

receptors_scores_abundant <- receptors_df %>% 
  
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
  
  # create score for SCbeta-cell melton2019 expression
  mutate(SCbeta_scRNA_melton2019_score = percent_rank(sc_beta_scRNA_Melton2019)) %>% 
  
  # create score for SCbeta-cell differential expression in melton 2019
  mutate(scRNA_up_in_SCbeta_melton2019_score = percent_rank(scRNA_up_in_SCbeta_melton2019)) %>% 
  
  # create score for SCbeta-cell differential expression in balboa 2022
  mutate(scRNA_up_in_SCbeta_balboa2022_score = percent_rank(scRNA_up_in_SCbeta_balboa2022)) %>% 
  
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
                                  SCbeta_scRNA_melton2019_score,
                                  scRNA_up_in_SCbeta_melton2019_score,
                                  scRNA_up_in_SCbeta_balboa2022_score,
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
                combined_sc_max_score,
                scRNA_up_in_eBC_score,
                SCbeta_scRNA_melton2019_score,
                scRNA_up_in_SCbeta_melton2019_score,
                scRNA_up_in_SCbeta_balboa2022_score,
                sc_islet_bulk_rna_score,
                bulkRNA_up_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_up_in_sc_score,
                aggregate_score)





############################ scarce rank function ####################################### 

# percent_rank() assigns a highest score of 1 and lowest score of 0. 
# The benefit of using percent_rank() over min_rank() is that min rank
# would have a different highest score for our receptor & ligands list, 
# because we have 422 ligands and 349 receptors

receptor_scores_scarce <- receptors_df %>% 
  
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
  
  # create scores for SC-islet bulk RNAseq log2FC vs H-islet of melton
  mutate(bulkRNA_down_in_SCislet_score = percent_rank(desc(bulkRNA_down_in_SCislet))) %>% 
  
  # create score for SCbeta-cell melton2019 expression
  mutate(SCbeta_scRNA_melton2019_score = percent_rank(desc(sc_beta_scRNA_Melton2019))) %>% 
  
  # create score for SCbeta-cell differential expression in melton 2019
  mutate(scRNA_down_in_SCbeta_melton2019_score = percent_rank(desc(scRNA_down_in_SCbeta_melton2019))) %>% 
  
  # create score for SCbeta-cell differential expression in balboa 2022
  mutate(scRNA_down_in_SCbeta_balboa2022_score = percent_rank(desc(scRNA_down_in_SCbeta_balboa2022))) %>% 
  
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
                                  SCbeta_scRNA_melton2019_score,
                                  scRNA_down_in_SCbeta_melton2019_score,
                                  scRNA_down_in_SCbeta_balboa2022_score,
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
                combined_sc_max_score,
                scRNA_down_in_eBC_score,
                SCbeta_scRNA_melton2019_score,
                scRNA_down_in_SCbeta_melton2019_score,
                scRNA_down_in_SCbeta_balboa2022_score,
                sc_islet_bulk_rna_score,
                bulkRNA_down_in_SCislet_score,
                all_sc_proteomics_score,
                proteomics_down_in_sc_score,
                aggregate_score)



################### save the ranks df's #####################################






#################### combine scarce and abundant df's for plotting #################

receptor_scores_plot <- left_join(receptors_scores_abundant,
                                  receptor_scores_scarce %>% 
                                    select(hgnc_symbol,
                                           scRNA_down_in_eBC_score,
                                           scRNA_down_in_SCbeta_melton2019_score,
                                           scRNA_down_in_SCbeta_balboa2022_score,
                                           bulkRNA_down_in_SCislet_score,
                                           proteomics_down_in_sc_score),
                                  by = "hgnc_symbol")



#################### heatmap using pheatmap #####################################


receptors_scores_rowname <- receptor_scores_plot
receptors_scores_rowname <- column_to_rownames(as.data.frame(receptor_scores_plot), 
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


# load column annotations table
col_annotations_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/heatmap_col_annotations.csv", 
                          sep = ","))  %>% 
  select(1:5)

col_annotations_rowname <- column_to_rownames(as.data.frame(col_annotations_df), 
                                              var = "column_names")

# heatmap of just ranks
receptors_scores_rowname %>% 
  dplyr::select(-aggregate_score,
                # reorder columns for plotting
                sorted_eBCs_score,
                combined_sc_max_score,
                SCbeta_scRNA_melton2019_score,
                sc_islet_bulk_rna_score,
                all_sc_proteomics_score,
                scRNA_up_in_eBC_score,
                scRNA_up_in_SCbeta_melton2019_score,
                scRNA_up_in_SCbeta_balboa2022_score,
                bulkRNA_up_in_SCislet_score,
                proteomics_up_in_sc_score,
                scRNA_down_in_eBC_score,
                scRNA_down_in_SCbeta_melton2019_score,
                scRNA_down_in_SCbeta_balboa2022_score,
                bulkRNA_down_in_SCislet_score,
                proteomics_down_in_sc_score) %>%
  relocate(sorted_eBCs_score,
          combined_sc_max_score,
          SCbeta_scRNA_melton2019_score,
          sc_islet_bulk_rna_score,
          all_sc_proteomics_score,
          scRNA_up_in_eBC_score,
          scRNA_up_in_SCbeta_melton2019_score,
          scRNA_up_in_SCbeta_balboa2022_score,
          bulkRNA_up_in_SCislet_score,
          proteomics_up_in_sc_score,
          scRNA_down_in_eBC_score,
          scRNA_down_in_SCbeta_melton2019_score,
          scRNA_down_in_SCbeta_balboa2022_score,
          bulkRNA_down_in_SCislet_score,
          proteomics_down_in_sc_score) %>% 
  # arrange(desc(aggregate_score)) %>% 
  # head(50) %>% 
  pheatmap(
    annotation_col = col_annotations_rowname,
    gaps_col = c(5,10),
    color = colorRampPalette(RColorBrewer::brewer.pal(7, "BuPu"))(50),
    border_color = "black",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize_col = 9,
    fontsize_row = 10,
    angle_col = 45) 



# labels_col = c("SCβ-cell scRNA-seq (Lynn)",
#                "SCβ-cell scRNA-seq (Melton)",
#                "Max SC-islet cell scRNA-seq (Lynn)",
#                "SC-islet bulk RNA-seq (Melton)",
#                "SC-islet proteomics (Johnson)",
#                "Up in SCβ-cell scRNA-seq (Lynn)",
#                "Up in SCβ-cell scRNA-seq (Melton)",
#                "Up in SCβ-cell scRNA-seq (Balboa)",
#                "Up in SC-islet bulk RNA-seq (Melton)",
#                "Up in SC-islet proteomics (Johnson)",
#                "Down in SCβ-cell scRNA-seq (Lynn)",
#                "Down in SCβ-cell scRNA-seq (Melton)",
#                "Down in SCβ-cell scRNA-seq (Balboa)",
#                "Down in SC-islet bulk RNA-seq (Melton)",
#                "Down in SC-islet proteomics (Johnson)"),




################# heatmap with complexheatmap ################################


# heatmap of just ranks
matrix1 <- receptors_scores_rowname %>% 
  dplyr::select(-aggregate_score,
                # reorder columns for plotting
                sorted_eBCs_score,
                combined_sc_max_score,
                SCbeta_scRNA_melton2019_score,
                sc_islet_bulk_rna_score,
                all_sc_proteomics_score,
                scRNA_up_in_eBC_score,
                scRNA_up_in_SCbeta_melton2019_score,
                scRNA_up_in_SCbeta_balboa2022_score,
                bulkRNA_up_in_SCislet_score,
                proteomics_up_in_sc_score,
                scRNA_down_in_eBC_score,
                scRNA_down_in_SCbeta_melton2019_score,
                scRNA_down_in_SCbeta_balboa2022_score,
                bulkRNA_down_in_SCislet_score,
                proteomics_down_in_sc_score) %>%
  relocate(sorted_eBCs_score,
           combined_sc_max_score,
           SCbeta_scRNA_melton2019_score,
           sc_islet_bulk_rna_score,
           all_sc_proteomics_score,
           scRNA_up_in_eBC_score,
           scRNA_up_in_SCbeta_melton2019_score,
           scRNA_up_in_SCbeta_balboa2022_score,
           bulkRNA_up_in_SCislet_score,
           proteomics_up_in_sc_score,
           scRNA_down_in_eBC_score,
           scRNA_down_in_SCbeta_melton2019_score,
           scRNA_down_in_SCbeta_balboa2022_score,
           bulkRNA_down_in_SCislet_score,
           proteomics_down_in_sc_score) %>% 
  data.matrix(rownames.force = TRUE)

# create column annotation object
col_annot_obj <- HeatmapAnnotation("Research Group" = col_annotations_df$Research.Group,
                                   "Method" = col_annotations_df$Technology,
                                   "Cell Type" = col_annotations_df$Cell.Group,
                                   show_annotation_name = FALSE,
                                   border = FALSE,
                                   col = list("Research Group" = c("Lynn/Johnson" = "#440D54",
                                                                   "Melton" = "#1F968B",
                                                                   "Otonkoski" = "#FDE725"),
                                              "Method" = c("bulk RNA-seq" = "gray40",
                                                           "proteomics" = "black",
                                                           "scRNA-seq"= "lightgray"),
                                              "Cell Type" = c("SC-islet" = "lightblue",
                                                              "SC?-cell" = "darkblue")),
                                   annotation_legend_param = list("Cell Type" = list(labels = c("SC-islets",
                                                                                                expression(paste("SC", beta, "-cells"))))))

heatmap_rec_all <- matrix1 %>% 
  Heatmap(cluster_rows = TRUE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Rank"),
          column_split = factor(col_annotations_df$Measurement,
                                levels = c("Expression", "Upregulation",
                                           "Downregulation"),
                                labels = c("Expression",
                                           "Upregulated",
                                           "Downregulated")),
          top_annotation = col_annot_obj,
          border = TRUE,
          height = unit(9, "cm"),
          row_title = "All receptors")
heatmap_rec_all



# heatmap of top 25 most abundant

# heatmap of just ranks
matrix2 <- receptors_scores_rowname %>% 
  arrange(desc(aggregate_score)) %>% 
  head(25) %>% 
  dplyr::select(-aggregate_score,
                # reorder columns for plotting
                sorted_eBCs_score,
                combined_sc_max_score,
                SCbeta_scRNA_melton2019_score,
                sc_islet_bulk_rna_score,
                all_sc_proteomics_score,
                scRNA_up_in_eBC_score,
                scRNA_up_in_SCbeta_melton2019_score,
                scRNA_up_in_SCbeta_balboa2022_score,
                bulkRNA_up_in_SCislet_score,
                proteomics_up_in_sc_score,
                scRNA_down_in_eBC_score,
                scRNA_down_in_SCbeta_melton2019_score,
                scRNA_down_in_SCbeta_balboa2022_score,
                bulkRNA_down_in_SCislet_score,
                proteomics_down_in_sc_score) %>%
  relocate(sorted_eBCs_score,
           combined_sc_max_score,
           SCbeta_scRNA_melton2019_score,
           sc_islet_bulk_rna_score,
           all_sc_proteomics_score,
           scRNA_up_in_eBC_score,
           scRNA_up_in_SCbeta_melton2019_score,
           scRNA_up_in_SCbeta_balboa2022_score,
           bulkRNA_up_in_SCislet_score,
           proteomics_up_in_sc_score,
           scRNA_down_in_eBC_score,
           scRNA_down_in_SCbeta_melton2019_score,
           scRNA_down_in_SCbeta_balboa2022_score,
           bulkRNA_down_in_SCislet_score,
           proteomics_down_in_sc_score) %>% 
  data.matrix(rownames.force = TRUE)

heatmap_rec_abund <- matrix2 %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Rank"),
          column_split = factor(col_annotations_df$Measurement,
                                levels = c("Expression", "Upregulation",
                                           "Downregulation")),
     #     top_annotation = col_annot_obj,
          border = TRUE,
          height = unit(9, "cm"),
          row_title = "Top 25 abundant receptors",
          show_heatmap_legend = FALSE)
heatmap_rec_abund



# plot top scarce receptors


# filter scarce df to just top 25
receptor_scores_scarce25 <- receptor_scores_scarce %>% 
  arrange(desc(aggregate_score)) %>% 
  head(25)

receptor_scores_plot_scarce25 <- receptor_scores_plot %>% 
  filter(hgnc_symbol %in% receptor_scores_scarce25$hgnc_symbol)

receptors_scores_scare25_rowname <- column_to_rownames(as.data.frame(receptor_scores_plot_scarce25), 
                                               var = "hgnc_symbol")
                          

matrix3 <- receptors_scores_scare25_rowname %>% 
  dplyr::select(-aggregate_score,
                # reorder columns for plotting
                sorted_eBCs_score,
                combined_sc_max_score,
                SCbeta_scRNA_melton2019_score,
                sc_islet_bulk_rna_score,
                all_sc_proteomics_score,
                scRNA_up_in_eBC_score,
                scRNA_up_in_SCbeta_melton2019_score,
                scRNA_up_in_SCbeta_balboa2022_score,
                bulkRNA_up_in_SCislet_score,
                proteomics_up_in_sc_score,
                scRNA_down_in_eBC_score,
                scRNA_down_in_SCbeta_melton2019_score,
                scRNA_down_in_SCbeta_balboa2022_score,
                bulkRNA_down_in_SCislet_score,
                proteomics_down_in_sc_score) %>%
  relocate(sorted_eBCs_score,
           combined_sc_max_score,
           SCbeta_scRNA_melton2019_score,
           sc_islet_bulk_rna_score,
           all_sc_proteomics_score,
           scRNA_up_in_eBC_score,
           scRNA_up_in_SCbeta_melton2019_score,
           scRNA_up_in_SCbeta_balboa2022_score,
           bulkRNA_up_in_SCislet_score,
           proteomics_up_in_sc_score,
           scRNA_down_in_eBC_score,
           scRNA_down_in_SCbeta_melton2019_score,
           scRNA_down_in_SCbeta_balboa2022_score,
           bulkRNA_down_in_SCislet_score,
           proteomics_down_in_sc_score) %>% 
  data.matrix(rownames.force = TRUE)

heatmap_rec_scarce <- matrix3 %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Rank"),
          column_split = factor(col_annotations_df$Measurement,
                                levels = c("Expression", "Upregulation",
                                           "Downregulation")),
          #     top_annotation = col_annot_obj,
          border = TRUE,
          height = unit(9, "cm"),
          row_title = "Top 25 scarce receptors",
          show_heatmap_legend = FALSE)

heatmap_rec_scarce







ht_list = heatmap_rec_all %v% heatmap_rec_abund %v% heatmap_rec_scarce 

ht_list

png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptor_scores/receptor_scores_complexheat.png",
    width = 6,
    height = 14,
    units = "in",
    res = 300)

draw(ht_list)

dev.off()


