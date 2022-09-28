# scoring the 349 receptors based on aggregate rankings of all data of interest
# that we have compiled. Receptor scores will then be combined with ligand
# scores based on known interactions from omnipath data. The combined final 
# ligand scores will be used to decide which ligands to order of in vitro 
# screening.

library(tidyverse)
library(cowplot)


############# load data ###############################


receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics_DE.tsv", 
                          sep = "\t"))




############### prep data for scoring ############################



# replace NAs with 0 in prot`eomics so score function doesn't break
# we will use pooled sc proteomics to capture all proteins of interest
receptors_df$all_sc_proteomics_0 <- 
  lapply(receptors_df$all_sc_proteomics,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))



# prep proteomics differential abundance data
# this returns p-values only for proteins upregulated in stemcell
# with a p value of less than 0.05
receptors_df <- receptors_df %>% 
# Fold change > 1 indicates higher abundancein stem cell  
mutate(proteomics_up_in_sc = ifelse(stemcell_NDislet_proteomics_FC > 1
                                    , ifelse(stemcell_NDislet_proteomics_adj_p_value < 0.05,
                                             stemcell_NDislet_proteomics_adj_p_value,
                                             NA), 
                                    NA))
# replace NAs with 1 in proteomics_up_in_sc score function doesn't break
receptors_df$proteomics_up_in_sc_0 <- 
  lapply(receptors_df$proteomics_up_in_sc,
         replace_na, replace = 1)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))


       

#################### scoring function ################################# 

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
  mutate(proteomics_up_in_sc_score = percent_rank(desc(proteomics_up_in_sc_0))) %>% 

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
         aggregate_score = sum(c(sorted_eBCs_score,
                               combined_sc_max_score,
                               all_sc_proteomics_score,
                               proteomics_up_in_sc_score))) %>% 

  ungroup() %>% 
  select(hgnc_symbol,
         sorted_eBCs_score,
         sorted_immature_beta_score,
         sorted_ins_gfp_gcg_sst_score,
         sorted_ins_gfp_sst_score,
         sorted_pancreatic_proj_score,
         unsorted_immature_beta_score,
         unsorted_ins_gfp_score,
         unsorted_ins_gfp_gcg_sst_score,
         unsorted_pancreatic_proj_score,
         sorted_sc_mean_score,
         unsorted_sc_mean_score,
         sorted_sc_max_score,
         unsorted_sc_max_score,
         combined_sc_max_score,
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


receptors_scores %>% 
  ggplot(aes(x = aggregate_score)) +
  geom_histogram(bins = 20,
                 fill = "#1F9274") +
  labs(x = "Aggregate Score",
       y = "Count") +
  theme_cowplot()


View(receptors_df %>% 
       filter(keep_in_list %in% c("Yes", "TBD")) %>% 
       filter(consensus_score > 7) %>% 
     select(hgnc_symbol,
              stemcell_NDislet_proteomics_FC,
              stemcell_NDislet_proteomics_adj_p_value,
              proteomics_up_in_sc_0))




