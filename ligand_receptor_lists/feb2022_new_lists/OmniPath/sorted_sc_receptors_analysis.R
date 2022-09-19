# to determine which receptors to target in the stem cell screening
# create an aggregate score of receptor expression based on all cell types 
# in the sorted stem cells dataset
# b/c it may be useful to target receptors in cell types other than the eBCs


library(tidyverse)


# load receptors data
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                          sep = "\t"))

receptors_sorted_sc_data <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         sorted_immature_beta,
         sorted_ins_gfp_gcg_sst,
         sorted_ins_gfp_sst,
         sorted_pancreatic_proj,
         sorted_sc_proteomics,
         sorted_sc_proteomics_rank) %>% 
  # create ranks
  mutate(eBCs_rank = dense_rank(desc(sorted_eBCs)),
         immature_beta_rank = dense_rank(desc(sorted_immature_beta)),
         ins_gfp_gcg_sst_rank = dense_rank(desc(sorted_ins_gfp_gcg_sst)),
         ins_gfp_sst_rank = dense_rank(desc(sorted_ins_gfp_sst)),
         pancreatic_proj_rank = dense_rank(desc(sorted_pancreatic_proj))) %>% 
  # create aggregate ranks, weighing all equally
  rowwise() %>% 
  mutate(sorted_sc_mean_rank = mean(c(eBCs_rank,
                                           immature_beta_rank,
                                           ins_gfp_gcg_sst_rank,
                                           ins_gfp_sst_rank,
                                           pancreatic_proj_rank)),
         sorted_sc_min_rank = min(c(eBCs_rank,
                                    immature_beta_rank,
                                    ins_gfp_gcg_sst_rank,
                                    ins_gfp_sst_rank,
                                    pancreatic_proj_rank))) %>% 
  ungroup()

receptors_sorted_sc_data %>% 
  select(hgnc_symbol,
         eBCs_rank,
         sorted_sc_mean_rank,
         sorted_sc_min_rank) %>% 
  View()



#### ensure that all receptors detected in proteomics are included in prioritized list

         

# replace NAs with 0 in sorted proteomics so score function doesn't break
receptors_sorted_sc_data$sorted_proteomics_0 <- 
  lapply(receptors_sorted_sc_data$sorted_sc_proteomics,
        replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_sorted_sc_data <- as.data.frame(lapply(receptors_sorted_sc_data, unlist))
  
# define function to score the proteomics to either 0 or 1000 (binary)
ScoreFn <- function(x){
    score <- 0
    if(x>1) {
      score <- 1000
    } else {
      score <- 0
    }
    
    return(score)
  }
  
# apply the scoring function to each row of the proteomics data
# using the column without NA's, as the score function can't handle NA's
receptors_sorted_sc_data$sorted_proteomics_score <- 
  lapply(receptors_sorted_sc_data$sorted_proteomics_0,
         ScoreFn)

# define score function 


