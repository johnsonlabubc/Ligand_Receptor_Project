# combining the receptors & ligands dataframes for convienence,
# starting with the dataframes with n=3 proteomics added


# install.packages("viridis")


library(tidyverse)
library(cowplot)
library(ggrepel) 
library(viridis)




######### combine ligand and receptor dataframes ###############################

# open receptors list with proteomics
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                          sep = "\t"))
# open ligands list with proteomics
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics.tsv", 
                        sep = "\t")) 

# join the dataframes, keeping all rows from both dataframes
genes_df <- receptors_df %>% 
  full_join(ligands_df,
            by = "hgnc_symbol")


# check how many receptors and ligands are in our filtered lists

# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  nrow()
# 422 ligands


# receptors
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  nrow()
# 349 receptors




########### plot human islet proteomics vs bulk RNAseq #############################


# ligands

# plot human islet bulk RNA-seq against proteomics
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(islet_tpm), 
             y = log2(all_islet_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 Human Islet RNA-seq",
       y = "Log2 Human Islet Proteomics") +
  ggtitle("Ligands") +
  scale_x_continuous(limits = c(-8,18)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(ligands_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4)), 
                  size = 3,
                  max.overlaps = 30,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



# receptors

# plot human islet bulk RNA-seq against proteomics
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(islet_tpm), 
             y = log2(all_islet_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 Human Islet RNA-seq",
       y = "Log2 Human Islet Proteomics") +
  ggtitle("Receptors") +
  scale_x_continuous(limits = c(-1,8)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(receptors_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 7)), 
                  size = 3,
                  max.overlaps = 20,
                  force = 1,
                  max.time = 2) +
  theme_cowplot()



# determine which genes have NA for one of the proteomics categories but not the other
# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  filter(!is.na(all_islet_proteomics)) %>%
  filter(is.na(sorted_sc_proteomics)) %>% 
  select(hgnc_symbol,
         sorted_sc_proteomics,
         unsorted_sc_proteomics,
         all_sc_proteomics,
         all_islet_proteomics) %>% 
  # save list of ligand proteins present in islet but missing in stem cell sorted
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_proteins_missing_in_stemcells.tsv")


# receptors
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  filter(!is.na(all_islet_proteomics)) %>%
  filter(is.na(sorted_sc_proteomics)) %>% 
  select(hgnc_symbol,
         sorted_sc_proteomics,
         unsorted_sc_proteomics,
         all_sc_proteomics,
         all_islet_proteomics) %>% 
  # save list of receptor proteins present in islet but missing in stem cell sorted
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_proteins_missing_in_stemcells.tsv")



# do the same for genes in stem cells but not in human islets
View(ligands_df %>% 
       filter(keep_in_list %in% c("Yes", "TBD")) %>% 
       filter(consensus_score > 7) %>% 
       filter(!is.na(all_sc_proteomics)) %>%
       filter(is.na(all_islet_proteomics)) %>% 
       select(hgnc_symbol,
              sorted_sc_proteomics,
              unsorted_sc_proteomics,
              all_sc_proteomics,
              all_islet_proteomics))
# no proteins detected for either ligands or receptors in sorted stem cell 
# & all stem cell that were not present in all human islet prep


########### plot all human islet vs sorted stem cell proteomics ##############

# create dataframe that replaces NA's in stem cell is pseudovalue for plotting
ligands_pseudo_for_NAs <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  # remove rows that are NA in all islet proteomics
  filter(!is.na(all_islet_proteomics)) %>%
  # convert NA's in sorted stem cell to psuedovalues to position on plot
  mutate_at(c("sorted_sc_proteomics"),
         ~replace(., is.na(.), 32.000)) # log2(32) = 5
  

# plot ligands with pseudovalues
ligands_pseudo_for_NAs %>% 
  ggplot(aes(x = log2(all_islet_proteomics), 
             y = log2(sorted_sc_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  # overwrite the colour for the pseudovalues
  geom_point(data = subset(ligands_pseudo_for_NAs, sorted_sc_proteomics == 32.0),
             colour = "darkorange1", size = 2) +
  labs(x = "Log2 Human Islet Proteomics",
       y = "Log2 Sorted SC Proteomics") +
  ggtitle("Ligands") +
  scale_y_continuous(limits = c(4,20)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(ligands_pseudo_for_NAs),
                  size = 3,
                  max.overlaps = 30,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



# receptors

# create dataframe that replaces NA's in stem cell is pseudovalue for plotting
receptors_pseudo_for_NAs <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  # remove rows that are NA in all islet proteomics
  filter(!is.na(all_islet_proteomics)) %>%
  # convert NA's in sorted stem cell to psuedovalues to position on plot
  mutate_at(c("sorted_sc_proteomics"),
            ~replace(., is.na(.), 32.000)) # log2(32) = 5


# plot receptors with pseudovalues
receptors_pseudo_for_NAs %>% 
  ggplot(aes(x = log2(all_islet_proteomics), 
             y = log2(sorted_sc_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  # overwrite the colour for the pseudovalues
  geom_point(data = subset(receptors_pseudo_for_NAs, sorted_sc_proteomics == 32.0),
             colour = "darkorange1", size = 2) +
  labs(x = "Log2 Human Islet Proteomics",
       y = "Log2 Sorted SC Proteomics") +
  ggtitle("Receptors") +
#  scale_y_continuous(limits = c(4,18)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(receptors_pseudo_for_NAs),
                  size = 3,
                  max.overlaps = 40,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



############ Plot consensus score histogram ################################

# plot histogram of consensus scores

#ligands
ligands_df  %>% 
  ggplot(aes(x = consensus_score)) +
  geom_histogram(bins = 15, fill = "#36226B") +
  geom_vline(xintercept = 5, color = 'black') +
  labs(x = "Omnipath consensus score") + 
  ggtitle("Ligands") +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_consensus_score.jpg")

# receptors
receptors_df  %>% 
  ggplot(aes(x = consensus_score)) +
  geom_histogram(bins = 15, fill = "#1F9274") +
  geom_vline(xintercept = 8, color = 'black') +
  labs(x = "Omnipath consensus score") + 
  ggtitle("Receptors") +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_consensus_score.jpg")


############# plot histogram of mean islet tpm ############################

# plot histogram of islet_tpm

#ligands
ligands_df  %>% 
  # display only our final 422 ligands
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(islet_tpm))) +
  geom_histogram(bins = 20, fill = "#36226B") +
  labs(x = "Log2 Human islet TPM") + 
  ggtitle("Ligands") +
 # scale_x_continuous(limits = c(0,24)) +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_islet_tpm_histo.jpg")


# receptors
receptors_df  %>% 
  # display only our final 422 ligands
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(islet_tpm))) +
  geom_histogram(bins = 20, fill = "#1F9274") +
  labs(x = "Log2 Human islet TPM") + 
  ggtitle("Receptors") +
  # scale_x_continuous(limits = c(0,24)) +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_islet_tpm_histo.jpg")


############# tables of top genes based on islet tpm ######################

ligands_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(Ligands = hgnc_symbol,
         islet_tpm) %>% 
  arrange(desc(islet_tpm)) %>% 
  head(25)

receptors_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(Ligands = hgnc_symbol,
         islet_tpm) %>% 
  arrange(desc(islet_tpm)) %>% 
  head(25)


############# plot eBCs vs Beta genes expression  ###############

#ligands
ligands_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(Beta), 
             y = log2(sorted_eBCs),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 Primary Beta cells counts",
       y = "Log2 sorted eBCs counts") +
  ggtitle("Ligands") +
  geom_text_repel(data=subset(ligands_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4),
                                log2(sorted_eBCs) > -10), # to label only certain points
                  size = 3,                  
                  max.overlaps = 12,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



# receptors
receptors_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(Beta), 
             y = log2(sorted_eBCs),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 Primary Beta cells counts",
       y = "Log2 sorted eBCs counts") +
  ggtitle("Receptors") +
  geom_text_repel(data=subset(receptors_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 7),
                              log2(sorted_eBCs) > -10), # to label only certain points
                  size = 3,                  
                  max.overlaps = 7,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()







############ tables of top genes in Beta & eBCs by expression & specificity

### ligands

# top 25 expressed in Beta
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         Beta) %>% 
  arrange(desc(Beta)) %>% 
  head(25)

# top 25 specific in Beta compared to GTEx
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(beta_gtex_specificity = Beta/gtex_mean) %>% 
  select(hgnc_symbol,
         beta_gtex_specificity) %>% 
  arrange(desc(beta_gtex_specificity)) %>% 
  head(25)


# top 25 expressed in eBCs
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         sorted_eBCs) %>% 
  arrange(desc(sorted_eBCs)) %>% 
  head(25)

# top 25 specific in eBCs compared to GTEx
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(sorted_eBCs_gtex_specificity = sorted_eBCs/gtex_mean) %>% 
  select(hgnc_symbol,
         sorted_eBCs_gtex_specificity) %>% 
  arrange(desc(sorted_eBCs_gtex_specificity)) %>% 
  head(25)



### receptors

# top 25 expressed in Beta
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         Beta) %>% 
  arrange(desc(Beta)) %>% 
  head(25)

# top 25 specific in Beta compared to GTEx
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(beta_gtex_specificity = Beta/gtex_mean) %>% 
  select(hgnc_symbol,
         beta_gtex_specificity) %>% 
  arrange(desc(beta_gtex_specificity)) %>% 
  head(25)


# top 25 expressed in eBCs
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         sorted_eBCs) %>% 
  arrange(desc(sorted_eBCs)) %>% 
  head(25)

# top 25 specific in eBCs compared to GTEx
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(sorted_eBCs_gtex_specificity = sorted_eBCs/gtex_mean) %>% 
  select(hgnc_symbol,
         sorted_eBCs_gtex_specificity) %>% 
  arrange(desc(sorted_eBCs_gtex_specificity)) %>% 
  head(25)



########## rank of ranks table #######################

### ligands


ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  
  mutate(sorted_eBCs_rank = dense_rank(desc(sorted_eBCs)),
         sorted_eBCs_specificity
         # had to add very small amounts to prevent infinity values
         # the values are so small that adding 0.01 screws up the data
         eBCs_human_beta_ratio = ((sorted_eBCs + 0.0000001) / (Beta + 0.0000001)),
         human_beta_eBCs_ratio = ((Beta + 0.0000001) / (sorted_eBCs + 0.0000001)))  
  
  
  select(islet_tpm_rank,
         islet_specificity_rank,
         beta_rank,
         beta_specificity_rank,
         sorted_eBCs_rank, 
         sorted_sc_proteomics_rank,
         all_islet_proteomics_rank,
         ) 

  
  
########### plot eBCs vs Beta difference in expression & GTEx specificity##############
  
#ligands
    
ligands_df_ratio_ranks <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(eBCs_human_beta_ratio = ((sorted_eBCs + 0.0000001) / (Beta + 0.0000001)),
         sorted_eBCs_gtex_specificity = ((sorted_eBCs + 0.0000001) / (gtex_mean + 0.0000001)),
         human_beta_gtex_specificity = ((Beta + 0.0000001) / (gtex_mean + 0.0000001)),
         eBCs_human_beta_gtex_specifictiy_ratio = ((sorted_eBCs_gtex_specificity + 0.0000001) / 
                                                     (human_beta_gtex_specificity + 0.0000001)),
         eBCs_human_beta_ratio_rank = dense_rank(desc(eBCs_human_beta_ratio)),
         eBCs_human_beta_gtex_specifictiy_ratio_rank = dense_rank(desc(eBCs_human_beta_gtex_specifictiy_ratio))) %>% 
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_human_beta_ratio,
         sorted_eBCs_gtex_specificity,
         human_beta_gtex_specificity,
         eBCs_human_beta_gtex_specifictiy_ratio,
         eBCs_human_beta_ratio_rank,
         eBCs_human_beta_gtex_specifictiy_ratio_rank) %>% 
  filter(sorted_eBCs != 0) %>%
  filter(Beta != 0)
  
  
# plot difference in specificity against difference in expression
ligands_df_ratio_ranks %>% 
      ggplot(aes(x = log2(eBCs_human_beta_ratio), 
               y = log2(eBCs_human_beta_gtex_specifictiy_ratio),
               label = hgnc_symbol)) +
    geom_point(colour = "#36226B", size = 2) +
    labs(x = "Log2 eBCs to beta-cell expression ratio",
         y = "Log2 eBCs to beta-cell specificity ratio") +
    ggtitle("Ligands") +
    geom_text_repel(size = 3,                  
                    max.overlaps = 30,
                    force = 1,
                    max.time = 10) +
    theme_cowplot()


# plot difference in expression against Beta expression
ligands_df_ratio_ranks %>% 
  ggplot(aes(x = log2(sorted_eBCs), 
             y = log2(eBCs_human_beta_ratio),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 beta cell expression",
       y = "Log2 eBCs to beta cell expression ratio") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()


# plot eBCs GTEx specificity vs beta-cell GTEx specificity
ligands_df_ratio_ranks %>% 
  ggplot(aes(x = log2(human_beta_gtex_specificity), 
             y = log2(sorted_eBCs_gtex_specificity),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 eBCs specificity",
       y = "Log2 beta cell specificity") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()




ligands_df_ratio_ranks <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(sorted_eBCs_rank = dense_rank(sorted_eBCs),
         beta_rank = dense_rank(Beta), 
         eBCs_human_beta_ratio = ((sorted_eBCs_rank) / (beta_rank)),
         sorted_eBCs_gtex_specificity = ((sorted_eBCs_rank) / (gtex_mean)),
         human_beta_gtex_specificity = ((beta_rank) / (gtex_mean)),
         eBCs_human_beta_gtex_specifictiy_ratio = ((sorted_eBCs_gtex_specificity) / 
                                                     (human_beta_gtex_specificity)),
         eBCs_human_beta_ratio_rank = dense_rank(desc(eBCs_human_beta_ratio)),
         eBCs_human_beta_gtex_specifictiy_ratio_rank = dense_rank(desc(eBCs_human_beta_gtex_specifictiy_ratio))) %>% 
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         beta_rank,
         sorted_eBCs_rank,
         eBCs_human_beta_ratio,
         sorted_eBCs_gtex_specificity,
         human_beta_gtex_specificity,
         eBCs_human_beta_gtex_specifictiy_ratio,
         eBCs_human_beta_ratio_rank,
         eBCs_human_beta_gtex_specifictiy_ratio_rank)



# plot difference in expression against Beta expression
ligands_df_ratio_ranks %>% 
  ggplot(aes(x = sorted_eBCs, 
             y = eBCs_human_beta_ratio,
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 beta cell expression",
       y = "Log2 eBCs to beta cell expression ratio") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()


# difference instead of ratio
ligands_df_difference <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_gtex_specif = ((sorted_eBCs + 0.0000001) / (gtex_mean + 0.0000001)),
         beta_gtex_specif = ((Beta + 0.0000001) / (gtex_mean + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_gtex_specif - beta_gtex_specif) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         eBCs_gtex_specif,
         beta_gtex_specif,
         eBCs_beta_specif_diff)

ligands_df_difference %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs beta specificity diff") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()

## receptors

# difference instead of ratio
receptors_df_difference <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_gtex_specif = ((sorted_eBCs + 0.0000001) / (gtex_mean + 0.0000001)),
         beta_gtex_specif = ((Beta + 0.0000001) / (gtex_mean + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_gtex_specif - beta_gtex_specif) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         eBCs_gtex_specif,
         beta_gtex_specif,
         eBCs_beta_specif_diff)

receptors_df_difference %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs beta specificity diff") +
  ggtitle("Receptors") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()

  
######### eBCs vs Beta difference in expression & in-sample specificity #########

# figures out which columns to average for sorted sample mean
ligands_df[48:52] %>% 
  colnames()

ligands_df$mean_sorted_scRNAseq <- rowMeans(ligands_df[48:52])

# difference instead of ratio
ligands_df_sample_diff <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_sample_specif = ((sorted_eBCs + 0.0000001) / (mean_sorted_scRNAseq + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_sample_specif - beta_specificity) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         mean_sorted_scRNAseq,
         mean_islet_scRNAseq,
         beta_specificity,
         eBCs_sample_specif,
         eBCs_beta_specif_diff)

# plot ligands
ligands_df_sample_diff %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs-primary-beta specificity diff") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



## repeat for receptors

# figures out which columns to average for sorted sample mean
receptors_df[48:52] %>% 
  colnames()

receptors_df$mean_sorted_scRNAseq <- rowMeans(receptors_df[48:52])

# difference instead of ratio
receptors_df_sample_diff <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_sample_specif = ((sorted_eBCs + 0.0000001) / (mean_sorted_scRNAseq + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_sample_specif - beta_specificity) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         mean_sorted_scRNAseq,
         mean_islet_scRNAseq,
         beta_specificity,
         eBCs_sample_specif,
         eBCs_beta_specif_diff)

# plot receptors
receptors_df_sample_diff %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs-primary-beta specificity diff") +
  ggtitle("Receptors") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()


