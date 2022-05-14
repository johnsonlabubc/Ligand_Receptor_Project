# add new proteomics data that includes n=3 technical replicates
# includes both the stem cell and human islet proteomics



# load libraries
library(tidyverse)


##### load dataframes ###########

# load genes list with all prior expression data and annotations
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_stemcell_updateddata.tsv", 
                      sep = "\t")


# open new proteomics data from Grace
proteomics_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/proteomics_n_3_from_grace.txt", 
                           sep = "\t")) %>% 
  mutate(uniprot_gn_id = Protein.Ids) %>% 
  select(-Protein.Ids)



# regenerate old table of uniprot ID's to hgnc symbol conversion
gene_uniprot_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_annot_manual.tsv", 
                      sep = "\t") %>% 
  rename(genecards_url = 17)

# create table with just hgnc symbols and uniprot ids
# will append expression data to this table
gene_uniprot_list <- gene_uniprot_list %>% 
  dplyr::select(hgnc_symbol,
                uniprot_gn_id)

# get max value for each hgnc gene
proteomics_uniprot <- left_join(gene_uniprot_list,
                                      proteomics_df,
                                      by = "uniprot_gn_id") %>% 
  # drop the uniprot gene id column
  dplyr::select(-uniprot_gn_id) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(sorted_sc_proteomics = max(sorted.avg, na.rm = TRUE),
            unsorted_sc_proteomics = max(unsorted.avg, na.rm = TRUE),
            all_sc_proteomics = max(stem.mean, na.rm = TRUE),
            ND_islet_proteomics = max(ND.mean, na.rm = TRUE),
            T2D_islet_proteomics = max(T2D.mean, na.rm = TRUE),
            all_islet_proteomics = max(human.mean, na.rm = TRUE)) %>% 
  # convert -inf values to NAs
  mutate_if(is.numeric, list(~na_if(., -Inf)))

# append proteomics data to our genes list using hgnc
gene_list_final <- left_join(gene_list,
                             proteomics_uniprot,
                             by = "hgnc_symbol") %>% 
  # calculate the ranks
  mutate(sorted_sc_proteomics_rank = dense_rank(desc(sorted_sc_proteomics)),
         unsorted_sc_proteomics_rank = dense_rank(desc(unsorted_sc_proteomics)),
         all_sc_proteomics_rank = dense_rank(desc(all_sc_proteomics)),
         ND_islet_proteomics_rank = dense_rank(desc(ND_islet_proteomics)),
         T2D_islet_proteomics_rank = dense_rank(desc(T2D_islet_proteomics)),
         all_islet_proteomics_rank = dense_rank(desc(all_islet_proteomics))) %>% 
  # reorder columns to put ranks next to abundance
    relocate(sorted_sc_proteomics,
         sorted_sc_proteomics_rank,
         unsorted_sc_proteomics,
         unsorted_sc_proteomics_rank,
         all_sc_proteomics,
         all_sc_proteomics_rank,
         ND_islet_proteomics,
         ND_islet_proteomics_rank,
         T2D_islet_proteomics,
         T2D_islet_proteomics_rank,
         all_islet_proteomics,
         all_islet_proteomics_rank,
         .after = sorted_pancreatic_proj) %>% 
   # fix minor prior error where beta_specificity column is in wrong spot
   relocate(beta_specificity,
            .after = delta_specificity)

# save full spreadsheet with all 12 new proteomics data columns
write_tsv(gene_list_final, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv")

