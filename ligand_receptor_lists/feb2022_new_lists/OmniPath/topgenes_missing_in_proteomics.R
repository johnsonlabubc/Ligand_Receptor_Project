# create list of top genes detected in humal islet RNA-seq that are
# missing in the proteomics. Will work with Jason in Foster lab to see if 
# possible to go through proteomics data again to try & detect these proteins
# make list for all genes, but also for top ligands/receptors


library(tidyverse)


# load data
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_bulkRNA.tsv", 
                          sep = "\t"))

ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_bulkRNA.tsv", 
                        sep = "\t"))




######################## ligands ########################################


# replace NAs with 0 in proteomics doesn't break
# we will use pooled sc proteomics to capture all proteins of interest
ligands_df$all_sc_proteomics_0 <- 
  lapply(ligands_df$all_sc_proteomics,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_df <- as.data.frame(lapply(ligands_df, unlist))


# select only columns of interest & filter for genes with "0" in proteomics
ligands_df_top <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         description,
         uniprot_gn_ids,
         Lund_h_islet_bulkRNA,
         Melton_h_islet_bulkRNA,
         sorted_eBCs,
         mean_islet_scRNAseq,
         all_sc_proteomics_0) %>% 
  rowwise() %>% 
  mutate(mean_rnaseq = mean(c(Lund_h_islet_bulkRNA,
                              Melton_h_islet_bulkRNA,
                              sorted_eBCs,
                              mean_islet_scRNAseq))) %>% 
  ungroup() %>% 
  filter(all_sc_proteomics_0 == 0) %>% 
  arrange(desc(Lund_h_islet_bulkRNA))




######################## receptors ########################################


# replace NAs with 0 in proteomics doesn't break
# we will use pooled sc proteomics to capture all proteins of interest
receptors_df$all_sc_proteomics_0 <- 
  lapply(receptors_df$all_sc_proteomics,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_df <- as.data.frame(lapply(receptors_df, unlist))


# select only columns of interest & filter for genes with "0" in proteomics
receptors_df_top <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         description,
         uniprot_gn_ids,
         Lund_h_islet_bulkRNA,
         Melton_h_islet_bulkRNA,
         sorted_eBCs,
         mean_islet_scRNAseq,
         all_sc_proteomics_0) %>% 
  rowwise() %>% 
  mutate(mean_rnaseq = mean(c(Lund_h_islet_bulkRNA,
                              Melton_h_islet_bulkRNA,
                              sorted_eBCs,
                              mean_islet_scRNAseq))) %>% 
  ungroup() %>% 
  filter(all_sc_proteomics_0 == 0) %>% 
  arrange(desc(Lund_h_islet_bulkRNA))





################ genome wide using the bulk rna on proteomics preps ####################

# for now just look at human islet data, but in future also look at stem cell data
# cause would be interesting to find proteins detected in stem cells but not h-islets

# load data from jelena with per donor bulk rna
bulk_rna_protsamples <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/bulk_rna_proteomics_samples/HumanIslets_OxfordStanfordCombined_GeneExp_ProteomicsSamples.txt", 
                          sep = "\t"))

# get mean of all donors
bulk_rna_protsamples %>% 
  mutate(mean_rnaseq = rowMeans(bulk_rna_protsamples[3:98])) %>% 
  select(1,2,99) %>% 
  View()
  

  

proteomics_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/bulk_rna_proteomics_samples/HumanIslets_JohnsonUBCCombined_Proteomics.csv"))

# load data from jelena with per donor proteomics


# save data
  


