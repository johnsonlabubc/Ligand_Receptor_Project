# explore the new proteomics data from Foster lab now that have all the 
# n's. Compare it to the bulk RNA-seq data



library(tidyverse)
library(cowplot)



# open new proteomics data from Grace
proteomics_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ave_proteomics_grace_mar2022.txt", 
                          sep = "\t")) %>% 
  mutate(uniprot_gn_ids = Protein.Group) 



# open receptors list
# this list has been filtered to exclude all receptors with a consensus score below 6
# and to remove all that were manually annotated as "No" to keep in list
# kept the "TBDs"
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_filtered_for_amino_acid_length_alldata.txt", 
                          sep = "\t")) %>% 
  select(hgnc_symbol,
         description,
         comments,
         ensembl_gene_ids,
         uniprot_gn_ids,
         consensus_score,
         islet_tpm,
         islet_tpm_rank,
         islet_proteomics,
         islet_proteomics_rank) %>% 
  # append new proteomics data
  left_join(proteomics_df,
            by = "uniprot_gn_ids")


# open ligands list
# this list has been filtered to exclude all ligands with a consensus score below 6
# and to remove all that were manually annoated as "No" to keep in list
# kept the "TBDs"
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_filtered_for_amino_acid_length_alldata.txt", 
                        sep = "\t")) %>% 
  select(hgnc_symbol,
         description,
         comments,
         ensembl_gene_ids,
         uniprot_gn_ids,
         consensus_score,
         islet_tpm,
         islet_tpm_rank,
         islet_proteomics,
         islet_proteomics_rank) %>% 
  # append new proteomics data
  left_join(proteomics_df,
            by = "uniprot_gn_ids")


# plot human islet bulk RNA-seq against proteomics
ggplot(receptors_df, aes(x = log2(islet_tpm), y = log2(proteomics_human_islets))) +
  geom_point(colour = "dodgerblue2", size = 2) +
  labs(x = "Log2 Human Islet RNA-seq",
       y = "Log2 Human Islet Proteomics") +
  ggtitle("Receptors") +
  scale_x_continuous(limits = c(-1,10)) +
  theme_cowplot()


# plot human islet bulk RNA-seq against proteomics
ggplot(ligands_df, aes(x = log2(islet_tpm), y = log2(proteomics_human_islets))) +
  geom_point(colour = "dodgerblue2", size = 2) +
  labs(x = "Log2 Human Islet RNA-seq",
       y = "Log2 Human Islet Proteomics") +
  ggtitle("Ligands") +
  scale_x_continuous(limits = c(-5,20)) +
  theme_cowplot()

# save tables with just the columns needed for the basics

#receptors
receptors_df_columns <- receptors_df %>% 
  select(hgnc_symbol,
         description,
         ensembl_gene_ids,
         uniprot_gn_ids,
         consensus_score,
         islet_tpm,
         proteomics_human_islets) %>% 
  arrange(desc(islet_tpm))

write_tsv(receptors_df_columns, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_islet_proteomics_rnaseq.tsv")


### save tables with just the columns needed for the basics

# ligands
ligands_df_columns <- ligands_df %>% 
  select(hgnc_symbol,
         description,
         ensembl_gene_ids,
         uniprot_gn_ids,
         consensus_score,
         islet_tpm,
         proteomics_human_islets) %>% 
  arrange(desc(islet_tpm))

write_tsv(ligands_df_columns, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_islet_proteomics_rnaseq.tsv")


