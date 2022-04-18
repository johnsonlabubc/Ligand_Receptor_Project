# create plot of distribution of the amino acid lengths of all ligands

library(biomaRt)
library(tidyverse)
library(cowplot)

# open ligands list
# this list has been filtered to exclude all ligands with a consensus score below 6
# and to remove all that were manually annoated as "No" to keep in list
# kept the "TBDs"
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_filtered_for_amino_acid_length.txt", 
                       sep = "\t"))


# downloaded from https://www.uniprot.org/uniprot/?query=proteome:UP000005640
uniprot_proteome <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/uniprot-proteome_UP000005640.tab", 
          sep = "\t")) %>% 
  mutate(uniprot = Entry)


ligands_df_lengths <- ligands_df %>% 
  left_join(uniprot_proteome,
            by = "uniprot") %>% 
  select(uniprot,
         genesymbol,
         consensus_score,
         Length)


ligands_df_lengths_filtered <- ligands_df %>% 
  select(hgnc_symbol,
         uniprot_gn_ids,
         description) %>% 
  mutate(genesymbol = hgnc_symbol) %>% 
  left_join(ligands_df_lengths,
            by = "genesymbol") %>% 
  select(hgnc_symbol,
         uniprot_gn_ids,
         description,
         uniprot_gn_ids,
         consensus_score,
         Length) %>% 
# filter out some additional that are wrong
  filter(!hgnc_symbol %in% c("RELN", "VWF", "C4B_2", "TG"))


summary(ligands_df_lengths_filtered$Length)

ggplot(ligands_df_lengths_filtered, aes(x = Length)) +
  geom_histogram(bins = 25, fill = "dodgerblue2") +
  geom_vline(aes(xintercept = min(Length)), color = 'red') +
  labs(x = "Amino Acid Length") + 
  theme_cowplot()

                                   
write_tsv(ligands_df_lengths_filtered, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_filtered_amino_acid_lengths.tsv")
