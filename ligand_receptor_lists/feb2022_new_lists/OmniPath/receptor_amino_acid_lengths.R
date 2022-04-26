# create plot of distribution of the amino acid lengths of all receptors

library(biomaRt)
library(tidyverse)
library(cowplot)

# open receptors list
# this list has been filtered to exclude all receptors with a consensus score below 6
# and to remove all that were manually annoated as "No" to keep in list
# kept the "TBDs"
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_filtered_for_amino_acid_length.txt", 
                        sep = "\t")) %>% 
  mutate(uniprot = uniprot_gn_ids)


# downloaded from https://www.uniprot.org/uniprot/?query=proteome:UP000005640
uniprot_proteome <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/uniprot-proteome_UP000005640.tab", 
                              sep = "\t")) %>% 
  mutate(uniprot = Entry)


receptors_df_lengths <- receptors_df %>% 
  left_join(uniprot_proteome,
            by = "uniprot") %>% 
  select(uniprot,
         consensus_score,
         Length,
         genesymbol = hgnc_symbol)


receptors_df_lengths_filtered <- receptors_df %>% 
  select(hgnc_symbol,
         uniprot_gn_ids,
         description) %>% 
  mutate(genesymbol = hgnc_symbol) %>% 
  left_join(receptors_df_lengths,
            by = "genesymbol") %>% 
  select(hgnc_symbol,
         uniprot_gn_ids,
         description,
         uniprot_gn_ids,
         consensus_score,
         Length) #%>% 
  # filter out some additional that are wrong
  #filter(!hgnc_symbol %in% c("RELN", "VWF", "C4B_2", "TG"))


summary(receptors_df_lengths_filtered$Length)

ggplot(receptors_df_lengths_filtered, aes(x = Length)) +
  geom_histogram(bins = 25, fill = "dodgerblue2") +
#  geom_vline(aes(xintercept = min(Length)), color = 'red') +
  labs(x = "Receptors Amino Acid Lengths") + 
  theme_cowplot()


write_tsv(receptors_df_lengths_filtered, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_filtered_amino_acid_lengths.tsv")
