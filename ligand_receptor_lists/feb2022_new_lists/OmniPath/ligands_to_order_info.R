# create list of all 422 ligands under consideration to be ordered
# purpose is to get quotes & availability info for each protein from suppliers
# include the uniprot protein ID and amino acid size info


library(tidyverse)
library(cowplot)


# load ligands data & filter
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics.tsv", 
                          sep = "\t")) %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4)

# old data for comparison
old_aa_data <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_filtered_amino_acid_lengths.tsv", 
                                           sep = "\t"))


# downloaded from https://www.uniprot.org/uniprot/?query=proteome:UP000005640
uniprot_proteome <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/uniprot-proteome_UP000005640.tab", 
                              sep = "\t")) %>% 
  mutate(uniprot_gn_ids = Entry)


ligands_df_lengths <- ligands_df %>% 
  left_join(uniprot_proteome,
            by = "uniprot_gn_ids") %>% 
  select(hgnc_symbol,
         description,
         uniprot_gn_ids,
         consensus_score,
         # clarify that the amino acid count is for precursor proteins
         precursor_aa_length = Length) %>% 
  # remove TG, which was supposed to already be removed from our ligands list
  filter(!hgnc_symbol %in% c("TG"))

# compute summary statistics
summary(ligands_df_lengths$precursor_aa_length)

# plot distribution
ggplot(ligands_df_lengths, aes(x = precursor_aa_length)) +
  geom_histogram(bins = 25, fill = "#36226B") +
  geom_vline(aes(xintercept = min(precursor_aa_length)), color = 'red') +
  labs(x = "Ligand precursor amino acid counts") + 
  theme_cowplot()


# save
write_tsv(ligands_df_lengths, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_precursor_aa_lengths.tsv")
