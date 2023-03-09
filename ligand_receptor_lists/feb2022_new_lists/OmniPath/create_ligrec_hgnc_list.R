library(tidyverse)


receptors <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_scores_3.tsv",
               sep = "\t") %>% 
  select(1) %>% 
  mutate(type = "receptor")


ligands <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligand_scores_3.tsv",
                      sep = "\t") %>% 
  select(1) %>% 
  mutate(type = "ligand")



ligrec_hgnc_symbols <- full_join(receptors, ligands)


# save
write_tsv(ligrec_hgnc_symbols,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_hgnc_symbols.tsv")
