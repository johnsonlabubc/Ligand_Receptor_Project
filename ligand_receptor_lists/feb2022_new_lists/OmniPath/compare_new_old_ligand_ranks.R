# compare new scarce lig abund rec ranks to old version





# load old ranking and order them by ranking
old <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/overall_ligand_ranks.tsv",
                sep = "\t") %>% 
  arrange(desc(overall_ligand_rank)) %>% 
  mutate(row_number = 1:422) %>% 
  select(ligand,
         old_rank = row_number)


new <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/scarce_lig_abund_rec_ranks.tsv",
                 sep = "\t") %>% 
  arrange(desc(scarce_lig_abund_rec_rank)) %>% 
  mutate(row_number = 1:395) %>% 
  select(ligand,
         new_rank = row_number)


# merge together
old_new <- left_join(old,
                     new) %>% 
  # check that all are equal
  mutate(equal = new_rank - old_rank)
