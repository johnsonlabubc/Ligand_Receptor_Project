# create list of interactions between our finalized ligands list & receptors list
# will get interactions data from omnipath
# and then will use this data for the final scoring of ligands for screening

library(OmnipathR)
library(tidyverse)
library(cowplot)



################# load interaction data from omnipath #########################

## The interactions are stored into a data frame.
interactions <- import_omnipath_interactions(resources = NULL,
                                             organism = 9606)


# for now, don't extra interactions that don't have literature references
# ligand-receptor interactions without literature reference
# interactions_extra <- import_ligrecextra_interactions(organism=9606)


# check how many duplicates for each gene at uniprot level
as.data.frame(interactions$source) %>% 
  distinct() %>% 
  nrow()
# 5270 unique uniprot ligands

# check how many duplicates for each gene at gene symbol level
as.data.frame(interactions$source_genesymbol) %>% 
  distinct() %>% 
  nrow()
# 5250 unique gene symbol ligands
# so only 20 were duplicates with different uniprot IDs


# check for duplicates and directionality
interactions_filtered <- interactions %>% 
  dplyr::filter(source_genesymbol == "INSR" | source_genesymbol == "INS") %>% 
  dplyr::filter(target_genesymbol == "INSR" | target_genesymbol == "INS")
# confirmed that there are receptors found in the "source column too"
# INS -> INSR and INSR -> INS are both present, but are not the same.
# they have different references




######## filter interactions to only include lig/rec's from our lists ##########

# load our ligand and receptor lists
ligand_list <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligand_scores_3.tsv", 
                           sep = "\t")) %>% 
  select(hgnc_symbol)

receptors_list <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_scores_3.tsv", 
                         sep = "\t")) %>% 
  select(hgnc_symbol)


# filter interactions list
interactions_filtered <- interactions %>% 
  dplyr::filter(source_genesymbol %in% ligand_list$hgnc_symbol) %>% 
  # 3398 interactions after filtering for our 422 ligands
  dplyr::filter(target_genesymbol %in% receptors_list$hgnc_symbol)
  # 1569 left after filtering for our receptors as targets as well


# check for duplicate gene symbols again
# check how many duplicates for each gene at uniprot level
as.data.frame(interactions_filtered$source) %>% 
  distinct() %>% 
  nrow()
# 396 unique uniprot ligands

# check how many duplicates for each gene at gene symbol level
as.data.frame(interactions_filtered$source_genesymbol) %>% 
  distinct() %>% 
  nrow()
# 395 unique gene symbol ligands
# so only 1 duplicate 

# create histogram of the number of resources which identified each interaction
interactions_filtered %>% 
  ggplot(aes(x = n_resources, fill = n_resources, group = n_resources)) +
  geom_histogram(bins = 32,
                 linewidth = .1,
                 color = "black") +
#  geom_vline(xintercept = 3, color = 'black') +
  scale_y_continuous(breaks = c(0, 100, 250)) +
  labs(x = "resources per interaction",
       y= "count") + 
  scale_fill_fermenter(palette = "YlGnBu",
                        direction = 1,
                        breaks = c(1,3,5,10,15)) +
  theme_cowplot()

# based on the histogram, may make sense to make a threshold of 2 or 3 for inclusion

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_resources_histo_5.png",
       scale = 1.7)



# save filtered interactions list
write_tsv(interactions_filtered,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_interactions_filtered.tsv")





