### repeat process that I did for the Omnipath ligands for the receptors
# this script is based off of omnipath_explore.R

# goal is to get full list of receptors from omnipath

library(OmnipathR)
library(tidyverse)


############## Prepare initial omnipath receptor list ##########################

## We check some of the different intercell categories
get_intercell_generic_categories()

# get list of all receptors from omnipath
intercell_receptors <- import_omnipath_intercell(categories = c("receptor"),
                                                   source = "composite") %>% 
                         filter(entity_type == "protein")

# check how many unique receptor gene symbols
as.data.frame(intercell_receptors$genesymbol) %>% 
  distinct() %>% 
  nrow()
# 4286 unique receptor gene symbols

# check what percent are secreted
intercell_receptors %>% 
  group_by(secreted) %>% 
  summarise(n())
# 555 of 4286 are secreted
# keep these all in our list for now


#confirm that all receptors are recievers 
intercell_receptors %>% 
  group_by(transmitter) %>% 
  summarise(n())

intercell_receptors %>% 
  group_by(receiver) %>% 
  summarise(n())

# confirmed that all receptors are receivers and are not transmitters


### lets continue for now with this receptors list
receptors_df <- intercell_receptors %>% 
  select(category,
         uniprot,
         genesymbol,
         consensus_score,
         secreted,
         plasma_membrane_transmembrane,
         plasma_membrane_peripheral) 



###  get table of all intercell signalling categories genes for reference
intercell_all <- import_omnipath_intercell(source = "composite") %>%
  filter(entity_type == "protein")
# gives us 53,000 genes

# create table of just the key categories from fig 2 in their paper
intercell_main <- import_omnipath_intercell(categories = c("adhesion",
                                                           "receptor",
                                                           "ligand",
                                                           "ecm",
                                                           "secreted_enzyme",
                                                           "cell_surface_enzyme",
                                                           "cell_surface_ligand",
                                                           "transporter"),
                                            source = "composite") %>% 
  filter(entity_type == "protein")
# gives us 8887 genes 


###################  append categories info receptor list ##################
### in 1 column

# first do with all categories
# widen our datatable using reshape2
intercell_all_cat <-intercell_all %>% 
  select(uniprot,
         category) %>% 
  reshape2::dcast(uniprot ~ category, function(x) any(length(x)))

# replace TRUE's with column names
intercell_final <- intercell_all_cat %>% 
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(omnipath_categories = paste0(key,collapse=", ")) 


## repeat with just the main categories

# widen our datatable using reshape2
intercell_main_cat <- intercell_main %>% 
  select(uniprot,
         category) %>% 
  reshape2::dcast(uniprot ~ category, function(x) any(length(x)))

# replace TRUE's with column names
intercell_main_final <- intercell_main_cat %>% 
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(omnipath_main_categories = paste0(key,collapse=", ")) 


## append both the new categories columns to our ligand list

# remove duplicates in ligand list
receptors_df <- receptors_df %>% 
  select(uniprot, genesymbol, consensus_score) %>% 
  # distinct based on uniprot only (to remove duplicates with diff consensus scores)
  distinct(uniprot, .keep_all = TRUE)
# brings us from 1036 ligands (with duplicates) to 942

# join tables
receptors_with_cat <- receptors_df %>% 
  left_join(intercell_main_final,
            by = "uniprot") %>% 
  left_join(intercell_final,
            by = "uniprot") %>% 
  rename(generic_categories = omnipath_main_categories) %>% 
  rename(all_categories = omnipath_categories) 




############### append known ligands to receptor list ##################
# also add the sources for each known ligand (omnipath has this data)

## We import the intercell data into a dataframe
interactions <- import_omnipath_interactions(resources = NULL,
                                             datasets = "omnipath",
                                             organism = 9606) # human


## get ligand uniprot symbols

# widen our datatable using reshape2
interactions_wide <- interactions %>% 
  select(uniprot = target,
         source) %>% 
  reshape2::dcast(uniprot ~ source, function(x) any(length(x))) %>% 
  # drop columns containing complexes
  select(-contains("COMPLEX"))

# replace TRUE's with column names
interactions_compact <- interactions_wide %>% 
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(ligand_uniprot = paste0(key,collapse=", ")) 


## repeat process with ligand gene symbols

# widen our datatable using reshape2
interactions_genesymbol <- interactions %>% 
  select(uniprot = target,
         source_genesymbol) %>% 
  reshape2::dcast(uniprot ~ source_genesymbol, function(x) any(length(x))) %>% 
  # drop columns containing
  select(-contains("COMPLEX")) %>% 
  # replace TRUE's with column names
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(ligand_symbol = paste0(key,collapse=", ")) 


# now get ligand references

# widen our datatable using reshape2
interactions_references <- interactions %>% 
  select(uniprot = target,
         source_genesymbol,
         references) %>%
  # remove rows where ligands symbol contains underscore (b/c means is a complex)
  dplyr::filter(str_detect(source_genesymbol, "_") != TRUE) %>% 
  # create new column containing ligand-reference pairs
  mutate(ligand_references = paste0(source_genesymbol, 
                                      " [",
                                      references,
                                      "]")) %>% 
  # select just the new column
  select(uniprot,
         ligand_references) %>% 
  group_by(uniprot) %>% 
  summarise(ligand_references = paste0(ligand_references, collapse=", ")) 


### append all our new receptor info columns to ligands datatable
receptors_with_annot <- receptors_with_cat %>% 
  left_join(interactions_compact,
            by = "uniprot") %>% 
  left_join(interactions_genesymbol,
            by = "uniprot") %>% 
  left_join(interactions_references,
            by = "uniprot")

# save ligands table
write_tsv(receptors_with_annot, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors.tsv")

