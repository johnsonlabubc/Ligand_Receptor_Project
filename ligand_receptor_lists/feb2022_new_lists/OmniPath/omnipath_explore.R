# initial exploration of OmniPathR package and data

# install omnipathR package from github b/c it's updated more than bioconductor package
# require(devtools)
# install_github('saezlab/OmnipathR')

library(OmnipathR)
library(tidyverse)


## We check some of the different interaction databases
interac_resoucres <- as.tibble(get_interaction_resources())

# save this list for future reference
write_tsv(interac_resoucres, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/interaction_resources.tsv")


## The interactions are stored into a data frame.
interactions <- import_omnipath_interactions(resources = NULL,
                                             organism = 9606)

View(head(interactions))

## We visualize the first interactions in the data frame.
print_interactions(head(interactions))


# ligand-receptor interactions without literature reference
interactions_extra <- import_ligrecextra_interactions(organism=9606)

View(interactions_extra)
# can see there are 3013 extra ligand-receptor interactions
# a lot of them are found from multiple references, including CellChat, connectomedb, etc


# save this list for future reference
write_tsv(interactions_extra, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrecextra_interactions.tsv")


## We check some of the different intercell categories
get_intercell_generic_categories()

get_intercell_categories()

## We import the intercell data into a dataframe
intercell <- import_omnipath_intercell() 


  
nrow(intercell)
# we see that there are 323572 intercellular communication 
# role records. This includes many duplicate records for each gene
# For each gene, there are multiple rows with each data source & annotation

# To check, lets see how many unique genes there are here
as.data.frame(intercell$genesymbol) %>% 
  distinct() %>% 
  nrow()
# there are 46207
# thats still a ridiculously high number

# now lets try with only ligands & receptors
intercell_ligrec <- import_omnipath_intercell(categories = c("ligand", 
                                                             "receptor"))

as.data.frame(intercell_ligrec$genesymbol) %>% 
  distinct() %>% 
  nrow()
# there are 8017... which is still A LOT


# lets try only going resource specific
get_intercell_resources()

intercell_ligrec <- import_omnipath_intercell(categories = c("ligand", 
                                                             "receptor"),
                                              resources = "OmniPath")
as.data.frame(intercell_ligrec$genesymbol) %>% 
  distinct() %>% 
  nrow()
# interestingly we still get 8017!
# so all of them are in the "OmniPath" resource

# lets try just getting ligands, and filter out complexes
intercell_ligands <- import_omnipath_intercell(categories = "ligand",
                                               source = "composite",
                                               secreted = TRUE) %>% 
                     filter(entity_type == "protein")
# 1479 ligands, including complexes
# tried changing to aspect = functional & scope = specific
# but both resulted in 0 genes being returned

nrow(intercell_ligands)
# 897 ligands, not including complexes

as.data.frame(intercell_ligands$uniprot) %>% 
  distinct() %>% 
  nrow()
# confirmed that all 897 ligands are unique
# NOTE: use uniprot for joins/merge's cause a few don't have genesymbol


intercell_ligands %>% 
  group_by(plasma_membrane_transmembrane) %>% 
  summarise(n())
# of the 897, 120 have a plasma memberane transmembrane
# will keep them cause they are also defined as secreted out of the cell

#confirm that all ligands are transmitters too
as.data.frame(intercell_ligands$transmitter) %>% 
  distinct()
# confirmed that all ligands are transmitters, are secreted, and are not receivers


# In the intercell database a `consensus_score` is available for each composite 
# record: the number of resources annotating the protein for a certain category

# noted that some ligands from ligrec_extra are not in the 897 ligand list
# check to see if they are in the cell_surface_ligand list
intercell_all_ligands <- import_omnipath_intercell(categories = c("ligand",
                                                                     "cell_surface_ligand"),
                                               source = "composite",
                                               secreted = TRUE) %>% 
  filter(entity_type == "protein")
# this added about 120 genes, of which 45 were unique, so total is now 1036
as.data.frame(intercell_all_ligands$uniprot) %>% 
  distinct() %>% 
  nrow()


# compare list to ligrec_extra
extra_ligands <- interactions_extra %>% 
  select(source) %>% 
  distinct()
# has 1651 unique ligands, though 30 of them are complexes

# filter out the complexes
extra_ligands <- extra_ligands %>% 
  filter(substr(source, 1, 3) != "COM")
# left with 1621 extra ligands

# check how many of the main list ligands are in the extra ligands
intercell_all_ligands$uniprot %in% extra_ligands$source %>% 
  summary()
# 866 are in extra list, 170 are not

# now check opposite (how many extra ligands are in our main list)
extra_ligands$source %in% intercell_all_ligands$uniprot %>% 
  summary()
# 775 are in it, 846 are not

### lets continue for now without the extra ligands
ligands_df <- intercell_all_ligands %>% 
  select(category,
         uniprot,
         genesymbol,
         consensus_score) #%>% 
#  filter(consensus_score > 2)
# the 30 cell surface ligands all have a consensus score <=4

### also get table of all intercell signalling categories genes for reference
intercell_all <- import_omnipath_intercell(source = "composite") %>%
                   filter(entity_type == "protein")
# gives us 53,000 genes

# check how many duplicate genes
as.data.frame(intercell_all$uniprot) %>% 
  distinct() %>% 
  nrow()
# 20,925 unique genes

# save table of all interaction categories genes
write_tsv(intercell_all, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/allcategories_intercell_genes.tsv")


# check all categories in our table
as.data.frame(intercell_all$category) %>% 
  distinct()

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

# check how many duplicate genes
as.data.frame(intercell_main$uniprot) %>% 
  distinct() %>% 
  nrow()
# 6271 unique genes here

# save table of main interaction categories genes
write_tsv(intercell_main, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/maincategories_intercell_genes.tsv")


### try to append all the categories that our ligand genes have to our ligands list
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
ligands_df <- ligands_df %>% 
  select(uniprot, genesymbol, consensus_score) %>% 
  # distinct based on uniprot only (to remove duplicates with diff consensus scores)
  distinct(uniprot, .keep_all = TRUE)
# brings us from 1036 ligands (with duplicates) to 942

# join tables
ligands_with_cat <- ligands_df %>% 
  left_join(intercell_main_final,
            by = "uniprot") %>% 
  left_join(intercell_final,
            by = "uniprot") %>% 
  rename(generic_categories = omnipath_main_categories) %>% 
  rename(all_categories = omnipath_categories) 

####### append known receptors to ligands list ###########
# also add the sources for each known receptor (omnipath has this data)

## We import the intercell data into a dataframe
interactions <- import_omnipath_interactions(resources = NULL,
                                             datasets = "omnipath",
                                             organism = 9606) # human

nrow(interactions)
# 40014 interactions
# check how many duplicates for each gene
as.data.frame(interactions$source) %>% 
  distinct() %>% 
  nrow()
# 5270 unique uniprot ligands

# widen our datatable using reshape2
interactions_wide <- interactions %>% 
  select(uniprot = source,
         target) %>% 
  reshape2::dcast(uniprot ~ target, function(x) any(length(x))) %>% 
  # drop columns containing complexes (6800 with complexes, 6398 without)
  select(-contains("COMPLEX"))

# replace TRUE's with column names
interactions_compact <- interactions_wide %>% 
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(receptor_uniprot = paste0(key,collapse=", ")) 
  # went from 5270 rows to 5094 cause some rows contained only complexes

## repeat process with receptor gene symbols

# widen our datatable using reshape2
interactions_genesymbol <- interactions %>% 
  select(uniprot = source,
         target_genesymbol) %>% 
  reshape2::dcast(uniprot ~ target_genesymbol, function(x) any(length(x))) %>% 
  # drop columns containing
  select(-contains("COMPLEX")) %>% 
# replace TRUE's with column names
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(receptor_symbol = paste0(key,collapse=", ")) 


## repeat process with receptor references

# free up RAM
rm(extra_ligands,
   interac_resoucres,
   intercell_filter,
   intercell_final)
gc()

# widen our datatable using reshape2
interactions_references <- interactions %>% 
  select(uniprot = source,
         target_genesymbol,
         references) %>%
  # create new column containing receptor-reference pairs
  mutate(receptor_references = paste0(target_genesymbol, 
                                      " (",
                                      references,
                                      "), ")) %>% 
  # select just the new column
  select(uniprot,
         receptor_references) %>% 
  reshape2::dcast(uniprot ~ receptor_references, function(x) any(length(x))) %>% 
  # drop columns containing
  select(-contains("_")) %>% 
  # replace TRUE's with column names
  gather(key = "key", value = "value", -uniprot) %>%  #Change to long format
  filter(value) %>% # Filter for value which are TRUE
  group_by(uniprot) %>% 
  summarise(receptor_references = paste0(key,collapse=", ")) 


# try another way

# widen our datatable using reshape2
interactions_references <- interactions %>% 
  select(uniprot = source,
         target_genesymbol,
         references) %>%
  # remove rows where receptor symbol contains underscore (b/c means is a complex)
  # wil remove 1837 rows
  dplyr::filter(str_detect(target_genesymbol, "_") != TRUE) %>% 
  # create new column containing receptor-reference pairs
  mutate(receptor_references = paste0(target_genesymbol, 
                                      " [",
                                      references,
                                      "]")) %>% 
  # select just the new column
  select(uniprot,
         receptor_references) %>% 
  group_by(uniprot) %>% 
  summarise(receptor_references = paste0(receptor_references, collapse=", ")) 


### append all our new receptor info columns to ligands datatable
ligands_with_annot <- ligands_with_cat %>% 
  left_join(interactions_compact,
            by = "uniprot") %>% 
  left_join(interactions_genesymbol,
            by = "uniprot") %>% 
  left_join(interactions_references,
            by = "uniprot")
  
# save ligands table
write_tsv(ligands_with_annot, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands.tsv")
View((interactions_genesymbol))

# merge all rows of the same ligand together



#### random stuff below still in progress

## now get annotations for our ligands
get_annotation_resources()

annotation_categories()

ligands_with_cat$uniprot %>% 
translate_ids(organism = 9606,
              uniprot_id = uniprot,
              esnt)

  

annot_df <- import_omnipath_annotations(ligands_with_annot$uniprot,
                            resources = "UniProt_family",
                            wide = TRUE)


View(intercell_ligands)
