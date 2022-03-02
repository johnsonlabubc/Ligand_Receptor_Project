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
interactions_extra <- import_ligrecextra_interactions(resources=c("iTALK",
                                                            "Baccin2019"), organism=9606)

View(interactions_extra)
# can see there are 3013 extra ligand-receptor interactions
# a lot of them are found from multiple references, including CellChat, connectomedb, etc


# save this list for future reference
write_tsv(interactions_extra, "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrecextra_interactions.tsv")


## We check some of the different intercell categories
get_intercell_generic_categories()


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
                                               resources = "OmniPath",
                                               secreted = TRUE) %>% 
                     filter(entity_type == "protein")
# 1479 ligands, including complexes

nrow(intercell_ligands)
# 897 lgands, not including complexes

intercell_ligands %>% 
  group_by(plasma_membrane_transmembrane) %>% 
  summarise(n())
# of the 897, 120 have a plasma memberane transmembrane
# will keep them cause they are also defined as secreted out of the cell

#confirm that all ligands are transmitters too
as.data.frame(intercell_ligands$transmitter) %>% 
  distinct()
# confirmed that all ligands are transmitters, are secreted, and are not receivers



intercell_filter <- intercell_consensus_filter(intercell,
                                               percentile = )

View(intercell_ligands)
