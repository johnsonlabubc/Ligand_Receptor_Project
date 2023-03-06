# try to get more annotations from omnipath to group ligands and receptors 
# based on signaling pathway family


library(OmnipathR)
library(tidyverse)


# load final ligands and receptors
ligands <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligand_scores_3.tsv", 
                     sep = "\t")) %>% 
  select(hgnc_symbol) %>% 
  mutate(type = "ligand")

receptors <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_scores_3.tsv", 
                       sep = "\t")) %>% 
  select(hgnc_symbol) %>% 
  mutate(type = "receptor")


# join ligands and receptors together
genes_df <- ligands %>% 
  full_join(receptors)


# check what annotation resources are available in omnipath
get_annotation_resources()


# get list of just HGNC symbols
# split into 2 because omnipath only allows 500 at a time
genes_list1 <- head(genes_df$hgnc_symbol,500)
genes_list2 <- tail(genes_df$hgnc_symbol,271)

# confirm that didnt miss any or have duplicates
identical(c(genes_list1, genes_list2), genes_df$hgnc_symbol)


# get annotations part 1
omnipath_annots1 <- import_omnipath_annotations(proteins = genes_list1,
                            resources = "CellChatDB",
                            wide = FALSE) %>%
  # keep only the useful annotations
  filter(label == "pathway")

# get annotations part 2
omnipath_annots2 <- import_omnipath_annotations(proteins = genes_list2,
                                                resources = "CellChatDB",
                                                wide = FALSE) %>%
  # keep only the useful annotations
  filter(label == "pathway")


# join annotations together
omnipath_annots <- omnipath_annots1 %>% 
  full_join(omnipath_annots2) %>% 
  rename(hgnc_symbol = genesymbol) %>% 
  select(hgnc_symbol,
         CellChat_Pathway = value)


# join back to genes_df
genes_df_annot <- genes_df %>% 
  left_join(omnipath_annots)



# summarize receptor annotations
rec_annots <- genes_df_annot %>% 
  filter(type == "receptor") %>% 
  group_by(CellChat_Pathway) %>% 
  summarise(n()) %>% 
  arrange(desc(`n()`))


# save annotations
write_tsv(genes_df_annot, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_cellchat_pathway_annot.tsv")







  
  