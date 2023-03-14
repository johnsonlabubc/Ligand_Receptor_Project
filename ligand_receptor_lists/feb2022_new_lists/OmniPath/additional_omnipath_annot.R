# try to get more annotations from omnipath to group ligands and receptors 
# based on signaling pathway family


library(OmnipathR)
library(tidyverse)
library(cowplot)


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


############# plot pathway data and consensus score in lists ##########################

# reopen
genes_df_annot <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_cellchat_pathway_annot.tsv", 
                     sep = "\t")) 

# load final ligands and receptors and get just consensus scores
ligands_pathway_scores <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_meltonbalboa.tsv", 
                     sep = "\t")) %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         consensus_score) %>% 
  # get mean consensus scores
  left_join(genes_df_annot,
             by = "hgnc_symbol") %>% 
  group_by(CellChat_Pathway) %>% 
  summarize(mean(consensus_score)) %>%
  drop_na(CellChat_Pathway) %>% 
  mutate(type = "ligand")
  

receptors_pathway_scores <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_meltonbalboa.tsv", 
                       sep = "\t")) %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         consensus_score) %>% 
  # get mean consensus scores
  left_join(genes_df_annot,
            by = "hgnc_symbol") %>% 
  group_by(CellChat_Pathway) %>% 
  summarize(mean(consensus_score)) %>% 
  drop_na(CellChat_Pathway) %>% 
  mutate(type = "receptor")


# join lists together
ligrec_pathway_scores <- full_join(ligands_pathway_scores,
                                   receptors_pathway_scores)
  

### get number of ligands/rec per pathway
  
# summarize receptor annotations
rec_count <- genes_df_annot %>% 
  filter(type == "receptor") %>% 
  group_by(CellChat_Pathway) %>% 
  summarise(n()) %>%
  drop_na(CellChat_Pathway) %>% 
  mutate(type = "receptor") %>% 
  arrange(desc(`n()`)) %>% 
  head(25)

# summarize receptor annotations
lig_count <- genes_df_annot %>% 
  filter(type == "ligand") %>% 
  group_by(CellChat_Pathway) %>% 
  summarise(n()) %>%
  drop_na(CellChat_Pathway) %>% 
  mutate(type = "ligand") %>% 
  arrange(desc(`n()`)) %>% 
  head(25)


# join lig and rec for plotting
ligrec_count <- full_join(rec_count,
                          lig_count)

# join counts and mean consensus score tables
pathways_df <- left_join(ligrec_count,
                         ligrec_pathway_scores,
                         by = c("CellChat_Pathway", "type"))


### plot with facet wrap

#ligands
pathways_df %>% 
  filter(type == "ligand") %>% 
  ggplot(aes(x = 1,y = fct_reorder(CellChat_Pathway,`n()`), group = type )) +
  geom_point(aes(size = `n()`,
                 color = `mean(consensus_score)`)) +
  facet_wrap("type",
             scales = "free") +
  scale_size_continuous(range = c(2, 7)) +
  scale_color_distiller(palette = "BuPu",
                        direction = 1) +
  theme_cowplot() +
  # reverse alphabetical orderS
#  scale_y_discrete(limits = rev) +
 # scale_x_discrete(labels = c("Johnson 2023", "Lund 2014")) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.box.spacing = unit(0, "cm"),
        axis.text.y = element_text(size = 9,
                                   colour = "black"),
        legend.text = element_text(size = 8,
                                   colour = "black"),
        legend.title = element_text(size = 10,
                                    colour = "black"),
        legend.key.size = unit(0.5, "cm")) +
  labs(colour = "Mean consensus score",
       size = "Number of ligands") +
  # order the 2 legends
  guides(color = guide_colorbar(order = 0),
         size = guide_legend(order = 1))
# save 
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligand_pathways_consensus_scores.png",
       scale = 1.7)


#receptrors
pathways_df %>% 
  filter(type == "receptor") %>% 
  ggplot(aes(x = 1,y = fct_reorder(CellChat_Pathway,`n()`), group = type )) +
  geom_point(aes(size = `n()`,
                 color = `mean(consensus_score)`)) +
  facet_wrap("type",
             scales = "free") +
  scale_size_continuous(range = c(2, 7)) +
  scale_color_distiller(palette = "BuGn",
                        direction = 1) +
  theme_cowplot() +
  # reverse alphabetical orderS
  #  scale_y_discrete(limits = rev) +
  # scale_x_discrete(labels = c("Johnson 2023", "Lund 2014")) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.box.spacing = unit(0, "cm"),
        axis.text.y = element_text(size = 9,
                                   colour = "black"),
        legend.text = element_text(size = 8,
                                   colour = "black"),
        legend.title = element_text(size = 10,
                                    colour = "black"),
        legend.key.size = unit(0.5, "cm")) +
  labs(colour = "Mean consensus score",
       size = "Number of receptors") +
  # order the 2 legends
  guides(color = guide_colorbar(order = 0),
         size = guide_legend(order = 1))
# save 
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptor_pathways_consensus_scores.png",
       scale = 1.7)


  