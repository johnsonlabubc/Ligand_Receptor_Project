# differential expression analysis of the islet & stem cell proteomics data
# data comes from grace in file proteomics_DE_analysis.xlsx
# start with pooled stem cell vs non-diabetic human islet data
# this stem cell data is n=18 with 3 technical replicates, including both 
# the unsorted and sorted stem cell groups

library(tidyverse)
library(cowplot)
library(ggrepel)
library(viridis)


##### load dataframes ###########

# load receptors data
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                          sep = "\t"))

# open new proteomics data from Grace
proteomics_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/proteomics_DE_stemcell_vs_NDislet.txt", 
                           sep = "\t")) %>% 
  rename(uniprot_gn_id = Protein.Ids)

# regenerate old table of uniprot ID's to hgnc symbol conversion
gene_uniprot_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_annot_manual.tsv", 
                              sep = "\t") %>% 
  rename(genecards_url = 17)

# create table with just hgnc symbols and uniprot ids
# will append expression data to this table
gene_uniprot_list <- gene_uniprot_list %>% 
  dplyr::select(hgnc_symbol,
                uniprot_gn_id)

# get max value for each hgnc gene
proteomics_uniprot <- left_join(gene_uniprot_list,
                                proteomics_df,
                                by = "uniprot_gn_id") %>% 
  # drop the uniprot gene id column
  dplyr::select(-uniprot_gn_id) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(stemcell_NDislet_proteomics_FC = max(stemcell_NDislet_proteomics_FC, na.rm = TRUE),
            stemcell_NDislet_proteomics_adj_p_value = max(stemcell_NDislet_proteomics_adj_p_value, na.rm = TRUE)) %>% 
  # convert -inf values to NAs
  mutate_if(is.numeric, list(~na_if(., -Inf)))

# append proteomics DE data to our genes list using hgnc
receptors_df_final <- left_join(receptors_df,
                             proteomics_uniprot,
                             by = "hgnc_symbol")

# save receptors list
write_tsv(receptors_df_final, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics_DE.tsv")

View(receptors_df_final %>% 
       filter(keep_in_list %in% c("Yes", "TBD")) %>% 
       filter(consensus_score > 7) %>% 
       select(hgnc_symbol,
              stemcell_NDislet_proteomics_FC,
              stemcell_NDislet_proteomics_adj_p_value))

###################### create receptors volcano plot #################################

# receptors
receptors_df_final %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 3) +
  labs(x = "stemcell/ND-islets Log2FC",
       y = "-log10(adjusted p-value)") +
  ggtitle("Receptors") +
  scale_x_continuous(limits = c(-1.5, 1.5)) +
  geom_text_repel(data=subset(receptors_df_final %>% 
                              filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                              filter(consensus_score > 7), 
                              log2(stemcell_NDislet_proteomics_FC) < -0.5),
                  size = 3) +
  geom_text_repel(data=subset(receptors_df_final %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 7), 
                              log2(stemcell_NDislet_proteomics_FC) > 0.5),
                  size = 3) +
  theme_cowplot()



############### repeat for ligands ########################################

# load ligands data
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics.tsv", 
                          sep = "\t"))

# open new proteomics data from Grace
proteomics_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/proteomics_DE_stemcell_vs_NDislet.txt", 
                           sep = "\t")) %>% 
  rename(uniprot_gn_id = Protein.Ids)

# regenerate old table of uniprot ID's to hgnc symbol conversion
gene_uniprot_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_annot_manual.tsv", 
                              sep = "\t") %>% 
  rename(genecards_url = 17)

# create table with just hgnc symbols and uniprot ids
# will append expression data to this table
gene_uniprot_list <- gene_uniprot_list %>% 
  dplyr::select(hgnc_symbol,
                uniprot_gn_id)

# get max value for each hgnc gene
proteomics_uniprot <- left_join(gene_uniprot_list,
                                proteomics_df,
                                by = "uniprot_gn_id") %>% 
  # drop the uniprot gene id column
  dplyr::select(-uniprot_gn_id) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(stemcell_NDislet_proteomics_FC = max(stemcell_NDislet_proteomics_FC, na.rm = TRUE),
            stemcell_NDislet_proteomics_adj_p_value = max(stemcell_NDislet_proteomics_adj_p_value, na.rm = TRUE)) %>% 
  # convert -inf values to NAs
  mutate_if(is.numeric, list(~na_if(., -Inf)))

# append proteomics DE data to our genes list using hgnc
ligands_df_final <- left_join(ligands_df,
                                proteomics_uniprot,
                                by = "hgnc_symbol")

# save ligands list
write_tsv(ligands_df_final, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics_DE.tsv")

View(ligands_df_final %>% 
       filter(keep_in_list %in% c("Yes", "TBD")) %>% 
       filter(consensus_score > 4)  %>% 
  select(hgnc_symbol,
         stemcell_NDislet_proteomics_FC,
         stemcell_NDislet_proteomics_adj_p_value))


################ ligands volcano plot ##############################

# ligands
ligands_df_final %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 3) +
  labs(x = "stemcell/ND-islets Log2FC",
       y = "-log10(adjusted p-value)") +
  ggtitle("Ligands") +
  scale_x_continuous(limits = c(-2.5, 2.6)) +
  geom_text_repel(data=subset(ligands_df_final %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4), 
                              log2(stemcell_NDislet_proteomics_FC) < -1),
                  size = 3) +
  geom_text_repel(data=subset(ligands_df_final %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4), 
                              log2(stemcell_NDislet_proteomics_FC) > 0.5),
                  size = 3) +
  theme_cowplot()





# ligands with only insulin labelled
ligands_df_final %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  # create column with log2FC threshold data for colouring the plot
  mutate(threshold = cut(stemcell_NDislet_proteomics_FC, 
                         breaks=c(-Inf,0.5,2,Inf),
                         labels=c("Up in human-islets", 
                                  "no change",
                                  "Up in SC-islets"))) %>%
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = threshold),
             size = 3.5,
             alpha = 1) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  # bracket notation for subscript
  labs(x = bquote(Log["2"]*FC("SC-islets/human-islets")),
       y = bquote(-Log["10"]("Adj. p-value"))) +
  ggtitle("Ligands") +
  scale_x_continuous(limits = c(-2.5, 2.6),
                     breaks = c(-2, 0, 2)) +
  scale_y_continuous(breaks = c(0,5,10)) +
  # add text labels to the colours
  annotate("text",
           x = -1.6, 
           y = 0.1, 
           label = "Up in human islets",
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 1.6, 
           y = 0.1, 
           label = "Up in SC-islets",
           size = 4,
           colour = "#DB5B52") +
  # add text label for insulin
  geom_text_repel(data=subset(ligands_df_final %>% 
                                filter(hgnc_symbol == "INS")),
                  size = 4) +
  theme_bw(base_size = 13) +
  theme(legend.position="none")



######### all proteomics proteins volcano plot #####################

proteomics_df %>% 
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value))) +
  geom_point(colour = "black", size = 2) +
  labs(x = bquote(Log["2"]*FC("SC-islets/human-islets")),
       y = bquote(Log["10"]("Adj. p-value"))) +
  scale_x_continuous(limits = c(-8,8)) +
  ggtitle("All proteins") +
  theme_cowplot()


# with fold change colour labels
proteomics_df %>% 
  # create column with log2FC threshold data for colouring the plot
  mutate(threshold = cut(stemcell_NDislet_proteomics_FC, 
                         breaks=c(-Inf,0.5,2,Inf),
                         labels=c("Up in human-islets", 
                                  "no change",
                                  "Up in SC-islets"))) %>%
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value))) +
  # colour points based on the threshold
  geom_point(aes(colour = threshold),
             size = 2,
             alpha = 0.8) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  labs(x = bquote(Log["2"]*FC("SC-islets/human-islets")),
       y = bquote(-Log["10"]("Adj. p-value"))) +
  scale_x_continuous(limits = c(-8,8)) +
  ggtitle("All proteins") +
  # add text labels to the colours
  annotate("text",
           x = -5.5, 
           y = -0.5, 
           label = "Up in human islets",
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 5.9, 
           y = -0.5, 
           label = "Up in SC-islets",
           size = 4,
           colour = "#DB5B52") +
  theme_bw(base_size = 13) +
  theme(legend.position="none")

