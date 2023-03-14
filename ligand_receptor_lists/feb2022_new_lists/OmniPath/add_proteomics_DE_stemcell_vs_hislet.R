# differential expression analysis of the islet & stem cell proteomics data
# data comes from grace in file proteomics_DE_analysis.xlsx
# start with pooled stem cell vs non-diabetic human islet data
# this stem cell data is n=18 with 3 technical replicates, including both 
# the unsorted and sorted stem cell groups

library(tidyverse)
library(cowplot)
library(ggrepel)


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


### receptors plot with colour to indicate proteins that significantly 
# differentially expressed, with both log2FC > +-0.5 & adj p value < 0.05

### using log2FC > +-0.5 for receptors, and > +-1 for ligands, because there 
# are very few hits for the receptors, and there's also an important DE'd receptor
# IGF2R that is just below the cutoff
### note - even tho i wrote this comment above, it looks like in code and figure
### the receptors and ligands are both still using log2FC > +/- 1
### and instead just manually labelled IGF2R


# first create new columns with data needed for threshold colouring of plot
receptors_df_final <- receptors_df_final %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  # create column with log2FC threshold data 
  mutate(log2FC_threshold = cut(stemcell_NDislet_proteomics_FC, 
                                # these values are before taking log2
                                # using 
                                breaks=c(-Inf,0.5,2,Inf), 
                                labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(stemcell_NDislet_proteomics_adj_p_value,
                                   breaks = c(-Inf,0.05,Inf),
                                   labels = c(1,0))) %>%
  # create column that combines both thresholds by multiplying
  # have to convert factor columns to numeric
  mutate(sig_change_threshold = as.numeric(as.character(log2FC_threshold)) 
         * as.numeric(as.character(adjpvalue_threshold)),
         sig_change_threshold = cut(sig_change_threshold,
                                    breaks=c(-Inf,-0.5,0.5,Inf), 
                                    labels=c("Up in human islets",
                                             "no change",
                                             "Up in SC-islets")))

# plot
receptors_df_final %>%   
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
             size = 3.5,
             alpha = 1) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  # bracket notation for subscript
  labs(x = bquote(Log["2"]*FC("SC-islets/human-islets")),
       y = bquote(-Log["10"]("Adj. p-value"))) +
  ggtitle("Receptors") +
  scale_x_continuous(limits = c(-2.5, 2.6),
                     breaks = c(-2, 0, 2)) +
  scale_y_continuous(breaks = c(0,5,10,15)) +
  # add text labels to the colours
  annotate("text",
           x = -1.6, 
           y = 0.1, 
           label = "Down in SC-islets",
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 1.6, 
           y = 0.1, 
           label = "Up in SC-islets",
           size = 4,
           colour = "#DB5B52") +
  geom_text_repel(data=receptors_df_final %>% 
                    filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                    filter(consensus_score > 4) %>% 
                    filter(sig_change_threshold %in% c("Up in human islets",
                                                       "Up in SC-islets")),
                  size = 3) +
  # add text label for IGF2R b/c it's interesting but just below log2FC threshold
  geom_text_repel(data=subset(receptors_df_final %>% 
                                filter(hgnc_symbol == "IGF2R")),
                  size = 3) +
  theme_bw(base_size = 13) +
  theme(legend.position="none") 

# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/proteomics_volcano_plots/receptors_proteomics_volcano_stemcell_vs_NDislets_3.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)





############### prep ligands data  ########################################

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



### ligands plot with colour to indicate proteins that significantly 
# differentially expressed, with both log2FC > +-1 & adj p value < 0.05


# first create new columns with data needed for threshold colouring of plot
ligands_df_final <- ligands_df_final %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  # create column with log2FC threshold data 
  mutate(log2FC_threshold = cut(stemcell_NDislet_proteomics_FC, 
                                breaks=c(-Inf,0.5,2,Inf), # these values are before taking log2
                                labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(stemcell_NDislet_proteomics_adj_p_value,
                                   breaks = c(-Inf,0.05,Inf),
                                   labels = c(1,0))) %>%
  # create column that combines both thresholds by multiplying
  # have to convert factor columns to numeric
  mutate(sig_change_threshold = as.numeric(as.character(log2FC_threshold)) 
         * as.numeric(as.character(adjpvalue_threshold)),
         sig_change_threshold = cut(sig_change_threshold,
                                    breaks=c(-Inf,-0.5,0.5,Inf), 
                                    labels=c("Up in human islets",
                                             "no change",
                                             "Up in SC-islets")))
  
  # plot
ligands_df_final %>%   
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
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
           label = "Down in SC-islets",
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 1.6, 
           y = 0.1, 
           label = "Up in SC-islets",
           size = 4,
           colour = "#DB5B52") +
  geom_text_repel(data=ligands_df_final %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4) %>% 
                                filter(sig_change_threshold %in% c("Up in human islets",
                                                                   "Up in SC-islets")),
                  size = 3) +
  theme_bw(base_size = 13) +
  theme(legend.position="none") 

# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_proteomics_volcano_stemcell_vs_NDislets_3.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)





############## ligands volcano with only insulin labelled ###############
ligands_df_final %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  # create column with log2FC threshold data 
  mutate(log2FC_threshold = cut(stemcell_NDislet_proteomics_FC, 
                         breaks=c(-Inf,0.5,2,Inf), # these values are before taking log2
                         labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(stemcell_NDislet_proteomics_adj_p_value,
                                   breaks = c(-Inf,0.05,Inf),
                                   labels = c(1,0))) %>%
  # create column that combines both thresholds by multiplying
  # have to convert factor columns to numeric
  mutate(sig_change_threshold = as.numeric(as.character(log2FC_threshold)) 
         * as.numeric(as.character(adjpvalue_threshold)),
         sig_change_threshold = cut(sig_change_threshold,
                                      breaks=c(-Inf,-0.5,0.5,Inf), 
                                      labels=c("Up in human islets",
                                               "no change",
                                               "Up in SC-islets"))) %>% 
  
  # plot
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
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

# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_proteomics_volcano_stemcell_vs_NDislets_INSlabel_3.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)



####################### all proteomics proteins volcano plot ################################

proteomics_df %>% 
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value))) +
  geom_point(colour = "black", size = 2) +
  labs(x = bquote(Log["2"]*FC("SC-islets/human-islets")),
       y = bquote(Log["10"]("Adj. p-value"))) +
  scale_x_continuous(limits = c(-8,8)) +
  ggtitle("All proteins") +
  theme_cowplot()


####### with fold change & p value colour labels
proteomics_df %>% 
  # create column with log2FC threshold data 
  mutate(log2FC_threshold = cut(stemcell_NDislet_proteomics_FC, 
                                breaks=c(-Inf,0.5,2,Inf), # these values are before taking log2
                                labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(stemcell_NDislet_proteomics_adj_p_value,
                                   breaks = c(-Inf,0.05,Inf),
                                   labels = c(1,0))) %>%
  # create column that combines both thresholds by multiplying
  # have to convert factor columns to numeric
  mutate(sig_change_threshold = as.numeric(as.character(log2FC_threshold)) 
         * as.numeric(as.character(adjpvalue_threshold)),
         sig_change_threshold = cut(sig_change_threshold,
                                    breaks=c(-Inf,-0.5,0.5,Inf), 
                                    labels=c("Up in human islets",
                                             "no change",
                                             "Up in SC-islets"))) %>% 
  # plot
  ggplot(aes(x = log2(stemcell_NDislet_proteomics_FC), 
             y = -log10(stemcell_NDislet_proteomics_adj_p_value),
             label = Genes)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
             size = 1.5,
             alpha = 0.9) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  labs(x = bquote(Log["2"]*FC("SC-islets/human-islets")),
       y = bquote(-Log["10"]("Adj. p-value"))) +
  scale_x_continuous(limits = c(-8,8)) +
  ggtitle("All proteins") +
  # add text labels to the colours
  annotate("text",
           x = -5.2, 
           y = -0.5, 
           label = "Down in SC-islets",
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 5.9, 
           y = -0.5, 
           label = "Up in SC-islets",
           size = 4,
           colour = "#DB5B52") +
  # add text label for IGF2BP2
  geom_text_repel(data=subset(proteomics_df %>% 
                                filter(Genes == "IGF2BP2")),
                  size = 3) +
  theme_bw(base_size = 13) +
  theme(legend.position="none")
  
  
  # save the plot
  ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/allproteins_proteomics_volcano_stemcell_vs_NDislets_IGF2BP_2.JPG",
         device = "jpg",
         width = 1500,
         height = 1500,
         units = "px",
         scale = 0.8)


  
#################  number of differentially abundant proteins #####################

# all proteins  
proteomics_df %>% 
    # create column with log2FC threshold data 
    mutate(log2FC_threshold = cut(stemcell_NDislet_proteomics_FC, 
                                  breaks=c(-Inf,0.5,2,Inf), # these values are before taking log2
                                  labels=c(-1,0,1))) %>%
    # create column with adj p value threshold
    mutate(adjpvalue_threshold = cut(stemcell_NDislet_proteomics_adj_p_value,
                                     breaks = c(-Inf,0.05,Inf),
                                     labels = c(1,0))) %>%
    # create column that combines both thresholds by multiplying
    # have to convert factor columns to numeric
    mutate(sig_change_threshold = as.numeric(as.character(log2FC_threshold)) 
           * as.numeric(as.character(adjpvalue_threshold)),
           sig_change_threshold = cut(sig_change_threshold,
                                      breaks=c(-Inf,-0.5,0.5,Inf), 
                                      labels=c("Up in human islets",
                                               "no change",
                                               "Up in SC-islets"))) %>% 
    group_by(sig_change_threshold) %>% 
    summarise(n())
  