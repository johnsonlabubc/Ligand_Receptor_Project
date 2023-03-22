# analyze seurat integration analysis data from Meltem
# compared the eBCs cell group from the lynn lab sorted stem cell dataset to
# the beta-cells group in the lynn lab human islet dataset


# load packages
library(tidyverse)
library(cowplot)
library(ggrepel) 


###################### load data ####################################



# open integration analysis data
integration_df <- (read.csv("single_cell_analysis/data/markers_humanislet_vs_enrichedbeta_minpct0_lfc0.txt", 
                          sep = "\t")) %>% 
  rename(hgnc_symbol = hgnc_gene_id) %>% 
  # also change the direction of the Fold change ratio for consistency
  # with the proteomics fold change, which is stem cell over human islet
  # data from meltem is human beta cells over eBCs
  # take the negative to invert the direction of the fold change
  mutate(inverse_avg_log2FC = -avg_log2FC) %>% 
  # replace absolute 0 adj p values with 1e-300 for graphical purposes
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj))


# open latest receptors list
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics_DE.tsv", 
                          sep = "\t"))
# open latest ligands list
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics_DE.tsv", 
                        sep = "\t")) 



######## append integration data to ligand/receptor lists ####################

# append to receptors
receptors_df <- receptors_df %>% 
  left_join(integration_df,
            by = "hgnc_symbol") %>% 
  # only keep log2FC and adj p value columns
  select(-avg_log2FC,
         -p_val,
         -pct.1,
         -pct.2) %>% 
  # rename columns for clarity in the main datatables
  rename(eBC_beta_scRNA_log2FC = inverse_avg_log2FC,
         eBC_beta_scRNA_adj_p_value = p_val_adj) %>% 
  # reorder columns for consistency with proteomics 
  relocate(eBC_beta_scRNA_adj_p_value,
           .after = eBC_beta_scRNA_log2FC)


# append to ligands
ligands_df <- ligands_df %>% 
  left_join(integration_df,
            by = "hgnc_symbol") %>% 
# only keep log2FC and adj p value columns
select(-avg_log2FC,
       -p_val,
       -pct.1,
       -pct.2) %>% 
  # rename columns for clarity in the main datatables
  rename(eBC_beta_scRNA_log2FC = inverse_avg_log2FC,
         eBC_beta_scRNA_adj_p_value = p_val_adj) %>% 
  # reorder columns for consistency with proteomics 
  relocate(eBC_beta_scRNA_adj_p_value,
           .after = eBC_beta_scRNA_log2FC)


# save new datatables
write_tsv(ligands_df,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_scRNA_DE.tsv")

write_tsv(receptors_df,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_scRNA_DE.tsv")


############## check number of genes with integration data #####################

# View data
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         Beta,
         eBC_beta_scRNA_log2FC,
         eBC_beta_scRNA_adj_p_value) %>%
  View()

# check number of ligands/receptors that arent NA in the integration
# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  filter(!is.na(eBC_beta_scRNA_log2FC)) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         Beta,
         eBC_beta_scRNA_log2FC,
         eBC_beta_scRNA_adj_p_value) %>%
  View()
# 202 ligands without NA's

# receptors
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  filter(!is.na(eBC_beta_scRNA_log2FC)) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         Beta,
         eBC_beta_scRNA_log2FC,
         eBC_beta_scRNA_adj_p_value) %>%
  View()
# 231 receptors



######################## ligands volcano plot ############################

# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = eBC_beta_scRNA_log2FC, 
             y = -log10(eBC_beta_scRNA_adj_p_value),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  ggtitle("Ligands") +
  scale_y_continuous(limits = c(0,200)) +
  scale_x_continuous(limits = c(-4,4)) +
  geom_text_repel(data=subset(ligands_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4) %>% 
                                filter(eBC_beta_scRNA_log2FC < -1)), 
                  size = 3) +
  theme_cowplot()






ligands_df_2 <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  # create column with log2FC threshold data 
  mutate(log2FC_threshold = cut(eBC_beta_scRNA_log2FC, 
                                #  based on log2 transformed breakpoints
                                breaks=c(-Inf,-1,1,Inf), 
                                labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(eBC_beta_scRNA_adj_p_value,
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

#figure out which ligands are the 2 others to label
ligands_df_2 %>% 
  select(hgnc_symbol,
         eBC_beta_scRNA_log2FC) %>% 
  View()

  # plot
ligands_df_2 %>% 
  ggplot(aes(x = eBC_beta_scRNA_log2FC, 
             y = -log10(eBC_beta_scRNA_adj_p_value),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
             size = 2.5,
             alpha = 0.7) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  labs(x = bquote(Log["2"]*FC(paste("SC", beta, "-cells/human ", beta, "-cells"))),
       y = bquote(-Log["10"]("Adj. p-value"))) +
  scale_x_continuous(limits = c(-4,4),
                     breaks = c(-3,0,3)) +
#  scale_y_continuous(limits = c(0,320)) +
  ggtitle("Ligands") +
  # add text labels to the colours
  annotate("text",
           x = -2.3, 
           y = -0.5, 
           label = expression(paste("Down in SC", beta, "-cells")),
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 2.5, 
           y = -0.5, 
           label = expression(paste("Up in SC", beta, "-cells")),
           size = 4,
           colour = "#DB5B52") +
  geom_text_repel(data=ligands_df_2 %>% 
                    filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                    filter(consensus_score > 4) %>% 
                    filter(sig_change_threshold %in% c("Up in human islets")),
                  size = 3) +
  geom_text_repel(data=ligands_df_2 %>% 
                    filter(hgnc_symbol %in% c("ADCYAP1",
                                              "GAST",
                                              "BMP5",
                                              "PYY",
                                              "DLK1")),
                  size = 3,
                  force_pull = 1,
                  max.overlaps = 5,
                  max.time = 3,
                  force = 5,
                  nudge_y = -1,
                  nudge_x = -.05) +
  theme_bw(base_size = 13) +
  theme(legend.position="none")


# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/scRNAseq_volcano_plots/ligands_scRNAseq_volcano_eBC_vs_betacell_2.png",
       scale = 1.5)





########################## receptors volcano plot ####################################
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = avg_log2FC, 
             y = -log10(p_val_adj),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  ggtitle("Receptors") +
  scale_y_continuous(limits = c(0,200)) +
  scale_x_continuous(limits = c(-2,2)) +
  geom_text_repel(data=subset(receptors_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 7)), 
                  size = 3,
                  max.overlaps = 30,
                  force = 30,
                  max.time = 10) +
  theme_cowplot()






receptors_df_2 <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  # create column with log2FC threshold data 
  mutate(log2FC_threshold = cut(eBC_beta_scRNA_log2FC, 
                                #  based on log2 transformed breakpoints
                                breaks=c(-Inf,-1,1,Inf), 
                                labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(eBC_beta_scRNA_adj_p_value,
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
receptors_df_2 %>% 
  ggplot(aes(x = eBC_beta_scRNA_log2FC, 
             y = -log10(eBC_beta_scRNA_adj_p_value),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
             size = 2.5,
             alpha = 0.7) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  labs(x = bquote(Log["2"]*FC(paste("SC", beta, "-cells/human ", beta, "-cells"))),
       y = bquote(-Log["10"]("Adj. p-value"))) +
  scale_x_continuous(limits = c(-2,2),
                     breaks = c(-2,0,2)) +
  scale_y_continuous(limits = c(0,320)) +
  ggtitle("Receptors") +
  # add text labels to the colours
  annotate("text",
           x = -1.2, 
           y = -0, 
           label = expression(paste("Down in SC", beta, "-cells")),
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 1.2, 
           y = -0, 
           label = expression(paste("Up in SC", beta, "-cells")),
           size = 4,
           colour = "#DB5B52") +
  geom_text_repel(data=receptors_df_2 %>% 
                    filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                    filter(consensus_score > 7) %>% 
                    filter(sig_change_threshold %in% c("Up in human islets",
                                                       "Up in SC-islets")),
                  size = 3,
                  force_pull = 1,
                  max.overlaps = 10,
                  max.time = 3,
                  force = 20) +
  geom_text_repel(data=receptors_df_2 %>% 
                    filter(hgnc_symbol %in% c("PLXNC1",
                                              "PLXNA2",
                                              "SSTR2",
                                              "NRXN1",
                                              "GLP1R",
                                              "GIPR",
                                              "FZD3")),
                  size = 3,
                  force_pull = 5,
                  max.overlaps = 10,
                  max.time = 3,
                  force = 20) +
  theme_bw(base_size = 13) +
  theme(legend.position="none")


# save the plot
# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/scRNAseq_volcano_plots/receptors_scRNAseq_volcano_eBC_vs_betacell_2.png",
       scale = 1.5)




# ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/scRNAseq_volcano_plots/receptors_scRNAseq_volcano_eBC_vs_betacell_log2fc0_5.JPG",
#        device = "jpg",
#        width = 1500,
#        height = 1500,
#        units = "px",
#        scale = 0.8)



########################## all genes volanco plot ####################################


integration_df %>% 
# create column with log2FC threshold data 
mutate(log2FC_threshold = cut(inverse_avg_log2FC, 
                              #  based on log2 transformed breakpoints
                              breaks=c(-Inf,-1,1,Inf), 
                              labels=c(-1,0,1))) %>%
  # create column with adj p value threshold
  mutate(adjpvalue_threshold = cut(p_val_adj,
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
  ggplot(aes(x = inverse_avg_log2FC, 
             y = -log10(p_val_adj),
             label = hgnc_symbol)) +
  # colour points based on the threshold
  geom_point(aes(colour = sig_change_threshold),
             size = 1.4,
             alpha = 0.6) +
  scale_colour_manual(values = c("#396FCB", "gray60", "#DB5B52")) +
  labs(x = bquote(Log["2"]*FC(paste("SC", beta, "-cells/human ", beta, "-cells"))),
       y = bquote(-Log["10"]("Adj. p-value"))) +
 # scale_x_continuous(limits = c(-4,6)) +
  ggtitle("All genes") +
  # add text labels to the colours
  annotate("text",
           x = -5.5, 
           y = -0.5, 
           label = expression(paste("Down in SC", beta, "-cells")),
           size = 4,
           colour = "#396FCB") +
  annotate("text",
           x = 3.9, 
           y = -0.5, 
           label = expression(paste("Up in SC", beta, "-cells")),
           size = 4,
           colour = "#DB5B52") +
  # add text label for genes that really stand out
  geom_text_repel(data=subset(integration_df %>% 
                                filter(hgnc_symbol %in% c("MALAT1",
                                                          "INS", 
                                                          "IAPP",
                                                          "eGFP",
                                                          "CHGA"))),
                  size = 3,
                  nudge_y = -3,
                  force = 1,
                  force_pull = 1) +
  theme_bw(base_size = 13) +
  theme(legend.position="none")

# save the plot
# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/scRNAseq_volcano_plots/allgenes_scRNAseq_volcano_eBC_vs_betacell.png",
       scale = 1.5)



# ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/scRNAseq_volcano_plots/allgenes_scRNAseq_volcano_eBC_vs_betacell.JPG",
#        device = "jpg",
#        width = 1500,
#        height = 1500,
#        units = "px",
#        scale = 0.8)
