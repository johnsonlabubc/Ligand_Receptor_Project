# analyze seurat integration analysis data from Meltem
# compared the eBCs cell group from the lynn lab sorted stem cell dataset to
# the beta-cells group in the lynn lab human islet dataset


# load packages
library(tidyverse)
library(cowplot)
library(ggrepel) 


###### join dataframes of new data and ligands/receptors lists #############

# open integration analysis data
integration_df <- (read.csv("single_cell_analysis/data/diff_exp_markers_based_on_SCT_assay_no_pct.txt", 
                          sep = "\t")) %>% 
  rename(hgnc_symbol = hgnc_gene_id)

# open latest receptors list
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                          sep = "\t"))
# open latest ligands list
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics.tsv", 
                        sep = "\t")) 


# append integration analysis data to main ligand/receptor datatables
receptors_df <- receptors_df %>% 
  left_join(integration_df,
            by = "hgnc_symbol")

ligands_df <- ligands_df %>% 
  left_join(integration_df,
            by = "hgnc_symbol")

############## check number of genes with integration data #####################

# View data
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         Beta,
         avg_log2FC,
         p_val_adj) %>%
  View()

# check number of ligands/receptors that arent NA in the integration
# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
 # filter(!is.na(avg_log2FC)) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         Beta,
         avg_log2FC,
         p_val_adj) %>%
  View()
# 28 ligands without NA's

# receptors
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  filter(!is.na(avg_log2FC)) %>% 
  select(hgnc_symbol,
         sorted_eBCs,
         Beta,
         avg_log2FC,
         p_val_adj) %>%
  View()
# 31 receptors

######################## create volcano plot ############################

# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = avg_log2FC, 
             y = -log10(p_val_adj),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  ggtitle("Ligands") +
  scale_y_continuous(limits = c(0,200)) +
  scale_x_continuous(limits = c(-4,4)) +
  geom_text_repel(data=subset(ligands_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4)), 
                  size = 3,
                  max.overlaps = 18,
                  force = 20,
                  max.time = 10) +
  theme_cowplot()


# receptors
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


# all genes in integration dataset
# this plot will establish baseline for the volcano plot
integration_df%>% 
  ggplot(aes(x = avg_log2FC, 
             y = -log10(p_val_adj),
             label = hgnc_symbol)) +
  geom_point(colour = "gray2", size = 2) +
  labs(x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  ggtitle("All genes") +
  geom_text_repel(data=subset(integration_df %>% 
                                filter(avg_log2FC > 2.5)), 
                  size = 3,
                  max.overlaps = 30,
                  force = 30,
                  max.time = 10) +
  geom_text_repel(data=subset(integration_df %>% 
                                filter(avg_log2FC < -2.5)), 
                  size = 3,
                  max.overlaps = 30,
                  force = 30,
                  max.time = 10) +
  theme_cowplot()


