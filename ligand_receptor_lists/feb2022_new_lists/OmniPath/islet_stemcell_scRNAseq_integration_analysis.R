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



######################## create volcano plot ############################

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
integration_df %>% 
  ggplot(aes(x = avg_log2FC, 
             y = -log10(p_val_adj),
             label = hgnc_symbol)) +
  geom_point(colour = "gray2", size = 1) +
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



