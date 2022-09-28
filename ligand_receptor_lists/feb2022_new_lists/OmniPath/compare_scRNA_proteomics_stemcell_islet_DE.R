# compare scRNAseq & proteomics differential expression/abundance between
# stem cell and human islet/beta cell


# load libraries
library(tidyverse)
library(cowplot)
library(ggrepel)


########################## load data ##################################

# open latest receptors list
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_scRNA_DE.tsv", 
                          sep = "\t"))
# open latest ligands list
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_scRNA_DE.tsv", 
                        sep = "\t"))


####### compare correlation of proteomics to scRNA-seq fold changes ##########

# ligands
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = eBC_beta_scRNA_log2FC, 
             y = log2(stemcell_NDislet_proteomics_FC),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 4) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_x_continuous(limits = c(0,1)) +
  labs(x = "eBCs/beta cell scRNA-seq Log2FC",
       y = "SCs/islet proteomics Log2FC") +
  ggtitle("Receptors") +
  geom_text_repel(size = 3) +
  theme_cowplot()

# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = eBC_beta_scRNA_log2FC, 
             y = log2(stemcell_NDislet_proteomics_FC),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 4) +
  labs(x = "eBCs/beta cell scRNA-seq Log2FC",
       y = "SCs/islet proteomics Log2FC") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3) +
  theme_cowplot()




#################### for all genes/proteins ##############################

## open scRNA-seq integration analysis data
scRNA_DE_df <- (read.csv("single_cell_analysis/data/diff_exp_markers_based_on_SCT_assay_no_pct.txt", 
                            sep = "\t")) %>% 
  rename(hgnc_symbol = hgnc_gene_id) %>% 
  # also change the direction of the Fold change ratio for consistency
  # with the proteomics fold change, which is stem cell over human islet
  # data from meltem is human beta cells over eBCs
  # take the negative to invert the direction of the fold change
  mutate(inverse_avg_log2FC = -avg_log2FC) %>% 
  dplyr::select(hgnc_symbol,
         scRNA_log2FC = inverse_avg_log2FC,
         scRNA_p_val_adj = p_val_adj)


# open proteomics integration data
proteomics_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/proteomics_DE_stemcell_vs_NDislet.txt", 
                           sep = "\t")) %>% 
  rename(uniprot_gn_id = Protein.Ids,
         hgnc_symbol = Genes)


# append proteomics DE & scRNAseq DE together
integration_df <- left_join(scRNA_DE_df,
                            proteomics_df,
                            by = "hgnc_symbol")

# all proteins plot    
integration_df %>% 
  ggplot(aes(x = scRNA_log2FC, 
             y = log2(stemcell_NDislet_proteomics_FC),
             label = hgnc_symbol)) +
  geom_point(colour = "black", size = 2) +
  scale_x_continuous(limits = c(-5,5)) +
  labs(x = bquote(Log["2"]*FC("eBCs/beta cell scRNA-seq")),
       y = "SCs/islet proteomics Log2FC") +
  ggtitle("All genes") +
  theme_cowplot()

