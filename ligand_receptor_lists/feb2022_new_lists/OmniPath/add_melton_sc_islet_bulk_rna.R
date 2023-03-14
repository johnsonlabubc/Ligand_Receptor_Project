# add bulk SC-islet RNA-seq data from melton 2020 cell stem cell paper
# to our ligand/receptor data

library(tidyverse)
library(biomaRt)



############# load data ###############################


receptors_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_scRNA_DE.tsv", 
                         sep = "\t")

ligands_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_scRNA_DE.tsv", 
                         sep = "\t")


sc_islet_bulk_rna_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/melton-sc-beta-bulk-RNA-seq/GSM4146881_SCbeta_RNA-seq_rep1_HTseq_byGencode.v19.htseq-count.out.txt", 
                              sep = "\t",
                              col.names = c("ensembl_gene_id",
                                            "sc_islet_bulk_rna"))

h_islet_bulk_rna_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/melton-sc-beta-bulk-RNA-seq/GSM4146884_Beta_RNA-seq_rep1_HTseq_byGencode.v19.htseq-count.out.txt", 
                                 sep = "\t",
                                 col.names = c("ensembl_gene_id",
                                               "h_islet_bulk_rna"))



######################### tidy data #################################

# trim version number from ensembl gene ids
sc_islet_bulk_rna_df$ensembl_gene_id <- gsub('\\..+$', '', sc_islet_bulk_rna_df$ensembl_gene_id)

h_islet_bulk_rna_df$ensembl_gene_id <- gsub('\\..+$', '', h_islet_bulk_rna_df$ensembl_gene_id)

# check for duplicate gene ids
sc_islet_bulk_rna_df %>% 
  group_by(ensembl_gene_id) %>% 
  summarise(n())
# all 57,819 ensembl gene ids are distinct


##### convert ensembl IDs

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes_df <- getBM(attributes = 
                    c("ensembl_gene_id",
                      "hgnc_symbol"),
                  filters=c('ensembl_gene_id'),
                  values =  sc_islet_bulk_rna_df[1],
                  mart = mart)

sc_islet_bulk_rna_df <- left_join(sc_islet_bulk_rna_df, 
                         genes_df, 
                         by = "ensembl_gene_id") %>% 
  left_join(h_islet_bulk_rna_df,
            by = "ensembl_gene_id")

# approx 500 of the 57819 ensembl gene IDs did not have a corresponding hgnc gene id
# seems like these were all rna genes, non-coding, etc, so its not a problem for us

#################### join to receptors & ligands lists ##############################

receptors_df2 <- left_join(receptors_df, 
                          sc_islet_bulk_rna_df[2:4], 
                               by = "hgnc_symbol") %>% 
  # for clarity, rename islet_tpm to lund_islet_tpm
  rename(Lund_h_islet_bulkRNA = islet_tpm,
         Melton_h_islet_bulkRNA = h_islet_bulk_rna)


ligands_df2 <- ligands_df %>% 
  # fix bug caused by gene with blank hgnc_symbol in ligands list
  # by removing that row from the dataframe
  filter(hgnc_symbol != "") %>% 
  # now join after fixing bug
  left_join(sc_islet_bulk_rna_df[2:4], 
                           by = "hgnc_symbol") %>% 
  # fix bug caused by gene with blank hgnc_symbol in ligands list
  # for clarity, rename islet_tpm to lund_islet_tpm
  rename(Lund_h_islet_bulkRNA = islet_tpm,
         Melton_h_islet_bulkRNA = h_islet_bulk_rna)


# check top ligands and receptors

receptors_df2 %>% 
  dplyr::select(hgnc_symbol,
         sc_islet_bulk_rna,
         Melton_h_islet_bulkRNA) %>% 
  View()


############ fold change vs HI-islet bulk RNA-seq ###################

# a decent amount of genes have very high log2fc's b/c they were only 
# detected in 1 of the 2 rna-seq experiments
# so worth doing a more direct comparison to HI-islet data from the same
# melton paper, even though it has fewer n's
# b/c melton comparison is same methods, processing, etc


# fold change compared to Lund human islet RNA-seq
# lowest value in data is 1, so using plus 0.1 as reasonable value to prevent -inf
receptors_df2 <- receptors_df2 %>% 
  mutate(SCislet_Hislet_bulkRNA_Lund_log2FC = log2((sc_islet_bulk_rna + 0.1) / 
                                                     (Lund_h_islet_bulkRNA + 0.1)),
         SCislet_Hislet_bulkRNA_Melton_log2FC = log2((sc_islet_bulk_rna + 0.1) / 
                                                       (Melton_h_islet_bulkRNA + 0.1)))
  
# take average of the 2 Lund & Melton log2FC's to potentialy use for scoring
receptors_df2 <- receptors_df2 %>% 
  rowwise() %>% 
  mutate(SCislet_Hislet_bulkRNA_Ave_log2FC = 
           mean(c(SCislet_Hislet_bulkRNA_Lund_log2FC, 
                  SCislet_Hislet_bulkRNA_Melton_log2FC
           ))) %>% 
  ungroup()

# view data
receptors_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  dplyr::select(hgnc_symbol,
                Lund_h_islet_bulkRNA,
                sc_islet_bulk_rna,
                Melton_h_islet_bulkRNA,
                SCislet_Hislet_bulkRNA_Lund_log2FC,
                SCislet_Hislet_bulkRNA_Melton_log2FC,
                SCislet_Hislet_bulkRNA_Ave_log2FC) %>% 
  View()



## repeat for ligands list

# lowest value in data is 1, so using plus 0.1 as reasonable value to prevent -inf
ligands_df2 <- ligands_df2 %>% 
  mutate(SCislet_Hislet_bulkRNA_Lund_log2FC = log2((sc_islet_bulk_rna + 0.1) / 
                                                     (Lund_h_islet_bulkRNA + 0.1)),
         SCislet_Hislet_bulkRNA_Melton_log2FC = log2((sc_islet_bulk_rna + 0.1) / 
                                                       (Melton_h_islet_bulkRNA + 0.1)))

# take average of the 2 Lund & Melton log2FC's to potentialy use for scoring
ligands_df2 <- ligands_df2 %>%
  rowwise() %>% 
  mutate(SCislet_Hislet_bulkRNA_Ave_log2FC = 
           mean(c(SCislet_Hislet_bulkRNA_Lund_log2FC, 
                  SCislet_Hislet_bulkRNA_Melton_log2FC
           ))) %>% 
  ungroup()

# view data
ligands_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  dplyr::select(hgnc_symbol,
                Lund_h_islet_bulkRNA,
                sc_islet_bulk_rna,
                Melton_h_islet_bulkRNA,
                SCislet_Hislet_bulkRNA_Lund_log2FC,
                SCislet_Hislet_bulkRNA_Melton_log2FC,
                SCislet_Hislet_bulkRNA_Ave_log2FC) %>% 
  View()



############## plot log2fc lund vs melton H-islet vs SC-islet ########################

# receptors

# log2FC comparision to Lund & Melton data are very discrepant
# try plotting association
receptors_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = SCislet_Hislet_bulkRNA_Lund_log2FC,
             y = SCislet_Hislet_bulkRNA_Melton_log2FC)) +
  geom_point(colour = "#1F9274",
             alpha = 0.8) +
  ggtitle("Receptors") +
  labs(x = bquote(Log["2"]*FC("SC-islet/Lund H-islet")),
       y = bquote(Log["2"]*FC("SC-islet/Melton H-islet"))) +
  theme_bw()
    

# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_bulkRNA_SCislet_melton_lund_Hislet_log2fc.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)




# ligands

# log2FC comparision to Lund & Melton data are very discrepant
# try plotting association
ligands_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = SCislet_Hislet_bulkRNA_Lund_log2FC,
             y = SCislet_Hislet_bulkRNA_Melton_log2FC)) +
  geom_point(colour = "#36226B",
             alpha = 0.8) +
  ggtitle("Ligands") +
  labs(x = bquote(Log["2"]*FC("SC-islet/Lund H-islet")),
       y = bquote(Log["2"]*FC("SC-islet/Melton H-islet"))) +
  theme_bw()


# save the plot
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_bulkRNA_SCislet_melton_lund_Hislet_log2fc.JPG",
       device = "jpg",
       width = 1500,
       height = 1500,
       units = "px",
       scale = 0.8)





# save new datatables
write_tsv(receptors_df2,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_bulkRNA.tsv")

write_tsv(ligands_df2,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_bulkRNA.tsv")


