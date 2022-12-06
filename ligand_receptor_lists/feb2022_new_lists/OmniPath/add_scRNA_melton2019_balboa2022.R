# add SCbeta scRNA-seq data from Veres/Melton 2019 & Balboa/Otonkoski 2022
# melton2019 hs both expression and differential expression in comparison to 
# adult human beta cells, while balboa2022 only has differential expression


library(tidyverse)



############ load data ###############################


receptors_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_bulkRNA.tsv", 
                         sep = "\t")

ligands_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_bulkRNA.tsv", 
                       sep = "\t")


melton2019_df <-  read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/melton2019_balboa2022_scRNA/41586_2019_1168_MOESM4_ESM_meltem2019.txt", 
                           sep = "\t") %>% 
  # inverse the log2fc b/c currently its human over stem cell
  mutate(SCbeta_Hbeta_scRNA_Melton2019_log2FC = -SCbeta_Hbeta_scRNA_Melton2019)

balboa2019_df <-  read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/melton2019_balboa2022_scRNA/41587_2022_1219_MOESM8_ESM_Balboa2022.txt", 
                           sep = "\t") %>% 
  # inverse the log2fc b/c currently its human over stem cell
  mutate(SCbeta_Hbeta_scRNA_Balboa2022_log2FC = -avg_log2FC,
         SCbeta_Hbeta_scRNA_Balboa2022_adj_p_value = p_val_adj) %>% 
  # rename gene symbols column to match
  rename(hgnc_symbol = Gene) %>% 
  select(hgnc_symbol,
         SCbeta_Hbeta_scRNA_Balboa2022_log2FC,
         SCbeta_Hbeta_scRNA_Balboa2022_adj_p_value)


#################### join to receptors & ligands lists ##############################

receptors_df2 <- receptors_df %>% 
  left_join(melton2019_df[-4],
                           by = "hgnc_symbol") %>% 
  left_join(balboa2019_df,
            by = "hgnc_symbol")


ligands_df2 <- ligands_df %>% 
  left_join(melton2019_df[-4],
            by = "hgnc_symbol") %>% 
  left_join(balboa2019_df,
            by = "hgnc_symbol")
                           
                           
                  
# view data
receptors_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  dplyr::select(hgnc_symbol,
                SCbeta_Hbeta_scRNA_Balboa2022_log2FC,
                SCbeta_Hbeta_scRNA_Balboa2022_adj_p_value,
                SCbeta_Hbeta_scRNA_Melton2019_log2FC) %>% 
  View()        


# view data
ligands_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  dplyr::select(hgnc_symbol,
                SCbeta_Hbeta_scRNA_Balboa2022_log2FC,
                SCbeta_Hbeta_scRNA_Balboa2022_adj_p_value,
                SCbeta_Hbeta_scRNA_Melton2019_log2FC) %>% 
  View()   
   

############## plot log2fc balboa vs melton  ########################

# receptors

# log2FC comparision to Lund & Melton data are very discrepant
# try plotting association
receptors_df2 %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = SCbeta_Hbeta_scRNA_Balboa2022_log2FC,
             y = SCbeta_Hbeta_scRNA_Melton2019_log2FC)) +
  geom_point(colour = "#1F9274",
             alpha = 0.8) +
  ggtitle("Receptors") +
  theme_bw()








# save new datatables
write_tsv(receptors_df2,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_meltonbalboa.tsv")

write_tsv(ligands_df2,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_meltonbalboa.tsv")


                           