# revist the bulk RNA-seq & human islet scRNA-seq to correc the data
# data in spreadsheet currently has 2 flaws:
# 1) NA's need to be replaced with 0's
# 2) when calculating ratios, added 0.01 to prevent infinite values.
# 0.01 is too large a value and instead should add 0.0000001 b/c some
# of the expression values are very small

# load libraries
library(tidyverse)


##### prepare starting dataframe ###########

# load genes list with all prior expression data and annotations
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_stemcell.tsv", 
                      sep = "\t") %>% 
  # drop columns with outdated info that we are replacing
  select(-islet_tpm_rank,
         -islet_proteomics,
         -islet_proteomics_rank,
         -islet_specificity,
         -islet_specificity_rank,
         -islet_proteomics,
         -delta_rank,
         -beta_rank,
         -alpha_rank,
         -mean_scRNAseq,
         -delta_specificity,
         -beta_specificity,
         -alpha_specificity,
         -delta_specificity_rank,
         -beta_specificity_rank,
         -delta_specificity,
         -alpha_specificity_rank)


# replace all NA's in human islet bulk & scRNA-seq data with zeros 
gene_list_no_na <- mutate_at(gene_list, 
                                    c("islet_tpm", 
                                      "Immune",
                                      "Duct",
                                      "Endothelial",
                                      "Pericytes",
                                      "Acinar",
                                      "INS.SST",
                                      "INS.GCG",
                                      "PPY",
                                      "Delta",
                                      "Beta",
                                      "Alpha"), 
                                    ~replace(., is.na(.), 0.000))

######### regenerate the ranks & specificities columns ########

gene_list_ranks <- gene_list_no_na %>% 
  mutate(islet_tpm_rank = dense_rank(desc(islet_tpm)))


########## Regenerate GTEx whole pancreas and all body tissues #########################

# load gtex median gene tpm data
gtex_data <- (read.csv("ligand_receptor_lists/initial_exploration_2021/data/GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv", 
                       sep = "\t")) 

#get mean gene expression accross all tissues
gtex_data$gtex_mean <- rowMeans(gtex_data[3:56])

# convert gencode IDs to ensebl gene ids by trimming off the version number
## AM NOT sure this is the correct way to do it, but gonna try it out at least
gtex_data$Name <- gsub('\\..+$', '', gtex_data$Name)

# change column names
gtex_data <- gtex_data %>% 
  dplyr::select(ensembl_gene_id = Name,
                gtex_mean)
  
# convert from ensembl ids to hgnc symbol
gene_ensembl_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_annot_manual.tsv", 
                      sep = "\t") %>% 
  rename(genecards_url = 17)

# remove all genes with consensus score below 5
gene_ensembl_list <- gene_ensembl_list %>% 
  filter(consensus_score > 4)


# create table with just hgnc symbols and ensembl ids
# will append expression data to this table
gene_ensembl_list <- gene_ensembl_list %>% 
  dplyr::select(hgnc_symbol,
                ensembl_gene_id)


# get max value for each hgnc gene
gtex_data_ensembl <- left_join(gene_ensembl_list,
                               gtex_data,
                               by = "ensembl_gene_id") %>% 
  # drop the ensembl gene id column
  dplyr::select(hgnc_symbol,
                gtex_mean) %>% 
  # drop NA's in tpm values
  drop_na(gtex_mean) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(gtex_mean = max(gtex_mean, na.rm = TRUE))


# append gtex data to our genes list using hgnc
gene_list_ranks <- left_join(gene_list_ranks,
                              gtex_data_ensembl,
                              by = "hgnc_symbol") %>% 
  # calculate the islet specificity                         
  mutate(islet_specificity = (islet_tpm + 0.0000001)/(gtex_mean + 0.0000001)) %>% 
  # calculate the ranks
  mutate(islet_specificity_rank = dense_rank(desc(islet_specificity)))


##### Regenerate human scRNA-seq ranks & specificities ####################


# calculate the ranks
gene_list_ranks <- gene_list_ranks %>% 
  mutate(delta_rank = dense_rank(desc(Delta)),
         beta_rank = dense_rank(desc(Beta)),
         alpha_rank = dense_rank(desc(Alpha)))


# calculate mean across all islet cell types in scRNA-seq
gene_list_ranks$mean_islet_scRNAseq <- rowMeans(gene_list_ranks[19:29])


# calculate cell type specifities
gene_list_ranks <- gene_list_ranks %>%   
  mutate(delta_specificity = (Delta + 0.0000001)/(mean_islet_scRNAseq + 0.0000001),
         beta_specificity = (Beta + 0.0000001)/(mean_islet_scRNAseq + 0.0000001),
         alpha_specificity = (Alpha + 0.0000001)/(mean_islet_scRNAseq + 0.0000001)) %>% 
  # calculate ranks of the specificities
  mutate(delta_specificity_rank = dense_rank(desc(delta_specificity)),
         beta_specificity_rank = dense_rank(desc(beta_specificity)),
         alpha_specificity_rank = dense_rank(desc(alpha_specificity)))


######## reorder columns and save ###########################################

gene_list_ranks_ordered <- gene_list_ranks %>% 
  # reorder columns to put everything back in original order
  relocate(islet_tpm_rank,
           gtex_mean,
           islet_specificity,
           islet_specificity_rank,
           .after = islet_tpm) %>% 
  relocate(delta_rank,
           beta_rank,
           alpha_rank,
           mean_islet_scRNAseq,
           delta_specificity,
           beta_specificity,
           alpha_specificity,
           delta_specificity_rank,
           beta_specificity_rank,
           alpha_specificity_rank,
           .after = Alpha)


# save full spreadsheet with all corrected data
write_tsv(gene_list_ranks_ordered, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_stemcell_updateddata.tsv")




