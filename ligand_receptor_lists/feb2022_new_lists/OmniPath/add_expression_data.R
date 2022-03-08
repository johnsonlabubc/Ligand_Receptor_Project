### add all bulk RNA-seq, scRNA-seq, and proteomics data to our annotated
### lists of ligands & receptors
### for each measurement, also convert the expression values to ranks



# load libraries
library(tidyverse)
library(biomaRt)
library(Seurat)

# load data
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_annot.tsv", 
                       sep = "\t")
                       



                       
########## Lund 2014 bulk human islet rna-seq ##########################

# load mean per gene TPM values data (already converted from per transcript)
lund_tpm <- read.csv("ligand_receptor_lists/initial_exploration_2021/data/islet_gene_tpm_annot.tsv", 
                     sep = "\t") %>% 
  dplyr::select(1,5:93) %>% 
  # calculate average expression per gene accross all n=89 donors
  transmute(ensembl_gene_id, #transmute drops all other columns
            islet_tpm = rowMeans(dplyr::select(., -ensembl_gene_id))) 


# append mean tpm values to our genes list
genes_expression <- left_join(gene_list,
                              lund_tpm,
                              by = "ensembl_gene_id") %>% 
  # calculate the ranks
  mutate(islet_tpm_rank = dense_rank(desc(islet_tpm)))





########## human islet proteomics data from jelena #########################

# These are all the proteins (and quantification) found in our human islet 
# pooled sample (pooled protein from 10 donors, 5 males, 5 females, wide range 
# of  BMI and age)

# load data
islet_proteomics_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/human_islet_proteomics_jelena.csv") %>% 
  # select only the columns we need and rename the columns
  dplyr::select(uniprot_gn_id = Protein.Ids,
         islet_proteomics = 6) 


# append proteomics data to our genes list using uniprot accessions
genes_expression <- left_join(genes_expression,
                              islet_proteomics_df,
                              by = "uniprot_gn_id") %>% 
  # calculate the ranks
  mutate(islet_proteomics_rank = dense_rank(desc(islet_proteomics)))




########## GTEx whole pancreas and all body tissues #########################

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


# append GTEx mean expression to genes list
genes_expression <- left_join(genes_expression,
                              gtex_data,
                              by = "ensembl_gene_id")

genes_expression <- View(genes_expression %>% 
  # calculate the islet specificity                         
  mutate(islet_specificity = (islet_tpm + 0.01)/(gtex_mean + 0.01)) %>% 
  # calculate the ranks
  mutate(islet_specificity_rank = dense_rank(desc(islet_specificity))) %>% 
  # drop GTEx median column
  dplyr::select(-gtex_mean))


########## stem cell derived cells bulk RNA-seq from Doug Melton lab ###########






########## human islet scRNA-seq from lynn lab ##############################

# open the 4.6 GB RDS file containing the already processed data
filename <- file.choose()
seurat_obj <- readRDS(filename)

AverageExpression(seurat_obj,
                  group.by = "ident")
  







########## stem cell derived cells scRNA-seq from lynn lab ####################w





