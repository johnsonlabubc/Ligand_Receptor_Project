### add all bulk RNA-seq, scRNA-seq, and proteomics data to our annotated
### lists of receptors
### for each measurement, also convert the expression values to ranks



# load libraries
library(tidyverse)
library(Seurat)

# load data that includes most recent manual annotations on the google sheet
# this include manual correction of 46 NA's for genes with consensus score above 4
# did not manually correct NA's with consensus score below 5 cause we can
# toss them out 
# plus some additional comments
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_annot_manual.tsv", 
                      sep = "\t") %>% 
  rename(genecards_url = 17)

# check how many NA's are starting with
sum(is.na(gene_list$hgnc_symbol))
# 338 NA's left


# remove all genes with consensus score below 5

gene_list <- gene_list %>% 
  filter(consensus_score > 4)

# confirm that this removed all the NA's
sum(is.na(gene_list$hgnc_symbol))
# 0 NA's left


########## Collapse genes down to the level of gene symbol ####################
# merge all rows that contain ensembl gene ID's that map to the same HGNC symbol
# when merging expression values, take the highest value

# check how many rows are starting with
gene_list$hgnc_symbol %>% 
  length()
# 6108 rows down from 10,000 before filtering by consensus score

# check how many unique hgnc symbols
gene_list$hgnc_symbol %>% 
  unique() %>% 
  length()
# 950 unique hgnc symbols

# check how many unique uniprot symbols
gene_list$uniprot_gn_id %>% 
  unique() %>% 
  length()
# 948 unique uniprot symbols



# create table with just hgnc symbols and ensembl ids
# will append expression data to this table
gene_ensembl_list <- gene_list %>% 
  dplyr::select(hgnc_symbol,
         ensembl_gene_id)

# create table with just hgnc symbols and uniprot ids
# will append expression data to this table
gene_uniprot_list <- gene_list %>% 
  dplyr::select(hgnc_symbol,
         uniprot_gn_id)


# some duplicates likely due to NA's
# also some weird duplicates where everything is the same except 1 column
# regardless, grouping by HGNC symbol will solve these issues


# create reduced table of annotations data
# this will be our main datatable
# for each column, paste together all unique values corresponding to a hgnc symbol
# all columns except 1 are strings so each to paste the strings together
gene_list_reduced <- gene_list %>% 
  group_by(hgnc_symbol) %>% 
  summarise(description = paste0(unique(description), collapse=", "),
            keep_in_list = paste0(unique(keep_in_list), collapse=", "),
            comments = paste0(unique(comments), collapse=", "),
            ensembl_gene_ids = paste0(unique(ensembl_gene_id), collapse=", "),
            uniprot_gn_ids = paste0(unique(uniprot_gn_id), collapse=", "),
            Protein.families = paste0(unique(Protein.families), collapse=", "),
            # consensus score is only numeric column here
            # take max value instead of pasting
            # only effects 1 gene, CALCA, which has 2 diff consensus scores
            consensus_score = max(unique(consensus_score)),
            generic_categories = paste0(unique(generic_categories), collapse=", "),
            all_categories = paste0(unique(all_categories), collapse=", "),
            ligand_uniprot = paste0(unique(ligand_uniprot), collapse=", "),
            ligand_symbol = paste0(unique(ligand_symbol), collapse=", "),
            ligand_references = paste0(unique(ligand_references), collapse=", "),            
            GO_cellular_component = paste0(unique(GO_cellular_component), collapse=", "),
            GO_molecular_function = paste0(unique(GO_molecular_function), collapse=", "),
            GO_biological_process = paste0(unique(GO_biological_process), collapse=", "),
            genecards_url = paste0(unique(genecards_url), collapse=", ")
  )




########## Lund 2014 bulk human islet rna-seq ##########################

# load mean per gene TPM values data (already converted from per transcript)
lund_tpm <- read.csv("ligand_receptor_lists/initial_exploration_2021/data/islet_gene_tpm_annot.tsv", 
                     sep = "\t") %>% 
  dplyr::select(1,5:93) %>% 
  # calculate average expression per gene accross all n=89 donors
  transmute(ensembl_gene_id, #transmute drops all other columns
            islet_tpm = rowMeans(dplyr::select(., -ensembl_gene_id))) 


# get max tpm value for each hgnc gene
lund_tpm_hgnc <- left_join(gene_ensembl_list,
                           lund_tpm,
                           by = "ensembl_gene_id") %>% 
  # drop the ensembl gene id column
  dplyr::select(hgnc_symbol,
         islet_tpm) %>% 
  # drop NA's in tpm values
  drop_na(islet_tpm) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(islet_tpm = max(islet_tpm, na.rm = TRUE)) %>% 
  # calculate the ranks
  mutate(islet_tpm_rank = dense_rank(desc(islet_tpm)))

# append max tpm values to our genes list
genes_expression <- left_join(gene_list_reduced,
                              lund_tpm_hgnc,
                              by = "hgnc_symbol")


########## human islet proteomics data from jelena #########################

# These are all the proteins (and quantification) found in our human islet 
# pooled sample (pooled protein from 10 donors, 5 males, 5 females, wide range 
# of  BMI and age)

# load data
islet_proteomics_df <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/human_islet_proteomics_jelena.csv") %>% 
  # select only the columns we need and rename the columns
  dplyr::select(uniprot_gn_id = Protein.Ids,
                islet_proteomics = 6) 

# get max value for each hgnc gene
islet_proteomics_uniprot <- left_join(gene_uniprot_list,
                                      islet_proteomics_df,
                                      by = "uniprot_gn_id") %>% 
  # drop the uniprot gene id column
  dplyr::select(hgnc_symbol,
         islet_proteomics) %>% 
  # drop NA's in tpm values
  drop_na(islet_proteomics) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(islet_proteomics = max(islet_proteomics, na.rm = TRUE)) %>% 
  # calculate the ranks
  mutate(islet_proteomics_rank = dense_rank(desc(islet_proteomics)))

# append proteomics data to our genes list using hgnc
genes_expression <- left_join(genes_expression,
                              islet_proteomics_uniprot,
                              by = "hgnc_symbol")






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


#### old way i calulated before grouping by hgnc symbols
# append GTEx mean expression to genes list
#genes_expression <- left_join(genes_expression,
#                              gtex_data,
#                              by = "ensembl_gene_id")
#
#genes_expression <- genes_expression %>% 
# calculate the islet specificity                         
#  mutate(islet_specificity = (islet_tpm + 0.01)/(gtex_mean + 0.01)) %>% 
# calculate the ranks
# mutate(islet_specificity_rank = dense_rank(desc(islet_specificity))) %>% 
# drop GTEx median column
#  dplyr::select(-gtex_mean)



# get max value for each hgnc gene
gtex_data_ensembl <- left_join(gene_ensembl_list,
                               gtex_data,
                               by = "ensembl_gene_id") %>% 
  # drop the uniprot gene id column
  dplyr::select(hgnc_symbol,
         gtex_mean) %>% 
  # drop NA's in tpm values
  drop_na(gtex_mean) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(hgnc_symbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(gtex_mean = max(gtex_mean, na.rm = TRUE))

# append gtex data to our genes list using hgnc
genes_expression <- left_join(genes_expression,
                              gtex_data_ensembl,
                              by = "hgnc_symbol") %>% 
  # calculate the islet specificity                         
  mutate(islet_specificity = (islet_tpm + 0.01)/(gtex_mean + 0.01)) %>% 
  # calculate the ranks
  mutate(islet_specificity_rank = dense_rank(desc(islet_specificity))) %>% 
  # drop GTEx median column
  dplyr::select(-gtex_mean)




########## stem cell derived cells bulk RNA-seq from Doug Melton lab ###########







########## human islet scRNA-seq from lynn lab ##############################
# this data is provided at the level of HGNC symbols
# this is common in scRNA-seq b/c abundances are so low

# open the 4.6 GB RDS file containing the already processed data
filename <- file.choose()
seurat_obj <- readRDS(filename)

# view different cell types
head(Idents(seurat_obj))

rm(gene_list)
rm(gtex_data)
rm(lund_tpm)
rm(islet_proteomics_df)


gc()
memory.limit(size=1000000)

# compute mean expression by cell type
#If slot is set to 'data', this function assumes that the data has been 
#log normalized and therefore feature values are exponentiated prior to averaging
#so that averaging is done in non-log space

#seurat_ave_expr <- AverageExpression(seurat_obj,
#                    group.by = "cell_type")

# output includes 3 separate assays, RNA, integrated, & SCT
# got Warning message:
# In PseudobulkExpression(object = object, pb.method = "average",  :
# Exponentiation yielded infinite values. `data` may not be log-normed.

View(seurat_ave_expr$RNA)
# we see that the RNA averages include infinte values
# got warning "`data` may not be log-normed."

# view the details of the 3 assays
seurat_obj@assays
# they are identical in the number of cells and features (genes)

#within the RNA assay, we see there is counts, data, scale.data, etc)
View(seurat_obj@assays$RNA@data@Dimnames[[1]])
View(seurat_obj@assays$RNA@counts)




seurat_ave_expr <- AverageExpression(seurat_obj,
                                     # group by the predetermined cell types
                                     group.by = "cell_type",
                                     # only do RNA, not SCT & integrated
                                     assays = "RNA", 
                                     # data is not log normalized
                                     # so don't do exponentiation
                                     slot = "counts")


View(seurat_ave_expr$RNA)
# data looks good now & does not appear to be log normalized
# does not have infinite values anymore & did not give a warning


# convert to dataframe 
seurat_ave_expr_df <- as.data.frame(seurat_ave_expr$RNA) %>% 
  # make rownames a column
  rownames_to_column(var = "hgnc_symbol")


# save average counts by cell type in scRNA-seq data
write_tsv(seurat_ave_expr_df, 
          "single_cell_analysis/Islet_sc_rnaseq_francis/ave_counts_cell_types.tsv")


# append human islet scRNA-seq ave counts to genes list
genes_expression <- left_join(genes_expression,
                              seurat_ave_expr_df,
                              by = "hgnc_symbol")


rm(seurat_obj)


# calculate the ranks
genes_expression <- genes_expression %>% 
  mutate(delta_rank = dense_rank(desc(Delta)),
         beta_rank = dense_rank(desc(Beta)),
         alpha_rank = dense_rank(desc(Alpha)))

# calculate mean across all islet cell types in scRNA-seq
genes_expression$mean_scRNAseq <- rowMeans(genes_expression[24:34])

# calculate cell type specifities
genes_expression <- genes_expression %>%   
  mutate(delta_specificity = (Delta + 0.01)/(mean_scRNAseq + 0.01),
         beta_specificity = (Beta + 0.01)/(mean_scRNAseq + 0.01),
         alpha_specificity = (Alpha + 0.01)/(mean_scRNAseq + 0.01)) %>% 
  # calculate ranks of the specificities
  mutate(delta_specificity_rank = dense_rank(desc(delta_specificity)),
         beta_specificity_rank = dense_rank(desc(beta_specificity)),
         alpha_specificity_rank = dense_rank(desc(alpha_specificity)))

########## stem cell derived cells scRNA-seq from lynn lab ####################





# save final file
write_tsv(genes_expression, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_annot_expr.tsv")




