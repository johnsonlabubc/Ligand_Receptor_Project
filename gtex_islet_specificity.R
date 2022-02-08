
# determine gene expression islet specificity by ratio of islet expression 
# to average expression accross all gtex tissues

library(tidyverse)


# load gtex median gene tpm data
gtex_data <- (read.csv("data/GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv", sep = "\t")) 

#get mean gene expression accross all tissues
gtex_data$gtex_mean <- rowMeans(gtex_data[3:56])

# convert gencode IDs to ensebl gene ids by trimming off the version number
## AM NOT sure this is the correct way to do it, but gonna try it out at least
gtex_data$Name <- gsub('\\..+$', '', gtex_data$Name)

# change column names
gtex_data <- gtex_data %>% 
  dplyr::select(ensembl_gene_id = Name,
                gtex_symbol = Description,
                gtex_mean)

# add log2 median TPM column
gtex_data <- gtex_data %>% 
  # add 0.1 to values so dont get -inf values
  mutate(log2_mean_gtex_tpm = log2(gtex_mean + 0.1))


#### append GTEx mean expression to manually annotated secreted factor list

# load secreted factors list
secreted_factors_df <- (read.csv("data/secreted_factors_annotated_manual.tsv", sep = "\t")) %>% 
  dplyr::filter(keep_in_list. == "Yes" ) %>% 
  dplyr::select(hgnc_symbol, description, mean_islet_tpm, ensembl_gene_id) %>% 
  mutate(log2_islet_tpm = log2(mean_islet_tpm + 0.1)) %>% 
  # append gtex data
  left_join(gtex_data, 
            by = "ensembl_gene_id")


# plot gtex mean expression against islet expression
secreted_factors_df  %>% 
  ggplot(aes(x = log2_islet_tpm, y = log2_mean_gtex_tpm)) +
  geom_point(col = "lightseagreen") + 
  ggtitle("Secreted Factors Pancreas vs Islet Expression")


# calculate islet/mean gtex ratio
secreted_factors_df <- secreted_factors_df %>% 
  mutate(islet_specificity = log2(mean_islet_tpm/gtex_mean)) %>% 
  arrange(desc(islet_specificity))


# look at list of top islet specific genes
secreted_factors_specific <- secreted_factors_df %>% 
  dplyr::select(hgnc_symbol,
                description,
                log2_islet_tpm,
                islet_specificity)

#save file
write_tsv(secreted_factors_specific, "data/secreted_factors_islet_specificity.tsv")



#### make same data for receptors list ##############################3

# append GTEx mean expression to manually annotated receptors list

# load receptors list
receptors__annot_df <- (read.csv("data/receptors_annotated_manual.tsv", sep = "\t")) %>% 
  dplyr::filter(keep_in_list. == "Yes") %>% 
  dplyr::select(hgnc_symbol, description, mean_islet_tpm, ensembl_gene_id) %>% 
  mutate(log2_islet_tpm = log2(mean_islet_tpm + 0.1)) %>% 
  # append gtex data
  left_join(gtex_data, 
            by = "ensembl_gene_id")


# plot gtex mean expression against islet expression
receptors__annot_df  %>% 
  ggplot(aes(x = log2_islet_tpm, y = log2_mean_gtex_tpm)) +
  geom_point(col = "lightseagreen") + 
  ggtitle("Secreted Factors Pancreas vs Islet Expression")


# calculate islet/mean gtex ratio
receptors__annot_df <- receptors__annot_df %>% 
  mutate(islet_specificity = log2(mean_islet_tpm/gtex_mean)) %>% 
  arrange(desc(islet_specificity))


# look at list of top islet specific genes
receptors__annot_df <- receptors__annot_df %>% 
  dplyr::select(hgnc_symbol,
                description,
                log2_islet_tpm,
                islet_specificity)

# save file
write_tsv(receptors__annot_df, "data/receptors_islet_specificity.tsv")
