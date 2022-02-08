library(tidyverse)


# load gtex median gene tpm data
gtex_data <- (read.csv("../../../RNA_Seq_Data/GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv", sep = "\t")) %>% 
   # keep only pancreas expression data 
   dplyr::select(Name,
                 Description,
                 Pancreas)  

# convert gencode IDs to ensebl gene ids by trimming off the version number
## AM NOT sure this is the correct way to do it, but gonna try it out at least
gtex_data$Name <- gsub('\\..+$', '', gtex_data$Name)

# change column names
gtex_data <- gtex_data %>% 
  dplyr::select(ensembl_gene_id = Name,
                gtex_symbol = Description,
                Pancreas_Median_TPM = Pancreas)

# add log2 median TPM column
gtex_data <- gtex_data %>% 
  mutate(log2_median_pancreas_tpm = log2(Pancreas_Median_TPM))

#### append GTEx pancreas expression to manually annotated secreted factor list

# load secreted factors list
secreted_factors_df <- (read.csv("secreted_factors_annotated_manual.tsv", sep = "\t")) %>% 
  dplyr::filter(keep_in_list. == "Yes" ) %>% 
  dplyr::select(hgnc_symbol, description, mean_islet_tpm, ensembl_gene_id) %>% 
  mutate(log2_islet_tpm = log2(mean_islet_tpm)) %>% 
  # append gtex data
  left_join(gtex_data, 
            by = "ensembl_gene_id")


# check number of NA's in log tpm values
sum(is.na(secreted_factors_df$log2_median_pancreas_tpm)) 
# 2 NA's in median pancreas tpm

sum(is.na(secreted_factors_df$log2_islet_tpm)) 
# 9 NA's in islet tpm

#so total has 11 NA's in secreted factors
# which is why plot is missing 11 geom_points

  
# plot pancreas expression against islet expression
secreted_factors_df  %>% 
  ggplot(aes(x = log2_islet_tpm, y = log2_median_pancreas_tpm)) +
  geom_point(col = "lightseagreen") + 
  ggtitle("Secreted Factors Pancreas vs Islet Expression")

ggsave("secreted_factors_pancreas_vs_islet_expression.png")





#### make same plot for receptors list ##############################3

# append GTEx pancreas expression to manually annotated receptors list

# load receptors list
receptors__annot_df <- (read.csv("receptors_annotated_manual.tsv", sep = "\t")) %>% 
  dplyr::filter(keep_in_list. == "Yes") %>% 
  dplyr::select(hgnc_symbol, description, mean_islet_tpm, ensembl_gene_id) %>% 
  mutate(log2_islet_tpm = log2(mean_islet_tpm)) %>% 
  # append gtex data
  left_join(gtex_data, 
            by = "ensembl_gene_id")

# check number of NA's in log tpm values
sum(is.na(receptors__annot_df$log2_median_pancreas_tpm)) 
# 99 NA's in median pancreas tpm

sum(is.na(receptors__annot_df$log2_islet_tpm)) 
# 9 NA's in islet tpm

#so total has 109 receptors with NAs
# which is why plot is missing 108 geom_points



# plot pancreas expression against islet expression
receptors__annot_df  %>% 
  ggplot(aes(x = log2_islet_tpm, y = log2_median_pancreas_tpm)) +
  geom_point(col = "lightseagreen") + 
  ggtitle("Receptors Pancreas vs Islet Expression")

ggsave("secreted_factors_pancreas_vs_islet_expression.png")

