library(tidyverse)
library(cowplot)

cowplot::theme_cowplot()

# open mean islet secreted factors tpm data

setwd("~/Johnson_Lab/Islet_Bulk_RNA_Seq_Data/Lund_2014_data_George/From_Lab_PC_Human_Islet_Bulk_RNA_Seq")
getwd()

data <- (read.csv("secreted_factors_annotated_manual.tsv", sep = "\t")) %>% 
  dplyr::filter(keep_in_list. == "Yes" ) %>% 
  select(hgnc_symbol, description, mean_islet_tpm, ensembl_gene_id) %>% 
  mutate(log2_islet_tpm = log2(mean_islet_tpm))

# plot histogram of log2 mean tpm values
data  %>% 
  ggplot(aes(x = log2_islet_tpm)) +
  geom_histogram(fill = "lightseagreen") +
  ggtitle("526 secreted factors") +
  theme_cowplot()

ggsave("secreted_factors__filtered_islet_tpm_histogram_526.jpg")


#### Make the same plot for the receptors

receptors_data <- (read.csv("receptors_annotated_manual.tsv", sep = "\t")) %>% 
  dplyr::filter(keep_in_list. == "Yes" ) %>% 
  dplyr::select(hgnc_symbol, description, mean_islet_tpm, ensembl_gene_id) %>% 
  mutate(log2_islet_tpm = log2(mean_islet_tpm))


# plot histogram of log2 mean tpm values
receptors_data  %>% 
  ggplot(aes(x = log2_islet_tpm)) +
  geom_histogram(fill = "lightseagreen") +
  ggtitle("321 Receptors") +
  theme_cowplot()


ggsave("receptors__filtered_islet_tpm_histogram_321.jpg")







### look at all 89 individual n's instead of just TPM mean

# open file
raw_tpm_df <- (read.csv("tpm_values.tsv", sep = "\t"))
