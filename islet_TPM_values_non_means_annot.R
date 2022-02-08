# using individual n (non-mean) tpm values from Lund 2014 islet rna-seq data
# from file tpm_values.tsv, will annotate gene names based on ensemble
# transcript IDS and convert from transcript level tpm to gene level tpm


# load libraries
library(tidyverse)
library(biomaRt)

# check working directory
getwd()


# open raw islet tpm values (created in ave_TPM_calculation.R)
# this file was created by looping over all n=89 folders from SRP029262
# for each n from the Lund 2014 project
transcripts_df <- (read.csv("tpm_values.tsv", sep = "\t")) %>%
  rename(ensembl_transcript_id_version = target_id)


# get gene annotations from biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes_df <- getBM(attributes = 
                    c("ensembl_transcript_id_version", 
                      "ensembl_gene_id",
                      "hgnc_symbol",
                      "gene_biotype", 
                      "description"),
                  filters=c('ensembl_transcript_id_version'),
                  values =  transcripts_df[1],
                  mart = mart)


# add raw tpm values back to our annotated genes info from biomart
df_with_tpm <- left_join(genes_df, 
                         transcripts_df, 
                         by = "ensembl_transcript_id_version")

# clean description column
df_with_tpm$description <- sub("\\[.*","", df_with_tpm$description)



# save datatable of transcript level tpm's with gene annotations
write_tsv(df_with_tpm, "data/islet_transcript_tpm_annot.tsv")


### sum transcripts from same gene to get gene level TPM counts

# check how many trancripts we have first
nrow(df_with_tpm)
# 181041 transcripts

# start with just the transcript
gene_tpm <- df_with_tpm %>%
  # remove all annotation columns except ensemble gene ID to calculate sums
  dplyr::select(-c(1,3,4,5)) %>%  
  group_by(ensembl_gene_id) %>%
  summarise_all(sum)

# now check number of rows to see how many reduced by going 
# from transcript to genes
nrow(gene_tpm)
# 39993

# so went from 181,041 transcripts to 39,993 genes

## add back annotations to gene level tpm list from biomart
gene_tpm_annot <- left_join(gene_tpm, 
                         df_with_tpm[2:5], 
                         by = "ensembl_gene_id") %>% 
                  # the join created duplicate rows, so lets remove them
                  distinct() %>% 
                  dplyr::select(1,91,93,92,2:90)


# save datatable of gene level tpm's with gene annotations
write_tsv(gene_tpm_annot, "data/islet_gene_tpm_annot.tsv")



