library(biomaRt)
library(UniprotR)
library(tidyverse)
library(rentrez)


getwd()
#open df with transcript IDs
transcripts_df <- (read.csv("mean_tpm_values.tsv", sep = "\t")) %>%
  rename(ensembl_transcript_id_version = target_id)
  

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

df_with_tpm <- left_join(genes_df[1:2], 
                         transcripts_df, 
                         by = "ensembl_transcript_id_version")

# sum transcripts from same gene to get gene level TPM counts
gene_level_tpm <- df_with_tpm[2:3] %>%
  group_by(ensembl_gene_id) %>%
  summarise_all(sum)
  
# open secreted factors list from sci signalling 2019 paper
secreted_factors <- read.csv("sci_signalling_secreted_factors.csv") %>%
  rename(ensembl_gene_id = Ensembl.gene.id)


# append islet rna-seq gene level tpm data to list of secreted factors
# 
secreted_factors_tpm <- left_join(secreted_factors, 
                                  gene_level_tpm, 
                                  by = "ensembl_gene_id")



# run biomart again to add back in annotations from biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes_df2 <- getBM(attributes = 
                    c("ensembl_gene_id",
                      "hgnc_symbol",
                      "gene_biotype", 
                      "description"),
                  filters=c('ensembl_gene_id'),
                  values =  secreted_factors[1],
                  mart = mart)



sec_fact_tpm_annotat <- left_join(secreted_factors_tpm,
                                  genes_df2,
                                  by = "ensembl_gene_id")
                       

colnames(secreted_factors_tpm)

#save file
write_tsv(sec_fact_tpm_annotat, "sci_sig_secreted_factors_tpm.tsv")

# double bracket converts the column from dataframe to vector
sec_fact_tpm_annotat[[3]]


# use uniprotR package to get subcellular location data
subcell_df <- GetSubcellular_location(sec_fact_tpm_annotat[[3]], 
                                           directorypath = NULL) %>% 
              rownames_to_column(var = "UniProt.accession")
#save file
write_tsv(subcell_df, "uniprot_subcell_loc_df.tsv")

             
# use uniprotR package to get protein family data
family_df <- GetFamily_Domains(sec_fact_tpm_annotat[[3]], 
                                      directorypath = NULL) %>% 
              rownames_to_column(var = "UniProt.accession")
#save file
write_tsv(family_df, "uniprot_family_df.tsv")
              

# use uniprotR package to get Gene Ontology (GO) data
geno_ontology_df <- GetProteinGOInfo(sec_fact_tpm_annotat[[3]], 
                                     directorypath = NULL) %>% 
                    rownames_to_column(var = "UniProt.accession") 
#save file
write_tsv(geno_ontology_df, "uniprot_GO_ontology_df.tsv")
                   

# clean description column
sec_fact_tpm_annotat$description <- sub("\\[.*","", sec_fact_tpm_annotat$description)



  
# add all 3 uniprot data frames back to original df with rest of data
sec_factors_final_df <- sec_fact_tpm_annotat %>% 
  left_join(subcell_df,
       by = "UniProt.accession") %>% 
  left_join(family_df,
       by = "UniProt.accession") %>% 
  left_join(geno_ontology_df,
       by = "UniProt.accession")


# reorder the columns and remove unnecessary ones
# there are conflicting select functions from different libraries so specifc dplyr
sec_factors_final_df_select <-  dplyr::select(sec_factors_final_df, 
              hgnc_symbol,
              description,
              mean_islet_tpm = mean_tpm,
              ensembl_gene_id,
              UniProt.accession,
              Protein.families,
              Annotated.category,
              Subcellular_location = Subcellular.location..CC.,
              GO_cellular_component = Gene.ontology..cellular.component.,
              GO_molecular_function = Gene.ontology..molecular.function.,
              GO_biological_process = Gene.ontology..biological.process.
              ) %>% 
  arrange(desc(mean_islet_tpm)) #sort descending based on islet expression

# clean subcellular loc column
sec_factors_final_df_select$Subcellular_location <- sub("SUBCELLULAR LOCATION: ","", sec_factors_final_df_select$Subcellular_location)


#save file
write_tsv(sec_factors_final_df_select, "secreted_factors_annotated.tsv")


### try to get Entrez Gene Summary data using rentrez


# first get entrez ids using biomart
entrez_ids <- getBM(attributes = 
                    c("ensembl_gene_id",
                      'entrezgene_id'),
                  filters=c('ensembl_gene_id'),
                  values =  sec_factors_final_df_select$ensembl_gene_id,
                  mart = mart)

# use rentrez to get summaries using entrez ids
entrez_summary_df <- entrez_ids

entrez_summary_df$entrezgene_id <- as.character(entrez_summary_df$entrezgene_id)

entrez_summary_df <- add_column(entrez_summary_df, gene_summary_full = "")



for (row in 1:nrow(entrez_summary_df)){
  entrez_summary_df$gene_summary_full[row] <-  entrez_summary(id = entrez_summary_df$entrezgene_id[row], 
                                                            db = "gene")[[17]]
}



