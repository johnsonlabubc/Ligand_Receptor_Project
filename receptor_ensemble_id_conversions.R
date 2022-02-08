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

# open receptors list from FANTOM 5
receptors_df <- read.csv("ExpressionLigRec_Fantom5.txt", sep = "\t") %>% 
  dplyr::select(ApprovedSymbol, Type) %>% 
  filter( Type == "receptor")


#get ensembl ids for receptors list
receptor_biomart <- getBM(attributes = 
                      c("ensembl_gene_id",
                      "hgnc_symbol",
                      "gene_biotype", 
                      "description"),
                  filters=c('hgnc_symbol'),
                  values =  receptors_df[1],
                  mart = mart)


#save just the list of ensembl ids
write_tsv(receptor_biomart[1], "receptors_ensebl_ids.tsv")

# open fantom 5 all genes list which has uniprot accession
fantom_5 <- (read.csv("Fantom5_all_genes.txt", sep = "\t")) %>% 
  dplyr::select(ApprovedSymbol, UniProtAC)

# append uniprot accession to fantom 5 list of receptors
receptors_uniprot <- left_join(receptors_df, 
                                   fantom_5, 
                                   by = "ApprovedSymbol") %>% 
                     rename(uniprot_gn_id = UniProtAC)


#get ensembl ids for receptors list using uniprot accession
receptor_biomart <- getBM(attributes = 
                            c("ensembl_gene_id",
                              "hgnc_symbol",
                              "gene_biotype", 
                              "description",
                              "uniprot_gn_id"),
                          filters=c('uniprot_gn_id'),
                          values =  receptors_uniprot[3],
                          mart = mart)

# append islet rna-seq gene level tpm data to list of receptors
receptors_tpm_annotat <- left_join(receptors_uniprot, 
                                   receptor_biomart, 
                                   by = "uniprot_gn_id")


# append islet rna-seq gene level tpm data to list of receptors
receptors_tpm_annotat <- left_join(receptors_tpm_annotat, 
                                   gene_level_tpm, 
                                   by = "ensembl_gene_id")


















# this list from biomart includes some duplicates
# cause some genes have more than 1 ensembl gene id


#save file
write_tsv(receptors_tpm_annotat, "receptors_islet_tpm.tsv")

# double bracket converts the column from dataframe to vector
receptors_tpm_annotat[[3]]


# use uniprotR package to get subcellular location data
subcell_df <- GetSubcellular_location(receptors_tpm_annotat[[3]], 
                                      directorypath = NULL) %>% 
  rownames_to_column(var = "uniprot_gn_id")
#save file
write_tsv(subcell_df, "uniprot_subcell_loc_receptors.tsv")


# use uniprotR package to get protein family data
family_df <- GetFamily_Domains(receptors_tpm_annotat[[3]], 
                               directorypath = NULL) %>% 
  rownames_to_column(var = "uniprot_gn_id")
#save file
write_tsv(family_df, "uniprot_family_receptors.tsv")


# use uniprotR package to get Gene Ontology (GO) data
geno_ontology_df <- GetProteinGOInfo(receptors_tpm_annotat[[3]], 
                                     directorypath = NULL) %>% 
  rownames_to_column(var = "uniprot_gn_id") 
#save file
write_tsv(geno_ontology_df, "uniprot_GO_ontology_receptors.tsv")


# clean description column
receptors_tpm_annotat$description <- sub("\\[.*","", receptors_tpm_annotat$description)




# add all 3 uniprot data frames back to original df with rest of data
receptors_final_df <- receptors_tpm_annotat %>% 
  left_join(subcell_df,
            by = "uniprot_gn_id") %>% 
  left_join(family_df,
            by = "uniprot_gn_id") %>% 
  left_join(geno_ontology_df,
            by = "uniprot_gn_id")


# reorder the columns and remove unnecessary ones
# there are conflicting select functions from different libraries so specifc dplyr
receptors_final_df <-  dplyr::select(receptors_final_df, 
                                              hgnc_symbol,
                                              description,
                                              mean_islet_tpm = mean_tpm,
                                              Protein.families,
                                              Subcellular_location = Subcellular.location..CC.,
                                              GO_cellular_component = Gene.ontology..cellular.component.,
                                              GO_molecular_function = Gene.ontology..molecular.function.,
                                              GO_biological_process = Gene.ontology..biological.process.,
                                              ensembl_gene_id,
                                              UniProt.accession = uniprot_gn_id,
                      ) %>% 
                      #sort descending based on islet expression
                      arrange(desc(mean_islet_tpm)) 

# clean subcellular loc column
receptors_final_df$Subcellular_location <- sub("SUBCELLULAR LOCATION: ","", receptors_final_df$Subcellular_location)


#save file
write_tsv(receptors_final_df, "receptors_annotated.tsv")


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



