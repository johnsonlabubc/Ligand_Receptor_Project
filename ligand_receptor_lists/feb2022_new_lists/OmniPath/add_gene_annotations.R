# add additional gene annotations to our ligand/receptor lists created in OmniPath
# will get annotations from gene ontology, uniprot, etc
# code below is mostly modified from receptor_ensemble_id_conversions.R


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biomaRt")

#install.packages("UniprotR")
#BiocManager::install("GenomicAlignments")

# load libraries
library(biomaRt)
library(UniprotR)
library(tidyverse)
library(rentrez)



gene_list <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands.tsv", 
                       sep = "\t"))

# get ensembl transcript IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes_df <- getBM(attributes = 
                    c("uniprot_gn_id",
                      "ensembl_gene_id",
                      "hgnc_symbol",
                      "description"),
                  filters= c('uniprot_gn_id'),
                  values =  gene_list[1],
                  mart = mart)


# clean description column
genes_df$description <- sub("\\[.*",
                            "", 
                            genes_df$description)

# use uniprotR package to get protein family data
family_df <- UniprotR::GetFamily_Domains(gene_list[[1]], 
                               directorypath = NULL) %>% 
  rownames_to_column(var = "uniprot_gn_id")

family_df <- family_df %>% 
  dplyr::select(uniprot_gn_id,
                    Protein.families)

# use uniprotR package to get Gene Ontology (GO) data
geno_ontology_df <- UniprotR::GetProteinGOInfo(gene_list[[1]], 
                                     directorypath = NULL) %>% 
  rownames_to_column(var = "uniprot_gn_id") 

geno_ontology_df <- geno_ontology_df %>% 
  dplyr::select(uniprot_gn_id,
                GO_cellular_component = Gene.ontology..cellular.component.,
                GO_molecular_function = Gene.ontology..molecular.function.,
                GO_biological_process = Gene.ontology..biological.process.)

# append manual annotations from google sheet
manual_annot <- (read.csv("ligand_receptor_lists/feb2022_new_lists/data/old_ligand_manual_annots_google_sheet.tsv", 
                       sep = "\t")) %>% 
  dplyr::select(uniprot_gn_id = UniProt.accession,
         keep_in_list = keep_in_list.,
         comments)


# add all data frames back to original df with rest of data
genes_annotated <- gene_list %>% 
  rename(uniprot_gn_id = uniprot) %>% 
  left_join(genes_df,
            by = "uniprot_gn_id") %>% 
  left_join(family_df,
            by = "uniprot_gn_id") %>% 
  left_join(geno_ontology_df,
            by = "uniprot_gn_id") %>% 
  left_join(manual_annot,
            by = "uniprot_gn_id") %>% 
  dplyr::select(hgnc_symbol,
         description,
         keep_in_list,
         comments,
         ensembl_gene_id,
         uniprot_gn_id,
         Protein.families,
         consensus_score,
         generic_categories,
         all_categories,
         receptor_uniprot,
         receptor_symbol,
         receptor_references,
         GO_cellular_component,
         GO_molecular_function,
         GO_biological_process)


genes_annotated %>% 
  group_by(keep_in_list) %>% 
  summarise(n())

#save file
write_tsv(genes_annotated, 
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_annot.tsv")


