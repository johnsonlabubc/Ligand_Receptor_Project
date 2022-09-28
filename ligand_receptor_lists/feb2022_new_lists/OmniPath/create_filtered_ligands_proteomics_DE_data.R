# create filtered tsv file of ligands with just the proteomics DE data for francis


# load ligands data
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics_DE.tsv", 
                        sep = "\t")) %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         description,
         Protein.families,
         all_sc_proteomics,
         ND_islet_proteomics,
         stemcell_NDislet_proteomics_FC,
         stemcell_NDislet_proteomics_adj_p_value)


write_tsv(ligands_df,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_proteomics_DE_forFrancis.tsv")
