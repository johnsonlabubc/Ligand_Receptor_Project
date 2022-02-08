# clean up Lund 2014 RNA-seq metadata

metadata_df <- (read.csv("data/SRP029262__Lund_metadata.txt", sep = ",")) %>% 
  dplyr::select(Run, 
                Age, 
                BMI, 
                Sex = gender,
                hba1c)

# save simplified metadata datatable
write_tsv(metadata_df, "data/SRP029262_Lund_metadata_donors.tsv")
