library(tidyverse)

#set working directory
setwd("~/Johnson_Lab/Islet_Bulk_RNA_Seq_Data/Lund_2014_data_George/From_Lab_PC_Human_Islet_Bulk_RNA_Seq/SRP029262/kallisto")

#get list of all folders
folder_list <- list.files()

#create initial df with just target_id
df <- (read.csv("SRR957884/abundance.tsv", sep = "\t")) %>%
  select(target_id)

df2 <- (read.csv("SRR957967/abundance.tsv", sep = "\t")) %>%
  select(target_id, tpm) %>%       #only load the rows we are interested in
  rename(!!ex := tpm)              #rename tpm column to SRR number from folder name

colnames((df2))



ex <- folder_list[1]


df_merged <- left_join(df, df2, by = "target_id")

# nrow(df)
# colnames(df)

# check that all rows are the same in each file
# df2 <- (read.csv("SRR957967/abundance.tsv", sep = "\t"))[1]
# df == df2


#loop over all folders and add the tpm values as new columns
for (i in folder_list){
  folder_name <- i
  file_path_string <- paste(i, "/abundance.tsv", sep = "")
  SRR_df <- (read.csv(file_path_string, sep = "\t")) %>%
  select(target_id, tpm) %>%       #only load the rows we are interested in
  rename(!!folder_name := tpm)
  df <- left_join(df, SRR_df, by = "target_id")
  }

nrow(df)
ncol(df)

#save output file with all values
write_tsv(df, "../../tpm_values.tsv")


#calculate average tpm and remove other rows
df_mean <- df %>%
  transmute(target_id,
            mean_tpm = rowMeans(select(., -target_id)))


#save output file with just mean tpm value
write_tsv(df_mean, "../../mean_tpm_values.tsv")

