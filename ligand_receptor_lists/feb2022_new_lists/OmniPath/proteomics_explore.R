# explore the new proteomics data from Foster lab now that have all the 
# n's. Compare it to the bulk RNA-seq data



library(tidyverse)
library(cowplot)


# open receptors list
# this list has been filtered to exclude all receptors with a consensus score below 6
# and to remove all that were manually annoated as "No" to keep in list
# kept the "TBDs"
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_filtered_for_amino_acid_length.txt", 
                          sep = "\t"))
