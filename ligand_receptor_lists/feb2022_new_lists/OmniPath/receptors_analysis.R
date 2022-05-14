# add several analysis columns & create figures to explore the combined datasets


install.packages("ggrepel")


# load libraries
library(tidyverse)
library(cowplot)
# for adding non-overlapping text labels to points on ggplot scatterplot
library(ggrepel) 

# load receptors list with all prior expression data and annotations
gene_list <- read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                      sep = "\t")
  

#### add ranks & analysis columns
gene_list <- gene_list %>% 
  mutate(sorted_eBCs_rank = dense_rank(desc(sorted_eBCs)),
         # had to add very small amounts to prevent infinity values
         # the values are so small that adding 0.01 screws up the data
         eBCs_human_beta_ratio = ((sorted_eBCs + 0.0000001) / (Beta + 0.0000001)),
         human_beta_eBCs_ratio = ((Beta + 0.0000001) / (sorted_eBCs + 0.0000001)))


## explore interesting filters and top genes

gene_list %>% 
  filter(keep_in_list == "Yes") %>% 
  select(hgnc_symbol,
         description,
         sorted_eBCs,
         Beta,
         human_beta_eBCs_ratio) %>% 
  arrange(desc(human_beta_eBCs_ratio)) %>% 
  View()

# plot human_beta_eBCs_ratio against huamn beta abundance
gene_list %>% 
  filter(keep_in_list == "Yes") %>% 
  ggplot(aes(x = log2(Beta), y = log2(human_beta_eBCs_ratio), label = hgnc_symbol)) +
  geom_point(colour = "dodgerblue2", size = 2) +
  labs(x = "Log2 beta cell RNA abundance",
       y = "Log2 beta cell / eBCs ratio") +
  # this uses the ggrepel to position
  geom_text_repel(data=subset(gene_list %>% filter(keep_in_list == "Yes"), 
                              log2(human_beta_eBCs_ratio) > 7),
                  size = 2.25,
                  max.overlaps = 15,
                  force = 0.5,
                  max.time = 2,
                  force_pull = 5) +
  theme_cowplot()
