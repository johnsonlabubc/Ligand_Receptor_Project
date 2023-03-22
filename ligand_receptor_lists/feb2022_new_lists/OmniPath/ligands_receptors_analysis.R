# combining the receptors & ligands dataframes for convienence,
# starting with the dataframes with n=3 proteomics added


# install.packages("viridis")


library(tidyverse)
library(cowplot)
library(ggrepel) 
library(viridis)
library(ggpubr)





######### combine ligand and receptor dataframes ###############################

# open receptors list with proteomics
receptors_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_with_proteomics.tsv", 
                          sep = "\t"))
# open ligands list with proteomics
ligands_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_with_proteomics.tsv", 
                        sep = "\t")) 

# join the dataframes, keeping all rows from both dataframes
genes_df <- receptors_df %>% 
  full_join(ligands_df,
            by = "hgnc_symbol")




# check how many receptors and ligands are in our filtered lists

# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  nrow()
# 422 ligands


# receptors
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  nrow()
# 349 receptors





# double check filtered lists for genes that dont belong
receptors_df %>% 
 # filter(keep_in_list %in% c("Yes", "TBD")) %>% 
 # filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         description,
         keep_in_list,
         consensus_score) %>% 
  View()



########### (OLD) plot human islet proteomics vs bulk RNAseq #############################


# ligands

# plot human islet bulk RNA-seq against proteomics
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(islet_tpm), 
             y = log2(all_islet_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 Human Islet RNA-seq",
       y = "Log2 Human Islet Proteomics") +
  ggtitle("Ligands") +
  scale_x_continuous(limits = c(-8,18)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(ligands_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4)), 
                  size = 3,
                  max.overlaps = 30,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



# receptors

# plot human islet bulk RNA-seq against proteomics
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(islet_tpm), 
             y = log2(all_islet_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 Human Islet RNA-seq",
       y = "Log2 Human Islet Proteomics") +
  ggtitle("Receptors") +
  scale_x_continuous(limits = c(-1,8)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(receptors_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 7)), 
                  size = 3,
                  max.overlaps = 20,
                  force = 1,
                  max.time = 2) +
  theme_cowplot()



# determine which genes have NA for one of the proteomics categories but not the other
# ligands
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  filter(!is.na(all_islet_proteomics)) %>%
  filter(is.na(sorted_sc_proteomics)) %>% 
  select(hgnc_symbol,
         sorted_sc_proteomics,
         unsorted_sc_proteomics,
         all_sc_proteomics,
         all_islet_proteomics) %>% 
  # save list of ligand proteins present in islet but missing in stem cell sorted
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_proteins_missing_in_stemcells.tsv")


# receptors
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  filter(!is.na(all_islet_proteomics)) %>%
  filter(is.na(sorted_sc_proteomics)) %>% 
  select(hgnc_symbol,
         sorted_sc_proteomics,
         unsorted_sc_proteomics,
         all_sc_proteomics,
         all_islet_proteomics) %>% 
  # save list of receptor proteins present in islet but missing in stem cell sorted
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_proteins_missing_in_stemcells.tsv")



# do the same for genes in stem cells but not in human islets
View(ligands_df %>% 
       filter(keep_in_list %in% c("Yes", "TBD")) %>% 
       filter(consensus_score > 7) %>% 
       filter(!is.na(all_sc_proteomics)) %>%
       filter(is.na(all_islet_proteomics)) %>% 
       select(hgnc_symbol,
              sorted_sc_proteomics,
              unsorted_sc_proteomics,
              all_sc_proteomics,
              all_islet_proteomics))
# no proteins detected for either ligands or receptors in sorted stem cell 
# & all stem cell that were not present in all human islet prep


########### plot all human islet vs sorted stem cell proteomics ##############

# create dataframe that replaces NA's in stem cell is pseudovalue for plotting
ligands_pseudo_for_NAs <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  # remove rows that are NA in all islet proteomics
  filter(!is.na(all_islet_proteomics)) %>%
  # convert NA's in sorted stem cell to psuedovalues to position on plot
  mutate_at(c("sorted_sc_proteomics"),
         ~replace(., is.na(.), 32.000)) # log2(32) = 5
  

# plot ligands with pseudovalues
ligands_pseudo_for_NAs %>% 
  ggplot(aes(x = log2(all_islet_proteomics), 
             y = log2(sorted_sc_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  # overwrite the colour for the pseudovalues
  geom_point(data = subset(ligands_pseudo_for_NAs, sorted_sc_proteomics == 32.0),
             colour = "darkorange1", size = 2) +
  labs(x = "Log2 Human Islet Proteomics",
       y = "Log2 Sorted SC Proteomics") +
  ggtitle("Ligands") +
  scale_y_continuous(limits = c(4,20)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(ligands_pseudo_for_NAs),
                  size = 3,
                  max.overlaps = 30,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



# receptors

# create dataframe that replaces NA's in stem cell is pseudovalue for plotting
receptors_pseudo_for_NAs <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  # remove rows that are NA in all islet proteomics
  filter(!is.na(all_islet_proteomics)) %>%
  # convert NA's in sorted stem cell to psuedovalues to position on plot
  mutate_at(c("sorted_sc_proteomics"),
            ~replace(., is.na(.), 32.000)) # log2(32) = 5


# plot receptors with pseudovalues
receptors_pseudo_for_NAs %>% 
  ggplot(aes(x = log2(all_islet_proteomics), 
             y = log2(sorted_sc_proteomics),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  # overwrite the colour for the pseudovalues
  geom_point(data = subset(receptors_pseudo_for_NAs, sorted_sc_proteomics == 32.0),
             colour = "darkorange1", size = 2) +
  labs(x = "Log2 Human Islet Proteomics",
       y = "Log2 Sorted SC Proteomics") +
  ggtitle("Receptors") +
#  scale_y_continuous(limits = c(4,18)) +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(data=subset(receptors_pseudo_for_NAs),
                  size = 3,
                  max.overlaps = 40,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



############ Plot consensus score histogram ################################

# regenerate receptors df b/c it's already been filtered to remove below
# consensus score of 5

# open old receptors list 
receptors_old <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors.tsv", 
                          sep = "\t")) %>% 
  # select only columns we need here
  select(genesymbol,
         consensus_score) %>% 
  drop_na(genesymbol) %>% 
  # group by hgnc symbol to merge duplicates
  group_by(genesymbol) %>% 
  # keep highest value for each hgnc symbol 
  summarise(consensus_score = max(consensus_score, na.rm = TRUE))



# Check how mayn receptors have a consensus score = 1
receptors_old %>% 
  filter(consensus_score == 1) %>% 
  nrow()




# plot histogram of consensus scores

#ligands
ligands_df  %>% 
  ggplot(aes(x = consensus_score)) +
  geom_histogram(bins = 15, fill = "#36226B") +
  geom_vline(xintercept = 5, color = 'black') +
  labs(x = "Omnipath consensus score") + 
  ggtitle("Ligands") +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_consensus_score.jpg")

# receptors
receptors_old  %>% 
  ggplot(aes(x = consensus_score)) +
  geom_histogram(bins = 15, fill = "#1F9274") +
  geom_vline(xintercept = 8, color = 'black') +
  labs(x = "Omnipath consensus score") + 
  ggtitle("Receptors") +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_consensus_score.jpg")


## facet wrap plot receptors and ligands histo side by side
# make combined dataframe
consensus_scores_lig <- ligands_df %>% 
  select(hgnc_symbol,
         consensus_score) %>% 
  mutate(type = "ligand") %>% 
  rename(genesymbol = hgnc_symbol)

consensus_scores_rec <- receptors_old %>% 
  select(genesymbol,
         consensus_score) %>% 
  mutate(type = "receptor")

consensus_scores_all <- full_join(consensus_scores_lig,
                                  consensus_scores_rec)


consensus_scores_all %>% 
  ggplot(aes(x = consensus_score, fill = type)) +
  geom_histogram(bins = 25) +
  facet_wrap("type",
             scales = "free_y") +
  geom_vline(data=filter(consensus_scores_all, type=="ligand"), aes(xintercept=4.3),
             color = 'black', linetype = "dashed") +
  geom_vline(data=filter(consensus_scores_all, type=="receptor"), aes(xintercept=7.3),
             color = 'black', linetype = "dashed") +
  # scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = "OmniPath consensus score",
       y= "count") + 
  scale_fill_manual(values = c("#36226B", "#1F9274")) +
  scale_x_continuous(breaks = c(0,10,20)) +
  # remove legend
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank(),axis.line=element_(color="black"))

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligrec_consensus_scores.png",
       scale = 1.5)



########### histogram of consensus scores showing manual filtering ##################


## facet wrap plot receptors and ligands histo side by side
# will stack included/excluded on top of each other in different colours

# make combined dataframe
consensus_scores_man_lig <- ligands_df %>% 
  select(hgnc_symbol,
         consensus_score,
         keep_in_list) %>% 
  mutate(type = "ligand") %>% 
  filter(consensus_score > 4) %>% 
  # add column to simplify whether was manually excluded or not
  mutate(manually_excluded = !(keep_in_list %in% c("Yes", "TBD")))

consensus_scores_man_rec <- receptors_df %>% 
  select(hgnc_symbol,
         consensus_score,
         keep_in_list) %>% 
  mutate(type = "receptor") %>% 
  filter(consensus_score > 7) %>% 
  # add column to simplify whether was manually excluded or not
  mutate(manually_excluded = !(keep_in_list %in% c("Yes", "TBD")))


consensus_scores_man_all <- consensus_scores_man_lig %>% 
  full_join(consensus_scores_man_rec)

consensus_scores_man_all %>% 
  ggplot(aes(x = consensus_score, fill = type)) +
  geom_histogram(bins = 17) +
  geom_histogram(data = subset(consensus_scores_man_all %>%
                                 filter(manually_excluded == "TRUE")),
                 fill = "white",
                 bins = 17) +
  geom_histogram(data = subset(consensus_scores_man_all %>%
                                 filter(manually_excluded == "TRUE")),
                 alpha = 0.3,
                 bins = 17,
                 colour = "black",
                 linewidth = 0.01) +
  facet_wrap("type",
             scales = "free") +
  # scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = "OmniPath consensus score",
       y= "count") + 
  scale_fill_manual(values = c("#36226B", "#1F9274")) +
 # scale_x_continuous(breaks = c(0,10,20)) +
  # remove legend
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank(),axis.line=element_line(color="black"))


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligrec_manual_exclusions_4.png",
       scale = 1.5)



############# plot histogram of mean islet tpm ################################

# plot histogram of islet_tpm

#ligands
ligands_df  %>% 
  # display only our final 422 ligands
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(islet_tpm))) +
  geom_histogram(bins = 20, fill = "#36226B") +
  labs(x = "Log2 Human islet TPM") + 
  ggtitle("Ligands") +
 # scale_x_continuous(limits = c(0,24)) +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_islet_tpm_histo.jpg")


# receptors
receptors_df  %>% 
  # display only our final 422 ligands
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(islet_tpm))) +
  geom_histogram(bins = 20, fill = "#1F9274") +
  labs(x = "Log2 Human islet TPM") + 
  ggtitle("Receptors") +
  # scale_x_continuous(limits = c(0,24)) +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_islet_tpm_histo.jpg")


## facet wrap plot receptors and ligands histo side by side

# make combined dataframe
islet_tpm_ligands <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         islet_tpm) %>% 
  mutate(type = "ligand")

islet_tpm_receptors <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         islet_tpm) %>% 
  mutate(type = "receptor")

islet_tpm_all <- islet_tpm_ligands %>% 
  full_join(islet_tpm_receptors)
  

islet_tpm_all %>% 
  ggplot(aes(x = log2(islet_tpm), fill = type)) +
  geom_histogram(bins = 20) +
  facet_wrap("type",
             scales = "free_y") +
  # scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = bquote(Log["2"]*("human islet TPM")),
       y= "count") + 
  scale_fill_manual(values = c("#36226B", "#1F9274")) +
  
  # scale_x_continuous(breaks = c(0,10,20)) +
  # remove legend
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank(),axis.line=element_line(color="black"))


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligrec_islet_lund_tpm_histo.png",
       scale = 1.7)




############### dot plot top lig/rec in bulk islet rna-seq #################

# open Johnson 2023 nulk islet RNAseq from howard
johnson_bulk_rna <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/RNA_islets_TPM_oneBatch_all_corrected.txt",
                                     sep = "\t", header = T) %>% 
  ## compute mean TPM
  # drop extra columns
  select(2,9:98) %>% 
  transmute(hgnc_symbol = external_gene_name,
            islet_tpm_johnson = rowMeans(select(., -external_gene_name))) 
  
  



#filter lists to just final lig/rec
ligands_filt <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         islet_tpm_lund = islet_tpm,
         gtex_mean,
         islet_specificity_lund = islet_specificity) %>% 
  # append johnson rna-seq
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") 

# replace NA's in johnson TPM with 0's
ligands_filt$islet_tpm_johnson <- 
  lapply(ligands_filt$islet_tpm_johnson,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_filt <- as.data.frame(lapply(ligands_filt, unlist))
  
  
# calculate johnson specificty
ligands_filt <- ligands_filt %>% 
  mutate(islet_specificity_johnson = (islet_tpm_johnson + 0.0000001)/(gtex_mean + 0.0000001)) %>% 
  # pivot longer
  select(hgnc_symbol,
         islet_tpm_lund,
         islet_tpm_johnson,
         islet_specificity_lund,
         islet_specificity_johnson) %>% 
  pivot_longer(cols = c(islet_tpm_lund,
                        islet_tpm_johnson,
                        islet_specificity_lund, 
                        islet_specificity_johnson),
               names_to = c("type", "group"),
               names_prefix = "islet_",
               names_sep = "_",
               values_to = "value") %>% 
  pivot_wider(id_cols = c("hgnc_symbol", "group"),
              names_from = "type") %>% 
  mutate(type = "ligand")


receptors_filt <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         islet_tpm_lund = islet_tpm,
         gtex_mean,
         islet_specificity_lund = islet_specificity) %>% 
  # append johnson rna-seq
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol")


# replace NA's in johnson TPM with 0's
receptors_filt$islet_tpm_johnson <- 
  lapply(receptors_filt$islet_tpm_johnson,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_filt <- as.data.frame(lapply(receptors_filt, unlist))

receptors_filt <- receptors_filt %>% 
  # calculate johnson specificty
  mutate(islet_specificity_johnson = (islet_tpm_johnson + 0.0000001)/(gtex_mean + 0.0000001)) %>% 
  # pivot longer
  select(hgnc_symbol,
         islet_tpm_lund,
         islet_tpm_johnson,
         islet_specificity_lund,
         islet_specificity_johnson) %>% 
  pivot_longer(cols = c(islet_tpm_lund,
                        islet_tpm_johnson,
                        islet_specificity_lund, 
                        islet_specificity_johnson),
               names_to = c("type", "group"),
               names_prefix = "islet_",
               names_sep = "_",
               values_to = "value") %>% 
  pivot_wider(id_cols = c("hgnc_symbol", "group"),
              names_from = "type") %>% 
  mutate(type = "receptor")


# get top ranked johnson ligands and receptors
johnson_bulk_rna_lig<- johnson_bulk_rna %>% 
  filter(hgnc_symbol %in% ligands_filt$hgnc_symbol) %>% 
  arrange(desc(islet_tpm_johnson)) %>% 
  head(25)

johnson_bulk_rna_rec<- johnson_bulk_rna %>% 
  filter(hgnc_symbol %in% receptors_filt$hgnc_symbol) %>% 
  arrange(desc(islet_tpm_johnson)) %>% 
  head(25)


# get top ranked lund ligands
lund_bulk_rna_lig <- ligands_filt %>% 
  filter(group == "lund") %>% 
  arrange(desc(tpm)) %>% 
  head(25)

# get top ranked lund receptors
lund_bulk_rna_rec <- receptors_filt %>% 
  filter(group == "lund") %>% 
  arrange(desc(tpm)) %>% 
  head(25)


# actually plot ligands and recptors seperately b/c 
# they need different colour and size scales

# plot ligands
ligands_filt %>% 
  filter(hgnc_symbol %in% johnson_bulk_rna_lig$hgnc_symbol) %>% 
  ggplot(aes(x = group,y = fct_reorder(hgnc_symbol, 1/tpm))) +
  geom_point(aes(size = log2(tpm),
                 color = log2(specificity))) +
  facet_wrap("type",
             scales = "free") +
  scale_size_continuous(range = c(0.01, 6)) +
  scale_color_distiller(palette = "BuPu",
                        direction = 1,
                        limits = c(-10,16)) +
  theme_cowplot() +
  # reverse alphabetical orders
  scale_y_discrete(limits = rev) +
  scale_x_discrete(labels = c("Johnson 2023", "Lund 2014")) +
  theme(axis.title = element_blank(),
        legend.box.spacing = unit(0, "cm"),
        axis.text.y = element_text(size = 9,
                                   colour = "black"),
        axis.text.x = element_text(size = 11,
                                   colour = "black",
                                   angle = 45,
                                   hjust = 1),
        legend.text = element_text(size = 8,
                                   colour = "black"),
        legend.title = element_text(size = 10,
                                    colour = "black"),
        legend.key.size = unit(0.5, "cm")) +
  labs(colour = bquote(Log["2"]*("islet specificty")),
       size = bquote(Log["2"]*("islet TPM"))) +
  # order the 2 legends
  guides(color = guide_colorbar(order = 0),
         size = guide_legend(order = 1))

# save 
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_top25injohnson_islet_tpm.png",
       scale = 1.7)

# receptors plot
receptors_filt %>% 
  filter(hgnc_symbol %in% johnson_bulk_rna_rec$hgnc_symbol) %>% 
  ggplot(aes(x = group,y = fct_reorder(hgnc_symbol,1/tpm))) +
  geom_point(aes(size = tpm,
                 color = log2(specificity))) +
  facet_wrap("type",
             scales = "free") +
  scale_size_continuous(range = c(.01, 6)) +
  scale_color_distiller(palette = "BuGn",
                        direction = 1) +
  theme_cowplot() +
  # reverse alphabetical orderS
  scale_y_discrete(limits = rev) +
  scale_x_discrete(labels = c("Johnson 2023", "Lund 2014")) +
  theme(axis.title = element_blank(),
        legend.box.spacing = unit(0, "cm"),
        axis.text.y = element_text(size = 9,
                                   colour = "black"),
        legend.text = element_text(size = 8,
                                   colour = "black"),
        axis.text.x = element_text(size = 11,
                                   colour = "black",
                                   angle = 45,
                                   hjust = 1),
        legend.title = element_text(size = 10,
                                    colour = "black"),
        legend.key.size = unit(0.5, "cm")) +
  labs(colour = bquote(Log["2"]*("islet specificty")),
       size = "islet TPM") +
  # order the 2 legends
  guides(color = guide_colorbar(order = 0),
         size = guide_legend(order = 1))



# save 
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_top25injohnson_islet_tpm2.png",
       scale = 1.7)

####  older plot of just lund
receptors_filt %>% 
  ggplot(aes(x = 1,y = fct_reorder(hgnc_symbol,1/islet_tpm))) +
  geom_point(aes(size = log2(islet_tpm),
                 color = log2(islet_specificity))) +
  facet_wrap("type",
             scales = "free") +
  scale_size_continuous(range = c(1, 5.5)) +
  scale_color_distiller(palette = "BuGn",
                        direction = 1) +
  theme_cowplot() +
  # reverse alphabetical orderS
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.box.spacing = unit(0, "cm"),
        axis.text.y = element_text(size = 9,
                                   colour = "black"),
        legend.text = element_text(size = 8,
                                   colour = "black"),
        legend.title = element_text(size = 10,
                                    colour = "black")) +
  labs(colour = bquote(Log["2"]*("islet TPM")),
       size = bquote(Log["2"]*("islet specificty")))

# save 
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptor_top25_islet_tpm1.png",
       scale = 1.3)


################### histo of johnson bulk rna #######################

# keep just ligands and receptors and group
johnson_bulk_rna_lig<- johnson_bulk_rna %>% 
  filter(hgnc_symbol %in% ligands_filt$hgnc_symbol) %>% 
  mutate(type = "ligand")
  
johnson_bulk_rna_rec<- johnson_bulk_rna %>% 
  filter(hgnc_symbol %in% receptors_filt$hgnc_symbol) %>% 
  mutate(type = "receptors")

johnson_bulk_rna_ligrec <- johnson_bulk_rna_lig %>% 
  full_join(johnson_bulk_rna_rec)

# plot
johnson_bulk_rna_ligrec %>% 
  ggplot(aes(x = log2(islet_tpm_johnson), fill = type)) +
  geom_histogram(bins = 20) +
  facet_wrap("type",
             scales = "free_y") +
  # scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = bquote(Log["2"]*("human islet TPM")),
       y= "count") + 
  scale_fill_manual(values = c("#36226B", "#1F9274")) +
  
  # scale_x_continuous(breaks = c(0,10,20)) +
  # remove legend
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank(),axis.line=element_line(color="black"))


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligrec_islet_johnson_tpm_histo.png",
       scale = 1.7)



######### scatter plot johnson vs lund tpm correlation ######################



ligands_lund_john <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  select(hgnc_symbol,
         islet_tpm_lund = islet_tpm) %>% 
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") %>% 
  mutate(type = "ligand") 


# replace NA's in johnson TPM with 0's
ligands_lund_john$islet_tpm_johnson <- 
  lapply(ligands_lund_john$islet_tpm_johnson,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
ligands_lund_john <- as.data.frame(lapply(ligands_lund_john, unlist))



receptors_lund_john <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  select(hgnc_symbol,
         islet_tpm_lund = islet_tpm) %>% 
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") %>% 
  mutate(type = "receptors") 


# replace NA's in johnson TPM with 0's
receptors_lund_john$islet_tpm_johnson <- 
  lapply(receptors_lund_john$islet_tpm_johnson,
         replace_na, replace = 0)  
# the lapply changed the format of the new column to a list which caused problems
# so convert all back to dataframe
receptors_lund_john <- as.data.frame(lapply(receptors_lund_john, unlist))


# plot lund vs johnson bulk rnaseq

# ligands 
ligands_lund_john %>% 
  ggplot(aes(x = log2(islet_tpm_lund), 
             y = log2(islet_tpm_johnson))) +
  geom_point(aes(colour = type), 
             size = 2, 
             alpha = 0.6) +
  labs(x = bquote(Log["2"]*("Lund 2014 TPM")),
       y = bquote(Log["2"]*("Johnson 2023 TPM"))) +
  # R coefficient
  stat_cor(method = "pearson", label.x.npc = 0.01, label.y.npc = 0.95, color = "red") +
  # trend line
  geom_smooth(method=lm, se = FALSE, color = "red", size = 1) +
  scale_color_manual(values = "#36226B") +
  facet_wrap("type",
             scales = "free") +
  guides(color = FALSE) +
  theme_cowplot()

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_lund_vs_johnson_scatter.png",
       scale = 2)

# receptors 
receptors_lund_john %>% 
  ggplot(aes(x = log2(islet_tpm_lund), 
             y = log2(islet_tpm_johnson))) +
  geom_point(aes(colour = type), 
             size = 2, 
             alpha = 0.6) +
  labs(x = bquote(Log["2"]*("Lund 2014 TPM")),
       y = bquote(Log["2"]*("Johnson 2023 TPM"))) +
  # R coefficient
  stat_cor(method = "pearson", label.x.npc = 0.01, label.y.npc = 0.95, color = "red") +
  # trend line
  geom_smooth(method=lm, se = FALSE, color = "red", size = 1) +
  scale_color_manual(values = "#1F9274") +
  facet_wrap("type",
             scales = "free") +
  guides(color = FALSE) +
  theme_cowplot()


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_lund_vs_johnson_scatter.png",
       scale = 2)




############# tables of top genes based on islet tpm ######################

ligands_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(Ligands = hgnc_symbol,
         islet_tpm) %>% 
  arrange(desc(islet_tpm)) %>% 
  head(25)

receptors_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(Ligands = hgnc_symbol,
         islet_tpm) %>% 
  arrange(desc(islet_tpm)) %>% 
  head(25)


############# plot eBCs vs Beta genes expression  ###############

#ligands
ligands_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  ggplot(aes(x = log2(Beta), 
             y = log2(sorted_eBCs),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 Primary Beta cells counts",
       y = "Log2 sorted eBCs counts") +
  ggtitle("Ligands") +
  geom_text_repel(data=subset(ligands_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 4),
                                log2(sorted_eBCs) > -10), # to label only certain points
                  size = 3,                  
                  max.overlaps = 12,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



# receptors
receptors_df %>%
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  ggplot(aes(x = log2(Beta), 
             y = log2(sorted_eBCs),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 Primary Beta cells counts",
       y = "Log2 sorted eBCs counts") +
  ggtitle("Receptors") +
  geom_text_repel(data=subset(receptors_df %>% 
                                filter(keep_in_list %in% c("Yes", "TBD")) %>% 
                                filter(consensus_score > 7),
                              log2(sorted_eBCs) > -10), # to label only certain points
                  size = 3,                  
                  max.overlaps = 7,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()







############ tables of top genes in Beta & eBCs by expression & specificity

### ligands

# top 25 expressed in Beta
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         Beta) %>% 
  arrange(desc(Beta)) %>% 
  head(25)

# top 25 specific in Beta compared to GTEx
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(beta_gtex_specificity = Beta/gtex_mean) %>% 
  select(hgnc_symbol,
         beta_gtex_specificity) %>% 
  arrange(desc(beta_gtex_specificity)) %>% 
  head(25)


# top 25 expressed in eBCs
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  select(hgnc_symbol,
         sorted_eBCs) %>% 
  arrange(desc(sorted_eBCs)) %>% 
  head(25)

# top 25 specific in eBCs compared to GTEx
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(sorted_eBCs_gtex_specificity = sorted_eBCs/gtex_mean) %>% 
  select(hgnc_symbol,
         sorted_eBCs_gtex_specificity) %>% 
  arrange(desc(sorted_eBCs_gtex_specificity)) %>% 
  head(25)



### receptors

# top 25 expressed in Beta
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         Beta) %>% 
  arrange(desc(Beta)) %>% 
  head(25)

# top 25 specific in Beta compared to GTEx
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(beta_gtex_specificity = Beta/gtex_mean) %>% 
  select(hgnc_symbol,
         beta_gtex_specificity) %>% 
  arrange(desc(beta_gtex_specificity)) %>% 
  head(25)


# top 25 expressed in eBCs
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>% 
  select(hgnc_symbol,
         sorted_eBCs) %>% 
  arrange(desc(sorted_eBCs)) %>% 
  head(25)

# top 25 specific in eBCs compared to GTEx
receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(sorted_eBCs_gtex_specificity = sorted_eBCs/gtex_mean) %>% 
  select(hgnc_symbol,
         sorted_eBCs_gtex_specificity) %>% 
  arrange(desc(sorted_eBCs_gtex_specificity)) %>% 
  head(25)



########## rank of ranks table #######################

### ligands


ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>% 
  
  mutate(sorted_eBCs_rank = dense_rank(desc(sorted_eBCs)),
         sorted_eBCs_specificity
         # had to add very small amounts to prevent infinity values
         # the values are so small that adding 0.01 screws up the data
         eBCs_human_beta_ratio = ((sorted_eBCs + 0.0000001) / (Beta + 0.0000001)),
         human_beta_eBCs_ratio = ((Beta + 0.0000001) / (sorted_eBCs + 0.0000001)))  
  
  
  select(islet_tpm_rank,
         islet_specificity_rank,
         beta_rank,
         beta_specificity_rank,
         sorted_eBCs_rank, 
         sorted_sc_proteomics_rank,
         all_islet_proteomics_rank,
         ) 

  
  
########### plot eBCs vs Beta difference in expression & GTEx specificity##############
  
#ligands
    
ligands_df_ratio_ranks <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(eBCs_human_beta_ratio = ((sorted_eBCs + 0.0000001) / (Beta + 0.0000001)),
         sorted_eBCs_gtex_specificity = ((sorted_eBCs + 0.0000001) / (gtex_mean + 0.0000001)),
         human_beta_gtex_specificity = ((Beta + 0.0000001) / (gtex_mean + 0.0000001)),
         eBCs_human_beta_gtex_specifictiy_ratio = ((sorted_eBCs_gtex_specificity + 0.0000001) / 
                                                     (human_beta_gtex_specificity + 0.0000001)),
         eBCs_human_beta_ratio_rank = dense_rank(desc(eBCs_human_beta_ratio)),
         eBCs_human_beta_gtex_specifictiy_ratio_rank = dense_rank(desc(eBCs_human_beta_gtex_specifictiy_ratio))) %>% 
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_human_beta_ratio,
         sorted_eBCs_gtex_specificity,
         human_beta_gtex_specificity,
         eBCs_human_beta_gtex_specifictiy_ratio,
         eBCs_human_beta_ratio_rank,
         eBCs_human_beta_gtex_specifictiy_ratio_rank) %>% 
  filter(sorted_eBCs != 0) %>%
  filter(Beta != 0)
  
  
# plot difference in specificity against difference in expression
ligands_df_ratio_ranks %>% 
      ggplot(aes(x = log2(eBCs_human_beta_ratio), 
               y = log2(eBCs_human_beta_gtex_specifictiy_ratio),
               label = hgnc_symbol)) +
    geom_point(colour = "#36226B", size = 2) +
    labs(x = "Log2 eBCs to beta-cell expression ratio",
         y = "Log2 eBCs to beta-cell specificity ratio") +
    ggtitle("Ligands") +
    geom_text_repel(size = 3,                  
                    max.overlaps = 30,
                    force = 1,
                    max.time = 10) +
    theme_cowplot()


# plot difference in expression against Beta expression
ligands_df_ratio_ranks %>% 
  ggplot(aes(x = log2(sorted_eBCs), 
             y = log2(eBCs_human_beta_ratio),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 beta cell expression",
       y = "Log2 eBCs to beta cell expression ratio") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()


# plot eBCs GTEx specificity vs beta-cell GTEx specificity
ligands_df_ratio_ranks %>% 
  ggplot(aes(x = log2(human_beta_gtex_specificity), 
             y = log2(sorted_eBCs_gtex_specificity),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 eBCs specificity",
       y = "Log2 beta cell specificity") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()




ligands_df_ratio_ranks <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(sorted_eBCs_rank = dense_rank(sorted_eBCs),
         beta_rank = dense_rank(Beta), 
         eBCs_human_beta_ratio = ((sorted_eBCs_rank) / (beta_rank)),
         sorted_eBCs_gtex_specificity = ((sorted_eBCs_rank) / (gtex_mean)),
         human_beta_gtex_specificity = ((beta_rank) / (gtex_mean)),
         eBCs_human_beta_gtex_specifictiy_ratio = ((sorted_eBCs_gtex_specificity) / 
                                                     (human_beta_gtex_specificity)),
         eBCs_human_beta_ratio_rank = dense_rank(desc(eBCs_human_beta_ratio)),
         eBCs_human_beta_gtex_specifictiy_ratio_rank = dense_rank(desc(eBCs_human_beta_gtex_specifictiy_ratio))) %>% 
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         beta_rank,
         sorted_eBCs_rank,
         eBCs_human_beta_ratio,
         sorted_eBCs_gtex_specificity,
         human_beta_gtex_specificity,
         eBCs_human_beta_gtex_specifictiy_ratio,
         eBCs_human_beta_ratio_rank,
         eBCs_human_beta_gtex_specifictiy_ratio_rank)



# plot difference in expression against Beta expression
ligands_df_ratio_ranks %>% 
  ggplot(aes(x = sorted_eBCs, 
             y = eBCs_human_beta_ratio,
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 beta cell expression",
       y = "Log2 eBCs to beta cell expression ratio") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()


# difference instead of ratio
ligands_df_difference <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_gtex_specif = ((sorted_eBCs + 0.0000001) / (gtex_mean + 0.0000001)),
         beta_gtex_specif = ((Beta + 0.0000001) / (gtex_mean + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_gtex_specif - beta_gtex_specif) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         eBCs_gtex_specif,
         beta_gtex_specif,
         eBCs_beta_specif_diff)

ligands_df_difference %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs beta specificity diff") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()

## receptors

# difference instead of ratio
receptors_df_difference <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_gtex_specif = ((sorted_eBCs + 0.0000001) / (gtex_mean + 0.0000001)),
         beta_gtex_specif = ((Beta + 0.0000001) / (gtex_mean + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_gtex_specif - beta_gtex_specif) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         eBCs_gtex_specif,
         beta_gtex_specif,
         eBCs_beta_specif_diff)

receptors_df_difference %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs beta specificity diff") +
  ggtitle("Receptors") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()

  
######### eBCs vs Beta difference in expression & in-sample specificity #########

# figures out which columns to average for sorted sample mean
ligands_df[48:52] %>% 
  colnames()

ligands_df$mean_sorted_scRNAseq <- rowMeans(ligands_df[48:52])

# difference instead of ratio
ligands_df_sample_diff <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_sample_specif = ((sorted_eBCs + 0.0000001) / (mean_sorted_scRNAseq + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_sample_specif - beta_specificity) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         mean_sorted_scRNAseq,
         mean_islet_scRNAseq,
         beta_specificity,
         eBCs_sample_specif,
         eBCs_beta_specif_diff)

# plot ligands
ligands_df_sample_diff %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#36226B", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs-primary-beta specificity diff") +
  ggtitle("Ligands") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()



## repeat for receptors

# figures out which columns to average for sorted sample mean
receptors_df[48:52] %>% 
  colnames()

receptors_df$mean_sorted_scRNAseq <- rowMeans(receptors_df[48:52])

# difference instead of ratio
receptors_df_sample_diff <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  mutate(eBCs_beta_diff = sorted_eBCs - Beta,
         eBCs_sample_specif = ((sorted_eBCs + 0.0000001) / (mean_sorted_scRNAseq + 0.0000001)),
         eBCs_beta_specif_diff = eBCs_sample_specif - beta_specificity) %>%
  select(hgnc_symbol,
         Beta,
         sorted_eBCs,
         eBCs_beta_diff,
         mean_sorted_scRNAseq,
         mean_islet_scRNAseq,
         beta_specificity,
         eBCs_sample_specif,
         eBCs_beta_specif_diff)

# plot receptors
receptors_df_sample_diff %>% 
  ggplot(aes(x = log2(eBCs_beta_diff),
             y = log2(eBCs_beta_specif_diff),
             label = hgnc_symbol)) +
  geom_point(colour = "#1F9274", size = 2) +
  labs(x = "Log2 eBCs-primary beta expr difference",
       y = "Log2 eBCs-primary-beta specificity diff") +
  ggtitle("Receptors") +
  geom_text_repel(size = 3,                  
                  max.overlaps = 10,
                  force = 1,
                  max.time = 10) +
  theme_cowplot()


############## Proteomics vs johnson bulk islet RNAseq(NEW HOWARD DATA) ##################

# open Johnson 2023 bulk islet RNAseq from howard
johnson_bulk_rna <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/RNA_islets_TPM_oneBatch_all_corrected.txt",
                               sep = "\t", header = T) %>% 
  ## compute mean TPM
  # drop extra columns
  select(2,9:98) %>% 
  transmute(hgnc_symbol = external_gene_name,
            islet_tpm_johnson = rowMeans(select(., -external_gene_name))) 


# get donor IDs from RNA to match the protein
rna_donor_ids <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/RNA_islets_TPM_oneBatch_all_corrected.txt",
           sep = "\t", header = T) %>%
  select(9:98) %>% 
  colnames() %>% 
  data.frame()


# # open Johnson 2023 n=90 human islet proteomics with AA length normalization
# proteomics_AAnorm <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/protein_islets_normLength.txt",
#                               sep = "\t", header = T) %>% 
#   select(hgnc_symbol = Genes,
#          AA_length = Length,
#          islet_proteomics_AAnorm = mean.signal.length)

# open johnson n=96 from grace that is donor matched but has 6 extra n's
# this data is better cause howard's had thresholded out some proteins 
proteomics_grace_96 <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/proteomics_matched_rnaseq_donors.txt",
                                     sep = "\t", header = T) %>% 
  select(-1) %>% 
  rename(hgnc_symbol = Genes)
                              

# keep only proteomics columns that are in RNA too
proteomics_grace_96_filt <- proteomics_grace_96 %>% 
  select(hgnc_symbol,
         # this selects column names that in the RNA IDs list
         any_of(rna_donor_ids$.)) %>% 
  # compute rowmeans
  transmute(hgnc_symbol = hgnc_symbol,
            proteomics_RNAmatched = rowMeans(select(., -hgnc_symbol),
                                             na.rm = TRUE)) 


# load AA length uniprot data
# use raw list b/c it has more proteins than howards filtered list
uniprot_size <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/uniprot_size.txt",
                           sep = "\t", header = T) %>% 
  select(hgnc_symbol = Genes,
         AA_length = Length)


# compute abundance normalized to length
proteomics_grace_96_filt <- proteomics_grace_96_filt %>% 
  left_join(uniprot_size) %>% 
  mutate(proteomics_RNAmatched_AAnorm = proteomics_RNAmatched/AA_length) %>% 
  # some genes have multiple AA size, so lets take the average
  select(hgnc_symbol,
         proteomics_RNAmatched_AAnorm) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(proteomics_RNAmatched_AAnorm_mean = mean(proteomics_RNAmatched_AAnorm))

ligands_prot <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  select(hgnc_symbol,
         ND_islet_proteomics,
         all_islet_proteomics) %>% 
  left_join(proteomics_grace_96_filt,
            by = "hgnc_symbol") %>% 
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") %>% 
  mutate(type = "ligand") 

receptors_prot <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  select(hgnc_symbol,
         ND_islet_proteomics,
         all_islet_proteomics) %>% 
  left_join(proteomics_grace_96_filt,
            by = "hgnc_symbol") %>% 
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") %>% 
  mutate(type = "receptor")

ligrec_prot <- ligands_prot %>% 
  full_join(receptors_prot)


# check number of lig/rec detected in proteomics
ligands_prot %>% 
  filter(is.na(proteomics_RNAmatched_AAnorm_mean) == FALSE) %>% 
  nrow()

receptors_prot %>% 
  filter(is.na(proteomics_RNAmatched_AAnorm_mean) == FALSE) %>% 
  nrow()

# plot human islet bulk RNA-seq against proteomics

# ligands 
ligands_prot %>% 
  ggplot(aes(x = log2(islet_tpm_johnson), 
             y = log2(proteomics_RNAmatched_AAnorm_mean),
             label = hgnc_symbol)) +
  geom_point(aes(colour = type), size = 2) +
  labs(x = bquote(Log["2"]*("RNA-seq TPM")),
       y = bquote(Log["2"]*("protein abundance"))) +
  scale_x_continuous(limits = c(-4,17)) +
  # R coefficient
  stat_cor(method = "pearson", label.x.npc = 0.01, label.y.npc = 0.95, color = "red") +
  # trend line
  geom_smooth(method=lm, se = FALSE, color = "red", size = 1) +
  scale_color_manual(values = "#36226B") +
  facet_wrap("type",
             scales = "free") +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(size = 3,
                  max.overlaps = 20,
                  force = 8,
                  max.time = 10,
                  segment.size = 0.2) +
  # remove legend
  guides(color = FALSE) +
  theme_cowplot()


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_protein_vs_rna_AAnorm_samedonors_2.png",
       scale = 2)




#receptors
receptors_prot %>% 
  ggplot(aes(x = log2(islet_tpm_johnson), 
             y = log2(proteomics_RNAmatched_AAnorm_mean),
             label = hgnc_symbol)) +
  geom_point(aes(colour = type), size = 2.5) +
  labs(x = bquote(Log["2"]*("RNA-seq TPM")),
       y = bquote(Log["2"]*("protein abundance"))) +
  scale_x_continuous(limits = c(-1.5,8)) +
  # R coefficient
  stat_cor(method = "pearson", label.x.npc = 0.00, label.y.npc = 0.93, color = "red") +
  # trend line
  geom_smooth(method=lm, se = FALSE, color = "red", size = 1) +
  scale_color_manual(values = "#1F9274") +
  facet_wrap("type",
             scales = "free") +
  geom_text_repel(size = 3,
                  max.overlaps = 10,
                  force = 8,
                  max.time = 10,
                  segment.size = 0.2) +
  # remove legend
  guides(color = FALSE) +
  theme_cowplot()

# save
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_protein_vs_rna_AAnorm_samedonors_2.png",
       scale = 2)



############ Proteomics vs Lund RNA-seq with new AA normalized proteomics #############

#
ligands_prot_lund <- ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  select(hgnc_symbol,
         ND_islet_proteomics,
         all_islet_proteomics,
         islet_tpm_lund = islet_tpm) %>% 
  left_join(proteomics_grace_96_filt,
            by = "hgnc_symbol") %>% 
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") %>% 
  mutate(type = "ligand") 

receptors_prot_lund <- receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  select(hgnc_symbol,
         ND_islet_proteomics,
         all_islet_proteomics,
         islet_tpm_lund = islet_tpm) %>% 
  left_join(proteomics_grace_96_filt,
            by = "hgnc_symbol") %>% 
  left_join(johnson_bulk_rna,
            by = "hgnc_symbol") %>% 
  mutate(type = "receptor")


# plot human islet bulk RNA-seq against proteomics

# ligands 
ligands_prot_lund %>% 
  ggplot(aes(x = log2(islet_tpm_lund), 
             y = log2(proteomics_RNAmatched_AAnorm_mean),
             label = hgnc_symbol)) +
  geom_point(aes(colour = type), size = 2) +
  labs(x = bquote(Log["2"]*("RNA-seq TPM")),
       y = bquote(Log["2"]*("protein abundance"))) +
  scale_x_continuous(limits = c(-4,17)) +
  # R coefficient
  stat_cor(method = "pearson", label.x.npc = 0.01, label.y.npc = 0.95, color = "red") +
  # trend line
  geom_smooth(method=lm, se = FALSE, color = "red", size = 1) +
  scale_color_manual(values = "#36226B") +
  facet_wrap("type",
             scales = "free") +
  # this uses the ggrepel to position
  #  geom_text_repel(data=subset(receptors_df %>% filter(keep_in_list == "Yes"), 
  #                              log2(islet_tpm) < 2), # to label only certain points
  #                  size = 3) +
  geom_text_repel(size = 3,
                  max.overlaps = 20,
                  force = 8,
                  max.time = 10,
                  segment.size = 0.2) +
  # remove legend
  guides(color = FALSE) +
  theme_cowplot()


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/ligands_protein_vs_rna_AAnorm_lund.png",
       scale = 2)




#receptors
receptors_prot_lund %>% 
  ggplot(aes(x = log2(islet_tpm_lund), 
             y = log2(proteomics_RNAmatched_AAnorm_mean),
             label = hgnc_symbol)) +
  geom_point(aes(colour = type), size = 2.5) +
  labs(x = bquote(Log["2"]*("RNA-seq TPM")),
       y = bquote(Log["2"]*("protein abundance"))) +
  scale_x_continuous(limits = c(-2.5,8)) +
  # R coefficient
  stat_cor(method = "pearson", label.x.npc = 0.00, label.y.npc = 0.93, color = "red") +
  # trend line
  geom_smooth(method=lm, se = FALSE, color = "red", size = 1) +
  scale_color_manual(values = "#1F9274") +
  facet_wrap("type",
             scales = "free") +
  geom_text_repel(size = 3,
                  max.overlaps = 10,
                  force = 8,
                  max.time = 10,
                  segment.size = 0.2) +
  # remove legend
  guides(color = FALSE) +
  theme_cowplot()

# save
ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/receptors_protein_vs_rna_AAnorm_lund.png",
       scale = 2)






#################### group by gene family #############################

ligands_df %>% 
  group_by(Protein.families) %>% 
  summarise(n()) %>% 
  arrange(desc(`n()`)) %>% 
  View()
        


##################### lig/rec with TPM > 1 #####################################

# lund 2014
ligands_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 4) %>%
  filter(islet_tpm > 1) %>% 
  nrow()

receptors_df %>% 
  filter(keep_in_list %in% c("Yes", "TBD")) %>% 
  filter(consensus_score > 7) %>%
  filter(islet_tpm > 1) %>% 
  nrow()

# check johnson 2023
ligands_prot %>% 
  filter(islet_tpm_johnson > 1) %>% 
  nrow()

receptors_prot %>% 
  filter(islet_tpm_johnson > 1) %>% 
  nrow()



