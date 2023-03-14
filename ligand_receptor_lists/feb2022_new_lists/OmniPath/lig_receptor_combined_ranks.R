# combine our filtered list of interactions with our ligand & receptor ranks
# to form an aggregated final rank for ligands to order


library(tidyverse)
library(pheatmap)
library(cowplot)
library(OmnipathR)

############# load data ###############################


ligand_scores <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligand_scores_3.tsv", 
                        sep = "\t")) %>% 
  rename(source_genesymbol = hgnc_symbol)

receptor_scores <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_scores_3.tsv", 
                            sep = "\t")) %>% 
  rename(target_genesymbol = hgnc_symbol)

interactions_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_interactions_filtered.tsv", 
                             sep = "\t"))




############# append tables and determine interaction scores ##############

# append ranks to interactions df
interactions_df <- interactions_df %>% 
  left_join(ligand_scores,
            by = "source_genesymbol") %>% 
  left_join(receptor_scores,
            by = "target_genesymbol",
            suffix = c(".lig", ".rec")) %>% # to differentiate the ligand from receptor ranks
  # create interaction rank by multiplying ligand and receptor ranks
  # taking the sum produces very similar results to the product
  mutate(interaction_score = aggregate_score.lig * aggregate_score.rec) %>% 
  # then compute the percent rank
  mutate(interaction_score = percent_rank(interaction_score))


# create column of combined interaction names
interactions_df <- interactions_df %>% 
  mutate(interaction_symbol = paste(source_genesymbol,
                                    target_genesymbol,
                                    sep = "-")) 

# save table with interaction references & duplicate genes
write_tsv(interactions_df,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/interaction_ranks_with_duplicates_references.tsv")


# there are some duplicate interactions with different uniprot IDS
# count number of duplicates
interactions_df %>%
  group_by(interaction_symbol) %>% 
  summarise(n())
# 1569 rows, but only 1552 are unique. So there are 17 duplicates to remove

# remove uniprot IDs to remove duplicate interactions  
interactions_df <- interactions_df %>% 
  # drop uniprot ID columns
  select(-source,
         -target) %>% 
    #     -is_directed,
    #     -is_stimulation,
    #     -is_inhibition,
    #     -consensus_direction,
    #     -consensus_stimulation,
    #     -consensus_inhibition,
    #     -sources,
    #     -references,
    #     -curation_effort,
    #     -n_references,
    #     -n_resources
  # move interaction name to first column
  relocate(interaction_symbol,
           .before = source_genesymbol) %>% 
  # remove rows that are completely duplicates of another
  # there were 7 rows that had same interaction symbol but were different 
  # in other omnipath columns. Keep just the first row in these cases
  distinct(interaction_symbol,
           .keep_all = TRUE)

interactions_df %>%
  group_by(interaction_symbol) %>% 
  summarise(n())


# save final processed interaction table without ranks and with duplicates removed
interactions_df %>% 
  select(1:14) %>% 
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_interactions_filtered_noduplicates_final.tsv")



# save interaction_df without duplicates & just the main columns
interactions_df_simple <- interactions_df %>% 
  select(interaction = interaction_symbol,
         ligand = source_genesymbol,
         receptor = target_genesymbol,
         ligand_rank = aggregate_score.lig,
         receptor_rank = aggregate_score.rec,
         interaction_rank = interaction_score) %>% 
  arrange(desc(interaction_rank))

write_tsv(interactions_df_simple,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/interaction_ranks.tsv")


#reopen
interactions_df_simple <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/interaction_ranks.tsv", 
sep = "\t"))

################# histograms of genome wide interactions data #######################

# plot histogram of the interaction score
interactions_df %>% 
  ggplot(aes(x = interaction_score)) +
  geom_histogram(bins = 20, fill = "cadetblue") +
  scale_y_continuous(breaks = c(0, 50, 100, 150)) +
  labs(x = "interaction scores",
       y= "count") + 
  theme_cowplot()

# plot scatterplot of interaction score vs number of resources
interactions_df %>% 
  ggplot(aes(x = n_resources, y = interaction_score)) +
  geom_point()

# noticed that our top ranked receptor, PLXNC1, is not in the interactions_df
# this is b/c it's ligands SEMA7A & SEMA6D are not on our ligands list at all
# likely there are other top genes with similar issues
# in this case appears is b/c they are membrane bound ligands

# check how many of our ligands are missing from the interactions_df
ligand_scores %>%
  select(source_genesymbol,
         aggregate_score) %>% 
  mutate(in_interactions = (source_genesymbol %in% interactions_df$source_genesymbol)) %>% 
  group_by(in_interactions) %>% 
  summarise(n())
 # 27 of the 422 ligands are missing from interactions_df
 # including 3 with very high rankings: PIP, CGB5, CXCL17

# check how many of our receptors are missing from the interactions_df
receptor_scores %>%
  select(target_genesymbol,
         aggregate_score) %>% 
  mutate(in_interactions = (target_genesymbol %in% interactions_df$target_genesymbol)) %>% 
 # View()
  group_by(in_interactions) %>% 
  summarise(n())
# 26 of the 349 receptors are missing from interactions_df
# including 8 with very high rankings: PLXNC1, PLXNB1, LINGO1, NCSTN, LRRC4C,
# RTN4RL1, GPR37, RTN4R 


  
# create histogram of number of interactions per ligand 
interactions_count <- interactions_df %>% 
  group_by(source_genesymbol) %>% 
  summarise(n()) %>% 
  # join with ligands_df to get the ligands with 0 interactions
  # for plotting purposes
  right_join(ligand_scores %>%
              select(source_genesymbol),
            by = "source_genesymbol") %>% 
  # replace NA's with 0
  mutate(`n()` = replace_na(`n()`, 0))

interactions_count %>% 
  ggplot(aes(x = `n()`)) +
  geom_histogram(bins = 15, fill = "#36226B") +
  #  geom_vline(xintercept = 3, color = 'black') +
  scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = "interactions per ligand",
       y= "count") + 
  theme_cowplot()
  

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interactions_per_ligand_histo_3.jpg",
       scale = 1.5)



# create histogram of number of interactions per receptor 
interactions_count_rec <- interactions_df %>% 
  group_by(target_genesymbol) %>% 
  summarise(n()) %>% 
  # join with receptors_df to get the receptors with 0 interactions
  # for plotting purposes
  right_join(receptor_scores %>%
               select(target_genesymbol),
             by = "target_genesymbol") %>% 
  # replace NA's with 0
  mutate(`n()` = replace_na(`n()`, 0))

interactions_count_rec %>% 
  ggplot(aes(x = `n()`)) +
  geom_histogram(bins = 15, fill = "#1F9274") +
  #  geom_vline(xintercept = 3, color = 'black') +
 # scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = "interactions per receptor",
       y= "count") + 
  theme_bw()


ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interactions_per_receptor_histo.jpg",
       scale = 1.5)



# histo of both receptors and ligands per interaction on same plot
interactions_count <-interactions_count %>% 
  mutate(type = "ligand") %>% 
  rename(genesymbol = source_genesymbol)

interactions_count_rec <-interactions_count_rec %>% 
  mutate(type = "receptor") %>% 
  rename(genesymbol = target_genesymbol)

interactions_count_all <- full_join(interactions_count,
                                    interactions_count_rec)



# check average number of interactions
interactions_count_all %>% 
  filter(type == "receptor") %>% 
  summarise(mean(`n()`))


interactions_count_all %>% 
  ggplot(aes(x = `n()`, fill = type)) +
  geom_histogram(bins = 20) +
  facet_wrap("type") +
  #  geom_vline(xintercept = 3, color = 'black') +
  # scale_y_continuous(breaks = c(0, 40, 80)) +
  labs(x = "interactions per ligand or receptor",
       y= "count") + 
  scale_fill_manual(values = c("#36226B", "#1F9274")) +
  scale_x_continuous(breaks = c(0,20,40)) +
  # remove legend
  guides(fill = FALSE) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank(),axis.line=element_line(color="black"))

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interactions_per_ligrec_3.png",
       scale = 1.5)




# bar plot of stimulation vs inhibition per reaction
## never mind this data is messy. Some are annotated as both, others neither.
## there is both is_ columns and consensus_ columns



# create histogram of the number of resources which identified each interaction
interactions_df %>% 
  ggplot(aes(x = n_resources)) +
  geom_histogram(bins = 20, fill = "gray15") +
  #  geom_vline(xintercept = 3, color = 'black') +
  scale_y_continuous(breaks = c(0, 150, 300)) +
  labs(x = "resources per interaction",
       y= "count") + 
  theme_cowplot()

# based on the histogram, may make sense to make a threshold of 2 or 3 for inclusion

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_resources_histo_2.jpg",
       scale = 1.5)



### create histogram of the number of references which identified each interaction

# first get number of unique references per interaction
### The interactions are stored into a data frame.
interactions <- import_omnipath_interactions(resources = NULL,
                                             organism = 9606,
                                             # this removes the duplicates
                                             references_by_resource = FALSE) %>% 
  select(source_genesymbol,
         target_genesymbol,
         references) 

interactions <- interactions %>% 
  mutate(n_references = str_split(references, pattern = ";")) %>% 
  rowwise() %>% 
  mutate(n_references_count = length(n_references)) %>% 
  ungroup() %>% 
  mutate(interaction_symbol = paste(source_genesymbol,
                                    target_genesymbol,
                                    sep = "-")) %>% 
  select(interaction_symbol,
         n_references_count)
  
interactions_df <- interactions_df %>% 
  left_join(interactions,
            by = "interaction_symbol") %>% 
  distinct(interaction_symbol,
           .keep_all = TRUE)


interactions_df %>% 
  ggplot(aes(x = n_references_count.x)) +
  geom_histogram(bins = 20, fill = "gray15") +
  #  geom_vline(xintercept = 3, color = 'black') +
  scale_y_continuous(breaks = c(0, 150, 300)) +
  labs(x = "resources per interaction",
       y= "count") + 
  theme_cowplot()
## ABORT THIS ISNT HELPFUL


interactions_df %>% 
  ggplot(aes(x = n_resources)) +
  geom_histogram(bins = 20, fill = "gray15") +
  #  geom_vline(xintercept = 3, color = 'black') +
  scale_y_continuous(breaks = c(0, 150, 300)) +
  labs(x = "resources per interaction",
       y= "count") + 
  theme_cowplot()

# based on the histogram, may make sense to make a threshold of 2 or 3 for inclusion

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_resources_histo_3.png",
       scale = 1.5)





################# line plot of interactions ##########################
interactions_df %>% 
  head(50) %>% 
  ggplot() +
  geom_text(aes(x = 1, y = interaction_symbol, label=source_genesymbol), nudge_x=-0.05, size=3, color="black") +
  geom_point(aes(x = 1, y = interaction_symbol), size=1) +
  geom_point(aes(x = 2, y = interaction_symbol), size=1) +
  geom_segment(aes(x=1, y=interaction_symbol, xend=2, yend=interaction_symbol), alpha = 0.2) +
  theme_void()




dat <- interactions_df %>% 
  select(1:3, n_resources) %>% 
  pivot_longer(cols = 2:3, names_to = c("type", "tmp"),
               names_sep = "_",
               values_to = "gene_symbol") %>% 
  select(-tmp)


ggplot(dat, aes(x = type, y = gene_symbol)) +
  geom_point() +
  geom_line(aes(group = interaction_symbol, colour = n_resources),
            alpha = 0.1) +
  theme_void()


####### max score for each ligand across all its interactions ##############

overall_ligand_ranks <- interactions_df %>% 
  select(ligand = source_genesymbol,
         ligand_rank = aggregate_score.lig,
         highest_receptor_rank = aggregate_score.rec,
         overall_ligand_rank = interaction_score) %>% 
  group_by(ligand) %>% 
  summarise_all(max)


# add the ligands that did not have any interactions with our receptors
# back to the table
overall_ligand_ranks <- overall_ligand_ranks %>% 
  right_join(ligand_scores %>% select(ligand = source_genesymbol,
                                         ligand_rank = aggregate_score),
            by = "ligand") %>% 
  # remove duplicate column & reorder
  select(ligand,
         ligand_rank = ligand_rank.y,
         highest_receptor_rank,
         overall_ligand_rank)

# save table
write_tsv(overall_ligand_ranks,
          "ligand_receptor_lists/feb2022_new_lists/OmniPath/data/overall_ligand_ranks.tsv")







############### heatmap to visualize interaction ranks #######################

interaction_ranks_rowname <- interactions_df %>% 
  select(interaction_symbol,
         ligand_rank = aggregate_score.lig,
         receptor_rank = aggregate_score.rec,
         interaction_rank = interaction_score)

interaction_ranks_rowname <- column_to_rownames(as.data.frame(interaction_ranks_rowname), 
                                             var = "interaction_symbol")


# heatmap of ranks
interaction_ranks_rowname %>%
  arrange(desc(interaction_rank)) %>% 
  #head(50) %>% 
  pheatmap(labels_col = c("Ligand rank",
                          "Receptor rank",
                          "Interaction rank"),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize_col = 12,
           fontsize_row = 7,
           angle_col = 90)








  
############### heatmap to visualize overall ligand ranks #######################

overall_ligand_rowname <- overall_ligand_ranks
overall_ligand_rowname <- column_to_rownames(as.data.frame(overall_ligand_ranks), 
                                               var = "ligand")


# heatmap of ranks
overall_ligand_rowname %>%
  arrange(desc(max_combined_rank)) %>% 
  #head(50) %>% 
  pheatmap(labels_col = c("Ligand rank",
                          "Max receptor rank",
                          "Overall ligand rank"),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize_col = 12,
           fontsize_row = 7,
           angle_col = 90)






