## plot ligand-receptor interactions using lines connecting dots
## using ggplot

library(tidyverse)
library(ggrepel)

theme_set(theme_classic())


# load interactions data with all 1552 interactions & pivot longer
# 2 rows for each interaction
interactions_df <- read.table("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_interactions_filtered_noduplicates_final.tsv",
                  header = T) %>% 
  select(1:3, n_resources) %>% 
  pivot_longer(cols = 2:3, names_to = c("type", "tmp"),
               names_sep = "_",
               values_to = "gene_symbol") %>% 
  select(-tmp) %>% 
  mutate(type = if_else(type == "source",
                        "ligands", "receptors"))
head(interactions_df)

### Plot all interactions
interactions_df <-  interactions_df %>% 
  group_by(type, gene_symbol) %>% 
  mutate(n_interactions = n())

ligands <- filter(interactions_df, type == "ligands") %>% 
  group_by(gene_symbol) %>% 
  slice(1) %>% 
  arrange(gene_symbol) 

recept <- filter(interactions_df, type == "receptors") %>% 
  group_by(gene_symbol) %>% 
  slice(1) %>% 
  arrange(gene_symbol)

# add row number into column to use for positioning on y-axis evenly
ligands <- tibble::rowid_to_column(ligands, "ID") %>%
  mutate(ID = (396-ID),
         ID = ID/396)
recept <- tibble::rowid_to_column(recept, "ID") %>% 
  mutate(ID = (323-ID),
         ID = ID/323)


# add column IDs to interactions_df
ligands_recept <- full_join(ligands,
                            recept)

interactions_df <- interactions_df %>% 
  left_join(ligands_recept %>% 
              select(gene_symbol,
                     ID),
            by = "gene_symbol")


# Plot
ggplot(interactions_df, aes(x = type, y = gene_symbol)) +
  geom_line(aes(group = interaction_symbol,
                colour = n_resources),
            alpha = 0.3, size = 0.3) +
  geom_point(data = ligands, aes(size = n_interactions),
             alpha = 0.6, colour = "#36226B") +
  geom_point(data = recept, aes(size = n_interactions),
             alpha = 0.6, colour = "#1F9274") +
  scale_color_fermenter(palette = "YlGnBu",
                        direction = 1,
                        breaks = c(1,3,5,10,15)) +
  scale_size_binned(range = c(0.1, 8),
                    breaks = c(5,10,15)) +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.box.spacing = unit(-2, "cm"),
        axis.text.x = element_text(size = 14, 
                                   colour = "black"),
        legend.text = element_text(size = 12, 
                                   colour = "black")) +
  guides(size = guide_legend(override.aes = list(color = "grey"))) +
  expand_limits(y = c(-15, 750)) +
  labs(colour = "Resources",
       size = "Interactions")

ggsave("interact_plot.png", width = 5, height = 10)

#################### plot with top lig/rec labelled ###################
labels_recept <- interactions_df %>%
  filter(type == "receptors",
         n_resources > 5) %>% 
  group_by(gene_symbol) %>% 
  summarise(n_interactions = unique(n_interactions),
            type = "receptors") %>% 
  left_join(ligands_recept %>% 
              select(gene_symbol,
                      ID)) %>% 
  arrange(desc(n_interactions)) %>%
  head(25)

labels_ligands <- interactions_df %>%
  filter(type == "ligands",
         n_resources > 5) %>% 
  group_by(gene_symbol) %>% 
  summarise(n_interactions = unique(n_interactions),
            type = "ligands") %>% 
  left_join(ligands_recept %>% 
              select(gene_symbol,
                     ID)) %>% 
  arrange(desc(n_interactions)) %>%
  head(25)




ggplot(interactions_df, aes(x = type, y = ID)) +
  geom_line(aes(group = interaction_symbol,
                colour = n_resources),
            alpha = 0.2, size = 0.8) +
  geom_point(data = ligands, aes(size = n_interactions),
             alpha = 0.6, colour = "#36226B") +
  geom_point(data = recept, aes(size = n_interactions),
             alpha = 0.6, colour = "#1F9274") +
  geom_text_repel(data = labels_recept, aes(label = gene_symbol),
                  nudge_x = 0.1, size = 5, hjust = 0,
                  direction = "y", box.padding = 0.04,
                  segment.size = 0.25) +
  geom_text_repel(data = labels_ligands, aes(label = gene_symbol),
                  nudge_x = -0.1, size = 5, hjust = 1,
                  direction = "y", box.padding = 0.04,
                  segment.size = 0.25) +
  # reverse y axis
#  scale_y_continuous(limits = rel) +
  scale_color_fermenter(palette = "YlGnBu",
                        direction = 1,
                        breaks = c(1,3,5,10,15)) +
  scale_size_binned(range = c(0.1, 13),
                    breaks = c(5,10,15)) +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.box.spacing = unit(-1, "cm"),
        axis.text.x = element_text(size = 18,
                                   colour = "black",
                                   face = "bold"),
        legend.text = element_text(size = 14,
                                   colour = "black"),
        legend.title = element_text(size = 16,
                                    colour = "black")) +
  guides(size = guide_legend(override.aes = list(color = "grey"))) +
 # expand_limits(y = c(-15, 750)) +
  scale_x_discrete(position = "top") +
  labs(colour = "Resources",
       size = "Interactions")

ggsave("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interactions_connections2.png", 
       width = 8, height = 16)

#basic:
ggplot(interactions_df, aes(x = type, y = gene_symbol)) +
  geom_point() +
  geom_line(aes(group = interaction_symbol),
            alpha = 0.2)  +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())

### Plot only interactions with >= 20 resources
interactions_df_20res <- interactions_df %>% 
  group_by(interaction_symbol) %>% 
  # arrange(desc(n_resources)) %>% 
  # filter(n_resources >= 20) %>% 
  group_by(type, gene_symbol) %>% 
  mutate(n_interactions = n()) # this is only number of interactions within top, not all interactions

# for labels:
interactions_df_20res_source <- filter(interactions_df_20res, type == "source") %>% 
  group_by(gene_symbol) %>% 
  slice(1)
interactions_df_20res_target <- filter(interactions_df_20res, type == "target") %>% 
  group_by(gene_symbol) %>% 
  slice(1)

# Plot
ggplot(interactions_df_20res, aes(x = type, y = gene_symbol)) +
  geom_point(aes(size = n_interactions,
                 colour = n_resources),
             alpha = 0.7) +
  geom_text(interactions_dfa = interactions_df_20res_source,
            aes(label = gene_symbol),
            nudge_x = -0.05, size = 1.4,
            hjust = 1) +
  geom_text(interactions_dfa = interactions_df_20res_target,
            aes(label = gene_symbol),
            nudge_x = 0.05, size = 1.4,
            hjust = 0) +
  geom_line(aes(group = interaction_symbol,
                colour = n_resources),
            alpha = 0.1) +
  scale_colour_gradient(low = "lightblue", high = "darkblue") +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
  
