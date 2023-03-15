# new combined ligand and receptor ranks 
# ranking for all 3 combinations


library(tidyverse)
library(ComplexHeatmap)
library(cowplot)



############### load data and keep just essential columns ###############################

#### ligands
ligand_scores_abundant <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_ranks_abundant.tsv", 
                           sep = "\t")) %>% 
  rename(ligand = hgnc_symbol) %>% 
  select(ligand, 
         ligand_abundant_rank = aggregate_score)

ligand_scores_scarce <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligands_ranks_scarce.tsv", 
                                    sep = "\t")) %>% 
  rename(ligand = hgnc_symbol) %>% 
  select(ligand, 
         ligand_scarce_rank = aggregate_score)


##### receptors
receptor_scores_abundant <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_ranks_abundant.tsv", 
                                    sep = "\t")) %>% 
  rename(receptor = hgnc_symbol) %>% 
  select(receptor, 
         receptor_abundant_rank = aggregate_score)

receptor_scores_scarce <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/receptors_ranks_scarce.tsv", 
                                  sep = "\t")) %>% 
  rename(receptor = hgnc_symbol) %>% 
  select(receptor, 
         receptor_scarce_rank = aggregate_score)



##### final filtered interactions with duplicates removed
interactions_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/ligrec_interactions_filtered_noduplicates_final.tsv", 
                             sep = "\t")) %>% 
  # keep just first 3 columns with basic info
  select(1:3) %>%
  rename(ligand = source_genesymbol,
         receptor = target_genesymbol)



################# create 1 dataframe with all interaction rank combos ####################

interaction_ranks <- interactions_df %>% 
  # add ligand ranks
  left_join(ligand_scores_scarce,
            by = "ligand") %>% 
  left_join(ligand_scores_abundant,
            by = "ligand") %>% 
  # add receptor ranks
  left_join(receptor_scores_scarce,
            by = "receptor") %>% 
  left_join(receptor_scores_abundant,
            by = "receptor") %>% 
  # create combined ligrec rank columns
  mutate(scarce_lig_abund_rec_rank = percent_rank(ligand_scarce_rank * receptor_abundant_rank),
         abund_lig_abund_rec_rank = percent_rank(ligand_abundant_rank * receptor_abundant_rank),
         abund_lig_scarce_rec_rank = percent_rank(ligand_abundant_rank * receptor_scarce_rank))


# save dataframe with all interaction ranks
interaction_ranks %>% 
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/interaction_ranks_all3ways.tsv")



# overall top scarce ligands with abund receptors

overal_scarce_lig_abund_rec_ranks <- interaction_ranks %>% 
  select(ligand,
         ligand_scarce_rank,
         receptor_abundant_rank,
         scarce_lig_abund_rec_rank) %>% 
  group_by(ligand) %>% 
  summarise_all(max) %>% 
  arrange(desc(scarce_lig_abund_rec_rank))

# save
overal_scarce_lig_abund_rec_ranks %>% 
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/scarce_lig_abund_rec_ranks.tsv")




# overall top abund ligands with abund receptors

overal_abund_lig_abund_rec_ranks <- interaction_ranks %>% 
  select(receptor,
         ligand_abundant_rank,
         receptor_abundant_rank,
         abund_lig_abund_rec_rank) %>% 
  group_by(receptor) %>% 
  summarise_all(max) %>% 
  arrange(desc(abund_lig_abund_rec_rank))

# save
overal_abund_lig_abund_rec_ranks %>% 
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/abund_lig_abund_rec_ranks.tsv")



# overall top abund ligands with scarce receptors

overal_abund_lig_scarce_rec_ranks <- interaction_ranks %>% 
  select(receptor,
         ligand_abundant_rank,
         receptor_scarce_rank,
         abund_lig_scarce_rec_rank) %>% 
  group_by(receptor) %>% 
  summarise_all(max) %>% 
  arrange(desc(abund_lig_scarce_rec_rank))

# save
overal_abund_lig_scarce_rec_ranks %>% 
  write_tsv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/abund_lig_scarce_rec_ranks.tsv")


################## heatmap prep ############################


## first try scarce lig with abund rec

interaction_ranks_rowname <- interaction_ranks
interaction_ranks_rowname <- column_to_rownames(as.data.frame(interaction_ranks), 
                                               var = "interaction_symbol")


############## heatmap scarce lig abund rec ###########################


heatmap_scarce_lig_interactions <- interaction_ranks_rowname %>% 
  arrange(desc(scarce_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(ligand_scarce_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuPu"),
          heatmap_legend_param = list(title = "Scarce ligand rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top interactions",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))


heatmap_scarce_lig_interactions

heatmap_scarce_lig_interactions_rec <- interaction_ranks_rowname %>% 
  arrange(desc(scarce_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(receptor_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Abundant receptor\nrank"),
          border = TRUE,
          height = unit(18, "cm"),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))


heatmap_scarce_lig_interactions_rec

ht1 = heatmap_scarce_lig_interactions + heatmap_scarce_lig_interactions_rec

ht1

png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_scores/top50_interactions_scarce_lig_abund_rec.png",
    width = 4,
    height = 8,
    units = "in",
    res = 300)

draw(ht1,
     ht_gap = unit(0.3, "cm"))

dev.off()



# overall top scarce ligands with abund receptors

overal_scarce_lig_abund_rec_ranks <- interaction_ranks %>% 
  select(ligand,
         ligand_scarce_rank,
         receptor_abundant_rank,
         scarce_lig_abund_rec_rank) %>% 
  group_by(ligand) %>% 
  summarise_all(max) %>% 
  arrange(desc(scarce_lig_abund_rec_rank))

overal_scarce_lig_abund_rec_ranks <- column_to_rownames(as.data.frame(overal_scarce_lig_abund_rec_ranks), 
                                                var = "ligand")

heatmap_overal_scarce_lig_abund <- overal_scarce_lig_abund_rec_ranks %>% 
  arrange(desc(scarce_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(ligand_scarce_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuPu"),
          heatmap_legend_param = list(title = "Scarce ligand rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top overall ligands",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))
heatmap_overal_scarce_lig_abund


heatmap_overal_scarce_lig_abund_rec <- overal_scarce_lig_abund_rec_ranks %>% 
  arrange(desc(scarce_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(receptor_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Highest abundant\nreceptor rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top overall ligands",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))

heatmap_overal_scarce_lig_abund_rec



ht2 = heatmap_overal_scarce_lig_abund + heatmap_overal_scarce_lig_abund_rec

ht2


png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_scores/top50_overall_scarce_lig_abund_rec.png",
    width = 4,
    height = 8,
    units = "in",
    res = 300)

draw(ht2,
     ht_gap = unit(0.3, "cm"))

dev.off()



################### heatmap abund lig abund rec ####################################


heatmap_abund_ligrec_interactions <- interaction_ranks_rowname %>% 
  arrange(desc(abund_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(ligand_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuPu"),
          heatmap_legend_param = list(title = "Abundant ligand rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top interactions",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))
heatmap_abund_ligrec_interactions


heatmap_abund_ligrec_interactions_rec <- interaction_ranks_rowname %>% 
  arrange(desc(abund_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(receptor_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Abundant receptor\nrank"),
          border = TRUE,
          height = unit(18, "cm"),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))


heatmap_scarce_lig_interactions_rec

ht3 = heatmap_abund_ligrec_interactions + heatmap_scarce_lig_interactions_rec

ht3

png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_scores/top50_interactions_abund_lig_abund_rec.png",
    width = 4,
    height = 8,
    units = "in",
    res = 300)

draw(ht3,
     ht_gap = unit(0.3, "cm"))

dev.off()


# overall top abund ligands with abund receptors

overal_abund_lig_abund_rec_ranks <- interaction_ranks %>% 
  select(receptor,
         ligand_abundant_rank,
         receptor_abundant_rank,
         abund_lig_abund_rec_rank) %>% 
  group_by(receptor) %>% 
  summarise_all(max) %>% 
  arrange(desc(abund_lig_abund_rec_rank))

overal_abund_lig_abund_rec_ranks <- column_to_rownames(as.data.frame(overal_abund_lig_abund_rec_ranks), 
                                                        var = "receptor")

heatmap_overal_abund_lig_abund <- overal_abund_lig_abund_rec_ranks %>% 
  arrange(desc(abund_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(ligand_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuPu"),
          heatmap_legend_param = list(title = "Highest abundant\nligand rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top overall ligands",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))
heatmap_overal_abund_lig_abund


heatmap_overal_abund_lig_abund_rec <- overal_abund_lig_abund_rec_ranks %>% 
  arrange(desc(abund_lig_abund_rec_rank)) %>% 
  head(50) %>% 
  select(receptor_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Abundant receptor\nrank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top overall ligands",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))

heatmap_overal_abund_lig_abund_rec



ht4 = heatmap_overal_abund_lig_abund + heatmap_overal_abund_lig_abund_rec

ht4


png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_scores/top50_overall_abund_lig_abund_rec.png",
    width = 4,
    height = 8,
    units = "in",
    res = 300)

draw(ht4,
     ht_gap = unit(0.3, "cm"))

dev.off()




################### heatmap abund lig scarce rec ####################################


heatmap_abund_lig_scarce_rec_interactions <- interaction_ranks_rowname %>% 
  arrange(desc(abund_lig_scarce_rec_rank)) %>% 
  head(50) %>% 
  select(ligand_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuPu"),
          heatmap_legend_param = list(title = "Abundant ligand rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top interactions",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))
heatmap_abund_lig_scarce_rec_interactions


heatmap_abund_lig_scarce_rec_interactions_rec <- interaction_ranks_rowname %>% 
  arrange(desc(abund_lig_scarce_rec_rank)) %>% 
  head(50) %>% 
  select(receptor_scarce_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Scarce receptor rank"),
          border = TRUE,
          height = unit(18, "cm"),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))


heatmap_abund_lig_scarce_rec_interactions_rec

ht5 = heatmap_abund_lig_scarce_rec_interactions + heatmap_abund_lig_scarce_rec_interactions_rec

ht5

png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_scores/top50_interactions_abund_lig_scarce_rec.png",
    width = 4,
    height = 8,
    units = "in",
    res = 300)

draw(ht5,
     ht_gap = unit(0.3, "cm"))

dev.off()


# overall top abund ligands with scarce receptors

overal_abund_lig_scarce_rec_ranks <- interaction_ranks %>% 
  select(receptor,
         ligand_abundant_rank,
         receptor_scarce_rank,
         abund_lig_scarce_rec_rank) %>% 
  group_by(receptor) %>% 
  summarise_all(max) %>% 
  arrange(desc(abund_lig_scarce_rec_rank))

overal_abund_lig_scarce_rec_ranks <- column_to_rownames(as.data.frame(overal_abund_lig_scarce_rec_ranks), 
                                                       var = "receptor")

heatmap_overal_abund_lig_scarce <- overal_abund_lig_scarce_rec_ranks %>% 
  arrange(desc(abund_lig_scarce_rec_rank)) %>% 
  head(50) %>% 
  select(ligand_abundant_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuPu"),
          heatmap_legend_param = list(title = "Highest abundant\nligand rank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top overall ligands",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))
heatmap_overal_abund_lig_scarce


heatmap_overal_abund_lig_scarce_rec <- overal_abund_lig_scarce_rec_ranks %>% 
  arrange(desc(abund_lig_scarce_rec_rank)) %>% 
  head(50) %>% 
  select(receptor_scarce_rank) %>% 
  data.matrix(rownames.force = TRUE) %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          show_column_names = FALSE,
          col = RColorBrewer::brewer.pal(7, "BuGn"),
          heatmap_legend_param = list(title = "Scarce receptor\nrank"),
          border = TRUE,
          height = unit(18, "cm"),
          # row_title = "Top overall ligands",
          row_title_gp = gpar(fontsize = 20),
          column_title_gp = gpar(font = 2,
                                 fontsize = 13),
          width = unit(0.8, "cm"))

heatmap_overal_abund_lig_scarce_rec



ht6 = heatmap_overal_abund_lig_scarce + heatmap_overal_abund_lig_scarce_rec

ht6


png("ligand_receptor_lists/feb2022_new_lists/OmniPath/figures/interaction_scores/top50_overall_abund_lig_scarce_rec.png",
    width = 4,
    height = 8,
    units = "in",
    res = 300)

draw(ht6,
     ht_gap = unit(0.3, "cm"))

dev.off()





