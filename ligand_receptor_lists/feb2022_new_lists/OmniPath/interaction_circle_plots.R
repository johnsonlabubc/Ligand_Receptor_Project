# create circle chord plots to visualize ligand receptor interactions


library(circlize)
library(tidyverse)


# load interaction data

interactions_df <- (read.csv("ligand_receptor_lists/feb2022_new_lists/OmniPath/data/interaction_ranks.tsv", 
                             sep = "\t"))



# Transform input data in a adjacency matrix
adjacencyData <- with(data, table(interactions_df$ligand, interactions_df$receptor))

# Make the circular plot
chordDiagram(adjacencyData, transparency = 0.5)

circos.link(interactions_df$ligand,
            interactions_df$receptor)


# Create an adjacency matrix: 
# a list of connections between 20 origin nodes, and 5 destination nodes:
numbers <- sample(c(1:1000), 100, replace = T)
data <- matrix( numbers, ncol=5)
rownames(data) <- paste0("orig-", seq(1,20))
colnames(data) <- paste0("dest-", seq(1,5))

# Make the circular plot
chordDiagram(data, transparency = 0.5)


# Create an edge list: a list of connections between 10 origin nodes, and 10 destination nodes:
origin <- paste0("orig ", sample(c(1:10), 20, replace = T))
destination <- paste0("dest ", sample(c(1:10), 20, replace = T))
data <- data.frame(origin, destination)

interactions_df_test <- interactions_df %>%
  select(ligand,
         receptor) %>% 
  head(100)

# Transform input data in a adjacency matrix
adjacencyData <- with(interactions_df_test, table(ligand, receptor))


chordDiagram(interactions_df_test, 
             transparency = 0.5)


# create a dataframe with connection between leaves (individuals)
all_leaves <- paste("subgroup", seq(1,100), sep="_")
connect <- rbind( 
  data.frame( from=sample(all_leaves, 100, replace=T) , to=sample(all_leaves, 100, replace=T)), 
  data.frame( from=sample(head(all_leaves), 30, replace=T) , to=sample( tail(all_leaves), 30, replace=T)), 
  data.frame( from=sample(all_leaves[25:30], 30, replace=T) , to=sample( all_leaves[55:60], 30, replace=T)), 
  data.frame( from=sample(all_leaves[75:80], 30, replace=T) , to=sample( all_leaves[55:60], 30, replace=T)) 
)


library(ggraph)
