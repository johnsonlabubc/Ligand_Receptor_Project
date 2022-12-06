# compare human islet vs sorted SC islet scRNAseq using cellchat
# start by loading the RDS files of the cellchat objects for each 
# they were both created using projectData onto the protein interaction network



# load datasets
cellchat_SC_islet <- readRDS(file = "single_cell_analysis/CellChat/data/cellchat_sorted_SCs.rds")
cellchat_H_islet <- readRDS(file = "single_cell_analysis/CellChat/data/cellchat_h_islet_projectData.rds")

########## compare SCsorted to human islet ######################

# join the datasets into 1 cellchat object
object.list <- list(SC_islet = cellchat_SC_islet, human_islet = cellchat_H_islet)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

compareInteractions(cellchat, show.legend = F, group = c(1,2))

computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

rankSimilarity(cellchat, type = "structural")



rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

