# CellChat exploration of sorted stem cell dataset from lynn lab



library(CellChat)
library(Seurat)
library(tidyverse)



######## Create a CellChat object ###################

filename <- file.choose()
seurat_obj <- readRDS(filename)

# got an error due to the seurat object being an old version, so updated it
seurat_obj_updated <- UpdateSeuratObject(seurat_obj)

# Got this error, so need to extract RNA assay instead of integrated
# Error in createCellChat(object = seurat_obj_updated, meta = meta, group.by = "labels") : 
# The data matrix contains negative values. Please ensure the normalized data matrix is used.
# In addition: Warning message:
# In createCellChat(object = seurat_obj_updated, meta = meta, group.by = "labels") :
# The data in the `integrated` assay is not suitable for CellChat analysis! Please use the `RNA` or `SCT` assay! 


# confirm the default assay is integrated
DefaultAssay(seurat_obj_updated)

# set to RNA
DefaultAssay(seurat_obj_updated) <- "RNA"

# confirm default assay was changed to RNA
DefaultAssay(seurat_obj_updated)

# don't specify meta because it takes the meta automatically from seurat object
cellchat_SC_islet <- createCellChat(object = seurat_obj_updated,
                           group.by = "clusters")




######## Set the ligand-receptor interaction database ###################

CellChatDB <- CellChatDB.human 

# CellChatDB in human contains 1,939 validated molecular interactions
# Here is the breakdown of them:
showDatabaseCategory(CellChatDB)

# Show the structure of the interaction database
View(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# set the used database in the object
cellchat_SC_islet@DB <- CellChatDB.use


################## Preprocessing the expression data ##################

# subset the expression data of signaling genes for saving computation cost
cellchat_SC_islet <- subsetData(cellchat_SC_islet) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

# open up memory space for next functions
rm(seurat_obj)
rm(seurat_obj_updated)
gc()

# Identify over-expressed signaling genes associated with each cell group
cellchat_SC_islet <- identifyOverExpressedGenes(cellchat_SC_islet)
# Identify over-expressed ligand-receptor interactions
cellchat_SC_islet <- identifyOverExpressedInteractions(cellchat_SC_islet)

### optional step
# project gene expression data onto PPI network (optional)
cellchat_SC_islet <- projectData(cellchat_SC_islet, PPI.human)




############## Inference of cell-cell communication network ################


# determine whether should set population.size = TRUE
# Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, 
# with the reason that abundant cell populations tend to send collectively 
# stronger signals than the rare cell populations.

# check number of each cell type in dataseet
cell_type_counts <- cellchat_SC_islet@meta %>% 
  as.tibble() %>% 
  group_by(clusters) %>% 
  summarise(n = n())

# save cell type counts table
write_tsv(cell_type_counts, "single_cell_analysis/stem_cell_scrnaseq_francis/lynn_sorted_SCs_celltype_counts.tsv")


# we should prob use population.size = TRUE

cellchat_SC_islet <- computeCommunProb(cellchat_SC_islet,
                              population.size = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_SC_islet <- filterCommunication(cellchat_SC_islet, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat_SC_islet <- computeCommunProbPathway(cellchat_SC_islet)

# Calculate the aggregated cell-cell communication network
cellchat_SC_islet <- aggregateNet(cellchat_SC_islet)

# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat_SC_islet@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat_SC_islet@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat_SC_islet@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")



# examine the signaling sent from each cell group
mat <- cellchat_SC_islet@net$weight
par(mfrow = c(3,4), xpd=TRUE)

# reset plot window to fix plotting glitch
dev.off()

# check that the rownames are the cell types
row.names(mat)

# comment out function that plots for all cell types
# cause can just do it for a couple cell types of interest manually
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0,
#                  nrow = nrow(mat),
#                  ncol = ncol(mat),
#                  dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2,
#                    vertex.weight = groupSize,
#                    weight.scale = T,
#                    edge.weight.max = max(mat),
#                    title.name = rownames(mat)[i])
#   }

# plot signals sent by beta cells
mat2 <- matrix(0,
               nrow = nrow(mat),
               ncol = ncol(mat),
               dimnames = dimnames(mat))
mat2["eBCs", ] <- mat["eBCs", ]
netVisual_circle(mat2,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 edge.weight.max = max(mat),
                 title.name = "eBCs")

# plot signals sent by alpha cells
mat2 <- matrix(0,
               nrow = nrow(mat),
               ncol = ncol(mat),
               dimnames = dimnames(mat))
mat2["Alpha", ] <- mat["Alpha", ]
netVisual_circle(mat2,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 edge.weight.max = max(mat),
                 title.name = "Alpha")



View(cellchat_SC_islet@LR$LRsig)

pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_SC_islet, 
                    signaling = pathways.show,  
                    vertex.receiver = vertex.receiver)



# Compute the network centrality scores
cellchat_SC_islet <- netAnalysis_computeCentrality(cellchat_SC_islet, 
                                          slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of 
# signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat_SC_islet, 
                                  signaling = pathways.show, 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)


# reset plot window to fix plotting glitch
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_SC_islet, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_SC_islet, pattern = "incoming")
ht1 + ht2

# combined both incoming and outgoing
netAnalysis_signalingRole_heatmap(cellchat_SC_islet, pattern = "all")


# save cellchat data file
saveRDS(cellchat_SC_islet, file = "single_cell_analysis/CellChat/data/cellchat_sorted_SCs.rds")


