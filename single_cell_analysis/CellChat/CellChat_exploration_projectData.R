# After going through CellChat_exploration.R with the Lynn lab scRNA-seq, had issue
# where seems only about 400 interactions were used from the total 1,939 interactions 
# in CellChatDB$interaction
# 
# May have to do with sequencing depth, so going to rerun the analysis here, this time 
# including optional step to project gene expression data onto PPI network
# using function projectData()
# 
# This function is useful when analyzing single-cell data with shallow sequencing depth 
# because the projection reduces the dropout effects of signaling genes, in particular 
# for possible zero expression of subunits of ligands/receptors


# Initial CellChat Exploration

# install a bunch of packages and their dependencies

# # install bioconductor
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# 
# # Biobase is needed for NMF
# BiocManager::install("Biobase")
# 
# # NFM is needed for CellChat
# install.packages('NMF')
# 
# 
# devtools::install_github("jokergoo/circlize")
# 
# devtools::install_github("jokergoo/ComplexHeatmap")
# 
# 
# # now finally install CellChat itself
# devtools::install_github("sqjin/CellChat")

# now everything is installed

library(CellChat)
library(Seurat)
library(tidyverse)



######## Create a CellChat object ###################



# open the 4.6 GB RDS file containing the already processed data
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


rm(seurat_obj)
gc()

# confirm the default assay is integrated
DefaultAssay(seurat_obj_updated)

# set to RNA
DefaultAssay(seurat_obj_updated) <- "RNA"

# confirm default assay was changed to RNA
DefaultAssay(seurat_obj_updated)

# don't specify meta because it takes the meta automatically from seurat object
cellchat <- createCellChat(object = seurat_obj_updated,
                           group.by = "ident")



######## Set the ligand-receptor interaction database ###################

CellChatDB <- CellChatDB.human 

# CellChatDB in human contains 1,939 validated molecular interactions
# Here is the breakdown of them:
showDatabaseCategory(CellChatDB)

# Show the structure of the interaction database
View(CellChatDB$interaction)

# save the full human database for later reference
write_tsv(CellChatDB$interaction, "single_cell_analysis/data/cellchat_human_interaction_db.tsv")

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use


################## Preprocessing the expression data ##################

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

# open up memory space for next functions
rm(seurat_obj)
rm(seurat_obj_updated)
gc()

# Identify over-expressed signaling genes associated with each cell group
cellchat <- identifyOverExpressedGenes(cellchat)
# Identify over-expressed ligand-receptor interactions
cellchat <- identifyOverExpressedInteractions(cellchat)


# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


############## Inference of cell-cell communication network ################


# determine whether should set population.size = TRUE
# Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, 
# with the reason that abundant cell populations tend to send collectively 
# stronger signals than the rare cell populations.

# check number of each cell type in dataseet
cell_type_counts <- cellchat@meta %>% 
  as.tibble() %>% 
  group_by(ident) %>% 
  summarise(n = n())

# save cell type counts table
write_tsv(cell_type_counts, "single_cell_analysis/data/lynn_human_scRNAseq_celltype_counts.tsv")

# we should prob use population.size = TRUE
# raw.use = FALSE to to use the projected data when analyzing single-cell data with
# shallow sequencing depth because the projected data could help to reduce the
# dropout effects of signaling genes, in particular for possible zero expression 
# of subunits of ligands/receptors.
cellchat <- computeCommunProb(cellchat,
                              raw.use = FALSE,
                              population.size = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


#### Extract the inferred cellular communication network as a data frame

# returns a data frame consisting of all the inferred cell-cell 
# communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)

# save datatable will all direct cell-cell interactions list
write_tsv(df.net, "single_cell_analysis/data/inferred_direct_cell_cell_comms.tsv")


# gives the inferred communications from and to specific cell groups
df.net <- subsetCommunication(cellchat, 
                              sources.use = c(4), 
                              targets.use = c(10))


# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")




# examine the signaling sent from each cell group
mat <- cellchat@net$weight
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
mat2["Beta", ] <- mat["Beta", ]
netVisual_circle(mat2,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 edge.weight.max = max(mat),
                 title.name = "Beta")

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



View(cellchat@LR$LRsig)

pathways.show <- c("WNT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, 
                    signaling = pathways.show,  
                    vertex.receiver = vertex.receiver)



# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of 
# signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = pathways.show, 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)

# reset plot window to fix plotting glitch
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# outgoing signals
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",
                                         font.size = 6)
# incoming signals
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                         font.size = 6)
ht1 + ht2

# both incoming and outgoing signals
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all",
                                  font.size = 6)


####### identify main pathways in alpha/beta signalling to pericytes ##########

#Rank ligand-receptor interactions for any pair of two cell groups
cellchat <- rankNetPairwise(cellchat)

# Identify all the significant interactions (L-R pairs) from some cell groups 
# to other cell groups
identifyEnrichedInteractions(cellchat,
                             from = "Alpha",
                             to = "Pericytes",
                             bidirection = FALSE,
                             pair.only = FALSE,
                             thresh = 0.05)


# attempt to fix "Error: not enough space for cells at track index '1'."

pdf(file ="single_cell_analysis/figures/alpha_pericyte_pathways.pdf",
    width = 10, 
    height = 10)

# show all the interactions received by Pericytes
netVisual_chord_gene(cellchat, 
                     sources.use = "Alpha", 
                     targets.use = "Pericytes", 
                     small.gap = 1,
                     big.gap = 10,
                     lab.cex = 1.3,
                     legend.pos.x = 25,
                     legend.pos.y = 45)

dev.off()



# save cellchat data file in case anything happens
saveRDS(cellchat, file = "single_cell_analysis/CellChat/data/cellchat_h_islet_projectData.rds")



#reopen data
cellchat <- readRDS(file = "single_cell_analysis/CellChat/data/cellchat_h_islet_projectData.rds")




############ explore beta-cell - duct cell signalling ########################


cellchat <- rankNetPairwise(cellchat)



# Identify all the significant interactions (L-R pairs) from some cell groups 
# to other cell groups
identifyEnrichedInteractions(cellchat,
                             from = "Beta",
                             to = "Duct",
                             bidirection = TRUE,
                             pair.only = FALSE,
                             thresh = 0.05)

# 
# png(file = "single_cell_analysis/CellChat/figures/human_islet/figures/beta_duct_pathways.png",
#     width = 5, 
#     height = 5,
#     units = "in",
#     res = 300)

# 
# netVisual_chord_gene(cellchat, 
#                      sources.use = "Beta", 
#                      targets.use = "Duct",
#                      slot.name = "netP",
#                      show.legend = FALSE)
#                      # small.gap = .1,
#                      # big.gap = 1,
#                      # lab.cex = 1.3,
#                      # legend.pos.x = 25,
#                      # legend.pos.y = 45)
# 
# 
# 
# dev.off()
# 
# CellChat::netAnalysis_dot(cellchat)
# 
# netVisual_chord_gene(cellchat, 
#                      sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

png(file = "single_cell_analysis/CellChat/figures/human_islet/ANGPT_pathway_signaling.jpg",
    width = 5,
    height = 3,
    units = "in",
    res = 300)


# show pathway by pathway the top ones between duct and beta
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = "ANGPT", 
                                  width = 8, 
                                  height = 2.5, 
                                  font.size = 10)


dev.off()


