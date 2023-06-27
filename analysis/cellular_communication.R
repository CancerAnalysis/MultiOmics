########### cell cell communication ########
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

seurat.data<-ascite.dmp

data.input <- GetAssayData(seurat.data, , slot = "counts") # normalized data matrix
labels <- Idents(seurat.data)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

#showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
### project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)
#cellchat <- projectData(cellchat, PPI.mouse)

### long time consume
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

####signaling role identification
# Compute the network centrality scores
#plan(sequential)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat<-rankNetPairwise(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

total_colors <- brewer.pal(10, "Set3")

par(mfrow = c(1,2), xpd=TRUE)

groupSize <- as.numeric(table(cellchat_sev@idents))
cellchat_sev@net$count<-cellchat_sev@net$count[,]
groupSize <- as.numeric(table(cellchat_sev@idents))
netVisual_circle(cellchat_sev@net$count, 
                  vertex.weight = groupSize, weight.scale = T, 
                  label.edge= F, title.name = "Number of interactions",
                  color.use= total_colors[-9],
                  vertex.label.cex = 1.3 )
