
#########Data integration##########

##standard QC
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
ascite[["percent.mt"]] <- PercentageFeatureSet(ascite, pattern = "^MT-")

#feature filtering
ascite <- subset(ascite, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization
ascite <- NormalizeData(ascite, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features
ascite <- FindVariableFeatures(ascite, selection.method = "vst", nfeatures = 2000)

### scailing (Dmension reduction)
all.genes <- rownames(ascite)
ascite <- ScaleData(ascite, features = all.genes)

#PCA
ascite <- RunPCA(ascite, features = VariableFeatures(object = ascite))

### Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# !!!!!!!!!! Time consumpion process !!!!!!!!!!
# computation time: 7 min !!!
ascite <- JackStraw(ascite, num.replicate = 100)
ascite <- ScoreJackStraw(ascite, dims = 1:20)

#replication time and pvalue
#JackStrawPlot(ascite, dims = 1:15)
#sd 
#ElbowPlot(ascite)

###Clustering
ascite <- FindNeighbors(ascite, dims = 1:10)
ascite <- FindClusters(ascite, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(ascite), 5)

# UMAP / tSNE
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
ascite <- RunUMAP(ascite, dims = 1:10)
ascite <- RunTSNE(ascite, dims = 1:10)

#######################demultiplexing#######################

############## After RunUMAP ##############
#####multiplet elimination#####
##ref:https://github.com/chris-mcginnis-ucsf/DoubletFinder
ascite<-readRDS("ascite.rds")

library(DoubletFinder)

## pK Identification (no ground-truth) 
sweep.res.list_ascite <- paramSweep_v3(ascite, PCs = 1:10, sct = FALSE)
sweep.stats_ascite <- summarizeSweep(sweep.res.list_ascite, GT = FALSE)
bcmvn_ascite <- find.pK(sweep.stats_ascite)

## Homotypic Doublet Proportion Estimate 
annotations <- Idents(ascite)
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.075*length(names(Idents(ascite))))  

## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
ascite <- doubletFinder_v3(ascite, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#ascite <- doubletFinder_v3(ascite, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_3803", sct = FALSE)

##3803, 3324 means the number of doublet 
summary(as.factor(ascite@meta.data[,8]))

## cell culturing ##
elim_cell<-rownames(ascite@meta.data[ascite@meta.data[,8]=="Doublet",])
length(elim_cell)

filt.cell<-rownames(ascite@meta.data[ascite@meta.data[,8]=="Singlet",])
length(filt.cell)

ascite.dmp<-subset(ascite, cells=filt.cell)
