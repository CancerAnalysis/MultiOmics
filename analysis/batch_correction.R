#### batch correction ###
library(rliger)

#ascite
fib.ascite.dmp<-readRDS("fib.pbmc.dmp.ascite.rds")
#primary
fib.deep.dmp<-readRDS("/data2/MetaCancer/scSeq/GSE167297/deep/fib.primary_dmp.rds")


ifnb_liger <- createLiger(list(ctrl = fib.deep.dmp, stim = fib.ascite.dmp))

ifnb_liger <- normalize(fib.combined)
ifnb_liger <- selectGenes(ifnb_liger)
ifnb_liger <- scaleNotCenter(ifnb_liger)


library(Seurat)
library(SeuratData)
library(patchwork)
# install dataset
InstallData("ifnb")
# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

