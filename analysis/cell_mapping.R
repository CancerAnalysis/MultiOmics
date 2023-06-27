########### SingleR cell mapping ##########

library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(RColorBrewer)

hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

ascite.dmp_sce <- as.SingleCellExperiment(ascite.dmp)


pred <- SingleR(sc_data=ascite.dmp_sce, ref=hpca.se, assay.type.test=1,labels = hpca.se$label.main)

clu<-data.frame(ascite.dmp@active.ident)
colnames(clu)<-c("cluster")

pred.d<-data.frame(pred)
pred.d<-cbind(pred.d,clu)

a<-table(pred.d$labels, pred.d$cluster)
b<-apply(a,2,max)

for ( i in 1:ncol(a)){
cell<-names(which(a[,i]==b[i]))
pt<-paste0("C",i-1,": ", cell)
cat(pt,"\n")

}

