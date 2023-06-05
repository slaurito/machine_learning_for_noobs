setwd("c:/Users/User/Lab/idisba/Lia")
library(readr)


#Pruebas de distintas tablas, Galaxy.BED la generé con la función Intersect de
#Galaxy entre EPIC_annotation y la tabla de coordenadas de genes Lia.
test2=read.table("EPIC_annotation_hg38.txt", header=TRUE, fill=TRUE)
test3=read.table("Galaxy.bed")
test4=read.table("data.txt",header=TRUE)
#limpieza de datos

clean_t3=test3[!duplicated(test3),]


#Masajeo de data de Beta values Lia y Lista de sondas GSEA

data_bvalue= test4[,1:15]
gsea_genes_prom= merge(data_bvalue, clean_t3, by.x=1, by.y=4)
gsea_genes_prom=gsea_genes_prom[!duplicated(gsea_genes_prom$TargetID),]
gsea_genes_prom=gsea_genes_prom[,1:15]
gsea_genes_prom_pca <- gsea_genes_prom
rownames(gsea_genes_prom_pca) = gsea_genes_prom_pca$TargetID
gsea_genes_prom_pca=gsea_genes_prom_pca[,-1]
gsea_genes_prom_pca = na.omit(gsea_genes_prom_pca)

gsea_genes_prom_pca2=data.frame(t(gsea_genes_prom_pca))

library(ChAMP)
rownames(data_bvalue)=data_bvalue$TargetID
data_bvalue<-data_bvalue[,-1]
batch_eff=as.data.frame(t(data_bvalue))



myNorm=batch_eff
myLoad=read.table("Sample.Table.txt",header=TRUE,fill=TRUE)
myLoad<-myLoad[,-c(1,3,10,11,12,13)]
colnames(myLoad)<-c("Sample_Name","Array","Slide", "Sample_Group1", "Sample_Group2","Sample_Group3","BatchArray")
myNormt = as.data.frame(t(myNorm))
myNormt = na.omit(myNormt)
colnames(myLoad)<-c("Sample_Name", "Array", "Slide","Mutation(8 or 13)","Mutated_or_Wild_type", "Heteroplasmy", "Batch_Array")

champ.SVD(beta=myNormt,pd=myLoad,RGEffect = TRUE) 


myCombat <- champ.runCombat(beta=myNormt,pd=myLoad,batchname="Batch_Array", variablename = "Heteroplasmy", logitTrans=FALSE)
batch_corr=as.data.frame(t(myCombat))
batch_corr2=as.data.frame(myCombat)
champ.SVD(beta=batch_corr2,pd=myLoad,RGEffect = TRUE)
library(ggplot2)
library(ggfortify)
myNorm2=as.data.frame(t(myNormt))
pca_before=prcomp(myNorm2, scale=TRUE)

myNorm3=myNorm2
myNorm3['Array']<-NA
myNorm3$Array=c("M","M","M","M","M","M","K","K","K","K","K","K","K","K")

autoplot(pca_before, data=myNorm3, colour='Array', label=TRUE, label.size=2, size=3)

myCombat2=as.data.frame(t(myCombat))
myCombat3=myCombat2
myCombat3['Array']<-NA
myCombat3$Array=c("M","M","M","M","M","M","K","K","K","K","K","K","K","K")
pca_after=prcomp(myCombat2, scale=TRUE)
autoplot(pca_after, data=myCombat3, colour='Array', label=TRUE, label.size=2, size=3)
