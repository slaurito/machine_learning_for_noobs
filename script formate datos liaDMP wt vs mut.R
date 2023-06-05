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

library(ggplot2)
library(ggfortify)

#PCA de las 850K de EPIC
raw_betas=data_bvalue
rownames(raw_betas) = raw_betas$TargetID
raw_betas=raw_betas[,-1]
raw_betas = na.omit(raw_betas)
raw_betas=as.data.frame(t(raw_betas))

pca_raw_betas=prcomp(raw_betas, scale=TRUE)
autoplot(pca_raw_betas, data=raw_betas,label=TRUE,colour='rownames',label.size=2,  size=3)


#PCA de las 348K sondas relacionadas con los pathways de interés(entodas las muestras)
pca_full_gsea=prcomp(gsea_genes_prom_pca2, scale=TRUE)

autoplot(pca_full_gsea,data=gsea_genes_prom_pca2, colour='rownames',label=TRUE,label.size=2, size=3)

#saco las WT para evitar el empuje
gsea_genes_prom_nowt=gsea_genes_prom_pca2[-c(1,2),]
pca_gsea_nowt=prcomp(gsea_genes_prom_nowt, scale=TRUE)
autoplot(pca_gsea_nowt,colour='rownames',label.size=2, size=3, data=gsea_genes_prom_nowt,label=TRUE)




#Creación de tablas separadas para compración por wilcoxon

dataWT = gsea_genes_prom_pca2[1:2,]
dataMut = gsea_genes_prom_pca2[-c(1:2),]  



meanWT = apply(dataWT, 2, mean)
meanMut = apply(dataMut, 2, mean)

fold = (meanMut-meanWT)

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : ncol(dataMut)) {
  x = dataMut[,i]
  y = dataWT[,i]
  
  t = wilcox.test(as.numeric(x), as.numeric(y), exact=FALSE) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}



combined  = cbind(pvalue, fold, abs(fold), data)
combined_t = as.data.frame(as.matrix(t(combined)))

# Filter per Z-Ratio and qvalue #
fold_cutoff = 0.1
pvalue_cutoff = 0.05

# fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
#dim(data[filter_by_fold, ]) 

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
#dim(data[filter_by_qvalue, ]) 

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue 

filtered = gsea_genes_prom_pca2[,filter_combined]
dim(filtered)
probes = as.data.frame(cbind("Probe" = names(fold),
                             "Fold" = fold,
                             "Abs fold" = abs(fold),
                             "P-value" = pvalue))
probes2 = probes[probes$`Abs fold` > fold_cutoff & probes$`P-value` < pvalue_cutoff,]

probes_with_changes = probes2


UpProbes = probes_with_changes[probes_with_changes$Fold >0, ]
DownProbes = probes_with_changes[probes_with_changes$Fold <0, ]

genes_lia= merge(probes2, clean_t3, by.x=1, by.y=4)

genes_lia2=genes_lia[!duplicated(genes_lia$Probe),]

symbols_lia=read.table("EPIC.hg38.manifest.gencode.v22.tsv",header=TRUE, fill=TRUE)
symbols_lia=symbols_lia[c(5,6)]
gsea_genes_wtvsmut=merge(symbols_lia, genes_lia2, by.x=1,by.y=1)
colnames(gsea_genes_wtvsmut)<-c("ProbeID", "GeneName", "Fold", "AbsFold", "P-value", "chr", "start", "end")
gsea_genes_wtvsmut= na.omit(gsea_genes_wtvsmut)

write_tsv(gsea_genes_wtvsmut, "DMP_genes_WT_vs_Mut.tsv", col_names=TRUE)
lia_probes=gsea_genes_wtvsmut[,1]

lia_probes=as.data.frame(lia_probes)
mut_vs_wt_probes_betas= merge(lia_probes, gsea_genes_prom, by.x=1, by.y=1)

row.names(mut_vs_wt_probes_betas) <- mut_vs_wt_probes_betas[,1]
mut_vs_wt_probes_betas <- mut_vs_wt_probes_betas[,-1]

mut_vs_wt_probes_betas_pca=data.frame(t(mut_vs_wt_probes_betas))


pca_wt_mut= prcomp(mut_vs_wt_probes_betas_pca, scale=TRUE)
pca_wt_mut2 = pca_wt_mut
pca_wt_mut2$x[11,1] = pca_wt_mut2$x[11,1] + 3
autoplot(pca_wt_mut2,colour='rownames',label.size=2, size=3,label=TRUE)

library(gplots)


#Heatmap preparation
mut_vs_wt_probes_betas_heatm<-as.data.frame(mut_vs_wt_probes_betas_pca)
mut_vs_wt_probes_betas_heatm['Sample_Type']<-NA
mut_vs_wt_probes_betas_heatm$Sample_Type= c("WildType", "WildType", "Mut13Low1", "Mut13Low2", "Mut13High1","Mut13High2","Mut13Low3", "Mut13High4","Mut8Low1","Mut8Low2","Mut8Low7", "Mut8High1", "Mut8High2", "Mut8High3")
mut_vs_wt_probes_betas_pca<- mut_vs_wt_probes_betas_pca[,-15267]
colv = as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(mut_vs_wt_probes_betas_pca))))))
rowv = as.dendrogram(hclust(as.dist(1-cor(as.matrix(mut_vs_wt_probes_betas_pca)))))  

heatmap.2(as.matrix(t(mut_vs_wt_probes_betas_pca)), Colv=colv,Rowv=rowv,key=TRUE,labCol=mut_vs_wt_probes_betas_heatm$Sample_Type, cexCol=0.7,col = rev(redblue(256)), scale = "row", trace=NULL, tracecol=NULL)
legend('bottomleft',legend = unique(mut_vs_wt_probes_betas_heatm$Sample_Type), border=T)
dev.off()

