#query and download de TCGA
library(TCGAbiolinks)
library(DT)
library(dplyr)
options(timeout=3600)

#busqueda de datos de tumores brca por subtipo
subtypes<-TCGAquery_subtype(tumor="BRCA")

#masajeo/filtro solo luminales A y basal

lumA_basclin=filter(subtypes,BRCA_Subtype_PAM50 == "LumA" |
                    BRCA_Subtype_PAM50== "Basal")

#download de datos de expresión
query_exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  file.type = "rna_seq.augmented_star_gene_counts.tsv",
  barcode=lumA_basclin$patient)
GDCdownload(query_exp)
exp_data<-GDCprepare(query_exp, summarizedExperiment = FALSE)

#masajeo de tabla datos de expresión
exp_data_shrink=filter(exp_data, gene_type=="protein_coding")
exp_data_shrink<-exp_data_shrink[,c(1:3, 856:1707)]
exp_data_shrink=as.data.frame(exp_data_shrink)
rownames(exp_data_shrink)<-exp_data_shrink$gene_id
exp_data_shrink<-exp_data_shrink[,-1]

#Filtro por pureza tumoral mayor a 60%
CPE_BRCA=read.table("c:/Users/User/Lab/idisba/FITR2/CPE_BCRA.txt",header=TRUE,fill=TRUE)
CPE_filt=filter(CPE_BRCA, CPE >= 0.6)

#Mergeo de  datos filtrados por CPE
lumA_basclin2<-lumA_basclin[,c(1,12)]
CPE_filt$TCGA_superred<-substr(CPE_filt$SampleID_Reduced, 1, 
                               nchar(CPE_filt$SampleID_Reduced) -3)

lumA_basclin2=merge(lumA_basclin2, CPE_filt, by.x=1, by.y=9)

exp_data_shrink2<-as.data.frame(t(exp_data_shrink))

exp_data_shrink2<-exp_data_shrink2[-c(1,2),]
exp_data_shrink2$TCGA_ID_red<-rownames(exp_data_shrink2)
exp_data_shrink2$TCGA_ID_red<-substr(exp_data_shrink2$TCGA_ID_red,16,
                                     nchar(exp_data_shrink2$TCGA_ID_red) -12)
exp_data_shrink2<-merge(exp_data_shrink2, lumA_basclin2, by.x=19963, by.y=3)

#tabla filtrada por CPE para DEA

exp_data_shrink2<-exp_data_shrink2[!duplicated(exp_data_shrink2$TCGA_ID_red),]

exp_data_shrink3<-exp_data_shrink2[,-c(19964:19972)]


rownames(exp_data_shrink3)<-exp_data_shrink3$TCGA_ID_red
exp_data_shrink3<-exp_data_shrink3[,-c(1)]


exp_data_shrink3=as.matrix(t(exp_data_shrink3))
#algunos datos son character, los paso a numéricos con la siguiente funcion
exp_data_shrink3 = data.frame(apply(exp_data_shrink3, 2, 
                   function(x) as.numeric(as.character(x))), 
                  row.names = rownames(exp_data_shrink3))

#DEA Analysis
library(DESeq2)

DEA_lumA_Bas <- DESeqDataSetFromMatrix(exp_data_shrink3, 
                                       colData=exp_data_shrink2, 
                                       design=~ BRCA_Subtype_PAM50)
exp_pre<-DESeq(DEA_lumA_Bas)
res_DEA <- results(exp_pre)
res_DEA2=as.data.frame(res_DEA)

genenames=exp_data_shrink[,-c(3:854)]
#macheo con Gene Symbol
res_DEA2=merge(res_DEA2,genenames, by='row.names')
#Filtro DEA por log2FoldChange y pvalue
res_DEA2=filter(res_DEA2, pvalue <=0.05 & (log2FoldChange >=1| log2FoldChange <= -1))

#Volcano plot (LumA vs Basal)
library(EnhancedVolcano)

EnhancedVolcano(res_DEA2,
                lab = res_DEA2$gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title="Volcano Plot",
                subtitle="Expresión diferencial Luminal A vs. Basal")





#con esta expresión se puede elegir que condición versus que condición hacemos.
res4<-results(exp_pre,contrast=c("BRCA_Subtype_PAM50","Basal","LumA"))
res5=as.data.frame(res4)



#arreglo de tablas para machine learning
machlear_data=merge( res_DEA2, exp_data_shrink3 , by.x=1, by.y='row.names')
machlear_data<-machlear_data[!duplicated(machlear_data$gene_name),]
rownames(machlear_data)<-machlear_data$gene_name
machlear_data=machlear_data[,-c(1:9)]
machlear_data= t(machlear_data)

#función para cmabiar de simbolo en rownames
rownames(machlear_data) = gsub(".", "-", rownames(machlear_data), fixed = TRUE)

lumA_basclin2<-lumA_basclin2[!duplicated(lumA_basclin2$patient),]
rownames(lumA_basclin2)<-lumA_basclin2$Sample_ID
machlear_data=merge(machlear_data, lumA_basclin2, by='row.names')
machlear_data=machlear_data[,-c(4991, 4993:5000)]



library(splitTools)

set.seed(123)

train_val_split <- partition(machlear_data$BRCA_Subtype_PAM50, p=c(train=0.7,valid=0.3), 
                             type="stratified",use_names = TRUE,
                             split_into_list = TRUE,shuffle=TRUE)
#con la función partition() obtengo 2 listas de las muestras 
#randomizadas por su gene mutation

#integro las listas de partition a maclear meta
train <- machlear_data[train_val_split$train, ]
valid <- machlear_data[train_val_split$valid, ]

#compruebo la proporción de casos
table(train$BRCA_Subtype_PAM50)
table(valid$BRCA_Subtype_PAM50)

train_dea=train[,-c(4991)]
rownames(train_dea)<-train_dea$Row.names
train_dea<- train_dea[,-c(1)]
train_dea<-as.matrix(t(train_dea))


#DEA set de training
library(DESeq2)

#DEA
DEA2 <- DESeqDataSetFromMatrix(train_dea, 
                               colData=train, 
                               design=~ BRCA_Subtype_PAM50, 
                               ignoreRank = TRUE)
exp_train2<-DESeq(DEA2)
train3e <- results(exp_train2)
train3e=as.data.frame(train3e)

#Filtro DEA por log2FoldChange y pvalue
train3e=filter(train3e, pvalue <=0.05 & (log2FoldChange >=1| log2FoldChange <= -1))

train_dea=merge(train_dea, train3e, by='row.names')
train_dea<-train_dea[,-c(458:463)]
rownames(train_dea)=train_dea$Row.names
train_dea<-train_dea[,-c(1)]
train_dea=t(train_dea)
train_dea=train_dea[-c(447:452),]
rownames(train)=train$Row.names
train<-train[,-c(1)]
trainrf<-train[,-c(1:4989)]
trainrf=as.factor(trainrf)

#Seleccion de variables para random forrest
library(varSelRF)

ranfor_train= varSelRF(train_dea, trainrf, vars.drop.frac = 0.3)
#Esto me da la lista de genes que constituyen una posible firma basado en OOB score
plot(ranfor_train)


#extraigo la lista de genes
sign=as.data.frame(ranfor_train$selected.vars)
colnames(sign)<-c("Predgenes")
rownames(sign)<-sign$Predgenes
train_meta=as.data.frame(train[,4990])
colnames(train_meta)<-"PAM50"
train2=t(train)
train2 = merge(train2,sign, by='row.names')

rownames(train2)=train2$Row.names
train2<-train2[,-c(1)]
train2<-t(train2)
train3<-merge(train2, train_meta, by='row.names')
rownames(train3)<-train3$Row.names
train3<-train3[,-c(1)]

ranfor2=randomForest(factor(PAM50) ~ ., data = train3, importance = TRUE, type = "classification")

votes=as.data.frame(ranfor2$votes)

library(pROC)
roc_train= roc(train3$PAM50, votes$LumA, direction="<")
plot(smooth(roc_train),asp=NA,col="red", print.auc=TRUE,auc.polygon=TRUE, legacy.axes=TRUE )

library(randomForest)
rownames(valid)=valid$Row.names
valid<-valid[,-c(1)]

p<-predict(ranfor2, valid, type="prob")

predvalid=as.data.frame(p)

roc_valid=roc(valid$BRCA_Subtype_PAM50, predvalid$LumA, direction= "<")

plot(roc_valid,asp=NA,col="red", print.auc=TRUE,auc.polygon=TRUE, legacy.axes=TRUE )

#Como el clasificador tiene 115 genes pruebo con una cantidad menor de los genes
#obtenidos en VarSelRF

allclasif=as.data.frame(ranfor_train$selec.history)


#prueba de signature de 5 genes
sign2=All_variables[30,]
sign2=t(sign2)
colnames(sign2)<-c("Predgenes")
sign2=as.data.frame(sign2)
sign2=sign2[!is.na(sign2$Predgenes),]
sign2=as.data.frame(sign2)
rownames(sign2)<-sign2$sign2
colnames(sign2)<-c("Predgenes")
train4=t(train)
train4 = merge(train4,sign2, by='row.names')

rownames(train4)=train4$Row.names
train4<-train4[,-c(1)]
train4<-t(train4)
train5<-merge(train4, train_meta, by='row.names')
rownames(train5)<-train5$Row.names
train5<-train5[,-c(1)]


ranfor3=randomForest(factor(PAM50) ~ ., data = train5, importance = TRUE, type = "classification")

votes=as.data.frame(ranfor3$votes)

library(pROC)
roc_train= roc(train5$PAM50, votes$LumA, direction="<")
plot(smooth(roc_train),asp=NA,col="red", print.auc=TRUE,auc.polygon=TRUE, legacy.axes=TRUE )

library(randomForest)


p2<-predict(ranfor3, valid, type="prob")

predvalid2=as.data.frame(p2)

roc_valid2=roc(valid$BRCA_Subtype_PAM50, predvalid$LumA, direction= "<")

plot(roc_valid2,asp=NA,col="red", print.auc=TRUE,auc.polygon=TRUE, legacy.axes=TRUE )







    
    
    
  