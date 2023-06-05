#Partición de expression data en cohorte de entrenamiento y cohorte
#de validación (70%-30%)


#
maclear=read.table("c:/Users/User/Lab/idisba/FITR2/maclear/exp_data_maclear.tsv")

#uso Janitor para ajustar la primera fila como colnames
library(janitor)
maclear<- maclear %>% row_to_names(row_number = 1)
maclear_meta=subset(maclear, select=c(1,2,3))
rownames(maclear)<-maclear$Expression_Barcode

#NO CORRER DESDE AQUI!!!
#Como los números son muy grandes en la tabla RNAseq aplico log2
maclear2=maclear[,-c(1,2,3)]
#FUNCION PARA CONVERTIR CHARACTER EN NUMERIC!!!!!! MUY IMPORTANTE!!!!!
maclear2=data.frame(apply(maclear2, 2, function(x) as.numeric(as.character(x))), row.names = rownames(maclear2))
maclear2=log2(maclear2+1)
maclear3=subset(maclear, select=c(1,3))
maclear=merge(maclear2, maclear3, by='row.names')
rownames(maclear)<-maclear$Row.names
maclear<-maclear[,-1]                             
#HASTA AQUI!!!!(esto es para achicar numeros con log 2, pero tener en cuenta que DESEQ2 funciona con integers 
#no normalizadas)

#esto es la library spliTools

#set.seed() es para que siempre largue el mismo randomized set
#ej set.seed(123), siempre que escriba ese codigo antes de generar un set random
#se generará el mismo set cada vez

library(splitTools)

set.seed(123)

train_val_split <- partition(maclear_meta$Gene, p=c(train=0.7,valid=0.3), 
                             type="stratified",use_names = TRUE,
                             split_into_list = TRUE,shuffle=TRUE)
#con la función partition() obtengo 2 listas de las muestras 
#randomizadas por su gene mutation

#integro las listas de partition a maclear meta
train <- maclear_meta[train_val_split$train, ]
valid <- maclear_meta[train_val_split$valid, ]

#compruebo la proporción de casos
table(train$Gene)
table(valid$Gene)

#fusiono maclear_meta con maclear (hice la partition sobre un df mas corto porque hacia overflow)
maclear_train=merge(maclear,train, by.x=60661, by.y=1)
maclear_valid=merge(maclear,valid, by.x=60661, by.y=1)

#DEA-Masajeo
maclear_train_DEA=maclear_train
rownames(maclear_train_DEA)<-maclear_train_DEA[,1]
maclear_train_DEA<-maclear_train_DEA[,-c(1,60662,60663,60664)]
rownames(maclear_train)<-maclear_train[,1]
colnames(maclear_train)[60662] <- "Gene"
maclear_train_DEA<-as.data.frame.matrix(t(maclear_train_DEA))
maclear_train_DEA<-data.matrix(maclear_train_DEA)

maclear_valid_DEA=maclear_valid
rownames(maclear_valid_DEA)<-maclear_valid_DEA[,1]
maclear_valid_DEA<-maclear_valid_DEA[,-c(1,60662,60663,60664)]
rownames(maclear_valid)<-maclear_valid[,1]
colnames(maclear_valid)[60662] <- "Gene"
maclear_valid_DEA<-as.data.frame.matrix(t(maclear_valid_DEA))


#DEA set de training
library(DESeq2)
#En caso de que DESEQ2 haga error porque los valores no son integer usar ROUND
maclear_train_DEA=round(maclear_train_DEA)

#DEA
dea_train <- DESeqDataSetFromMatrix(maclear_train_DEA, colData=maclear_train, design=~ Gene, ignoreRank = TRUE)
exp_train<-DESeq(dea_train)
res_exp_train <- results(exp_train)


#filtro los genes p menor a 0.05
library(dplyr)
res_exp_train_mtx= as.data.frame.matrix(res_exp_train)

#hago un volcano plot para observar aprox cuantos genes tengo
library(EnhancedVolcano)
EnhancedVolcano(res_exp_train_mtx,
                lab = rownames(res_exp_train_mtx),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff=0.5,
                pCutoff=0.05,
                ylim = c(0, -log10(10e-12)))
res_exp_train_mtx<-filter(res_exp_train_mtx, pvalue <= 0.05 & (log2FoldChange >= 1.5 | log2FoldChange<=-0.5))


#masajeo-obtengo lista de counts de mRNAS de genes mas significativos para el train set
train_data_exp=merge(maclear_train_DEA,res_exp_train_mtx, by='row.names')
train_data_exp<-train_data_exp[,-c(233,234,235,236,237,238,238)]
rownames(train_data_exp)<-train_data_exp[,1]
train_data_exp<-train_data_exp[,-1]
class_train<-subset(maclear_train,selec=c(60662))
class_train<-as.factor(t(class_train))
train_data_exp<-t(train_data_exp)

#Seleccion de variables para random forrest
library(varSelRF)

ranfor_train= varSelRF(train_data_exp, class_train, vars.drop.frac = 0.3)
#Esto me da la lista de genes que constituyen una posible firma basado en OOB score
plot(ranfor_train)
#extraigo la lista de genes
genes16=as.data.frame(ranfor_train$selected.vars)
colnames(genes16)<-c("Predgenes")
rownames(genes16)<-genes16$Predgenes
maclear_train_trans= as.data.frame(t(maclear_train))
maclear_train_trans=merge(maclear_train_trans,genes16, by='row.names')
maclear_train_16=as.data.frame(t(maclear_train_trans))
maclear_train_16<- maclear_train_16 %>% row_to_names(row_number = 1)
rownames(maclear_meta)<-maclear_meta$Expression_Barcode
maclear_meta_16=merge(maclear_train_16,maclear_meta, by='row.names')
maclear_train_16<-maclear_train_16[-232,]
rownames(maclear_meta_16)<-maclear_meta_16$Row.names
maclear_meta_16<-subset(maclear_meta_16, select=c( 19, 20))

#Contrasto la lista de genes obtenidos con la cohorte de entrenamiento para comprobar su eficiencia

ranfor2=randomForest(factor(Gene) ~ ., data = dd, importance = TRUE, type = "classification")
dd = merge(maclear_train_16, maclear_meta_16, by=0)
votes=as.data.frame(ranfor2$votes)

library(pROC)
roc_train= roc(dd$Gene, votes$PIK3CA, direction=">")
plot(smooth(roc_train),asp=NA,col="red", print.auc=TRUE,auc.polygon=TRUE, legacy.axes=TRUE )

library(randomForest)
maclear_valid_DEA=as.data.frame(maclear_valid_DEA)

p<-predict(ranfor2, maclear_valid, type="prob")

pi3kca_v=as.data.frame(p)

roc_valid=roc(maclear_valid$Gene, pi3kca_v$PIK3CA, direction= ">")

plot(smooth(roc_valid),asp=NA,col="red", print.auc=TRUE,auc.polygon=TRUE, legacy.axes=TRUE )

