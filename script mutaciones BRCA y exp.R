#query and download de TCGA
library(TCGAbiolinks)
library(DT)
library(dplyr)
options(timeout=3600)
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE,
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)


GDCdownload(query)
prep<-GDCprepare(query)

prep<-prep[,-1]

#query and download clinical data
query2 <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)
GDCdownload(query2)
clinical.BCRtab.all <- GDCprepare(query2)
names(clinical.BCRtab.all)
clin=tibble::tibble(sort(colnames(clinical.BCRtab.all$clinical_patient_brca)))
#masajeo
prep<-subset(prep, !is.na(HGVSc))
#genes mas mutados en BRCA(por mutacion específica)
#table() cuenta la frecuencia
prep$MutGenType <- paste(prep$Hugo_Symbol,prep$HGVSc)
freq=table(prep$MutGenType)
freq=data.frame(freq)
freq=as.data.frame.matrix(freq)
colnames(freq) <- c("Gene_Mut", "Freq")
freq2=filter(freq, Freq >= 6)
freq2 <- freq2[order(freq2$Freq,decreasing = TRUE),]

#Para filtrar mutaciones por ER status
plyr::count(clinical.BCRtab.all$clinical_patient_brca$er_status_by_ihc)
clin_and_mut_ER <- grep("^er",colnames(clinical.BCRtab.all$clinical_patient_brca))
er_stat=clinical.BCRtab.all$clinical_patient_brca[,c(2,clin_and_mut_ER)] 

er_stat<-er_stat[-c(1,2),]
#masajeo y merge de tablas para asignar ER status a las muestras(mutaciones)#
#relocate() funciona solo en library dplyr#
prep$TCGA_ID_part<-NA
prep<-prep %>% relocate(TCGA_ID_part)
#con substr() recorto parte del TCGA code para poder mergear con la tabla clinica#
prep$TCGA_ID_part<-substr(prep$Tumor_Sample_Barcode,1,12)

mut_brca_er_stat=merge(prep, er_stat, by.x=1,by.y=1)

#para acomodar las columnas de ER mas cerca del principio
mut_brca_er_stat<-mut_brca_er_stat %>% relocate(er_status_by_ihc, .after=TCGA_ID_part)
mut_brca_er_stat<-mut_brca_er_stat %>% relocate(er_status_ihc_Percent_Positive, .after=er_status_by_ihc)
#filtrar mutaciones erpos y erneg
mut_brca_er_stat_pos=mut_brca_er_stat[mut_brca_er_stat$`er_status_by_ihc`=="Positive",]
mut_brca_er_stat_neg=mut_brca_er_stat[mut_brca_er_stat$er_status_by_ihc=="Negative",]
#frecuencia de mutaciones en er positivos
freq_mut_er_pos=as.data.frame(table(mut_brca_er_stat_pos$MutGenType))
freq_mut_er_neg=as.data.frame(table(mut_brca_er_stat_neg$MutGenType))

#DIFFERENTIAL GENE EXPRESSION EN TUMORES ER POS P53 VS PIK3CA


#es al pedo pero si quiero hacer subsets es así
mut_brca_er_stat_neg_p53= subset(mut_brca_er_stat_neg, (Hugo_Symbol=="TP53"))
mut_brca_er_stat_neg_pi3= subset(mut_brca_er_stat_neg, (Hugo_Symbol=="PIK3CA"))

#filtro los que tienen los genes p53 o pik3ca, el operador "|"(tecla altgr+1) es "or"
mut_brca_er_stat_pos_p53_pi3= subset(mut_brca_er_stat_pos, (Hugo_Symbol=="TP53"| Hugo_Symbol=="PIK3CA"))
#saco los que están duplicados
duplo=as.data.frame(table(mut_brca_er_stat_pos_p53_pi3$Tumor_Sample_Barcode))
mut_brca_er_stat_pos_p53_pi3=mut_brca_er_stat_pos_p53_pi3[!duplicated(mut_brca_er_stat_pos_p53_pi3$Tumor_Sample_Barcode),]



#PLOTS GENERALES
#Barplot-argumento "las" cambia orientación de las labels#
barplot(height=freq2$Freq, names=freq2$Gene_Mut, las=2,cex.axis=1, cex.names=0.5, ylab="Frecuencia", main= "Mutaciones más frecuentes en BCRA")

#plot generales

library(maftools)

prep2 = read.maf(maf = prep)

datatable(getSampleSummary(prep2),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = prep2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

oncoplot(maf =prep2, top = 10, removeNonMutated = TRUE)
titv = titv(maf = prep2, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)

#lollipop plot

lollipopPlot(
  maf = prep2,
  gene = 'TP53',
  showMutationRate = TRUE,)
#Somatic interactions
somint=somaticInteractions(maf = prep2, top = 25, pvalue = c(0.05, 0.1))
library(readr)
write_tsv(somint, "C:/Users/User/Lab/idisba/FITR2/somatic_interactions.tsv")

#descarga de datos de expression

query_exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  file.type = "rna_seq.augmented_star_gene_counts.tsv",)
GDCdownload(query_exp)
exp_data<-GDCprepare(query_exp)

#traspngo la tabla para machar las muestras de mutación
library(SummarizedExperiment)
BRCA_exp_Matrix <- assay(exp_data,"unstranded")
BRCA_exp_Matrix<-as.data.frame.matrix(t(BRCA_exp_Matrix))

#MASAJEO para machear las ER positivas mutadas con los datos de expresion
#relocate funciona solo en library dplyr
BRCA_exp_Matrix$partialID<-NA
BRCA_exp_Matrix$completeID<-NA
BRCA_exp_Matrix<-BRCA_exp_Matrix %>% relocate(partialID)
BRCA_exp_Matrix<-BRCA_exp_Matrix %>% relocate(completeID)
BRCA_exp_Matrix$completeID<-rownames(BRCA_exp_Matrix)
#hago una columna de TCGA barcode hasta el character 15, que corresponde a muestra
BRCA_exp_Matrix$partialID<-substr(BRCA_exp_Matrix$completeID,1,15)
mut_brca_er_stat_pos_p53_pi3$TCGA_ID_part2<-substr(mut_brca_er_stat_pos_p53_pi3$Tumor_Sample_Barcode,1,15)
#idem comentario anterior. Además creo una tabla solo con los barcodes parciales de mutacion
mut_brca_er_stat_pos_p53_pi3<-mut_brca_er_stat_pos_p53_pi3 %>% relocate(TCGA_ID_part2)
er_pos_samples<-data.frame(mut_brca_er_stat_pos_p53_pi3$TCGA_ID_part2,
mut_brca_er_stat_pos_p53_pi3$er_status_by_ihc,mut_brca_er_stat_pos_p53_pi3$Tumor_Sample_Barcode)
er_pos_samples<-er_pos_samples[!duplicated(er_pos_samples$mut_brca_er_stat_pos_p53_pi3.TCGA_ID_part2),]

#Merge de tablas
BRCA_exp_er_pos=merge(BRCA_exp_Matrix,er_pos_samples,by.x=2,by.y=1)
#ajuste de tabla de expresión de muestras pi3 y p53 mutadas con er positivo
rownames(BRCA_exp_er_pos)<-BRCA_exp_er_pos[,2]
#Filtro por pureza tumoral mayor a 60%
CPE_BRCA=read.table("c:/Users/User/Lab/idisba/FITR2/CPE_BCRA.txt",header=TRUE,fill=TRUE)
CPE_filt=filter(CPE_BRCA, CPE >= 0.6)
BRCA_exp_er_pos=merge(BRCA_exp_er_pos,CPE_filt, by.x=1,by.y=2)
BRCA_exp_er_pos=BRCA_exp_er_pos[!duplicated(BRCA_exp_er_pos$completeID),]
rownames(BRCA_exp_er_pos)<-BRCA_exp_er_pos[,2]
#lista final de muestras ER positivas con p53 y pi3 mutado
BRCA_exp_er_pos_Clean=BRCA_exp_er_pos[,-c(1,2,60663,60664,60665,60666,60667,60668,60669,60670,60671)]
BRCA_exp_er_pos_Clean2=as.data.frame(t(BRCA_exp_er_pos_Clean))
samplelist<-c(rownames(BRCA_exp_er_pos_Clean))

#data query&download de las muestras filtradas

query_exp2 <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  barcode=samplelist,)
GDCdownload(query_exp2)
exp_data_POSTA<-GDCprepare(query_exp2)
#conversión de matriz
BRCAMatrix <- assay(exp_data_POSTA,"unstranded")
#análisis de outliers
BRCA.RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(exp_data_POSTA)

metadata_exp=data.frame(BRCA_exp_er_pos$completeID,BRCA_exp_er_pos$mut_brca_er_stat_pos_p53_pi3.Tumor_Sample_Barcode)
metadata_mut=data.frame(mut_brca_er_stat_pos_p53_pi3$Hugo_Symbol,mut_brca_er_stat_pos_p53_pi3$Tumor_Sample_Barcode)
colnames(metadata_exp)<-c("Expression_Barcode", "Mutation_Barcode")
colnames(metadata_mut)<-c("Gene", "Mutation_Barcode")
exp_metadata=merge(metadata_exp,metadata_mut,by.x=2,by.y=2)


#Differential Expression Analysis TP53 vs PIK3CA

library(DESeq2)

ddsSE2 <- DESeqDataSetFromMatrix(BRCAMatrix, colData=exp_metadata, design=~ Gene)
exp_pre<-DESeq(ddsSE2)
res <- results(exp_pre)

#macheo de datos ENSEMBL con Gene symbol
res2=as.data.frame(res)
symb=read.table("c:/Users/User/Lab/idisba/FITR2/ensembl_to_hugo.txt",header=TRUE, fill=TRUE)
symb<-symb[,c(2,5)]
res2$Ensembl<-NA
res2$Ensembl<-rownames(res2)
res3=merge(res2,symb, by.x=7, by.y=1)
res3<-res3[!duplicated(res3$Ensembl),]
#con esta expresión se puede elegir que condición versus que condición hacemos.
res4<-results(exp_pre,contrast=c("Gene","PIK3CA","TP53"))

#masajeo para plot
res3<-res3 %>% relocate(stable.1)
res3<-res3[!duplicated(res3$stable.1),]
rownames(res3)<-res3$stable.1
res3<-na.omit(res3)
res3<-res3[,-c(1,2)]

library(EnhancedVolcano)

EnhancedVolcano(res3,
                lab = rownames(res3),
                x = 'log2FoldChange',
                y = 'padj')


#preparo tabla para empezar machine learning
data_maclear=merge(metadata_exp,metadata_mut, by.x=2,by.y=2)
data_maclear2=merge(data_maclear,BRCA_exp_er_pos_Clean, by.x=2,by.y='row.names')
library(readr)
write_tsv(data_maclear2, "c:/Users/User/Lab/idisba/FITR2/maclear/exp_data_maclear.tsv")
#Adios