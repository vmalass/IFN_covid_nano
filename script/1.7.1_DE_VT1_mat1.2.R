# 1-library -----
library(openxlsx)
library(umap)
library(factoextra)
library(S4Vectors)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(scater)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(clustree)
library(fpc)
library(scran)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(forcats)
library(scales)
library(pheatmap)
library(apeglm)

# 2-ouverture des fichier----
rm(list = ls())
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE R vs NR en VT1 all genes -----
## 3.1-Creation de la matrice----
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] 
mat <- mat_VT1$REPONSE %in% "RP"
mat_VT1 <- mat_VT1[mat == F,]
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

# mat<- mat_VT1[mat_VT1$numero_patient==2,]                
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==13,])     
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==50,])     
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==62,])     
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==49,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==31,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==16,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==8,]) #creation d'une lmat avec R et NR choisi


coldata_mat<-as.data.frame(mat_VT1$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT1)
coldata_num<-mat_VT1$numero_patient
mat_VT1[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene----
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT1), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap


##### NR
dds_mat$condition <- relevel(dds_mat$condition, ref = "NR")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "NR", "R"))

# NR vs R
res_NR_R <- lfcShrink(dds_mat, coef = "condition_R_vs_NR", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

NR_R <- res_NR_R[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
NR_R_all_gene <- as.data.frame(NR_R)

data_NR_R_VT1 <- assay(rld)[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_NR_R_VT1)

heatmap_NR_R_VT1 <- pheatmap(data_NR_R_VT1, 
                             scale="row", 
                             fontsize_row=12, 
                             fontsize_col=15, 
                             annotation_col = annC_VT2_DE,
                             annotation_colors = list(condition = c( NR = "#7570BE",
                                                                     R="#F15854")),
                             color = my_palette, 
                             cutree_rows = 2,
                             main = "DE à partir de tout les gènes et des prélèvements VT1 avec NRvsR", 
                             labels_col = coldata_num,
                             cluster_cols = T,
                             cluster_rows = T,)

# 4-DE R vs NRen VT1 HVG----
## 4.1-Creation de la matrice---

# mat<- mat_VT1[mat_VT1$numero_patient==2,]                
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==13,])     
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==50,])     
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==62,])     
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==49,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==31,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==16,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==8,]) #creation d'une lmat avec R et NR choisi

mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] 
mat <- mat_VT1$REPONSE %in% "RP"
mat_VT1 <- mat_VT1[mat == F,]
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

coldata_num<- mat_VT1$numero_patient

top.hvgs<-as.data.frame(top.hvgs)    ##creation d'un DF
test <- colnames(mat_VT1) %in% top.hvgs$V1  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_HVG<- mat_VT1[,test==T]   #ne garde en colone que les True

all(rownames(mat_HVG) == rownames(mat_VT1)) #verif ordre identique
all(top.hvgs$V1 %in% colnames(mat_HVG)) #verif contient bien les meme nom de genes

coldata_mat<-as.data.frame(mat_VT1$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT1)

## 4.2-Gene DE DESeq2 HVG----
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_HVG), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap


##### R
dds_mat$condition <- relevel(dds_mat$condition, ref = "R")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "NR", "R"))
res <- as.data.frame(res)

# NR vs R
res_NR_R <- lfcShrink(dds_mat, coef = "condition_NR_vs_R", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

NR_R <- res_NR_R[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
NR_R_HVG <- as.data.frame(NR_R)

data_NR_R_VT1 <- assay(rld)[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_NR_R_VT1)

heatmap_NR_R_VT1_HVG <- pheatmap(data_NR_R_VT1, 
                                 scale="row", 
                                 fontsize_row=12,
                                 fontsize_col=15,
                                 annotation_col = annC_VT2_DE, 
                                 annotation_colors = list(condition = c( NR = "#7570BE",
                                                                         R="#F15854")),
                                 color = my_palette, 
                                 cutree_rows = 2, 
                                 main = "DE à partir des HVG et des prélèvements VT1 avec NRvsR", 
                                 labels_col = coldata_num,
                                 cluster_cols = T,
                                 cluster_rows = T,)

name_NR_R_HVG<- rownames(NR_R_HVG)
name_NR_R_all_gene<- rownames(NR_R_all_gene)
setdiff(name_NR_R_HVG, name_NR_R_all_gene)
setdiff(name_NR_R_all_gene,name_NR_R_HVG)
union(name_NR_R_HVG, name_NR_R_all_gene)

# 5-heatmap VT1 NRvsR all gene T-----
## 5.1-creation de la matrice----

# mat<- mat_VT1[mat_VT1$numero_patient==2,]
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==13,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==50,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==60,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==62,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==49,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==9,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==31,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==16,])
# mat<-rbind(mat,mat_VT1[mat_VT1$numero_patient==8,]) #creation d'une lmat avec R et NR choisi


mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T
# mat<-rbind(mat, mat_T)
# 
# mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),]
# mat_RP <- mat_VT1[which(mat_VT1$REPONSE %in% "RP"),]
# mat<-rbind(mat, mat_RP)
# 
# mat_R <- mat_VT1[which(mat_VT1$REPONSE %in% "R"),]
# mat_R<- mat_R[mat_R$numero_patient!=31,]
# mat_R<- mat_R[mat_R$numero_patient!=16,]
# mat_R<- mat_R[mat_R$numero_patient!=8,]
# mat_R<- mat_R[mat_R$numero_patient!=9,]
# mat_R<- mat_R[mat_R$numero_patient!=60,]
# mat<-rbind(mat, mat_R)

mat<-rbind(mat_VT1,mat_T)

## 5.2-Gene DE DESeq2 HVG----

mat_meta<-mat[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat$numero_patient

name_NR_R_all_gene<-as.data.frame(name_NR_R_all_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_NR_R_all_gene$name_NR_R_all_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_NR_R_HVG_all, 
         scale="row", 
         fontsize_row=10, 
         fontsize_col = 10,
         annotation_colors = list(mat_meta.REPONSE = c(NR = "#7570BE",
                                                       R="#F15854",
                                                       RP = "#882255",
                                                       "T"= "#117733" )),
         annotation_col = annC_VT4,
         color = my_palette, 
         cutree_rows = 2,
         main = "DE à partir de tout les gènes et des prélèvements VT1 avec NRvsR", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

# 6-heatmap VT1 NRvsR HVG avec T-----
## 6.1-Creation de la matrice-----
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat<-rbind(mat_VT1,mat_T)

mat_meta<-mat[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat$numero_patient

name_NR_R_HVG<-as.data.frame(name_NR_R_HVG)    ##creation d'un DF
test <- colnames(mat) %in% name_NR_R_HVG$name_NR_R_HVG  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_NR_R_HVG_all,
         scale="row",
         fontsize_row=10, 
         fontsize_col = 10,
         annotation_colors = list(mat_meta.REPONSE = c(NR = "#7570BE",
                                                       R="#F15854",
                                                       RP = "#882255",
                                                       "T"= "#117733" )),
         annotation_col = annC_VT4, 
         color = my_palette, 
         cutree_rows = 2,
         main = "DE à partir des HVG et des prélèvements VT1 avec NRvsR",
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

# pheatmap(mat_DE_NR_R_HVG_all, scale="row", fontsize_row=10, fontsize_col = 10,annotation_colors = list(mat_meta.REPONSE = c(REA = "#0000FF", NR = "#7570BE", "T"= "#117733", R="#F15854")),annotation_col = annC_VT4, color = my_palette, cutree_rows = 1, main = "DE à partir des HVG et des prélèvements VT1 avec NRvsR", labels_col = coldata_num,cluster_cols = F,cluster_rows = T,)


# 7-heatmap VT1 NRvsR HVG et all gene T-----
## 7.1-Creation de la matrice----
name_DE_NR_R_HVG_all<- union(name_NR_R_HVG$name_NR_R_HVG,name_NR_R_all_gene$name_NR_R_all_gene)

mat_meta<-mat[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat$numero_patient

name_DE_NR_R_HVG_all<-as.data.frame(name_DE_NR_R_HVG_all)    ##creation d'un DF
test <- colnames(mat) %in% name_DE_NR_R_HVG_all$name_DE_NR_R_HVG_all  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap
annC_VT4<-rename(annC_VT4, Groupe = mat_meta.REPONSE)

pheatmap(mat_DE_NR_R_HVG_all,
         scale="row",
         fontsize_row=12,
         fontsize_col = 15,
         annotation_colors = list(Groupe = c(NR = "#7570BE",
                                                       R="#F15854",
                                                       # RP = "#882255",
                                                       "T"= "#117733" )),
         annotation_col = annC_VT4, 
         color = my_palette, 
         cutree_rows = 2, 
         main = "Heatmap à partir des gènes uniques du DE de tout les gènes et des HVG en VT1", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)
# pheatmap(mat_DE_NR_R_HVG_all, scale="row", fontsize_row=7, fontsize_col = 10,annotation_colors = list(mat_meta.REPONSE = c(REA = "#0000FF", NR = "#7570BE", "T"= "#117733", R="#F15854")),annotation_col = annC_VT4, color = my_palette, cutree_rows = 1, main = "Heatmap à partir des gènes uniques du DE de tout les gènes et des HVG", labels_col = coldata_num,cluster_cols = F,cluster_rows = T,)

# 8- Save data-----
write.table(name_DE_NR_R_HVG_all, file = "data/name_DE_NR_R_HVG_all_VT1_mat1.2.txt")
a <- read.table("data/name_DE_NR_R_HVG_all_VT1_mat1.2.txt")
a

