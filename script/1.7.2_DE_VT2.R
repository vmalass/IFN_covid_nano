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
setwd("~/Documents/JM/NanoString/NanoString_Covid/nanostring_covid/data") #folder data
rm(list = ls())
load("1.3_mat_pat_clean_final.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE R vs RP en VT2 all genes -----
## 3.1-Creation de la matrice----
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] 

mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat


coldata_mat<-as.data.frame(mat_VT2$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT2)
coldata_num<-mat_VT2$numero_patient
mat_VT2[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene----
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT2), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap


##### RP
dds_mat$condition <- relevel(dds_mat$condition, ref = "RP")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "R", "RP"))

# RP vs R
res_NR_R <- lfcShrink(dds_mat, coef = "condition_R_vs_RP", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

NR_R <- res_NR_R[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
NR_R_all_gene <- as.data.frame(NR_R)

data_NR_R_VT2 <- assay(rld)[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_NR_R_VT2)

heatmap_NR_R_VT2 <- pheatmap(data_NR_R_VT2, 
                             scale="row", 
                             fontsize_row=10, 
                             annotation_col = annC_VT2_DE,
                             annotation_colors = list(condition = c( NR = "#7570BE",
                                                                     R="#F15854",
                                                                     RP = "#882255")),
                             color = my_palette, 
                             cutree_rows = 2,
                             main = "DE à partir de tout les gènes et des prélèvements VT2 avec NRvsR", 
                             labels_col = coldata_num,
                             cluster_cols = T,
                             cluster_rows = T,)


# 4-DE R vs RP en VT2 HVG----
## 4.1-Creation de la matrice---

mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] 

mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat

coldata_num<- mat_VT2$numero_patient

top.hvgs<-as.data.frame(top.hvgs)    ##creation d'un DF
test <- colnames(mat_VT2) %in% top.hvgs$V1  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_HVG<- mat_VT2[,test==T]   #ne garde en colone que les True

all(rownames(mat_HVG) == rownames(mat_VT2)) #verif ordre identique
all(top.hvgs$V1 %in% colnames(mat_HVG)) #verif contient bien les meme nom de genes

coldata_mat<-as.data.frame(mat_VT2$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT2)


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
res <- results(dds_mat, contrast=c("condition", "RP", "R"))
res <- as.data.frame(res)

# RP vs R
res_NR_R <- lfcShrink(dds_mat, coef = "condition_RP_vs_R", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

NR_R <- res_NR_R[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
NR_R_HVG <- as.data.frame(NR_R)

data_NR_R_VT2 <- assay(rld)[res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange > 1 | res_NR_R$svalue < 0.05 & !is.na(res_NR_R$svalue) & res_NR_R$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_NR_R_VT2)

heatmap_NR_R_VT2_HVG <- pheatmap(data_NR_R_VT2, 
                                 scale="row", 
                                 fontsize_row=10, 
                                 annotation_col = annC_VT2_DE, 
                                 annotation_colors = list(condition = c(  NR = "#7570BE",
                                                                          R="#F15854",
                                                                          RP = "#882255")),
                                 color = my_palette, 
                                 cutree_rows = 2, 
                                 main = "DE à partir des HVG et des prélèvements VT2 avec NRvsR", 
                                 labels_col = coldata_num,
                                 cluster_cols = T,
                                 cluster_rows = T,)

name_NR_R_HVG<- rownames(NR_R_HVG)
name_NR_R_all_gene<- rownames(NR_R_all_gene)
setdiff(name_NR_R_HVG, name_NR_R_all_gene)
setdiff(name_NR_R_all_gene,name_NR_R_HVG)

# 5-heatmap VT2 RvsRP all gene T-----
## 5.1-creation de la matrice----
name_DE_NR_R_HVG_all<- name_NR_R_all_gene

mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] 

mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat

mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T
mat_VT2<-rbind(mat_VT2, mat_T)

## 5.2-Gene DE DESeq2 HVG----

mat_meta<-mat_VT2[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat_VT2$numero_patient

name_DE_NR_R_HVG_all<-as.data.frame(name_DE_NR_R_HVG_all)    ##creation d'un DF
test <- colnames(mat_VT2) %in% name_DE_NR_R_HVG_all$name_DE_NR_R_HVG_all  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat_VT2[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_NR_R_HVG_all, 
         scale="row", 
         fontsize_row=10, 
         fontsize_col = 10,
         annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
                                                        R="#F15854",
                                                        RP = "#882255",
                                                       "T"= "#117733" )),
         annotation_col = annC_VT4,
         color = my_palette, 
         cutree_rows = 1,
         main = "DE à partir de tout les gènes et des prélèvements VT2 avec NRvsR", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

# 6-heatmap VT2 RvsRP HVG avec T-----
## 6.1-Creation de la matrice-----
name_DE_NR_R_HVG_all<- name_NR_R_HVG

mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] 

mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat

mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T
mat_VT2<-rbind(mat_VT2, mat_T)

mat_meta<-mat_VT2[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat_VT2$numero_patient

name_DE_NR_R_HVG_all<-as.data.frame(name_DE_NR_R_HVG_all)    ##creation d'un DF
test <- colnames(mat_VT2) %in% name_DE_NR_R_HVG_all$name_DE_NR_R_HVG_all  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat_VT2[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_NR_R_HVG_all,
         scale="row",
         fontsize_row=10, 
         fontsize_col = 10,
         annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
                                                        R="#F15854",
                                                        RP = "#882255",
                                                       "T"= "#117733" )),
         annotation_col = annC_VT4, 
         color = my_palette, 
         cutree_rows = 1,
         main = "DE à partir des HVG et des prélèvements VT2 avec NRvsR",
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

# pheatmap(mat_DE_NR_R_HVG_all, scale="row", fontsize_row=10, fontsize_col = 10,annotation_colors = list(mat_meta.REPONSE = c(REA = "#0000FF", RP = "#7570BE", "T"= "#117733", R="#F15854")),annotation_col = annC_VT4, color = my_palette, cutree_rows = 1, main = "DE à partir des HVG et des prélèvements VT2 avec NRvsR", labels_col = coldata_num,cluster_cols = F,cluster_rows = T,)


# 7-heatmap VT2 RvsRP HVG et all gene T-----
## 7.1-Creation de la matrice----
name_DE_NR_R_HVG_all<- union(name_NR_R_HVG,name_NR_R_all_gene)

mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] 

mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat

mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T
mat_VT2<-rbind(mat_VT2, mat_T)

mat_meta<-mat_VT2[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat_VT2$numero_patient

name_DE_NR_R_HVG_all<-as.data.frame(name_DE_NR_R_HVG_all)    ##creation d'un DF
test <- colnames(mat_VT2) %in% name_DE_NR_R_HVG_all$name_DE_NR_R_HVG_all  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat_VT2[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_NR_R_HVG_all,
         scale="row",
         fontsize_row=10,
         fontsize_col = 17,
         annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
                                                        R="#F15854",
                                                        RP = "#882255",
                                                       "T"= "#117733" )),
         annotation_col = annC_VT4, 
         color = my_palette, 
         cutree_rows = 2, 
         main = "Heatmap à partir des gènes uniques du DE de tout les gènes et des HVG", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)
# pheatmap(mat_DE_NR_R_HVG_all, scale="row", fontsize_row=7, fontsize_col = 10,annotation_colors = list(mat_meta.REPONSE = c(REA = "#0000FF", RP = "#7570BE", "T"= "#117733", R="#F15854")),annotation_col = annC_VT4, color = my_palette, cutree_rows = 1, main = "Heatmap à partir des gènes uniques du DE de tout les gènes et des HVG", labels_col = coldata_num,cluster_cols = F,cluster_rows = T,)

# 8-heatmap VT2/VT3/VT4 à partir des gènes unique du DE de tout les gènes et des HVG
## 8.1-Creation de la matrice
name_DE_NR_R_HVG_all<- union(name_NR_R_HVG,name_NR_R_all_gene)

mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] 
mat_VT2<- arrange(mat_VT2,REPONSE) #ordo 
mat_VT3 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT3"),] 
mat_VT3<- arrange(mat_VT3,REPONSE) #ordo 
mat_VT4 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT4"),] 
mat_VT4<- arrange(mat_VT4,REPONSE) #ordo 
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat_VT<-rbind(mat_VT2, mat_VT3)
mat_VT<-rbind( mat_VT,mat_VT4)
mat_VT<-rbind( mat_VT,mat_T)

mat_meta<-mat_VT[,737:743]
coldata_VT4<-as.data.frame(mat_meta$REPONSE)
row.names(coldata_VT4)<-row.names(mat_meta)

coldata_num<- mat_VT$numero_patient

name_DE_NR_R_HVG_all<-as.data.frame(name_DE_NR_R_HVG_all)    ##creation d'un DF
test <- colnames(mat_VT) %in% name_DE_NR_R_HVG_all$name_DE_NR_R_HVG_all  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_NR_R_HVG_all<- mat_VT[,test==T]   #ne garde en colone que les True

mat_DE_NR_R_HVG_all<-log1p(t(mat_DE_NR_R_HVG_all)) # passage en log 1 p pour les heatmap
annC_VT4 <- data.frame(condition= coldata_VT4) # pour les annotation en col pour la heatmap
all(rownames(annC_VT4) == colnames(mat_DE_NR_R_HVG_all))  # pour les annotation en col pour la heatmap
annC_VT4$real_time_point <- mat_meta$real_time_point  # ajout real_time_point

pheatmap(mat_DE_NR_R_HVG_all,
         scale="row",
         fontsize_row=12,
         fontsize_col = 10,
         annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
                                                        R="#F15854",
                                                        RP = "#882255",
                                                        "T"= "#117733" ),
                                  real_time_point = c(VT2="gold3",
                                                      VT3="azure3",
                                                      VT4="chocolate3",
                                                      "T"="#117733")),
         annotation_col = annC_VT4, 
         color = my_palette, 
         cutree_rows = 2, 
         main = "Heatmap à partir des gènes uniques du DE de tout les gènes et des HVG", 
         labels_col = coldata_num,
         cluster_cols = F,
         cluster_rows = T,)

pheatmap(mat_DE_NR_R_HVG_all,
         scale="row",
         fontsize_row=12,
         fontsize_col = 10,
         annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
                                                        R="#F15854",
                                                        RP = "#882255",
                                                        "T"= "#117733" ),
                                  real_time_point = c(VT2="gold3",
                                                      VT3="azure3",
                                                      VT4="chocolate3",
                                                      "T"="#117733")),
         annotation_col = annC_VT4, 
         color = my_palette, 
         cutree_rows = 2, 
         main = "Heatmap à partir des gènes uniques du DE de tout les gènes et des HVG", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

