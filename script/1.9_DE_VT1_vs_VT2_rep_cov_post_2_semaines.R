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
library(ggvenn)
# 2-ouverture des fichier----
rm(list = ls())
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE T vs R+RP en VT1 all genes -----
## 3.1-Creation de la matrice----
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] #DF avce tout les VT1
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T
mat_VT1<-rbind(mat_VT1, mat_T) # DF avec VT1 & T
mat_NR<-mat_VT1[1:5,]#conservation des NR
mat_VT1<-mat_VT1[6:42,] #suppression des NR pour pouvoir faire le DE avec la colonne condition biologique (covid vs T)

mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat


coldata_mat<-as.data.frame(mat_VT1$condition_biologique)
row.names(coldata_mat)<-row.names(mat_VT1)
coldata_num<-mat_VT1$numero_patient
mat_VT1[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene----
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT1), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

##### T
dds_mat$condition <- relevel(dds_mat$condition, ref = "T")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "T", "Covid"))

# covid vs T
res_cov_T <- lfcShrink(dds_mat, coef = "condition_Covid_vs_T", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

cov_T <- res_cov_T[res_cov_T$svalue < 0.05 & !is.na(res_cov_T$svalue) & res_cov_T$log2FoldChange > 1 | res_cov_T$svalue < 0.05 & !is.na(res_cov_T$svalue) & res_cov_T$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
cov_T_all_gene <- as.data.frame(cov_T)

data_cov_T_VT2 <- assay(rld)[res_cov_T$svalue < 0.05 & !is.na(res_cov_T$svalue) & res_cov_T$log2FoldChange > 1 | res_cov_T$svalue < 0.05 & !is.na(res_cov_T$svalue) & res_cov_T$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_cov_T_VT2)

heatmap_NR_R_VT2 <- pheatmap(data_cov_T_VT2, 
                             scale="row", 
                             fontsize_row=10, 
                             annotation_col = annC_VT2_DE,
                             color = my_palette, 
                             cutree_rows = 2,
                             main = "DE à partir de tout les gènes et des prélèvements VT2 avec NRvsR", 
                             labels_col = coldata_num,
                             cluster_cols = T,
                             cluster_rows = T,)

# 4-DE T vs RP en VT2 all gene-----
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] #DF avce tout les VT1
mat_VT2 <-rbind(mat_VT2, mat_T) # DF avec VT1 & T

mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat


coldata_mat<-as.data.frame(mat_VT2$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT2)
coldata_num<-mat_VT2$numero_patient
mat_VT2[,737:744]<-NULL

## 4.2-Gene DE DESeq2 all gene----
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT2), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

##### T
dds_mat$condition <- relevel(dds_mat$condition, ref = "T")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "T", "RP"))

# RP vs T
res_RP_T <- lfcShrink(dds_mat, coef = "condition_RP_vs_T", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

RP_T <- res_RP_T[res_RP_T$svalue < 0.05 & !is.na(res_RP_T$svalue) & res_RP_T$log2FoldChange > 1 | res_RP_T$svalue < 0.05 & !is.na(res_cov_T$svalue) & res_cov_T$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
RP_T_all_gene <- as.data.frame(RP_T)

data_RP_T_VT2 <- assay(rld)[res_RP_T$svalue < 0.05 & !is.na(res_RP_T$svalue) & res_RP_T$log2FoldChange > 1 | res_RP_T$svalue < 0.05 & !is.na(res_cov_T$svalue) & res_cov_T$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_RP_T_VT2)

heatmap_NR_R_VT2 <- pheatmap(data_RP_T_VT2, 
                             scale="row", 
                             fontsize_row=10, 
                             annotation_col = annC_VT2_DE,
                             annotation_colors = list(mat_meta.REPONSE = c(NR = "darkorange",
                                                                           `NR-` = "#DDCC77",
                                                                           R = "cornflowerblue",
                                                                           `RP-` = "#882255" ,
                                                                           RP = "brown3", 
                                                                           `T` = "chartreuse4", 
                                                                           A = "#BBBBBB")),
                             color = my_palette, 
                             cutree_rows = 2,
                             main = "DE à partir de tout les gènes et des prélèvements VT2 avec RPvsT", 
                             labels_col = coldata_num,
                             cluster_cols = T,
                             cluster_rows = T,)


name_RP_T_all_gene<-rownames(RP_T_all_gene)
name_cov_T_all_gene<-rownames(cov_T_all_gene)
name_inter_cov_RP_T_gene<-intersect(name_RP_T_all_gene,name_cov_T_all_gene) #intersection entre les deux vecteurs
name_unique_RP_T_gene<-setdiff(name_RP_T_all_gene,name_cov_T_all_gene)
name_unique_cov_T_gene<-setdiff(name_cov_T_all_gene,name_RP_T_all_gene)

set.seed(12)
liste_venn<-list(VT1=name_cov_T_all_gene ,
                 VT2=name_RP_T_all_gene)

ggvenn(liste_venn, fill_color = c("#CC6677", "#999933"), )+
  ggtitle("Diagramme de Venn")+
  theme(plot.title = element_text(size=25, hjust = 0.5, vjust = -5))

# 5-Resume----
## 5.1-18genes VT1+VT2-----

mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] #DF avce tout les VT1
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] #DF avce tout les VT1
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat<-rbind(mat_VT1,mat_VT2)
mat<-rbind(mat,mat_T)
mat_orig<-mat

name_unique_cov_T_gene<-as.data.frame(name_unique_cov_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_unique_cov_T_gene$name_unique_cov_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat<- mat[,test==T]   #ne garde en colone que les True

coldata<-mat_orig[,737:743] # creation du coldata pour les heatmap
coldata_num<- coldata$numero_patient  # ordo
coldata[,1:6]<-NULL # sup des data non utile

mat<-t(mat)  # translation de la matrice
all(rownames(coldata) == colnames(mat)) #verif col et row identique
all(rownames(coldata) %in% colnames(mat)) #attention de verif pas l'ordre mais seulement présence
mat_log<-log1p(mat) # passage en log 1 p pour les heatmap
annC <- data.frame(condition= coldata) # pour les annotation en col pour la heatmap
rownames(annC) <- colnames(mat_log)  # pour les annotation en col pour la heatmap
all(rownames(annC) == colnames(mat)) #vérif du même ordre / position
annC$real_time_point <- mat_orig$real_time_point  # ajout real_time_point

pheatmap(mat_log,
         scale="row",
         fontsize_row=15,
         fontsize_col = 12,
         annotation_col = annC,
         annotation_colors = list(mat_meta.REPONSE = c(NR = "darkorange",
                                                        `NR-` = "#DDCC77",
                                                        R = "cornflowerblue",
                                                        `RP-` = "#882255" ,
                                                        RP = "brown3", 
                                                        `T` = "chartreuse4", 
                                                        A = "#BBBBBB"), 
                                  real_time_point = c(VT1="deepskyblue",
                                                      VT2="gold3",
                                                      "T"="#117733")),
         color = my_palette, 
         cutree_rows = 1, 
         main = "18 gènes avec VT1+VT2", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

## 5.2-65 genes VT1+VT2-----

mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] #DF avce tout les VT1
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] #DF avce tout les VT1
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat<-rbind(mat_VT1,mat_VT2)
mat<-rbind(mat,mat_T)
mat_orig<-mat


name_inter_cov_RP_T_gene<-as.data.frame(name_inter_cov_RP_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_inter_cov_RP_T_gene$name_inter_cov_RP_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat<- mat[,test==T]   #ne garde en colone que les True

coldata<-mat_orig[,737:743] # creation du coldata pour les heatmap
coldata_num<- coldata$numero_patient  # ordo
coldata[,1:6]<-NULL # sup des data non utile

mat<-t(mat)  # translation de la matrice
all(rownames(coldata) == colnames(mat)) #verif col et row identique
all(rownames(coldata) %in% colnames(mat)) #attention de verif pas l'ordre mais seulement présence
mat_log<-log1p(mat) # passage en log 1 p pour les heatmap
annC <- data.frame(condition= coldata) # pour les annotation en col pour la heatmap
rownames(annC) <- colnames(mat_log)  # pour les annotation en col pour la heatmap
all(rownames(annC) == colnames(mat)) #vérif du même ordre / position
annC$real_time_point <- mat_orig$real_time_point  # ajout real_time_point

pheatmap(mat_log,
         scale="row",
         fontsize_row=11,
         fontsize_col = 12,
         annotation_col = annC,
         annotation_colors = list(mat_meta.REPONSE = c(NR = "darkorange",
                                                        `NR-` = "#DDCC77",
                                                        R = "cornflowerblue",
                                                        `RP-` = "#882255" ,
                                                        RP = "brown3", 
                                                        `T` = "chartreuse4", 
                                                        A = "#BBBBBB"), 
                                   real_time_point = c(VT1="deepskyblue",
                                                       VT2="gold3",
                                                       "T"="#117733")),
         color = my_palette, 
         color = my_palette, 
         cutree_rows = 1, 
         main = "65 gènes avec VT1+VT2", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

## 5.3-2 genes VT1+VT2-----

mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] #DF avce tout les VT1
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] #DF avce tout les VT1
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat<-rbind(mat_VT1,mat_VT2)
mat<-rbind(mat,mat_T)
mat_orig<-mat


name_unique_RP_T_gene<-as.data.frame(name_unique_RP_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_unique_RP_T_gene$name_unique_RP_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat<- mat[,test==T]   #ne garde en colone que les True

coldata<-mat_orig[,737:743] # creation du coldata pour les heatmap
coldata_num<- coldata$numero_patient  # ordo
coldata[,1:6]<-NULL # sup des data non utile

mat<-t(mat)  # translation de la matrice
all(rownames(coldata) == colnames(mat)) #verif col et row identique
all(rownames(coldata) %in% colnames(mat)) #attention de verif pas l'ordre mais seulement présence
mat_log<-log1p(mat) # passage en log 1 p pour les heatmap
annC <- data.frame(condition= coldata) # pour les annotation en col pour la heatmap
rownames(annC) <- colnames(mat_log)  # pour les annotation en col pour la heatmap
all(rownames(annC) == colnames(mat)) #vérif du même ordre / position
annC$real_time_point <- mat_orig$real_time_point  # ajout real_time_point

pheatmap(mat_log,
         scale="row",
         fontsize_row=11,
         fontsize_col = 12,
         annotation_col = annC,
         annotation_colors = list(REPONSE = c(NR = "#7570BE",  
                                              R="#F15854",
                                              RP="#882255",
                                              "T"= "#117733"), 
                                  real_time_point = c(VT1="deepskyblue",
                                                      VT2="gold3",
                                                      "T"="#117733")),
         color = my_palette, 
         cutree_rows = 1, 
         main = "2 gènes avec VT1+VT2", 
         labels_col = coldata_num,
         cluster_cols = F,
         cluster_rows = T,)

## 5.4-18+65 genes VT1+VT2- -----
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] #DF avce tout les VT1
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] #DF avce tout les VT1
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat<-rbind(mat_VT1,mat_VT2)
mat<-rbind(mat,mat_T)
mat_orig<-mat

name_unique_cov_T_gene<-as.data.frame(name_unique_cov_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_unique_cov_T_gene$name_unique_cov_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_18<- mat[,test==T]   #ne garde en colone que les True


name_inter_cov_RP_T_gene<-as.data.frame(name_inter_cov_RP_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_inter_cov_RP_T_gene$name_inter_cov_RP_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_65<- mat[,test==T]   #ne garde en colone que les True

mat<-cbind(mat_18,mat_65)

coldata<-mat_orig[,737:743] # creation du coldata pour les heatmap
coldata_num<- coldata$numero_patient  # ordo
coldata[,1:6]<-NULL # sup des data non utile

mat<-t(mat)  # translation de la matrice
all(rownames(coldata) == colnames(mat)) #verif col et row identique
all(rownames(coldata) %in% colnames(mat)) #attention de verif pas l'ordre mais seulement présence
mat_log<-log1p(mat) # passage en log 1 p pour les heatmap
annC <- data.frame(condition= coldata) # pour les annotation en col pour la heatmap
rownames(annC) <- colnames(mat_log)  # pour les annotation en col pour la heatmap
all(rownames(annC) == colnames(mat)) #vérif du même ordre / position
annC$real_time_point <- mat_orig$real_time_point  # ajout real_time_point

pheatmap(mat_log,
         scale="row",
         fontsize_row=10,
         fontsize_col = 12,
         annotation_col = annC,
         annotation_colors = list(REPONSE = c(NR = "#7570BE",  
                                              R="#F15854",
                                              RP="#882255",
                                              "T"= "#117733"), 
                                  real_time_point = c(VT1="deepskyblue",
                                                      VT2="gold3",
                                                      "T"="#117733")),
         color = my_palette, 
         cutree_rows = 1, 
         main = "Gènes de la réponse covid à VT1 avec clusterisation", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,)

pheatmap(mat_log,
         scale="row",
         fontsize_row=10,
         fontsize_col = 12,
         annotation_col = annC,
         annotation_colors = list(REPONSE = c(NR = "#7570BE",  
                                              R="#F15854",
                                              RP="#882255",
                                              "T"= "#117733"), 
                                  real_time_point = c(VT1="deepskyblue",
                                                      VT2="gold3",
                                                      "T"="#117733")),
         color = my_palette, 
         cutree_rows = 1, 
         main = "Gènes de la réponse covid à VT1", 
         labels_col = coldata_num,
         cluster_cols = F,
         cluster_rows = T,)


## 5.5-2+18+65 genes VT1+VT2- -----
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] #DF avce tout les VT1
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),] #DF avce tout les VT1
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T

mat<-rbind(mat_VT1,mat_VT2)
mat<-rbind(mat,mat_T)
mat_orig<-mat


name_unique_cov_T_gene<-as.data.frame(name_unique_cov_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_unique_cov_T_gene$name_unique_cov_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_18<- mat[,test==T]   #ne garde en colone que les True


name_inter_cov_RP_T_gene<-as.data.frame(name_inter_cov_RP_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_inter_cov_RP_T_gene$name_inter_cov_RP_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_65<- mat[,test==T]   #ne garde en colone que les True

name_unique_RP_T_gene<-as.data.frame(name_unique_RP_T_gene)    ##creation d'un DF
test <- colnames(mat) %in% name_unique_RP_T_gene$name_unique_RP_T_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_2<- mat[,test==T]   #ne garde en colone que les True

mat<-cbind(mat_18,mat_65)
mat<-cbind(mat,mat_2)

coldata<-mat_orig[,737:743] # creation du coldata pour les heatmap
coldata_num<- coldata$numero_patient  # ordo
coldata[,1:6]<-NULL # sup des data non utile

mat<-t(mat)  # translation de la matrice
all(rownames(coldata) == colnames(mat)) #verif col et row identique
all(rownames(coldata) %in% colnames(mat)) #attention de verif pas l'ordre mais seulement présence
mat_log<-log1p(mat) # passage en log 1 p pour les heatmap
annC <- data.frame(condition= coldata) # pour les annotation en col pour la heatmap
rownames(annC) <- colnames(mat_log)  # pour les annotation en col pour la heatmap
all(rownames(annC) == colnames(mat)) #vérif du même ordre / position
annC$real_time_point <- mat_orig$real_time_point  # ajout real_time_point

pheatmap(mat_log,
         scale="row",
         fontsize_row=10,
         fontsize_col = 12,
         annotation_col = annC,
         annotation_colors = list(REPONSE = c(NR = "#7570BE",  
                                              R="#F15854",
                                              RP="#882255",
                                              "T"= "#117733"), 
                                  real_time_point = c(VT1="deepskyblue",
                                                      VT2="gold3",
                                                      "T"="#117733")),
         color = my_palette, 
         cutree_rows = 1, 
         main = "2 + 18 + 65 gènes avec VT1+VT2", 
         labels_col = coldata_num,
         cluster_cols = F,
         cluster_rows = T,)

