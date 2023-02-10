# 1-library --------------------------------------------------------------------
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

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
# load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE T vs VT1 all genes ------------------------------------------------------
## 3.1-Creation de la matrice---------------------------------------------------
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

coldata_mat<-as.data.frame(mat_VT1$real_time_point)
row.names(coldata_mat)<-row.names(mat_VT1)
coldata_num<-mat_VT1$numero_patient
mat_VT1[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene--------------------------------------------------
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT1), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

##### TvsVT1
dds_mat$condition <- relevel(dds_mat$condition, ref = "T")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "VT1", "T"))

res_VT1_T <- lfcShrink(dds_mat, coef = "condition_VT1_vs_T", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

VT1_T <- res_VT1_T[res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange > 1 | res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
VT1_T_all_gene <- as.data.frame(VT1_T)

data_VT1_T <- assay(rld)[res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange > 1 | res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT1_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT1_DE) <- colnames(data_VT1_T)

name_VT1_T_all_gene <- rownames(VT1_T_all_gene)

# 6-heatmap VT1 T_R all gene----------------------------------------------------
## 6.1-creation de la matrice---------------------------------------------------
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

## 6.2-Gene DE avec liste name_VT1_T_all_gene------------------------------------
mat_meta<-mat_VT1[,737:743]
coldata_VT1<-as.data.frame(mat_meta$REPONSE)
colnames(coldata_VT1) <- "REPONSE"
row.names(coldata_VT1)<-row.names(mat_meta)

coldata_num<- mat_VT1$numero_patient

name_VT1_T_all_gene<-as.data.frame(name_VT1_T_all_gene)    ##creation d'un DF
mat <- colnames(mat_VT1) %in% name_VT1_T_all_gene$name_VT1_T_all_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_VT1_T_all<- mat_VT1[,mat==T]   #ne garde en colone que les True

mat_DE_VT1_T_all<-log1p(t(mat_DE_VT1_T_all)) # passage en log 1 p pour les heatmap
annC_VT1 <- data.frame(condition= coldata_VT1) # pour les annotation en col pour la heatmap
all(rownames(annC_VT1) == colnames(mat_DE_VT1_T_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_VT1_T_all, 
         scale="row", 
         fontsize_row=10, 
         fontsize_col = 15, 
         fontsize = 15, 
         # annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
         #                                                R="#F15854",
         #                                                RP = "#882255",
         #                                                "T"= "#117733" )),
         annotation_colors = list(REPONSE = c(NR = "#7570BE",
                                              `NR-` = "darkorange",
                                              R = "cornflowerblue",
                                              `RP-` = "#882255" ,
                                              RP = "brown3",
                                              `T` = "chartreuse4",
                                              A = "#BBBBBB")),
         annotation_col = annC_VT1,
         color = my_palette, 
         cutree_rows = 2,
         cutree_cols = 4,
         main = "DE à partir de tous les gènes avec VT1vsT en VT1", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,
         border_color = "gray50",
         angle_col = 315)


# 6-heatmap VT2 T_R all gene----------------------------------------------------
## 6.1-creation de la matrice---------------------------------------------------
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

## 6.2-Gene DE avec liste name_VT1_T_all_gene------------------------------------
mat_meta<-mat_VT1[,737:743]
coldata_VT1<-as.data.frame(mat_meta$REPONSE)
colnames(coldata_VT1) <- "REPONSE"
row.names(coldata_VT1)<-row.names(mat_meta)

coldata_num<- mat_VT1$numero_patient

name_VT1_T_all_gene<-as.data.frame(name_VT1_T_all_gene)    ##creation d'un DF
mat <- colnames(mat_VT1) %in% name_VT1_T_all_gene$name_VT1_T_all_gene  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_DE_VT1_T_all<- mat_VT1[,mat==T]   #ne garde en colone que les True

mat_DE_VT1_T_all<-log1p(t(mat_DE_VT1_T_all)) # passage en log 1 p pour les heatmap
annC_VT1 <- data.frame(condition= coldata_VT1) # pour les annotation en col pour la heatmap
all(rownames(annC_VT1) == colnames(mat_DE_VT1_T_all))  # pour les annotation en col pour la heatmap

pheatmap(mat_DE_VT1_T_all, 
         scale="row", 
         fontsize_row=10, 
         fontsize_col = 15, 
         fontsize = 15, 
         # annotation_colors = list(mat_meta.REPONSE = c( NR = "#7570BE",
         #                                                R="#F15854",
         #                                                RP = "#882255",
         #                                                "T"= "#117733" )),
         annotation_colors = list(REPONSE = c(NR = "#7570BE",
                                              `NR-` = "darkorange",
                                              R = "cornflowerblue",
                                              `RP-` = "#882255" ,
                                              RP = "brown3",
                                              `T` = "chartreuse4",
                                              A = "#BBBBBB")),
         annotation_col = annC_VT1,
         color = my_palette, 
         cutree_rows = 2,
         cutree_cols = 2,
         main = "DE à partir de tous les gènes avec VT1vsT en VT2", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,
         border_color = "gray50",
         angle_col = 315)



# 9- Save data------------------------------------------------------------------
name_VT1_T_all_gene
write.table(name_VT1_T_all_gene, file = "data/gene_DE_VT1_T_all_gene_mat1.3.txt")




mat1.2 <- read.table("data/gene_DE_VT1_T_all_gene_mat1.2.txt")
mat1.3 <- read.table("data/gene_DE_VT1_T_all_gene_mat1.3.txt")

setdiff(mat1.2,mat1.3)
setdiff(mat1.3,mat1.2)
intersect(mat1.3,mat1.2)



