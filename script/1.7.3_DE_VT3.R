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
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE R vs RP en VT3 all genes -----
## 3.1-Creation de la matrice----
mat_VT3 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT3"),] 
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] #DF avec tout les T
mat_VT3<- rbind(mat_VT3, mat_T)

mat_VT3<- arrange(mat_VT3,numero_patient) #ordo par n° pat


coldata_mat<-as.data.frame(mat_VT3$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT3)
coldata_num<-mat_VT3$numero_patient
mat_VT3[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene----
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT3), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

##### RP
dds_mat$condition <- relevel(dds_mat$condition, ref = "RP")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "R", "RP"))

# RP vs R
res_R_RP <- lfcShrink(dds_mat, coef = "condition_R_vs_RP", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

R_RP <- res_R_RP[res_R_RP$svalue < 0.05 & !is.na(res_R_RP$svalue) & res_R_RP$log2FoldChange > 1 | res_R_RP$svalue < 0.05 & !is.na(res_R_RP$svalue) & res_R_RP$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
R_RP_all_gene <- as.data.frame(NR_R)

data_R_RP_VT3 <- assay(rld)[res_R_RP$svalue < 0.05 & !is.na(res_R_RP$svalue) & res_R_RP$log2FoldChange > 1 | res_R_RP$svalue < 0.05 & !is.na(res_R_RP$svalue) & res_R_RP$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT3_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT3_DE) <- colnames(data_R_RP_VT3)

heatmap_R_RP_VT3 <- pheatmap(data_R_RP_VT3, 
                             scale="row", 
                             fontsize = 15,
                             fontsize_row=15, 
                             fontsize_col = 15,
                             annotation_col = annC_VT3_DE,
                             annotation_colors = list(condition = c( NR = "#7570BE",
                                                                     R="#F15854",
                                                                     RP = "#882255",
                                                                     "T"= "#117733" )),
                             color = my_palette, 
                             cutree_rows = 2,
                             main = "DE à partir de tout les gènes et des prélèvements VT3 avec RvsRP", 
                             labels_col = coldata_num,
                             cluster_cols = T,
                             cluster_rows = T,) 





