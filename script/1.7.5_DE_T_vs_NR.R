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
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE R vs NR en VT1 all genes ------------------------------------------------
## 3.1-Creation de la matrice---------------------------------------------------
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),] 
mat_T <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),] 
mat_VT1 <- rbind(mat_VT1,mat_T)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

coldata_mat<-as.data.frame(mat_VT1$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT1)
coldata_num<-mat_VT1$numero_patient
mat_VT1[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene--------------------------------------------------
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT1), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

##### NR
dds_mat$condition <- relevel(dds_mat$condition, ref = "NR")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "T", "NR"))

# NR vs T
res_T_NR <- lfcShrink(dds_mat, coef = "condition_T_vs_NR", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

T_NR <- res_T_NR[res_T_NR$svalue < 0.05 & !is.na(res_T_NR$svalue) & res_T_NR$log2FoldChange > 1 | res_T_NR$svalue < 0.05 & !is.na(res_T_NR$svalue) & res_T_NR$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
T_NR_all_gene <- as.data.frame(T_NR)

T_NR_all_gene  ### 1 gène DE : LILRA3 (49-50 en NR)

