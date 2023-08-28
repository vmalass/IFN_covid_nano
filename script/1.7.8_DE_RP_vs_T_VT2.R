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

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

# 3-DE T vs RP all genes -------------------------------------------------------
## 3.1-Creation de la matrice---------------------------------------------------
mat_VT2 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),]
mat_VT2 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT2)
mat_VT2<- arrange(mat_VT2,numero_patient) #ordo par n° pat

coldata_mat<-as.data.frame(mat_VT2$REPONSE)
row.names(coldata_mat)<-row.names(mat_VT2)
coldata_num<-mat_VT2$numero_patient
mat_VT2[,737:744]<-NULL

## 3.2-Gene DE DESeq2 all gene--------------------------------------------------
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT2), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

##### TvsVT2
dds_mat$condition <- relevel(dds_mat$condition, ref = "T")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "RP", "T"))

res_VT2_T <- lfcShrink(dds_mat, coef = "condition_RP_vs_T", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

VT2_T <- res_VT2_T[res_VT2_T$svalue < 0.05 & !is.na(res_VT2_T$svalue) & res_VT2_T$log2FoldChange > 1 | res_VT2_T$svalue < 0.05 & !is.na(res_VT2_T$svalue) & res_VT2_T$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
VT2_T_all_gene <- as.data.frame(VT2_T)

data_VT2_T <- assay(rld)[res_VT2_T$svalue < 0.05 & !is.na(res_VT2_T$svalue) & res_VT2_T$log2FoldChange > 1 | res_VT2_T$svalue < 0.05 & !is.na(res_VT2_T$svalue) & res_VT2_T$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT2_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT2_DE) <- colnames(data_VT2_T)

name_VT2_T_all_gene <- rownames(VT2_T_all_gene)

name_VT2_T_all_gene

name_VT1_T_all_gene <- read.table("data/gene_DE_VT1_T_all_gene_mat1.3.txt")

name_VT1_T_all_gene$name_VT1_T_all_gene

intersect(name_VT2_T_all_gene, name_VT1_T_all_gene$name_VT1_T_all_gene)
setdiff(name_VT2_T_all_gene, name_VT1_T_all_gene$name_VT1_T_all_gene)
setdiff(name_VT1_T_all_gene$name_VT1_T_all_gene, name_VT2_T_all_gene)

a <- as.data.frame(colnames(mat_VT2))

T_F <- colnames(mat_pat_clean) %in% "FCGR1A/B"
mat_macro <- mat_pat_clean[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean[,736 : 743])

T_F <- mat_macro$REPONSE %in% c("R", "RP", "T")
mat_macro <- mat_macro[T_F == T, ]

ggplot(mat_macro, aes(x = REPONSE, y = mat_macro, fill = real_time_point)) +
  geom_violin()+
  labs(title = "FCGR1A/B",
       y = "count",
       x= "")


