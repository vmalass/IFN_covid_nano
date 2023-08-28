Formatage pour GSEA

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
library(ggrepel)
if (!require('org.Hs.eg.db')) BiocManager::install('org.Hs.eg.db'); library('org.Hs.eg.db')
if (!require('AnnotationDbi')) BiocManager::install('AnnotationDbi'); library('AnnotationDbi')


# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg
load("data/HVG_scran.rds") #ouverture de la svg

# 3 DE T vs VT1 all genes ------------------------------------------------------
### matrix ###
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
T_F <- mat_VT1$REPONSE %in% c("Other", "NR")
mat_VT1 <- mat_VT1[T_F==F,]
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

coldata_mat<-as.data.frame(mat_VT1$real_time_point)
row.names(coldata_mat)<-row.names(mat_VT1)
meta_VT1 <- mat_VT1[,737:743]
mat_VT1[,737:743]<-NULL

### DESeq2 all gene ###
colnames(coldata_mat)<-"condition"

dds_mat <- DESeqDataSetFromMatrix(countData = t(mat_VT1), colData = coldata_mat,
                                  design = ~ condition) #creation de l'obj deseq2

dds_mat <- dds_mat[rowSums(counts(dds_mat)) >= 10,] #pre-filtrage sup les genes inf ou egale a 10

rld <- rlogTransformation(dds_mat, blind=FALSE) #transforme le dds en log pour utilisation en heatmap

### DE ###
dds_mat$condition <- relevel(dds_mat$condition, ref = "T")
dds_mat <- DESeq(dds_mat)
resultsNames(dds_mat)
res <- results(dds_mat, contrast=c("condition", "VT1", "T"))

res_VT1_T <- lfcShrink(dds_mat, coef = "condition_VT1_vs_T", type = "apeglm", lfcThreshold = 1) #resultat avec le calcule de la pval ajuste a partir de 1 et pas 0 (par rapport au threshold) : https://support.bioconductor.org/p/113664/

data_VT1_T <- assay(rld)[res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange > 1 | res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT1_DE <- data.frame(condition= coldata_mat)

all(colnames(data_VT1_T) == rownames(mat_VT1))


# 4 Formatage pour GSEA---------------------------------------------------------
data_VT1_T  #table sortie deseq2
meta_VT1    #table des meta data de data_VT1_T

t_meta_VT1 <- t(meta_VT1)
all(colnames(data_VT1_T) == colnames(t_meta_VT1))

nom_col <- colnames(data_VT1_T)

mat_ENSG = as.data.frame(mapIds(org.Hs.eg.db,
                                keys= rownames(data_VT1_T),
                                column="ENSEMBL",
                                keytype="SYMBOL",
                                multiVals="first"))

data_gsea <- cbind(mat_ENSG, data_VT1_T)

# write.xlsx(data_gsea, file = "/Users/victor/Documents/JM/NanoString/Papier/Supplementary/VT1_tab_gsea.xlsx", rowNames = T)
write.xlsx(meta_VT1, file = "/Users/victor/Documents/JM/NanoString/Papier/Supplementary/meta_tab_gsea.xlsx", rowNames = T)



# table_suple <- t(mat_DE_VT1_T_all)
# table_suple <- cbind(annC_VT1, table_suple)
# table_suple <- cbind(mat_VT1$numero_patient, table_suple)
# table_suple %>% 
#   dplyr::rename(numero_patient = `mat_VT1$numero_patient`) %>% 
#   dplyr::rename(Groups = REPONSE) -> table_suple
# table_suple %>%
#   group_by(Groups) %>%
#   summarise(nb = n())
# 
# DE_VT1_T_all_gene <- as.data.frame(res_VT1_T)
# DE_VT1_T_all_gene %>% 
#   arrange(log2FoldChange) -> DE_VT1_T_all_gene 
# 
# VT1_T_all_gene %>% 
#   arrange(log2FoldChange) -> VT1_T_all_gene

