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
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg
# load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)


# 3 HEATMAP-----------------------------------------------------------------------
## 3.1 DE T vs VT1 all genes ------------------------------------------------------
### matrix ###
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

coldata_mat<-as.data.frame(mat_VT1$real_time_point)
row.names(coldata_mat)<-row.names(mat_VT1)
coldata_num<-mat_VT1$numero_patient
mat_VT1[,737:744]<-NULL

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

VT1_T <- res_VT1_T[res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange > 1 | res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange < -1 , ]   #tris des genes avec sval<5% et L2FC <-1 & >1
VT1_T_all_gene <- as.data.frame(VT1_T)

data_VT1_T <- assay(rld)[res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange > 1 | res_VT1_T$svalue < 0.05 & !is.na(res_VT1_T$svalue) & res_VT1_T$log2FoldChange < -1 , ]    #tris des genes avec sval<5% et L2FC <-1 & >1 + utilisation de assay(rld) pour passer en log

annC_VT1_DE <- data.frame(condition= coldata_mat)
rownames(annC_VT1_DE) <- colnames(data_VT1_T)

name_VT1_T_all_gene <- rownames(VT1_T_all_gene)

## 3.2 VT1 Visualisation------------------------------------------------------------
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT1"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

T_F <- mat_VT1$REPONSE %in% "Other"
mat_VT1 <- mat_VT1[T_F==F,]

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
         annotation_colors = list(REPONSE = c(NR = "#7570BE",
                                              R = "#F15854",
                                              RP = "#882255",
                                              `T` = "chartreuse4")),  #, Other = "#BBBBBB"
         annotation_col = annC_VT1,
         color = my_palette, 
         cutree_rows = 2,
         cutree_cols = 2,
         main = "DE from all genes with VT1vsT to VT1", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,
         border_color = "gray50",
         angle_col = 315)


## 3.3 VT2 Visualisation------------------------------------------------------------
mat_VT1 <- mat_pat_clean[which(mat_pat_clean$real_time_point %in% "VT2"),]
mat_VT1 <- rbind(mat_pat_clean[which(mat_pat_clean$real_time_point %in% "T"),], mat_VT1)
mat_VT1<- arrange(mat_VT1,numero_patient) #ordo par n° pat

T_F <- mat_VT1$REPONSE %in% "Other"
mat_VT1 <- mat_VT1[T_F==F,]

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
         annotation_colors = list(REPONSE = c(NR = "#7570BE",
                                              R = "#F15854",
                                              RP = "#882255",
                                              `T` = "chartreuse4")),  # , Other = "#BBBBBB"
         annotation_col = annC_VT1,
         color = my_palette, 
         cutree_rows = 2,
         cutree_cols = 2,
         main = "DE from all genes with VT1vsT to VT2", 
         labels_col = coldata_num,
         cluster_cols = T,
         cluster_rows = T,
         border_color = "gray50",
         angle_col = 315)








# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg
gene_DE_VT1_T <- read.table("data/gene_DE_VT1_T_all_gene_mat1.3.txt")
Metadata <- mat_pat_clean[737:743]

mat_DE_VT1 <- mat_pat_clean[,match(gene_DE_VT1_T$name_VT1_T_all_gene, colnames(mat_pat_clean))] 
mat_DE_VT1 <- cbind(Metadata, mat_DE_VT1)



# 4 PCA DE VT1vsT---------------------------------------------------------------
### matrix ###
ma<-mat_DE_VT1$real_time_point %in% "REA"
mat<- mat_DE_VT1[ma==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:72]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))


ma <- mergeMetaAPC$REPONSE %in% "Other"
mergeMetaAPC <- mergeMetaAPC[ma == F,]

label<- as.character(mergeMetaAPC$numero_patient)
## 4.1 Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC avec tout les times points sur les genes DE VT1vsT", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 5),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE, shape = real_time_point, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = 20, 
                  size = 5)	+
  scale_color_manual(breaks = c("R", "RP", "T", "Other"),
                     values = c("cornflowerblue", "brown3", "chartreuse4", "white")) +
  labs(title="PCA à partir de tous les prelevements et des genes DE VT1vsT") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme_bw()

## 4.1 Extraction PC1 & visualisation-------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
# write.table(PC1_VT1, "data/PC1_VT_gene_DE_VT1vsT_mat1.3.txt")

ma <- PC1_VT1$REPONSE %in% "Other"
PC1_VT1 <- PC1_VT1[ma == F,]

ggplot(PC1_VT1, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE, shape = real_time_point, label = label)) +
  geom_point(size = 3) + 
  # geom_boxplot() +
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf, 
                  size = 5)	+
  scale_color_manual(breaks = c("R", "RP", "T", "Other"),
                     values = c("cornflowerblue", "brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous les prelevements et des genes DE VT1vsT en fonction des jours de prélèvement",
       x = "jours après début des symptômes",
       y = "PC1") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15)) +
  theme_bw()




















# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg
T_F <- mat_pat_clean$REPONSE %in% "REA"
mat_pat_clean_sans_REA <- mat_pat_clean[T_F == F,]

T_F <- mat_pat_clean_sans_REA$REPONSE %in% "A"
mat_pat_clean_sans_REA <- mat_pat_clean_sans_REA[T_F == F,]
a <- data.frame(colnames(mat_pat_clean_sans_REA))


T_F <- colnames(mat_pat_clean_sans_REA) %in% c("CD163", "CD68", "CD84", "MS4A4A")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]


ggplot(mat_macro, aes(x = REPONSE, y = CD163, fill = real_time_point)) +
  geom_violin()+

ggplot(mat_macro, aes(x = REPONSE, y = CD68, fill = real_time_point)) +
  geom_violin()+

ggplot(mat_macro, aes(x = REPONSE, y = CD84, fill = real_time_point)) +
  geom_violin()+

ggplot(mat_macro, aes(x = REPONSE, y = MS4A4A, fill = real_time_point)) +
  geom_violin()




T_F <- colnames(mat_pat_clean_sans_REA) %in% c("CD14", "CX3CR1", "CD38", "CD86")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = CD14, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = CX3CR1, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = CD38, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = CD86, fill = real_time_point)) +
  geom_violin()




T_F <- colnames(mat_pat_clean_sans_REA) %in% c("CSF1R", "FAS", "HLA-DRA", "HLA-DRB")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = CSF1R, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = FAS, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = `HLA-DRA`, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = `HLA-DRB`, fill = real_time_point)) +
  geom_violin()



T_F <- colnames(mat_pat_clean_sans_REA) %in% c("GZMA", "GZMB", "GZMH", "XCR1")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = GZMA, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = GZMB, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = GZMH, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = XCR1, fill = real_time_point)) +
  geom_violin()


T_F <- colnames(mat_pat_clean_sans_REA) %in% c("XCL1/2", "IL1B", "IL1R1", "IL1R2")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = `XCL1/2`, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = IL1B, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = IL1R1, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = IL1R2, fill = real_time_point)) +
  geom_violin()


T_F <- colnames(mat_pat_clean_sans_REA) %in% c("TLR1", "TLR2", "TLR3", "TLR4")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = TLR1, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = TLR2, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = TLR3, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = TLR4, fill = real_time_point)) +
  geom_violin()


T_F <- colnames(mat_pat_clean_sans_REA) %in% c("TLR5", "TLR6", "TLR7", "TLR8", "TLR9")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = TLR5, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = TLR6, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = TLR7, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = TLR8, fill = real_time_point)) +
  geom_violin()


ggplot(mat_macro, aes(x = REPONSE, y = TLR9, fill = real_time_point)) +
  geom_violin()




T_F <- colnames(mat_pat_clean_sans_REA) %in% c("CCL13", "CCL2")
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = CCL13, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = CCL2, fill = real_time_point)) +
  geom_violin()








# libraries:
library(ggplot2)
library(gganimate)

mat_macro <- arrange(mat_macro, real_time_point)

a <- ggplot(mat_macro, aes(x=REPONSE, y=CD163, fill=REPONSE)) + 
  geom_bar(stat='identity') +
  theme_bw() +
  transition_states(
    real_time_point,
    transition_length = 2,
    state_length = 1) +
  ease_aes('sine-in-out')

animate(a, duration = 5, fps = 20, width = 200, height = 200, renderer = gifski_renderer())
anim_save("output.gif")

