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

# 2-ouverture des fichier ------------------------------------------------------
rm(list = ls())
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]

# 3-PCA SANS REA & T------------------------------------------------------------
MataData_PCA<-mat_pat_clean_sans_R_T[,737:743]
mat_PCA<-mat_pat_clean_sans_R_T[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             # select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous les gènes et tous les prélèvements") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 3.1-Extraction --------------------------------------------------------------
PC1_VT <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT) <- rownames(SelectPCA)
all(rownames(PC1_VT) == rownames(MataData_PCA))
PC1_VT <- cbind(PC1_VT, MataData_PCA)

PC2_VT <- as.data.frame(SelectPCA$PC2)
row.names(PC2_VT) <- rownames(SelectPCA)
all(rownames(PC2_VT) == rownames(MataData_PCA))
PC2_VT <- cbind(PC2_VT, MataData_PCA)

### 3.1.1-Visualisation---------------------------------------------------------
ggplot(PC1_VT, aes(x = jours_prelevement, y = `SelectPCA$PC1` , color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  # geom_boxplot() +
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.5, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PC1 issu de la PCA de tous les gènes et prélèvement en fonction du jours de prélèvement") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

ggplot(PC1_VT, aes(x = jours_prelevement, y = `SelectPCA$PC1`)) +
  geom_boxplot() +
  geom_point(aes(color = REPONSE, shape = real_time_point), 
             size = 3) +
  geom_text(aes(color = REPONSE),
            label = label,
            nudge_x=0.5,
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB"))
  
ggplot(PC2_VT, aes(x = jours_prelevement, y = `SelectPCA$PC2` , color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.5, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PC2 issu de la PCA de tous les gènes et prélèvement en fonction du jours de prélèvement") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

# 4-PCA VT1---------------------------------------------------------------------
mat <- mat_pat_clean_sans_R_T$real_time_point %in% "VT1"
mat_VT1 <- mat_pat_clean_sans_R_T[mat == T,]

MataData_PCA<-mat_VT1[,737:743]
mat_PCA<-mat_VT1[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             # select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous genes à VT1") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 4.1-Extraction --------------------------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)

PC2_VT1 <- as.data.frame(SelectPCA$PC2)
row.names(PC2_VT1) <- rownames(SelectPCA)
all(rownames(PC2_VT1) == rownames(MataData_PCA))
PC2_VT1 <- cbind(PC2_VT1, MataData_PCA)

### 4.1.1-Visualisation---------------------------------------------------------
ggplot(PC1_VT1, aes(x = jours_prelevement, y = `SelectPCA$PC1` , color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous genes à VT2") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

ggplot(PC2_VT1, aes(x = jours_prelevement, y = `SelectPCA$PC2` , color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous genes à VT2") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))



# 5-PCA VT2---------------------------------------------------------------------
mat <- mat_pat_clean_sans_R_T$real_time_point %in% "VT2"
mat_VT2 <- mat_pat_clean_sans_R_T[mat == T,]

MataData_PCA<-mat_VT2[,737:743]
mat_PCA<-mat_VT2[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             # select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous genes à VT2") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 5.1-Extraction PC------------------------------------------------------------
PC1_VT2 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT2) <- rownames(SelectPCA)
all(rownames(PC1_VT2) == rownames(MataData_PCA))
PC1_VT2 <- cbind(PC1_VT2, MataData_PCA)

PC2_VT2 <- as.data.frame(SelectPCA$PC2)
row.names(PC2_VT2) <- rownames(SelectPCA)
all(rownames(PC2_VT2) == rownames(MataData_PCA))
PC2_VT2 <- cbind(PC2_VT2, MataData_PCA)

### 5.1.1-Visualisation---------------------------------------------------------
ggplot(PC1_VT2, aes(x = jours_prelevement, y = `SelectPCA$PC1` , color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous genes à VT2") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

ggplot(PC2_VT2, aes(x = jours_prelevement, y = `SelectPCA$PC2` , color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape() +
  geom_text(label = label,
            nudge_x=0.8, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous genes à VT2") +
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))








         


