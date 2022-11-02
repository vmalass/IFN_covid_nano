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
load("data/1.1_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]

# 3-PCA SANS REA & T----
MataData_PCA<-mat_pat_clean_sans_R_T[,737:742]
mat_PCA<-mat_pat_clean_sans_R_T[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "PCA sans REA et témoin", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  scale_color_manual(breaks = c("V1","V2","V3","V4"),
                     values = c("cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point")

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = real_time_point)) +
  geom_point() + 
  scale_color_manual(breaks = c("VT1","VT2","VT3","VT4"),
                     values = c("cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA real time point")

# 3-PCA AVEC REA & T----
MataData_PCA<-mat_pat_clean[,737:742]
mat_PCA<-mat_pat_clean[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "PCA avec REA et témoin", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  scale_color_manual(breaks = c("REA","T","V1","V2","V3","V4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point")

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = real_time_point)) +
  geom_point() + 
  scale_color_manual(breaks = c("REA","T","VT1","VT2","VT3","VT4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  # scale_color_brewer(palette="Dark2")+
  labs(title="PCA real time point")

# 4- Etude par real time point en visualisation PCA-----

## 4.1-PCA AVEC VT1-----
mat_VT1<-mat_pat_clean$real_time_point %in% "VT1"
mat_VT1<-mat_pat_clean[mat_VT1==T,]

MataData_PCA <- mat_VT1[,737:742]
mat_PCA<- mat_VT1[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA VT1", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","V1","V2","V3","V4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point VT1")


## 4.2-PCA AVEC VT2-----
mat_VT2<-mat_pat_clean$real_time_point %in% "VT2"
mat_VT2<-mat_pat_clean[mat_VT2==T,]

MataData_PCA <- mat_VT2[,737:742]
mat_PCA<- mat_VT2[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA VT2", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","V1","V2","V3","V4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point VT2")

## 4.3-PCA AVEC VT3-----
mat_VT3<-mat_pat_clean$real_time_point %in% "VT3"
mat_VT3<-mat_pat_clean[mat_VT3==T,]

MataData_PCA <- mat_VT3[,737:742]
mat_PCA<- mat_VT3[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA VT3", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","V1","V2","V3","V4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point VT3")

## 4.4-PCA AVEC VT4-----
mat_VT4<-mat_pat_clean$real_time_point %in% "VT4"
mat_VT4<-mat_pat_clean[mat_VT4==T,]

MataData_PCA <- mat_VT4[,737:742]
mat_PCA<- mat_VT4[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA VT4", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","V1","V2","V3","V4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point VT4")

# 5- Etude par time point en visualisation PCA-----

## 5.1-PCA AVEC V1-----
mat_VT1<-mat_pat_clean$time_point %in% "V1"
mat_VT1<-mat_pat_clean[mat_VT1==T,]

MataData_PCA <- mat_VT1[,737:742]
mat_PCA<- mat_VT1[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA V1", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = real_time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","VT1","VT2","VT3","VT4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point V1")


# 5-PCA AVEC V2-----
mat_VT2<-mat_pat_clean$time_point %in% "V2"
mat_VT2<-mat_pat_clean[mat_VT2==T,]

MataData_PCA <- mat_VT2[,737:742]
mat_PCA<- mat_VT2[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA V2", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = real_time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","VT1","VT2","VT3","VT4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point V2")

# 6-PCA AVEC V3-----
mat_VT3<-mat_pat_clean$time_point %in% "V3"
mat_VT3<-mat_pat_clean[mat_VT3==T,]

MataData_PCA <- mat_VT3[,737:742]
mat_PCA<- mat_VT3[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA V3", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = real_time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","VT1","VT2","VT3","VT4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point V3")

# 7-PCA AVEC V4-----
mat_VT4<-mat_pat_clean$time_point %in% "V4"
mat_VT4<-mat_pat_clean[mat_VT4==T,]

MataData_PCA <- mat_VT4[,737:742]
mat_PCA<- mat_VT4[,1:736]
mat_PCA<-scale(mat_PCA)
PCA<- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main="PCA V4", addlabels = TRUE)
SelectPCA<- as.data.frame(PCA$x)
mergeMetaAPC<- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = real_time_point)) +
  geom_point() + 
  geom_text(label = label,
            nudge_x=1.5, nudge_y=0.1,
            check_overlap=F)+
  scale_color_manual(breaks = c("REA","T","VT1","VT2","VT3","VT4"),
                     values = c("darkorange","chartreuse4","cornflowerblue","darkorchid2","black","brown3"))+
  labs(title="PCA time point V4")





# 8-Figure papier master----
mat_T<-mat_pat_clean$real_time_point %in% "T"
mat_T<-mat_pat_clean[mat_T==T,]
mat_pat_clean_sans_R<- rbind(mat_T, mat_pat_clean_sans_R_T)
mat_pat_clean_sans_R<-rename(mat_pat_clean_sans_R, Groupe = real_time_point)
MataData_PCA<-mat_pat_clean_sans_R[,737:742]
mat_PCA<-mat_pat_clean_sans_R[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "PCA sans REA et témoin", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = time_point)) +
  geom_point() + 
  scale_color_manual(breaks = c("T","V1","V2","V3","V4"),
                     values = c("#059748","cornflowerblue","darkorchid2","#909495","brown3"))+
  labs(title="PCA time point")

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, col = Groupe)) +
  geom_point() + 
  scale_color_manual(breaks = c("T","VT1","VT2","VT3","VT4"),
                     values = c("chartreuse4","cornflowerblue","darkorchid2","#909495","brown3"))+
  labs(x="PC1 : 27,9%", y="PC2 : 15%")+
  theme(legend.text = element_text(size = 13),legend.title = element_text(size = 15)) 

