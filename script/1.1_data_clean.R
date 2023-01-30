# 1-library ------
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

# 2-Preparation des data -----
#setwd("~/Documents/JM/NanoString/NanoString_Covid/nanostring_covid/data") #folder data
rm(list = ls())
seuratObj<-readRDS(file = "data/seuratObj.rds")
MetaData_sans_Hcov <- read.xlsx("data/METADATA_sans_Hcov.xlsx", colNames = T, rowNames = T)
mat_pat_gene_count <- t(seuratObj@assays[["RNA"]]@data) #sélection des data count dans le seurat obj
mat_pat_gene_count<- as.data.frame(mat_pat_gene_count) #transformation en data frame
mat_pat_gene_count<- arrange(mat_pat_gene_count, row.names(mat_pat_gene_count)) #ordo en fonction du nom des ligne
MetaData_sans_Hcov<- arrange(MetaData_sans_Hcov, row.names(MetaData_sans_Hcov)) #ordo en fct du nom des lignes
all(row.names(MetaData_sans_Hcov) == row.names(mat_pat_gene_count))  # verif que les lignes sont dans le meme ordre pour les deux data frame
mat_pat_gene_count <- cbind(mat_pat_gene_count, MetaData_sans_Hcov) # création d'un fichier contenant toute les info

mat_pat_gene_count$numero_patient <- as.factor(mat_pat_gene_count$numero_patient) #transfo en factor 

mat_pat_gene_count$jours_prelevement <- as.numeric(as.character(mat_pat_gene_count$jours_prelevement )) #transfo en character puis numeric pour ordo les jours

mat_pat_gene_count$numero_patient <- as.numeric(as.character(mat_pat_gene_count$numero_patient ))  #transfo en character puis numeric pour ordo les n° pat
mat_pat_gene_count$numero_patient <- as.factor(mat_pat_gene_count$numero_patient) #pour transfo les num patient en factor avec lvl

mat_pat_gene_count<- arrange(mat_pat_gene_count,jours_prelevement) #ordo par jours

mat_pat_gene_count<-mat_pat_gene_count[mat_pat_gene_count$numero_patient !=28,]
mat_pat_gene_count<-mat_pat_gene_count[mat_pat_gene_count$numero_patient !=56,]
mat_pat_gene_count<-mat_pat_gene_count[mat_pat_gene_count$numero_patient !=64,]
mat_pat_gene_count<-mat_pat_gene_count[mat_pat_gene_count$numero_patient !=66,] # sup les pat 28, 56, 64, 66 qui n'ont qu'un plvmt
# mat_pat_gene_count<-mat_pat_gene_count[mat_pat_gene_count$numero_patient !=24,] # sup le pat 24 qui n'a pas de V1

mat_pat_gene_count<- arrange(mat_pat_gene_count,real_time_point) #ordo par real time point

mat_pat_gene_count_R_T<- mat_pat_gene_count[1:19,] # conserve les REA & Temoins dans un obj

# 3- Data filtre ----
mat_pat_clean<-mat_pat_gene_count

mat_pat_clean[,1:736] <- sapply(mat_pat_clean[,1:736], as.numeric) #converti les col count avec les gene en num
mat_pat_clean$numero_patient <- as.factor(mat_pat_clean$numero_patient) # transfo les num pat en factor lvl

length(unique(mat_pat_clean$numero_patient)) #nombre de patient

mat_pat_clean$jours_prelevement <- as.factor(mat_pat_clean$jours_prelevement) #transforme les jours de plvt en factor lvl
mat_pat_clean <- arrange(mat_pat_clean, jours_prelevement) #tris en fct des jours de plvt
# mat_pat_clean$real_time_point<- as.numeric(mat_pat_clean$real_time_point)
mat_pat_clean$real_time_point <- as.factor(mat_pat_clean$real_time_point)  #transforme les real time point en factor lvl
mat_pat_clean <- arrange(mat_pat_clean, real_time_point) # tris en fct des real time point

# 4- Svg fichier ----
save(mat_pat_clean, file = "data/1.1_mat_pat_clean.rds") #svg du data
load("data/1.1_mat_pat_clean.rds") #ouverture de la svg
