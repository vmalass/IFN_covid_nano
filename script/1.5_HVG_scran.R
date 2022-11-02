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
load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]

# 3- HVG avec Scran----
mat_pat_clean_count<-mat_pat_clean[,1:736] #DF avec les counts
mat_pat_clean_count<-log1p(mat_pat_clean_count) #DF count scale avec log1p

HVGscran <- SingleCellExperiment(list(logcounts=t(mat_pat_clean_count)))

dec <- modelGeneVar(HVGscran,parametric=TRUE)

top.hvgs<- as.matrix(getTopHVGs(dec, var.threshold=0)) # Selection de tout les genes pos = gene significatif (rep bio)
top20scran <- head(top.hvgs, 20)

# 4-Visualisation----

plot1 <- plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance", main = " highly Variable Genes Scran")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

# 5-Svg fichier ----

save(top.hvgs, file = "data/HVG_scran.rds") # svg des HVG
load("data/HVG_scran.rds") #ouverture de la svg

