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

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
load("data/1.6_mat_pat_clean_FIGURE_RP.rds") #ouverture de la svg
T_F <- mat_pat_clean$REPONSE %in% "REA"
mat_pat_clean_sans_REA <- mat_pat_clean[T_F == F,]

T_F <- mat_pat_clean$REPONSE %in% "T"
mat_T <- mat_pat_clean[T_F == T,]

T_F <- mat_pat_clean_sans_REA$REPONSE %in% "T"
mat_pat_clean_sans_REA <- mat_pat_clean_sans_REA[T_F == F,]

T_F <- colnames(mat_T) %in% c("IFI27", "RSAD2", "IFIT1", "CXCL10")
mat_T <- mat_T[,T_F == T]

mat_pat_clean_sans_REA[,1:736] <- log10(mat_pat_clean_sans_REA[,1:736])
mat_T<- log10(mat_T)

T_F <- colnames(mat_pat_clean_sans_REA) %in% c("IFI27", "RSAD2", "IFIT1", "CXCL10")
mat <- mat_pat_clean_sans_REA[,T_F == T]
mat <- cbind(mat, mat_pat_clean_sans_REA[737:743])
T_F <- mat$REPONSE %in% "RP"
mat <- mat[T_F == T,]
# 3-Moyenne des temoins---------------------------------------------------------
summary(mat_T)
mean_RSAD2 <- mean(mat_T$RSAD2)
mean_CXCL10 <- mean(mat_T$CXCL10)
mean_IFI27 <- mean(mat_T$IFI27)
mean_IFIT1 <- mean(mat_T$IFIT1)

# visualisation cinetique-------------------------------------------------------
ggplot(data = mat_pat_clean_sans_REA, 
       aes_string(x="jours_prelevement", 
                  y= "RSAD2", 
                  color = "REPONSE", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = mean_RSAD2, linetype="dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(breaks = c("RP", "R"),
                     values = c("#CB2027", "gray70"),
                     labels=c("persistent response", "response"))+   
  labs(title="RSAD2",
       color = "Groups",
       x = "days after onset of symptoms",
       y = "log10 normalized counts") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.85)) 


ggplot(data = mat_pat_clean_sans_REA, 
       aes_string(x="jours_prelevement", 
                  y= "IFIT1", 
                  color = "REPONSE", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = mean_IFIT1, linetype="dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(breaks = c("RP", "R"),
                     values = c("#CB2027", "gray70"),
                     labels=c("persistent response", "response"))+   
  labs(title="IFIT1",
       color = "Groups",
       x = "days after onset of symptoms",
       y = "log10 normalized counts") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.85)) 


ggplot(data = mat_pat_clean_sans_REA, 
       aes_string(x="jours_prelevement", 
                  y= "CXCL10", 
                  color = "REPONSE", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = mean_CXCL10, linetype="dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(breaks = c("RP", "R"),
                     values = c("#CB2027", "gray70"),
                     labels=c("persistent response", "response"))+   
  labs(title="CXCL10",
       color = "Groups",
       x = "days after onset of symptoms",
       y = "log10 normalized counts") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.85)) 

ggplot(data = mat_pat_clean_sans_REA, 
       aes_string(x="jours_prelevement", 
                  y= "IFI27", 
                  color = "REPONSE", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = mean_IFI27, linetype="dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(breaks = c("RP", "R"),
                     values = c("#CB2027", "gray70"),
                     labels=c("persistent response", "response"))+   
  labs(title="IFI27",
       color = "Groups",
       x = "days after onset of symptoms",
       y = "log10 normalized counts") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.85)) 

