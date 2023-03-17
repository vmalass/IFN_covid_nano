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
library(cowplot)


# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg
T_F <- mat_pat_clean$REPONSE %in% "REA"
mat_pat_clean_sans_REA <- mat_pat_clean[T_F == F,]

T_F <- mat_pat_clean_sans_REA$REPONSE %in% "A"
mat_pat_clean_sans_REA <- mat_pat_clean_sans_REA[T_F == F,]

T_F <- mat_pat_clean_sans_REA$REPONSE %in% "Other"
mat_pat_clean_sans_REA <- mat_pat_clean_sans_REA[T_F == F,]

dico_gene <- data.frame(colnames(mat_pat_clean_sans_REA))

mat_pat_clean_sans_REA$REPONSE <- as.factor(mat_pat_clean_sans_REA$REPONSE)

## B cells----------------------------------------------------------------------

gene <- c("BLK",
  "CD19",
  "FAM30A",
  "FCRL2",
  "MS4A1",
  "PNOC",
  "SPIB",
  "TCL1A",
  "TNFRSF17")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = BLK, fill = real_time_point)) +
  geom_violin()
  
b <- ggplot(mat_macro, aes(x = REPONSE, y = CD19, fill = real_time_point)) +
  geom_violin()
  
c <- ggplot(mat_macro, aes(x = REPONSE, y = FAM30A, fill = real_time_point)) +
  geom_violin()
  
d <- ggplot(mat_macro, aes(x = REPONSE, y = FCRL2, fill = real_time_point)) +
  geom_violin()

e <- ggplot(mat_macro, aes(x = REPONSE, y = MS4A1, fill = real_time_point)) +
  geom_violin()
  
f <- ggplot(mat_macro, aes(x = REPONSE, y = PNOC, fill = real_time_point)) +
  geom_violin()

g <- ggplot(mat_macro, aes(x = REPONSE, y = SPIB, fill = real_time_point)) +
  geom_violin()

h <- ggplot(mat_macro, aes(x = REPONSE, y = TCL1A, fill = real_time_point)) +
  geom_violin()

i <- ggplot(mat_macro, aes(x = REPONSE, y = TNFRSF17, fill = real_time_point)) +
  geom_violin()

plot_grid(a, b, c, d, e, f, g, h, i, 
          labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
          ncol = 3, nrow = 3)  # Plot avec 3 figures


## CD 45------------------------------------------------------------------------

gene <- c("PTPRC")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = mat_macro, fill = real_time_point)) +
  labs(y = "PTPRC")+
  geom_violin()

plot_grid(a,
          labels=c("A"), 
          ncol = 1, nrow = 1)  # Plot avec 3 figures


## CD8 T cells----------------------------------------------------------------------

gene <- c("CD8A",
          "CD8B")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CD8A, fill = real_time_point)) +
  geom_violin()

b <- ggplot(mat_macro, aes(x = REPONSE, y = CD8B, fill = real_time_point)) +
  geom_violin()


plot_grid(a, b,
          labels=c("A", "B"), 
          ncol = 2, nrow = 1)  # Plot avec 3 figures


## Cytotoxic cells----------------------------------------------------------------------

gene <- c("CTSW",
          "GNLY",
          "GZMA",
          "GZMB",
          "GZMH",
          "KLRB1",
          "KLRD1",
          "KLRK1",
          "NKG7",
          "PRF1")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CTSW, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = GNLY, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = GZMA, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

d <- ggplot(mat_macro, aes(x = REPONSE, y = GZMB, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

e <- ggplot(mat_macro, aes(x = REPONSE, y = GZMH, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

f <- ggplot(mat_macro, aes(x = REPONSE, y = KLRB1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

g <- ggplot(mat_macro, aes(x = REPONSE, y = KLRD1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

h <- ggplot(mat_macro, aes(x = REPONSE, y = KLRK1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

i <- ggplot(mat_macro, aes(x = REPONSE, y = NKG7, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

j <- ggplot(mat_macro, aes(x = REPONSE, y = PRF1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

plot_grid(a, b, c, d, e, f, g, h, i, j,
          labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), 
          ncol = 5, nrow = 2)  # Plot avec 3 figures



## DC---------------------------------------------------------------------------

gene <- c("CCL13",
          "CD209",
          "HSD11B1")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = mat_macro, fill = real_time_point)) +
  labs(y = "CCL13")+
  geom_violin()+
  theme(legend.position = "none")

plot_grid(a,
          labels=c("A"), 
          ncol = 1, nrow = 1)  # Plot avec 3 figures



## Exhausted CD8----------------------------------------------------------------------

gene <- c("CD244",
          "EOMES",
          "LAG3",
          "PTGER4")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CD244, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = EOMES, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = LAG3, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

d <- ggplot(mat_macro, aes(x = REPONSE, y = PTGER4, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b, c, d,
          labels=c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2)  # Plot avec 3 figures


## Macrophages----------------------------------------------------------------------

gene <- c("CD163",
          "CD68",
          "CD84",
          "MS4A4A")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CD163, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = CD68, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = CD84, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

d <- ggplot(mat_macro, aes(x = REPONSE, y = MS4A4A, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b, c, d,
          labels=c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2)  # Plot avec 3 figures



## Mast cells----------------------------------------------------------------------

gene <- c("CPA3",
          "HDC",
          "MS4A2",
          "TPSAB1/B2")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CPA3, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = HDC, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = MS4A2, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

d <- ggplot(mat_macro, aes(x = REPONSE, y = `TPSAB1/B2`, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b, c, d,
          labels=c("A", "B", "C", "D"), 
          ncol = 2, nrow = 2)  # Plot avec 3 figures





## Neutrophils----------------------------------------------------------------------

gene <- c("CEACAM3",
          "CSF3R",
          "FCAR",
          "FCGR3A/B",
          "FPR1",
          "S100A12",
          "SIGLEC5")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CEACAM3, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = CSF3R, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = FCAR, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

d <- ggplot(mat_macro, aes(x = REPONSE, y = `FCGR3A/B`, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

e <- ggplot(mat_macro, aes(x = REPONSE, y = FPR1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

f <- ggplot(mat_macro, aes(x = REPONSE, y = S100A12, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

g <- ggplot(mat_macro, aes(x = REPONSE, y = SIGLEC5, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

plot_grid(a, b, c, d, e, f, g,
          labels=c("A", "B", "C", "D", "E", "F", "G"), 
          ncol = 4, nrow = 2)  # Plot avec 3 figures



## NK CD56dim cells----------------------------------------------------------------------

gene <- c("IL21R",
          "KIR2DL3",
          "KIR3DL1/2")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = IL21R, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = KIR2DL3, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = `KIR3DL1/2`, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b, c,
          labels=c("A", "B", "C"), 
          ncol = 2, nrow = 2)  # Plot avec 3 figures



## NK cells----------------------------------------------------------------------

gene <- c("NCR1",
          "XCL1/2")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = NCR1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = `XCL1/2`, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b,
          labels=c("A", "B"), 
          ncol = 2, nrow = 1)  # Plot avec 3 figures




## T-cells----------------------------------------------------------------------

gene <- c("CD3D",
          "CD3E",
          "CD3G",
          "CD6",
          "SH2D1A",
          "TRAT1")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CD3D, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = CD3E, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = CD3G, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

d <- ggplot(mat_macro, aes(x = REPONSE, y = CD6, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

e <- ggplot(mat_macro, aes(x = REPONSE, y = SH2D1A, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

f <- ggplot(mat_macro, aes(x = REPONSE, y = TRAT1, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b, c, d, e, f,
          labels=c("A", "B", "C", "D", "E", "F"), 
          ncol = 3, nrow = 2)  # Plot avec 3 figures



## Th1 cells & Treg----------------------------------------------------------------------

gene <- c("TBX21",
          "FOXP3")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = TBX21, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = FOXP3, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")


plot_grid(a, b,
          labels=c("A : Th1 cells", "B : Treg"), 
          ncol = 2, nrow = 1)  # Plot avec 3 figures





## Sophia marker----------------------------------------------------------------------

gene <- c("CXCR3",
          "HLA-DR",
          "CD38",
          # "CD163",
          "SELL",
          # "S100A12",
          "CD11C")

T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,737 : 743])


a <- ggplot(mat_macro, aes(x = REPONSE, y = CXCR3, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

b <- ggplot(mat_macro, aes(x = REPONSE, y = CD38, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")

c <- ggplot(mat_macro, aes(x = REPONSE, y = SELL, fill = real_time_point)) +
  geom_violin()+
  theme(legend.position = "none")



plot_grid(a, b, c,
          labels=c("A", "B", "C"), 
          ncol = 2, nrow = 2)  # Plot avec 3 figures





