rm(list = ls())
gene_DE_R_RP_HVG_all_VT2 <- read.table("data/gene_DE_R_RP_HVG_all_VT2_bis.txt")
gene_DE_NR_R_HVG_all_VT1 <- read.table("data/gene_DE_NR_R_HVG_all_VT1_bis.txt")
a <- setdiff(gene_DE_NR_R_HVG_all_VT1$x, gene_DE_R_RP_HVG_all_VT2$x) # unique gene_DE_NR_R_HVG_all_VT1
z <- setdiff(gene_DE_R_RP_HVG_all_VT2$x, gene_DE_NR_R_HVG_all_VT1$x) # unique gene_DE_R_RP_HVG_all_VT2
e <- intersect(gene_DE_NR_R_HVG_all_VT1$x, gene_DE_R_RP_HVG_all_VT2$x) # commun entre les deux
r <- union(gene_DE_NR_R_HVG_all_VT1, gene_DE_R_RP_HVG_all_VT2) # unique au deux 
a
z
e
