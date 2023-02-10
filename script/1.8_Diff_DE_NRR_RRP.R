rm(list = ls())
gene_DE_R_RP_HVG_all_VT2 <- read.table("data/gene_DE_R_RP_HVG_all_VT2_mat1.3.txt")
gene_DE_NR_R_HVG_all_VT1 <- read.table("data/gene_DE_NR_R_HVG_all_VT1_mat1.3.txt")
a <- setdiff(gene_DE_NR_R_HVG_all_VT1$commun_DE_all_HVG, gene_DE_R_RP_HVG_all_VT2$commun_DE_all_HVG) # unique gene_DE_NR_R_HVG_all_VT1
z <- setdiff(gene_DE_R_RP_HVG_all_VT2$commun_DE_all_HVG, gene_DE_NR_R_HVG_all_VT1$commun_DE_all_HVG) # unique gene_DE_R_RP_HVG_all_VT2
e <- intersect(gene_DE_NR_R_HVG_all_VT1$commun_DE_all_HVG, gene_DE_R_RP_HVG_all_VT2$commun_DE_all_HVG) # commun entre les deux
r <- union(gene_DE_NR_R_HVG_all_VT1, gene_DE_R_RP_HVG_all_VT2) # unique au deux 
a
z
e
r