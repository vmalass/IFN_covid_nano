PCA à partir des gènes DE TvsR

# 1-library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
gene_DE_T_R <- read.table("data/gene_DE_T_R_all_gene_VT1_mat1.2.txt")
Metadata <- mat_pat_clean[737:743]

unique_DE_T_R <- setdiff(gene_DE_T_R$x, gene_DE_VT1_T$name_VT1_T_all_gene) #gènes unique à DE sur les gènes HVG
unique_DE_VT1_T <- setdiff(gene_DE_VT1_T$name_VT1_T_all_gene,gene_DE_T_R$x) #gènes unique à DE sur tous les gènes
intersect_DE <- intersect(gene_DE_T_R$x, gene_DE_VT1_T$name_VT1_T_all_gene) #gènes commun entre les deux DE

# 3-creation du data------------------------------------------------------------
mat_DE_VT1<- mat_pat_clean[,match(gene_DE_T_R$x, colnames(mat_pat_clean))] 
mat_DE_VT1 <- cbind(Metadata, mat_DE_VT1)

# 4-PCA VT1 gene DE TvsR-------------------------------------------------------
## 4.1-Dataframe----------------------------------------------------------------
mat_VT1 <- mat_DE_VT1[which(mat_DE_VT1$real_time_point %in% "VT1"),]
mat <- rbind(mat_DE_VT1[which(mat_DE_VT1$real_time_point %in% "T"),], mat_VT1)
mat<- arrange(mat,numero_patient) #ordo par n° pat

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:74]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 4.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC à VT1 sur le gene DE TvsR", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 9),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point(size = 3) + 
  geom_text(label = label,
            nudge_x=0.4, 
            nudge_y=0.2,
            check_overlap=F,
            size = 6)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des prelevement VT1 et des genes DE TvsR") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 4.3-Extraction PC1 VT1-------------------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
write.table(PC1_VT1, "data/PC1_VT1_gene_DE_TvsR_mat1.2.txt")

ggplot(PC1_VT1, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE)) +
  geom_point(size = 3) + 
  geom_text(label = label,
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des prelevement VT1 et des genes DE TvsR en fonction des jours de prélèvement",
       x = "jours après début des symptômes",
       y = "PC1") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))




