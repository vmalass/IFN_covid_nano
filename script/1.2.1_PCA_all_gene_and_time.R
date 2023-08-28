PCA tous les gènes et times points 

# 1-library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
# load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
gene_DE_VT1_T <- read.table("data/gene_DE_VT1_T_all_gene_mat1.2.txt")


Metadata <- mat_pat_clean[737:743]
ma<-mat_pat_clean$real_time_point %in% "REA"
mat<- mat_pat_clean[ma==F,]

MataData_PCA<-mat[,737:743]
mat_PCA<-mat[,1:736]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 6.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC avec tout les times points sur tous les genes", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 15),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = real_time_point, shape = real_time_point)) + # , shape = REPONSE, label = label
  geom_point(size = 1.5 ) + 
  scale_shape()+
  scale_color_manual(name = "Time point",
                     breaks = c("VT1","VT2","VT3","VT4", "T"),
                     values = c("#d1ab75","#758bd1", "#8fd175" ,"#BBBBBB", "#CC3311"))+
  scale_shape_manual(name = "Time point",
                     breaks = c("VT1","VT2","VT3","VT4", "T"),
                     values = c(19,19, 19, 19, 17))+
  labs(title="PCA :  all genes and all time-points",
       col = "Time point",
       x = "PC1 : 27,8%",
       y = "PC2 : 14,9%") +
  theme_classic() +
    theme(legend.position = c(0.9,0.2))  # export PDF 5x7.5 inche



