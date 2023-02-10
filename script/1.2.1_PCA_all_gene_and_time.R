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

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = real_time_point, shape = REPONSE, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf, 
                  size = 5)	+
  scale_color_manual(breaks = c("VT1","VT2","VT3","VT4", "T"),
                     values = c("darkorange","cornflowerblue", "#DDCC77" ,"brown3", "chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous les prelevements et de tous les genes",
       x = "PC1 : 27,8%",
       y = "PC2 : 14,9%") + 
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))