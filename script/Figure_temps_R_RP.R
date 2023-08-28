PCA à partir des gènes DE VT1vsT

# 1-library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
# load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
gene_DE_VT1_T <- read.table("data/gene_DE_VT1_T_all_gene_mat1.3.txt")
Metadata <- mat_pat_clean[737:743]

# 3-creation du data------------------------------------------------------------
mat_DE_VT1 <- mat_pat_clean[,match(gene_DE_VT1_T$name_VT1_T_all_gene, colnames(mat_pat_clean))] 
mat_DE_VT1 <- cbind(Metadata, mat_DE_VT1)

# 6-PCA gene DE VT1vsT (VT1,2,3,4)-----------------------------------------------
## 6.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT1$real_time_point %in% "REA"
mat<- mat_DE_VT1[ma==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:72]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 6.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC avec tout les times points sur les genes DE VT1vsT", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 5),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE, shape = real_time_point, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = 20, 
                  size = 5)	+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("#7570BE", "darkorange","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous les prelevements et des genes DE VT1vsT") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 6.3-Extraction PC1 & visualisation-------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
# write.table(PC1_VT1, "data/PC1_VT_gene_DE_VT1vsT_mat1.3.txt")

T_F <- PC1_VT1$REPONSE %in% c("R", "RP")
mat <- PC1_VT1[T_F == T,]

ggplot(mat, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE)) +
  geom_point(size = 1.5) + 
  scale_color_manual(breaks = c("R","RP"),
                     values = c("gray60","#CB2027"))+
  labs(title="PCA all time points with gene DE VT1vsT, gene accordinf to the day after onset of symptoms",
       x = "days after onset of symptoms",
       y = "PC1 DE VT1vsT") + 
  theme_classic() +  
  theme(legend.position = c(0.85,0.8)) 
