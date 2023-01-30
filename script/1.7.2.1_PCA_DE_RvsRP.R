# 1-library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
# load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
gene_DE_R_RP_HVG_all_VT2 <- read.table("data/gene_DE_R_RP_HVG_all_VT2_bis.txt")
Metadata <- mat_pat_clean[737:743]

# 3-creation du data------------------------------------------------------------
mat_DE_VT2<- mat_pat_clean[,match(gene_DE_R_RP_HVG_all_VT2$commun_DE_all_HVG, colnames(mat_pat_clean))] 
mat_DE_VT2 <- cbind(Metadata, mat_DE_VT2)

# 4-PCA VT2 gene DE RvsRP-------------------------------------------------------
## 4.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT2$real_time_point %in% "VT2"
mat<- mat_DE_VT2[ma==T,]
ma<-mat_DE_VT2$real_time_point %in% "T"
mat_T<-mat_DE_VT2[ma==T,]
mat <- rbind(mat_T, mat)

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:79]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 4.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC à VT2 sur le gene DE VT2 RvsRP", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 5),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point(size = 3) + 
  geom_text(label = label,
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des prelevement VT2 et des genes DE VT2 RvsRP") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 4.3-Extraction PC1 VT2-------------------------------------------------------
PC1_VT2 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT2) <- rownames(SelectPCA)
all(rownames(PC1_VT2) == rownames(MataData_PCA))
PC1_VT2 <- cbind(PC1_VT2, MataData_PCA)
write.table(PC1_VT2, "data/PC1_VT2_gene_DE_RvsRP.txt")

ggplot(PC1_VT2, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE)) +
  geom_point(size = 3) + 
  geom_text(label = label,
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des prelevement VT2 et des genes DE VT2 RvsRP en fonction du jours de prélèvement") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

# 5-PCA gene DE RPvsR-----------------------------------------------------------
## 5.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT2$real_time_point %in% "REA"
mat<- mat_DE_VT2[ma==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:79]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 5.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC avec tous les times points sur les genes DE VT2 RvsRP", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 5),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text(label = label,
            nudge_x=0.5, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous les times points sur les genes DE VT2 RvsRP") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 5.3-Extraction PC1 & visualisation-------------------------------------------
PC1_VT2 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT2) <- rownames(SelectPCA)
all(rownames(PC1_VT2) == rownames(MataData_PCA))
PC1_VT2 <- cbind(PC1_VT2, MataData_PCA)
write.table(PC1_VT2, "data/PC1_VT_gene_DE_RvsRP.txt")

ggplot(PC1_VT2, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  # geom_boxplot() +
  scale_shape()+
  geom_text(label = label,
            nudge_x=0.5, 
            nudge_y=0.3,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PC1 issu des genes DE VT2 RvsRP en fonction du jour de prélèvement") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

### Boxplot ###
ggplot(PC1_VT2, aes(x = jours_prelevement, y = `SelectPCA$PC1`)) +
  geom_boxplot() +
  geom_point(aes(color = REPONSE, shape = real_time_point), 
             size = 3) +
  geom_text(aes(color = REPONSE),
            label = label,
            nudge_x=0.5,
            nudge_y=0.3,
            check_overlap=F,
            size = 6) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB"))
  