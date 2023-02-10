# 1-library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')

# 2-ouverture des fichier-------------------------------------------------------
rm(list = ls())
# load("data/1.2_mat_pat_clean.rds") #ouverture de la svg
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
gene_DE_NR_R_HVG_all_VT1 <- read.table("data/gene_DE_NR_R_HVG_all_VT1_mat1.3.txt")
Metadata <- mat_pat_clean[737:743]

# 3-creation du data------------------------------------------------------------
mat_DE_VT1<- mat_pat_clean[,match(gene_DE_NR_R_HVG_all_VT1$commun_DE_all_HVG, colnames(mat_pat_clean))] 
mat_DE_VT1 <- cbind(Metadata, mat_DE_VT1)

# 4-PCA VT1 gene DE NRvsR-------------------------------------------------------
## 4.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT1$real_time_point %in% "VT1"
mat<- mat_DE_VT1[ma==T,]
ma<-mat_DE_VT1$real_time_point %in% "T"
mat_T<-mat_DE_VT1[ma==T,]
mat <- rbind(mat_T, mat)

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:79]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 4.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC à VT1 sur le gene DE VT1 NRvsR", addlabels = TRUE)
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
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des prelevement VT1 et des genes DE VT1 NRvsR") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 4.3-Extraction PC1 VT1-------------------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
write.table(PC1_VT1, "data/PC1_VT1_gene_DE_NRvsR_mat1.3.txt")

ggplot(PC1_VT1, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE)) +
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
  labs(title="PCA à partir des prelevement VT1 et des genes DE VT1 NRvsR en fonction des jours de prélèvement") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


# 5-PCA VT2 gene DE NRvsR-------------------------------------------------------
## 5.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT1$real_time_point %in% "VT2"
mat<- mat_DE_VT1[ma==T,]
ma<-mat_DE_VT1$real_time_point %in% "T"
mat_T<-mat_DE_VT1[ma==T,]
mat <- rbind(mat_T, mat)

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:79]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 5.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC à VT2 sur le gene DE VT1 NRvsR", addlabels = TRUE)
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
  labs(title="PCA à partir des prelevement VT2 et des genes DE VT1 NRvsR") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 5.3-Extraction PC1 VT2-------------------------------------------------------
PC1_VT2 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT2) <- rownames(SelectPCA)
all(rownames(PC1_VT2) == rownames(MataData_PCA))
PC1_VT2 <- cbind(PC1_VT2, MataData_PCA)

ggplot(PC1_VT2, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE)) +
  geom_point(size = 3) + 
  geom_text(label = label,
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  labs(title="PCA à partir des prelevement VT2 et des genes DE VT1 NRvsR en fonction du jours de prélèvement") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


# 6-PCA gene DE NRvsR (VT1,2,3,4)-----------------------------------------------
## 6.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT1$real_time_point %in% "REA"
mat<- mat_DE_VT1[ma==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:79]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 6.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC avec tout les times points sur lse genes DE VT1 NRvsR", addlabels = TRUE)
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
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir de tous les prelevements et des genes DE VT1 NRvsR") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 6.3-Extraction PC1 & visualisation-------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
write.table(PC1_VT1, "data/PC1_VT_gene_DE_NRvsR_mat1.3.txt")

ggplot(PC1_VT1, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  # geom_boxplot() +
  scale_shape()+
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
  labs(title="PC1 issu des genes DE VT1 NRvsR en fonction du jour de prélèvement") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

### Boxplot ###
ggplot(PC1_VT1, aes(x = jours_prelevement, y = `SelectPCA$PC1`)) +
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

# 7-PCA gene DE NRvsR (V1)------------------------------------------------------
## 7.1-Dataframe----------------------------------------------------------------
ma<-mat_DE_VT1$time_point %in% "V1"
mat<- mat_DE_VT1[ma==T,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:79]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

## 7.2-Visualisation------------------------------------------------------------
fviz_eig(PCA, main = "Variance par PC avec tout les times points sur lse genes DE VT1 NRvsR", addlabels = TRUE)
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
            nudge_x=0.2, 
            nudge_y=0.15,
            check_overlap=F,
            size = 6)+
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des prélèvement V1 et des gènes DE NRvsR (en VT1)") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 6.3-Extraction PC1 & visualisation-------------------------------------------
PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
write.table(PC1_VT1, "data/PC1_V1_gene_DE_NRvsR_mat1.3.txt")

ggplot(PC1_VT1, aes(x =jours_prelevement , y = `SelectPCA$PC1`, color = REPONSE, shape = real_time_point)) +
  geom_point(size = 3) + 
  geom_boxplot() +
  scale_shape()+
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
  labs(title="PC1 issu des genes DE VT1 NRvsR en fonction du jour de prélèvement") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

### Boxplot ###
ggplot(PC1_VT1, aes(x = jours_prelevement, y = `SelectPCA$PC1`)) +
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
