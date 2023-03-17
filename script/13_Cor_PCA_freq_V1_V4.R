# 1 Library---------------------------------------------------------------------
library(readxl)
library(corrplot)
library(factoextra)
library(ggrepel)

# 2 import data-----------------------------------------------------------------
rm(list = ls())
set.seed(123)

Freq_V1_V4_R_RP <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# 3 Groupe----------------------------------------------------------------------
# nom <- Freq_V1_V4_R_RP$numero_patient_victor
# rownames(Freq_V1_V4_R_RP) <- c("1",  "6",  "7",  "8", "11", "25", "12", "58",
#                                "14", "16", "63", "17",  "4", "48", "21", "10",
#                                "15")

T_F <- Freq_V1_V4_R_RP$Reponse %in% "R"
data_R <- Freq_V1_V4_R_RP[T_F == T,]
data_RP <- Freq_V1_V4_R_RP[T_F == F,]

nom_R <- data_R$numero_patient_victor
rownames(data_R) <- c("1",  "6",  "7",  "8", "11", "25", "12", "58", "14", "16",
                      "63", "17")
nom_RP <- data_RP$numero_patient_victor
rownames(data_RP) <- c("4", "48", "21", "10", "15")

# 3 Correlation-----------------------------------------------------------------
# 3.1 Correlation groupe R-------------------------------------------------------
Meta_data_R <- data_R[,1:5]
data_R[,1:5] <- NULL

### Variable ###
data <- scale(data_R)  # normalisation

data_cor<- cor(data, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 0.6, #reduit la taille de l'ecriture des variables
         number.cex = 0.2, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
         
)


a <- data_cor > 0.5 | data_cor < -0.5
b <- apply(a, 1 , sum) >= 1
data_cor[!a] <- NA

write.csv(data_cor, file = "data/cor_V1_V4_R.csv")

### individus ###
rownames(data) <- c("1",  "6",  "7",  "8", "11", "25", "12", "58", "14", "16", "63", "17")

data_t<-t(data)  # transposition pour cor en fct des lignes (patients)

data_cor<- cor(data_t, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 2, #reduit la taille de l'ecriture des variables
         number.cex = 1, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
)     



# 3.2 Correlation groupe RP-------------------------------------------------------
Meta_data_RP <- data_RP[,1:5]
data_RP[,1:5] <- NULL

### Variable ###
data <- scale(data_RP)  # normalisation

data_cor<- cor(data, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 0.6, #reduit la taille de l'ecriture des variables
         number.cex = 0.2, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
         
)


a <- data_cor > 0.5 | data_cor < -0.5
b <- apply(a, 1 , sum) >= 1
data_cor[!a] <- NA

write.csv(data_cor, file = "result/cor_V1_V4_RP.csv")

### individus ###
rownames(data) <- c("4", "48", "21", "10", "15")

data_t<-t(data)  # transposition pour cor en fct des lignes (patients)

data_cor<- cor(data_t, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 2, #reduit la taille de l'ecriture des variables
         number.cex = 1, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
)     



# 3.3 Correlation groupe R&RP---------------------------------------------------
Meta_data <- Freq_V1_V4_R_RP[,1:5]
Freq_V1_V4_R_RP[,1:5] <- NULL

### Variable ###
data <- scale(Freq_V1_V4_R_RP)  # normalisation

data_cor<- cor(data, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 0.6, #reduit la taille de l'ecriture des variables
         number.cex = 0.2, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
         
)

a <- data_cor > 0.5 | data_cor < -0.5
b <- apply(a, 1 , sum) >= 1
data_cor[!a] <- NA

write.csv(data_cor, file = "data/cor_V1_V4_R_RP.csv")

### individus ###
rownames(data) <- c("1",  "6",  "7",  "8", "11", "25", "12", "58", "14", "16", "63", "17", "4", "48", "21", "10", "15")

data_t<-t(data)  # transposition pour cor en fct des lignes (patients)

data_cor<- cor(data_t, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 2, #reduit la taille de l'ecriture des variables
         number.cex = 1, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
)     











# 4 PCA-------------------------------------------------------------------------
## 4.1 Frequence V1+V4----------------------------------------------------------
mat_PCA <- scale(Freq_V1_V4_R_RP)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, Meta_data))
label<- as.character(mergeMetaAPC$numero_patient_victor)

### Visu ###
fviz_eig(PCA, main = "Variance cluster V1 et V4", addlabels = TRUE)

### PC1 & PC2 ###
fviz_pca_var(PCA, 
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()


fviz_contrib(PCA, choice = "var", axes = 1) + coord_flip()+
  labs(title="Comtribution des patients dans la PC1 (19.3%)",
       x = "Cluster") +
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

fviz_contrib(PCA, choice = "var", axes = 2) + coord_flip()+
  labs(title="Comtribution des patients dans la PC2 (12.5%)",
       x = "Cluster") +
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

fviz_contrib(PCA, choice = "var", axes = c(1:2)) + coord_flip()+
  labs(title="Comtribution des patients dans la PC1 & 2 (31.8%)",
       x = "Cluster") +
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = Responder, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = 5,
                  size = 5)	+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("#7570BE", "darkorange","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des clusters non supervisé de V1 et V4",
       x = "PC1 = 19.3%",
       y = "PC2 = 12.5%") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


### PC3 & PC4 ###
fviz_pca_var(PCA, 
             axes = c(3,4),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC3 et PC4") +
  theme_minimal()

fviz_contrib(PCA, choice = "var", axes = 3) + coord_flip()+
  labs(title="Comtribution des patients dans la PC3 (12.2%)",
       x = "Cluster") +
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

fviz_contrib(PCA, choice = "var", axes = 4) + coord_flip()+
  labs(title="Comtribution des patients dans la PC4 (10.9%)",
       x = "Cluster") +
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

fviz_contrib(PCA, choice = "var", axes = c(3:4)) + coord_flip()+
  labs(title="Comtribution des patients dans la PC3 & 4 (23.8%)",
       x = "Cluster") +
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

ggplot(mergeMetaAPC, aes(x = PC3, y = PC4, color = Responder, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = 5,
                  size = 5)	+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("#7570BE", "darkorange","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des clusters non supervisé de V1 et V4",
       x = "PC3 = 12.2%",
       y = "PC4 = 10.9%") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))



## 4.2 Frequence V1-------------------------------------------------------------
Freq_V1_V4_R_RP <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx")
data_V1 <- Freq_V1_V4_R_RP[,43:75]

mat_PCA <- scale(data_V1)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, Meta_data))
label<- as.character(mergeMetaAPC$numero_patient_victor)

### Visu ###
fviz_eig(PCA, main = "Variance cluster V1", addlabels = TRUE)

### PC1 & PC2 ###
fviz_pca_var(PCA, 
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = Responder, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = 5,
                  size = 5)	+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("#7570BE", "darkorange","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des clusters non supervisé de V1",
       x = "PC1 = 24%",
       y = "PC2 = 14.2%") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))




## 4.3 Frequence V4-------------------------------------------------------------
data_V4 <- Freq_V1_V4_R_RP[,6:42]

mat_PCA <- scale(data_V4)
PCA <- prcomp(mat_PCA, scale. = F)
SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, Meta_data))
label<- as.character(mergeMetaAPC$numero_patient_victor)

### Visu ###
fviz_eig(PCA, main = "Variance cluster V4", addlabels = TRUE)

### PC1 & PC2 ###
fviz_pca_var(PCA, 
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 20),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = Responder, label = label)) +
  geom_point(size = 3) + 
  scale_shape()+
  geom_text_repel(show.legend = F,
                  max.overlaps  = 5,
                  size = 5)	+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("#7570BE", "darkorange","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA à partir des clusters non supervisé de V4",
       x = "PC1 = 21.6%",
       y = "PC2 = 19.9%") + 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

