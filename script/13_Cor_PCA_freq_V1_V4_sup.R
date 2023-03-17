# 1 Library---------------------------------------------------------------------
library(readxl)
library(corrplot)
library(factoextra)
library(ggrepel)
library(dplyr)

# 2 import data-----------------------------------------------------------------
rm(list = ls())
set.seed(123)

Freq_V1_V4_R_RP <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx", 3)

Freq_V1_V4_R_RP <- arrange(Freq_V1_V4_R_RP, REPONSE)



col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Separation Groupe-------------------------------------------------------------
T_F <- Freq_V1_V4_R_RP$REPONSE %in% "R"
data_R <- Freq_V1_V4_R_RP[T_F == T,]
data_RP <- Freq_V1_V4_R_RP[T_F == F,]


# Groupe R RP-------------------------------------------------------------------
Meta_data <- Freq_V1_V4_R_RP[,1:7]
Freq_V1_V4_R_RP[,1:7] <- NULL

# Freq_V1 <- Freq_V1_V4_R_RP[,1:49]
# Freq_V4 <- Freq_V1_V4_R_RP[,50:104]
# 
# data_V1 <- scale(Freq_V1)  # normalisation
# data_V4 <- scale(Freq_V4)  # normalisation
# 
# data <- cbind(data_V1, data_V4)

### Variable ###
data <- scale(Freq_V1_V4_R_RP)  # normalisation

data_cor<- cor(data, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 0.4, #reduit la taille de l'ecriture des variables
         number.cex = 0.2, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
         )



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


## Cor avec R ------------------------------------------------------------------
Meta_data_R <- data_R[,1:7]
data_R[,1:7] <- NULL

### Variable ###
data <- scale(data_R)  # normalisation

data_cor<- cor(data, 
               method = "pearson")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 0.4, #reduit la taille de l'ecriture des variables
         number.cex = 0.2, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
         )


## Cor avec RP ------------------------------------------------------------------
Meta_data_RP <- data_RP[,1:7]
data_RP[,1:7] <- NULL

### Variable ###
data <- scale(data_RP)  # normalisation

data_cor<- cor(data, 
               method = "spearman")

corrplot(data_cor, method="color", 
         col=col(200),  
         tl.cex = 0.4, #reduit la taille de l'ecriture des variables
         number.cex = 0.2, #reduit la taille de l'ecriture des chiffres
         type="lower",  
         addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         diag=T # Cacher les coefficients de corrélation sur la diagonale
         )



















