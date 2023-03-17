# 1 Library---------------------------------------------------------------------
library(readxl)
library(corrplot)
library(factoextra)
library(ggrepel)
library(dplyr)
library(Hmisc)

# SUPERVISED ###################################################################
# 2 import data-----------------------------------------------------------------
rm(list = ls())
set.seed(123)

Freq_V1_V4_R_RP <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx", 4)

Freq_V1_V4_R_RP <- arrange(Freq_V1_V4_R_RP, REPONSE)



col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

# Separation Groupe-------------------------------------------------------------
T_F <- Freq_V1_V4_R_RP$REPONSE %in% "R"
data_R <- Freq_V1_V4_R_RP[T_F == T,]
data_RP <- Freq_V1_V4_R_RP[T_F == F,]


## Cor avec R ------------------------------------------------------------------
Meta_data_R <- data_R[,1:7]
data_R[,1:7] <- NULL
data_R_V1 <- data_R[,1:36]
data_R_V4 <- data_R[,37:81]

### V1--------------------------------------------------------------------------
### Variable ###
data <- scale(data_R_V1)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = F)] <- 0 ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = F)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

# dat <- dat[-19,-19]   # supression des var non corrélé
# dat_pval <- dat_pval[-19,-19]  # supression des var non corrélé


corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.7, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1"
)



### V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_R_V4)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

dat[upper.tri(dat, diag = T)] <- 0  ## passer en matrice triangle en remplacent le reste par des 0
dat_pval[upper.tri(dat_pval, diag = T)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.7, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V4"
)

### V1&V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_R)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = T)] <- 0  ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = T)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne


corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.7, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1 et V4"
)

dat <- dat[-c(1:36),-c(37:81)]
dat_pval <- dat_pval[-c(1:36),-c(37:81)]

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.8, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         # type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=90, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1 et V4"
)




## Cor avec RP ------------------------------------------------------------------
Meta_data_RP <- data_RP[,1:7]
data_RP[,1:7] <- NULL
data_RP_V1 <- data_RP[,1:36]
data_RP_V4 <- data_RP[,37:81]

### V1--------------------------------------------------------------------------
### Variable ###
data <- scale(data_RP_V1)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = F)] <- 0  ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = F)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(6,15,29),-c(6,15,29)]
dat_pval <- dat_pval[-c(6,15,29),-c(6,15,29)]

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne


corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.7, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe RP à V1"
)


### V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_RP_V4)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(37),-c(37)]
dat_pval <- dat_pval[-c(37),-c(37)]

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne


corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.7, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe RP à V4"
)


### V1&V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_RP)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = T)] <- 0  ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = T)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(6,15,29), -c(6,15,29)]
dat_pval <- dat_pval[-c(6,15,29), -c(6,15,29)]

## "Early_NK|count" & "CD3+_TCRgd+|count" -1 et 1 = 0... 
colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.7, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe RP à V1 et V4"
)


dat <- dat[-c(1:36),-c(37:81)]
dat_pval <- dat_pval[-c(1:36),-c(37:81)]

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.8, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         # type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=90, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1 et V4"
)








# UNSUPERVISED #################################################################
# 2 import data-----------------------------------------------------------------
rm(list = ls())
set.seed(123)

Freq_V1_V4_R_RP <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx", 1)

Freq_V1_V4_R_RP <- arrange(Freq_V1_V4_R_RP, Reponse)

col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

# Separation Groupe-------------------------------------------------------------
T_F <- Freq_V1_V4_R_RP$Reponse %in% "R"
data_R <- Freq_V1_V4_R_RP[T_F == T,]
data_RP <- Freq_V1_V4_R_RP[T_F == F,]

## Cor avec R ------------------------------------------------------------------
Meta_data_R <- data_R[,1:5]
data_R[,1:5] <- NULL
data_R_V1 <- data_R[,38:70]
data_R_V4 <- data_R[,1:37]

### V1--------------------------------------------------------------------------
### Variable ###
data <- scale(data_R_V1)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = F)] <- 0 ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = F)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(10,16),-c(10,16)]   # supression des var non corrélé
dat_pval <- dat_pval[-c(10,16),-c(10,16)]  # supression des var non corrélé

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 1, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1"
)



### V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_R_V4)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(14,19),-c(14,19)]   # supression des var non corrélé
dat_pval <- dat_pval[-c(14,19),-c(14,19)]  # supression des var non corrélé

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 1, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V4"
)

### V1&V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_R)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-19,-19]   # supression des var non corrélé
dat_pval <- dat_pval[-19,-19]  # supression des var non corrélé

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.9, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1 et V4"
)

dat <- dat[-c(1:36),-c(37:69)]
dat_pval <- dat_pval[-c(1:36),-c(37:69)]

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 1, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         # type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=90, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1 et V4"
)
         

## Cor avec RP ------------------------------------------------------------------
Meta_data_RP <- data_RP[,1:5]
data_RP[,1:5] <- NULL
data_RP_V1 <- data_RP[,38:70]
data_RP_V4 <- data_RP[,1:37]

### V1--------------------------------------------------------------------------
### Variable ###
data <- scale(data_RP_V1)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = F)] <- 0  ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = F)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(11,17,30),-c(11,17,30)]
dat_pval <- dat_pval[-c(11,17,30),-c(11,17,30)]

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne


corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 1, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe RP à V1"
)


### V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_RP_V4)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(7,8,17,28),-c(7,8,17,28)]
dat_pval <- dat_pval[-c(7,8,17,28),-c(7,8,17,28)]

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne


corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 1, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe RP à V4"
)


### V1&V4--------------------------------------------------------------------------
### Variable ###
data <- scale(data_RP)  # normalisation

p_cor <- rcorr(data, type = "spearman")

dat <- as.matrix(p_cor$r)
dat_pval <- as.matrix(p_cor$P)

b <- dat_pval < 0.05
dat[!b] <- 0

# dat[upper.tri(dat, diag = T)] <- 0  ## passer en matrice triangle en remplacent le reste par des 0
# dat_pval[upper.tri(dat_pval, diag = T)] <- 0

colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

dat <- dat[-c(7,54), -c(7,54)]
dat_pval <- dat_pval[-c(7,54), -c(7,54)]

## "Early_NK|count" & "CD3+_TCRgd+|count" -1 et 1 = 0... 
colnames(dat)[which(colSums(dat)==1)]  # a supprimé par colonne
rownames(dat)[which(rowSums(dat)==1)]  # a supprimé par ligne

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 0.9, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=45, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe RP à V1 et V4"
)


dat <- dat[-c(1:36),-c(37:68)]
dat_pval <- dat_pval[-c(1:36),-c(37:68)]

corrplot(dat,
         method = 'square', 
         col=col(200),  
         tl.cex = 1, #reduit la taille de l'ecriture des variables
         # number.cex = 0.4, #reduit la taille de l'ecriture des chiffres
         # type="lower",
         # addCoef.col = "black",# Ajout du coefficient de corrélation
         tl.col="black", 
         tl.srt=90, #Rotation des etiquettes de textes
         # addrect = 2,
         diag=F, # Cacher les coefficients de corrélation sur la diagonale
         # order = "hclust",
         p.mat = dat_pval,
         sig.level = c(0.001, 0.01, 0.05),
         pch.cex = 0.5,
         pch.col = 'green', 
         insig = 'label_sig',
         title = "Correlation groupe R à V1 et V4"
)



