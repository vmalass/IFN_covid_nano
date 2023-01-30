rm(list = ls())

gene_DE_NR_R_HVG_all_VT1 <- read.table("data/gene_DE_NR_R_HVG_all_VT1_bis.txt")
gene_DE_R_RP_HVG_all_VT2 <- read.table("data/gene_DE_R_RP_HVG_all_VT2_bis.txt")
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
Metadata <- mat_pat_clean[737:743]
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("data/HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('org.Hs.eg.db')) BiocManager::install('org.Hs.eg.db'); library('org.Hs.eg.db')
if (!require('AnnotationDbi')) BiocManager::install('AnnotationDbi'); library('AnnotationDbi')
if (!require('enrichplot')) BiocManager::install('enrichplot'); library('enrichplot')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler'); library('clusterProfiler')
if (!require('ReactomePA')) BiocManager::install('ReactomePA'); library('ReactomePA')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')


# Gene set IFN------------------------------------------------------------------
### REACTOME_INTERFERON_SIGNALING ###
a <- read_tsv('data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.tsv')
members <- a[19,2]
members <- str_split(members, ',')
members <- as.data.frame(members[[1]])


inter <- intersect(members$`members[[1]]`,colnames(mat_pat_clean[1:736]))
mat_IFN <- mat_pat_clean[,match(inter, colnames(mat_pat_clean))] 
mat_IFN <- cbind(Metadata, mat_IFN)
print(paste0("Sur les ", length(t(members)), " gènes du gene set IFN, ", length(inter), " sont dans le jeu de donnée"))

## PCA VT1 gene IFN-------------------------------------------------------------
ma<-mat_IFN$real_time_point %in% "VT1"
mat<- mat_IFN[ma==T,]
ma<-mat_IFN$real_time_point %in% "T"
mat_T<-mat_IFN[ma==T,]
mat <- rbind(mat_T, mat)


MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:78]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC à VT1 sur le gene set IFN", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, # Avoid text overlapping
             select.var = list(contrib = 100),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()


SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point(size = 3) + 
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
  labs(title="PCA à partir de VT1 sur le gene set IFN") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

PC1_VT1 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT1) <- rownames(SelectPCA)
all(rownames(PC1_VT1) == rownames(MataData_PCA))
PC1_VT1 <- cbind(PC1_VT1, MataData_PCA)
write.table(PC1_VT1, "/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT1_IFN_geneset.txt")


## PCA VT2 gene IFN-------------------------------------------------------------
ma<-mat_IFN$real_time_point %in% "VT2"
mat<- mat_IFN[ma==T,]
ma<-mat_IFN$real_time_point %in% "T"
mat_T<-mat_IFN[ma==T,]
mat <- rbind(mat_T, mat)

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:78]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC à VT2 sur le gene set IFN", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T, # Avoid text overlapping
             select.var = list(contrib = 100),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point(size = 3) + 
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
  labs(title="PCA à partir de VT2 sur le gene set IFN") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

PC1_VT2 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT2) <- rownames(SelectPCA)
all(rownames(PC1_VT2) == rownames(MataData_PCA))
PC1_VT2 <- cbind(PC1_VT2, MataData_PCA)
write.table(PC1_VT2, "/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT2_IFN_geneset.txt")

## PCA VT gene IFN-------------------------------------------------------------
ma<-mat_IFN$real_time_point %in% "REA"
mat<- mat_IFN[ma==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:78]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC à VT2 sur le gene set IFN", addlabels = TRUE)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T, # Avoid text overlapping
             select.var = list(contrib = 100),
             alpha.var="contrib") +
  ggtitle(label = "Contribution des variables dans les PC1 et PC2") +
  theme_minimal()

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))
label<- as.character(mergeMetaAPC$numero_patient)

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
  labs(title="PCA à partir de tous les times points sur le gene set IFN") + 
  theme(plot.title = element_text(size=20),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


PC1_VT2 <- as.data.frame(SelectPCA$PC1)
row.names(PC1_VT2) <- rownames(SelectPCA)
all(rownames(PC1_VT2) == rownames(MataData_PCA))
PC1_VT2 <- cbind(PC1_VT2, MataData_PCA)
write.table(PC1_VT2, "/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT_IFN_geneset.txt")























### WP_HOSTPATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_INTERFERON_INDUCTION ### -----
a <- read_tsv('data/WP_HOSTPATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_INTERFERON_INDUCTION.v2022.1.Hs.tsv')
members <- a[19,2]
members <- str_split(members, ',')
members <- as.data.frame(members[[1]])


inter <- intersect(members$`members[[1]]`,colnames(mat_pat_clean[1:736]))
mat_IFN <- mat_pat_clean[,match(inter, colnames(mat_pat_clean))] 
mat_IFN <- cbind(Metadata, mat_IFN)
print(paste0("Sur les ", length(t(members)), " gènes du gene set IFN, ", length(inter), " sont dans le jeu de donnée"))


## PCA VT1 gene IFN-------------------------------------------------------------

ma<-mat_IFN$real_time_point %in% "VT1"
mat<- mat_IFN[ma==T,]
mata<-mat$real_time_point %in% "REA"
mat<-mat[mata==F,]
# mat_R<-mat$real_time_point %in% "T"
# mat<-mat[mat_R==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:38]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point() + 
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA sur VT1 gène IFN")

## PCA VT2 gene IFN-------------------------------------------------------------
ma<-mat_IFN$real_time_point %in% "VT2"
mat<- mat_IFN[ma==T,]
mata<-mat$real_time_point %in% "REA"
mat<-mat[mata==F,]
mat_R<-mat$real_time_point %in% "T"
mat<-mat[mat_R==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:38]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point() + 
  # scale_color_manual(breaks = c("NR","R","RP","T"),
  #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  labs(title="PCA sur VT2 gène IFN")






