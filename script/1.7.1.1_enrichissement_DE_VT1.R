rm(list = ls())

name_DE_NR_R_HVG_all<- read.table("data/name_DE_NR_R_HVG_all.txt")
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
library(dplyr)



# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Hs.eg.db, 
                        keys=name_DE_NR_R_HVG_all[,1], 
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene          = geneListEntrez,
                   organism      = "hsa",
                   pAdjustMethod = "BH", 
                   minGSSize     = 3,
                   pvalueCutoff  = 0.05, # slider
                   qvalueCutoff  = 0.05) # slider

# Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=20, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
dotplot(Kres, showCategory=20) + ggtitle("Dotplot of functional enrichment by KEGG")

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Hs.eg.db, 
                        keys=name_DE_NR_R_HVG_all[,1], 
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH", 
  minGSSize     = 3,
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05, # slider
  readable      = TRUE)

# enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=20, order=T) + ggtitle("Barplot of functional enrichment by GO")
dotplot(enrichedRes_GO, showCategory=20) + ggtitle("Dotplot of functional enrichment by GO")

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Hs.eg.db, 
                        keys=name_DE_NR_R_HVG_all[,1], 
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GS <- enrichPathway(gene          = geneListEntrez_PA,
                                organism      = 'human',
                                pAdjustMethod = "BH", 
                                minGSSize     = 3,
                                pvalueCutoff  = 0.05, # slider
                                qvalueCutoff  = 0.05) # slider

barplot(enrichedRes_GS, showCategory=20, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
dotplot(enrichedRes_GS, showCategory=20) + ggtitle("Dotplot of functional enrichment by REACTOM_PA")

# GSEA--------------------------------------------------------------------------


# Gene set IFN------------------------------------------------------------------
par(mfrow=c(2,1))
### REACTOME_INTERFERON_SIGNALING ###
a <- read_tsv('data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.tsv')
members <- a[19,2]
members <- str_split(members, ',')
members <- as.data.frame(members[[1]])


inter <- intersect(members$`members[[1]]`,colnames(mat_pat_clean[1:736]))
mat_IFN <- mat_pat_clean[,match(inter, colnames(mat_pat_clean))] 
mat_IFN <- cbind(Metadata, mat_IFN)

## PCA VT1 gene IFN-------------------------------------------------------------

ma<-mat_IFN$real_time_point %in% "VT1"
mat<- mat_IFN[ma==T,]
mata<-mat$real_time_point %in% "REA"
mat<-mat[mata==F,]
mat_R<-mat$real_time_point %in% "T"
mat<-mat[mat_R==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:78]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point() + 
  scale_color_manual(breaks = c("NR","R","RP"),
                     values = c("darkorange","chartreuse4","cornflowerblue"))+
  labs(title="PCA sur VT1 gène IFN")

## PCA VT2 gene IFN-------------------------------------------------------------
ma<-mat_IFN$real_time_point %in% "VT2"
mat<- mat_IFN[ma==T,]
mata<-mat$real_time_point %in% "REA"
mat<-mat[mata==F,]
mat_R<-mat$real_time_point %in% "T"
mat<-mat[mat_R==F,]

MataData_PCA<-mat[,1:7]
mat_PCA<-mat[,8:78]
mat_PCA<-scale(mat_PCA)
PCA <- prcomp(mat_PCA, scale. = F)
fviz_eig(PCA, main = "Variance par PC", addlabels = TRUE)

SelectPCA <- as.data.frame(PCA$x)
mergeMetaAPC <- as.data.frame(cbind(SelectPCA, MataData_PCA))

ggplot(mergeMetaAPC, aes(x = PC1, y = PC2, color = REPONSE)) +
  geom_point() + 
  scale_color_manual(breaks = c("NR","R","RP"),
                     values = c("darkorange","chartreuse4","cornflowerblue"))+
  labs(title="PCA sur VT2 gène IFN")

### WP_HOSTPATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_INTERFERON_INDUCTION ###
a <- read_tsv('data/WP_HOSTPATHOGEN_INTERACTION_OF_HUMAN_CORONAVIRUSES_INTERFERON_INDUCTION.v2022.1.Hs.tsv')
members <- a[19,2]
members <- str_split(members, ',')
members <- as.data.frame(members[[1]])


inter <- intersect(members$`members[[1]]`,colnames(mat_pat_clean[1:736]))
mat_IFN <- mat_pat_clean[,match(inter, colnames(mat_pat_clean))] 
mat_IFN <- cbind(Metadata, mat_IFN)

## PCA VT1 gene IFN-------------------------------------------------------------

ma<-mat_IFN$real_time_point %in% "VT1"
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
  scale_color_manual(breaks = c("NR","R","RP"),
                     values = c("darkorange","chartreuse4","cornflowerblue"))+
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
  scale_color_manual(breaks = c("NR","R","RP"),
                     values = c("darkorange","chartreuse4","cornflowerblue"))+
  labs(title="PCA sur VT2 gène IFN")







