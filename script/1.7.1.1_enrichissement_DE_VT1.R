rm(list = ls())

name_DE_NR_R_HVG_all_VT1<- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/name_DE_NR_R_HVG_all_VT1.txt")
name_DE_R_RP_HVG_all_VT2<- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/name_DE_R_RP_HVG_all_VT2.txt")
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
#https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/

inteferon_gmt <- read.gmt("data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.gmt")

#dowload http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2022.1.Hs/c2.all.v2022.1.Hs.symbols.gmt (create count)
msigdb_C2_gmt <- read.gmt("data/c2.all.v2022.1.Hs.symbols.gmt")

enrichedRes_MSigDBC2 <- enricher(name_DE_NR_R_HVG_all_VT1$name_DE_NR_R_HVG_all,
                                 TERM2GENE = msigdb_C2_gmt,
                                 universe = colnames(mat_pat_clean[,1:736]))

enrichedRes_MSigDBC2["REACTOME_INTERFERON_SIGNALING",]
summary(enrichedRes_MSigDBC2)

barplot(enrichedRes_MSigDBC2, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by MSigDBC2")
dotplot(enrichedRes_MSigDBC2, showCategory=20) + ggtitle("Dotplot of functional enrichment by MSigDBC2")

# enrichedRes_GO_mat = as.matrix(enrichedRes_MSigDBC2[])


# geneorderlist = 153:1
# names(geneorderlist) = name_DE_NR_R_HVG_all_VT1$name_DE_NR_R_HVG_all
# 
# GSEA_MSigDBC2 = GSEA(geneorderlist,  TERM2GENE = msigdb_C2_gmt, pvalueCutoff = 1)
# 
# gseaplot(GSEA_MSigDBC2, "REACTOME_INTERFERON_SIGNALING")
# 
# GSEA_MSigDBC2_df = as.data.frame(GSEA_MSigDBC2)
# 
# GSEA_MSigDBC2_df[GSEA_MSigDBC2_df$Description == "REACTOME_INTERFERON_SIGNALING",] 




