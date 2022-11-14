rm(list = ls())

name_DE_NR_R_HVG_all<- read.table("data/name_DE_NR_R_HVG_all.txt")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('org.Hs.eg.db')) BiocManager::install('org.Hs.eg.db'); library('org.Hs.eg.db')
if (!require('AnnotationDbi')) BiocManager::install('AnnotationDbi'); library('AnnotationDbi')
if (!require('enrichplot')) BiocManager::install('enrichplot'); library('enrichplot')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler'); library('clusterProfiler')
if (!require('ReactomePA')) BiocManager::install('ReactomePA'); library('ReactomePA')



# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Hs.eg.db, 
                        keys=name_DE_NR_R_HVG_all[,1], 
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene         = geneListEntrez,
                   organism     = "hsa",
                   pvalueCutoff = 0.05)

Kresmat = as.matrix(Kres[])

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
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05,# slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=20, order=T) + ggtitle("Barplot of functional enrichment by GO")
dotplot(enrichedRes_GO, showCategory=20) + ggtitle("Dotplot of functional enrichment by GO")

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Hs.eg.db, 
                        keys=name_DE_NR_R_HVG_all[,1], 
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GS <- enrichPathway(gene = geneListEntrez_PA,
                                organism = 'human',
                                pvalueCutoff = 0.05)

barplot(enrichedRes_GS, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
dotplot(enrichedRes_GS, showCategory=20) + ggtitle("Dotplot of functional enrichment by REACTOM_PA")

# GSEA--------------------------------------------------------------------------
a <- read_tsv('data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.tsv')
members <- a[18,2]
members <- str_split(members, ',')
members <- as.data.frame(members[[1]])


x <- enricher(name_DE_NR_R_HVG_all, TERM2GENE = members)


