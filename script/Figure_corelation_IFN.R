
Cor entre les PC1

library(ggrepel)
library(ggplot2)

rm(list = ls())
PC1_VT_gene_DE_VT1vsT <- read.table("data/PC1_VT_gene_DE_VT1vsT_mat1.2.txt")
PC1_VT_IFN_geneset <- read.table("data/PC1_VT_IFN_geneset_mat1.2.txt")

##################################
data_fi <- inner_join(PC1_VT_IFN_geneset, PC1_VT_gene_DE_VT1vsT, by = c("jours_prelevement" = "jours_prelevement",
                                                                        "time_point"="time_point", 
                                                                        "real_time_point"="real_time_point", 
                                                                        "condition_biologique"="condition_biologique",
                                                                        "REA"="REA",
                                                                        "numero_patient"="numero_patient",
                                                                        "REPONSE"="REPONSE")) %>%
  mutate(label_patient = as.character(numero_patient))

data_fi <- data_fi[data_fi$REPONSE != "NR",]

test_cor <- cor.test(data_fi$SelectPCA.PC1.x, data_fi$SelectPCA.PC1.y, method= "pearson")
options( "digits"=3)
ggplot(data_fi, aes(x = SelectPCA.PC1.x, y = SelectPCA.PC1.y)) + #, label = label_patient
  geom_point(aes(color = REPONSE),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=0, 
           y=22, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation",
       x = "PC1 geneset IFN reactom",
       y = "PC1 gene DE VT1vsT") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 
  
