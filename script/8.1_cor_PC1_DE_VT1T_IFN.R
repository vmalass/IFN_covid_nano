Cor entre les PC1

library(ggrepel)
library(ggplot2)

rm(list = ls())
gene_DE_T_R <- read.table("data/PC1_VT1_gene_DE_TvsR_mat1.2.txt")
gene_DE_VT1_T <- read.table("data/PC1_VT1_gene_DE_VT1vsT_mat1.2.txt")
PC1_VT1_IFN_geneset <- read.table("data/PC1_VT1_IFN_geneset_mat1.2.txt")

gene_DE_VT1_T_VT2 <- read.table("data/PC1_VT2_gene_DE_VT1vsT_mat1.2.txt")
PC1_VT2_IFN_geneset <- read.table("data/PC1_VT2_IFN_geneset_mat1.2.txt")

PC1_VT_gene_DE_VT1vsT <- read.table("data/PC1_VT_gene_DE_VT1vsT_mat1.2.txt")
PC1_VT_IFN_geneset <- read.table("data/PC1_VT_IFN_geneset_mat1.2.txt")


##################################
data_fi <- inner_join(gene_DE_T_R, gene_DE_VT1_T, by = c("jours_prelevement" = "jours_prelevement",
                                                                         "time_point"="time_point", 
                                                                         "real_time_point"="real_time_point", 
                                                                         "condition_biologique"="condition_biologique",
                                                                         "REA"="REA",
                                                                         "numero_patient"="numero_patient",
                                                                         "REPONSE"="REPONSE")) %>%
  mutate(label_patient = as.character(numero_patient))

test_cor <- cor.test(data_fi$SelectPCA.PC1.x, data_fi$SelectPCA.PC1.y, method= "pearson")
options( "digits"=3)
ggplot(data_fi, aes(x = SelectPCA.PC1.x, y = SelectPCA.PC1.y, label = label_patient)) + 
  geom_point(aes(color = REPONSE, shape = real_time_point), 
             size = 3)+
  geom_text_repel(# nudge_x=0.6,
    # nudge_y=0.15,
    show.legend = F,
    max.overlaps  = Inf)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=-10, 
           y=-5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  labs(title = "PC1_gene_DE_T_R / PC1_gene_DE_VT1_T",
       x = "PC1_gene_DE_T_R",
       y = "PC1_gene_DE_VT1_T") + 
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size=24),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))



##################################
data_fi <- inner_join(PC1_VT1_IFN_geneset, gene_DE_VT1_T, by = c("jours_prelevement" = "jours_prelevement",
                                                                 "time_point"="time_point", 
                                                                 "real_time_point"="real_time_point", 
                                                                 "condition_biologique"="condition_biologique",
                                                                 "REA"="REA",
                                                                 "numero_patient"="numero_patient",
                                                                 "REPONSE"="REPONSE")) %>%
  mutate(label_patient = as.character(numero_patient))

test_cor <- cor.test(data_fi$SelectPCA.PC1.x, data_fi$SelectPCA.PC1.y, method= "pearson")
options( "digits"=3)
ggplot(data_fi, aes(x = SelectPCA.PC1.x, y = SelectPCA.PC1.y, label = label_patient)) + 
  geom_point(aes(color = REPONSE, shape = real_time_point), 
             size = 3)+
  geom_text_repel(# nudge_x=0.6,
    # nudge_y=0.15,
    show.legend = F,
    max.overlaps  = Inf)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=5, 
           y=-5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  labs(title = "PC1_VT1_IFN_geneset / PC1_gene_DE_VT1_T",
       x = "PC1_VT1_IFN_geneset",
       y = "PC1_gene_DE_VT1_T") + 
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size=24),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))


##################################
data_fi <- inner_join(PC1_VT2_IFN_geneset, gene_DE_VT1_T_VT2, by = c("jours_prelevement" = "jours_prelevement",
                                                                 "time_point"="time_point", 
                                                                 "real_time_point"="real_time_point", 
                                                                 "condition_biologique"="condition_biologique",
                                                                 "REA"="REA",
                                                                 "numero_patient"="numero_patient",
                                                                 "REPONSE"="REPONSE")) %>%
  mutate(label_patient = as.character(numero_patient))

test_cor <- cor.test(data_fi$SelectPCA.PC1.x, data_fi$SelectPCA.PC1.y, method= "pearson")
options( "digits"=3)
ggplot(data_fi, aes(x = SelectPCA.PC1.x, y = SelectPCA.PC1.y, label = label_patient)) + 
  geom_point(aes(color = REPONSE, shape = real_time_point), 
             size = 3)+
  geom_text_repel(# nudge_x=0.6,
    # nudge_y=0.15,
    show.legend = F,
    max.overlaps  = Inf)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=0, 
           y=-30, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  labs(title = "PC1_VT2_IFN_geneset / PC1_gene_DE_VT1_T_VT2",
       x = "PC1_VT2_IFN_geneset",
       y = "PC1_gene_DE_VT1_T_VT2") + 
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size=24),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))


##################################
data_fi <- inner_join(PC1_VT_IFN_geneset, PC1_VT_gene_DE_VT1vsT, by = c("jours_prelevement" = "jours_prelevement",
                                                                     "time_point"="time_point", 
                                                                     "real_time_point"="real_time_point", 
                                                                     "condition_biologique"="condition_biologique",
                                                                     "REA"="REA",
                                                                     "numero_patient"="numero_patient",
                                                                     "REPONSE"="REPONSE")) %>%
  mutate(label_patient = as.character(numero_patient))

test_cor <- cor.test(data_fi$SelectPCA.PC1.x, data_fi$SelectPCA.PC1.y, method= "pearson")
options( "digits"=3)
ggplot(data_fi, aes(x = SelectPCA.PC1.x, y = SelectPCA.PC1.y, label = label_patient)) + 
  geom_point(aes(color = REPONSE, shape = real_time_point), 
             size = 3)+
  geom_text_repel(# nudge_x=0.6,
    # nudge_y=0.15,
    show.legend = F,
    max.overlaps  = Inf)+
  scale_color_manual(breaks = c("NR","R","RP","T"),
                     values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
  # scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
  #                    values = c("darkorange", "#DDCC77","cornflowerblue",
  #                               "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=10, 
           y=0, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  labs(title = "PC1_VT_IFN_geneset / PC1_VT_gene_DE_VT1vsT",
       x = "PC1_VT_IFN_geneset",
       y = "PC1_VT_gene_DE_VT1vsT") + 
  theme(plot.title = element_text(size=30),
        axis.title = element_text(size=24),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))




