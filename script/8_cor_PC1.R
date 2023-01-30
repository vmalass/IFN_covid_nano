library(ggrepel)
library(ggplot2)

rm(list = ls())

PC1_VT_gene_DE_NRvsR <- read.table("data/PC1_VT_gene_DE_NRvsR.txt")
PC1_VT_gene_DE_RvsRP <- read.table("data/PC1_VT_gene_DE_RvsRP.txt")
PC1_VT_IFN_geneset <- read.table("data/PC1_VT_IFN_geneset.txt")


PC1_VT1_gene_DE_NRvsR <- read.table("data/PC1_VT1_gene_DE_NRvsR.txt")
PC1_VT1_IFN_geneset <- read.table("data/PC1_VT1_IFN_geneset.txt")

PC1_VT2_IFN_geneset <- read.table("data/PC1_VT2_IFN_geneset.txt")
PC1_VT2_gene_DE_RvsRP <- read.table("data/PC1_VT2_gene_DE_RvsRP.txt")


######################### VT #################################################
data_fi <- inner_join(PC1_VT_gene_DE_NRvsR, PC1_VT_gene_DE_RvsRP, by = c("numero_patient" = "numero_patient",
                                                                         "REPONSE"="REPONSE", 
                                                                         "jours_prelevement"="jours_prelevement", 
                                                                         "time_point"="time_point",
                                                                         "real_time_point"="real_time_point",
                                                                         "condition_biologique"="condition_biologique",
                                                                         "REA"="REA")) %>%
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
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) + 
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=10, 
           y=-5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  ggtitle("PC1_VT_gene_DE_NRvsR / PC1_VT_gene_DE_RvsRP") + 
  theme(plot.title = element_text(size=23),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))

# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
data_fi <- inner_join(PC1_VT_gene_DE_NRvsR, PC1_VT_IFN_geneset, by = c("numero_patient" = "numero_patient",
                                                                         "REPONSE"="REPONSE", 
                                                                         "jours_prelevement"="jours_prelevement", 
                                                                         "time_point"="time_point",
                                                                         "real_time_point"="real_time_point",
                                                                         "condition_biologique"="condition_biologique",
                                                                         "REA"="REA")) %>%
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
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) + 
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=10, 
           y=-5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  ggtitle("PC1_VT_gene_DE_NRvsR / PC1_VT_IFN_geneset") + 
  theme(plot.title = element_text(size=23),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))

# ********************************************************************************
# ********************************************************************************
# ********************************************************************************
data_fi <- inner_join(PC1_VT_gene_DE_RvsRP, PC1_VT_IFN_geneset, by = c("numero_patient" = "numero_patient",
                                                                       "REPONSE"="REPONSE", 
                                                                       "jours_prelevement"="jours_prelevement", 
                                                                       "time_point"="time_point",
                                                                       "real_time_point"="real_time_point",
                                                                       "condition_biologique"="condition_biologique",
                                                                       "REA"="REA")) %>%
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
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) + 
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=10, 
           y=-5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  ggtitle("PC1_VT_gene_DE_RvsRP / PC1_VT_IFN_geneset") + 
  theme(plot.title = element_text(size=23),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))

######################### VT1 #################################################
data_fi <- inner_join(PC1_VT1_gene_DE_NRvsR, PC1_VT1_IFN_geneset, by = c("numero_patient" = "numero_patient",
                                                                         "REPONSE"="REPONSE", 
                                                                         "jours_prelevement"="jours_prelevement", 
                                                                         "time_point"="time_point",
                                                                         "real_time_point"="real_time_point",
                                                                         "condition_biologique"="condition_biologique",
                                                                         "REA"="REA")) %>%
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
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) + 
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=7, 
           y=-5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  ggtitle("PC1_VT1_gene_DE_NRvsR / PC1_VT1_IFN_geneset") + 
  theme(plot.title = element_text(size=23),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))

######################### VT2 ##################################################
data_fi <- inner_join(PC1_VT2_IFN_geneset, PC1_VT2_gene_DE_RvsRP, by = c("numero_patient" = "numero_patient",
                                                 "REPONSE"="REPONSE", 
                                                 "jours_prelevement"="jours_prelevement", 
                                                 "time_point"="time_point",
                                                 "real_time_point"="real_time_point",
                                                 "condition_biologique"="condition_biologique",
                                                 "REA"="REA")) %>%
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
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) + 
  # geom_smooth(method = lm, se = FALSE) +
  geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=20, 
           y=0, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="red",
           size = 13) +
  ggtitle("PC1_VT2_IFN_geneset / PC1_VT2_gene_DE_RvsRP") + 
  theme(plot.title = element_text(size=23),
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15))






