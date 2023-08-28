Correlation IFN gamma vs PCA1(VT2)/IFAS1(MBGP-NC)/CD57(V1-4)
Uniquement avec les patoent IFN gamma à VT2  = 4R ET 3RP




# import datat------------------------------------------------------------------
rm(list = ls())
dosage_IFN <- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/Dosage_IFN.xlsx", sheet = 3)
dosage_IFA <- read.xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/titre_virale_IFI27_Ac.xlsx", sheet = 8)
cyto <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 1)

PC1_VT_gene_DE_VT1vsT <- read.table("data/PC1_VT_gene_DE_VT1vsT_mat1.2.txt")

# Modification table------------------------------------------------------------
IFN_gamma <- dosage_IFN[, c("IFN_gamma", "Reponse", "temps", "sample ID", 
                            "numero_patient_victor", "Age", "Sexe", 
                            "severité de COVID", "new_time_point")]
T_F <- IFN_gamma$new_time_point %in% "VT2"
IFN_gamma_VT2 <- IFN_gamma[T_F == T,]

IFA_mb_nc <- dosage_IFA[, c("numero_patient", "Groupe", "ID", "IFA.MBGP", 
                            "IFA.NC", "time_point")]
T_F <- IFA_mb_nc$time_point %in% "S1"
IFA_mb_nc_S1 <- IFA_mb_nc[T_F == T,]

CD57_pos <- cyto[, c("ID", "Groupe", "Time", "CD8_TEMRA_CD57+")]
CD57_pos <- na.omit(CD57_pos)
T_F <- CD57_pos$Time %in% "V1"
CD57_pos_V1 <- CD57_pos[T_F == T,]
CD57_pos_V4 <- CD57_pos[T_F == F,]

CD57_neg <- cyto[, c("ID", "Groupe", "Time", "CD8_TEMRA_CD57-")]
CD57_neg <- na.omit(CD57_neg)
T_F <- CD57_neg$Time %in% "V1"
CD57_neg_V1 <- CD57_neg[T_F == T,]
CD57_neg_V4 <- CD57_neg[T_F == F,]

# Correlation IFNgamma - PC1----------------------------------------------------
data <- inner_join(IFN_gamma_VT2, PC1_VT_gene_DE_VT1vsT, by = c("numero_patient_victor" = "numero_patient",
                                                       # "numero_patient_victor" = "numero_patient",
                                                       # "temps" = "time_point",
                                                       "new_time_point" = "real_time_point"))

test_cor <- cor.test(data$IFN_gamma, data$SelectPCA.PC1, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = SelectPCA.PC1)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 

# Correlation IFNgamma - IFA----------------------------------------------------
data <- inner_join(IFN_gamma_VT2, IFA_mb_nc_S1, by = c("Reponse" = "Groupe",
                                                # "numero_patient_victor" = "numero_patient",
                                                # "temps" = "time_point",
                                                "sample ID" = "ID"))

test_cor <- cor.test(data$IFN_gamma, data$IFA.MBGP, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = IFA.MBGP)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=5, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 



test_cor <- cor.test(data$IFN_gamma, data$IFA.NC, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = IFA.NC)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=3, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 


# Correlation IFNgamma - CD57pos V1----------------------------------------------------
data <- inner_join(IFN_gamma_VT2, CD57_pos_V1, by = c("Reponse" = "Groupe",
                                                       # "numero_patient_victor" = "numero_patient",
                                                       # "temps" = "time_point",
                                                       "sample ID" = "ID"))

test_cor <- cor.test(data$IFN_gamma, data$`CD8_TEMRA_CD57+`, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = `CD8_TEMRA_CD57+`)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=3, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 

# Correlation IFNgamma - CD57pos V4----------------------------------------------------
data <- inner_join(IFN_gamma_VT2, CD57_pos_V4, by = c("Reponse" = "Groupe",
                                                      # "numero_patient_victor" = "numero_patient",
                                                      # "temps" = "time_point",
                                                      "sample ID" = "ID"))

test_cor <- cor.test(data$IFN_gamma, data$`CD8_TEMRA_CD57+`, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = `CD8_TEMRA_CD57+`)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=3, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 

# Correlation IFNgamma - CD57neg VA----------------------------------------------------
data <- inner_join(IFN_gamma_VT2, CD57_neg_V1, by = c("Reponse" = "Groupe",
                                                      # "numero_patient_victor" = "numero_patient",
                                                      # "temps" = "time_point",
                                                      "sample ID" = "ID"))

test_cor <- cor.test(data$IFN_gamma, data$`CD8_TEMRA_CD57-`, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = `CD8_TEMRA_CD57-`)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=3, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 

# Correlation IFNgamma - CD57neg V4----------------------------------------------------
data <- inner_join(IFN_gamma_VT2, CD57_neg_V4, by = c("Reponse" = "Groupe",
                                                      # "numero_patient_victor" = "numero_patient",
                                                      # "temps" = "time_point",
                                                      "sample ID" = "ID"))

test_cor <- cor.test(data$IFN_gamma, data$`CD8_TEMRA_CD57-`, method= "pearson")
options( "digits"=3)
ggplot(data, aes(x = IFN_gamma, y = `CD8_TEMRA_CD57-`)) + #, label = label_patient
  geom_point(aes(color = Reponse),  #, shape = real_time_point
             size = 1.5)+
  scale_color_manual(breaks = c("R","RP","T"),
                     values = c("gray60","#CB2027","chartreuse4"))+
  # geom_rug(aes(color = REPONSE)) +
  annotate(geom="text", 
           x=100, 
           y=3, 
           label= paste0("cor ",round(test_cor$estimate, 2), " pvalue ", signif(test_cor$p.value, 3)),
           color="gray40",
           size = 6) +
  labs(title = "Correlation") + 
  theme_classic() + 
  theme(legend.position = c(0.85,0.2)) 
