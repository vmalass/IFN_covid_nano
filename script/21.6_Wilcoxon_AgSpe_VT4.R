Test stat frequence AG spé chez tout les patients à VT4







# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 1)
antigene <- data[,c("Groupe", "Ag_spe", "Time")]

T_F <- antigene$Time %in% c("VT4")
antigene <- antigene[T_F == T,]

dplyr::count(antigene, Time, Groupe)

ggplot(antigene , aes(x = factor(Time, levels = c("VT4")), 
                   y = Ag_spe, 
                   color =  Groupe)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(1), 
              xmin = c(0.8), 
              xmax = c(1.2),
              annotation = c("0.2")) + 
  labs(title = "Quantity of all specific antigens  at VT4",
       subtitle = "VT4 (13 R / 5 RP)",
       y =  "Frequence of Ag-specific CD8+ cells IFNg+",
       x = "")

wilcox_test(antigene,formula(paste0( "Ag_spe", "~ Groupe")))





antigene <- data[,c("Groupe", "IFNg_sous_clust_3_CD8_AG_SPE", "IFNg_sous_clust_4_CD8_AG_SPE", "IFNg_sous_clust_5_CD8_AG_SPE", "Time")]

T_F <- antigene$Time %in% c("VT4")
antigene <- antigene[T_F == T,]
antigene <- na.omit(antigene)
antigene[,2:4] <- lapply(antigene[,2:4], as.numeric)

dplyr::count(antigene, Time, Groupe)

ggplot(antigene , aes(x = factor(Time, levels = c("VT4")), 
                      y = IFNg_sous_clust_3_CD8_AG_SPE, 
                      color =  Groupe)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(35), 
              xmin = c(0.8), 
              xmax = c(1.2),
              annotation = c("**")) + 
  labs(title = "Quantity of antigen specific IFN gamma clust 3 at VT4",
       subtitle = "VT4 (11 R / 5 RP)",
       y =  "Frequence of Ag-specific CD8+ cells IFNg+ cluster 3",
       x = "")

wilcox_test(antigene,formula(paste0( "IFNg_sous_clust_3_CD8_AG_SPE", "~ Groupe")))



ggplot(antigene , aes(x = factor(Time, levels = c("VT4")), 
                      y = IFNg_sous_clust_4_CD8_AG_SPE, 
                      color =  Groupe)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(115), 
              xmin = c(0.8), 
              xmax = c(1.2),
              annotation = c("*")) + 
  labs(title = "Quantity of antigen specific IFN gamma clust 4 at VT4",
       subtitle = "VT4 (11 R / 5 RP)",
       y =  "Frequence of Ag-specific CD8+ cells IFNg+ cluster 4",
       x = "")

wilcox_test(antigene,formula(paste0( "IFNg_sous_clust_4_CD8_AG_SPE", "~ Groupe")))



ggplot(antigene , aes(x = factor(Time, levels = c("VT4")), 
                      y = IFNg_sous_clust_5_CD8_AG_SPE, 
                      color =  Groupe)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(40), 
              xmin = c(0.8), 
              xmax = c(1.2),
              annotation = c("**")) + 
  labs(title = "Quantity of antigen specific IFN gamma clust 5 at VT4",
       subtitle = "VT4 (11 R / 5 RP)",
       y =  "Frequence of Ag-specific CD8+ cells IFNg+ cluster 5",
       x = "")

wilcox_test(antigene,formula(paste0( "IFNg_sous_clust_5_CD8_AG_SPE", "~ Groupe")))
