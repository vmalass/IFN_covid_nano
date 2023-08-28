



# Import data-------------------------------------------------------------------
rm(list = ls())
load("data/1.1_mat_pat_clean.rds") #ouverture de la svg
data_GWAS <- read_xlsx("data/49_gene_id_GWAS_nature.xlsx", sheet = 2)

# commun gene-------------------------------------------------------------------

gene_nano <- colnames(mat_pat_clean[1:736])
gene_GWAS <- data_GWAS$`Nearest gene`

commun <- intersect(gene_nano, gene_GWAS)

# Visualization-----------------------------------------------------------------
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg
T_F <- mat_pat_clean$REPONSE %in% "REA"
mat_pat_clean_sans_REA <- mat_pat_clean[T_F == F,]

T_F <- mat_pat_clean_sans_REA$REPONSE %in% "A"
mat_pat_clean_sans_REA <- mat_pat_clean_sans_REA[T_F == F,]
a <- data.frame(colnames(mat_pat_clean_sans_REA))


T_F <- colnames(mat_pat_clean_sans_REA) %in% commun
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]


ggplot(mat_macro, aes(x = REPONSE, y = IFNAR2, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = IL10RB, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = JAK1, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = OAS1, fill = real_time_point)) +
  geom_violin() +
  
  ggplot(mat_macro, aes(x = REPONSE, y = TYK2, fill = real_time_point)) +
  geom_violin()


## Test stat---------------------------------------------------------------------
T_F <- mat_macro$REPONSE %in% "NR"
mat_macro <- mat_macro[T_F == F, ]
summary(mat_macro)
for (gene in commun) {
  print(mat_macro %>% 
          group_by(real_time_point) %>%
          get_summary_stats(gene, type = "common"))    #### problème avec titre col ####
  
  # Test kruskal
  res.kruskal <- mat_macro %>% kruskal_test(formula(paste0( gene," ~ real_time_point")))
  print(res.kruskal)
  
  # Taille de l’effet
  print(mat_macro %>% kruskal_effsize(formula(paste0(gene," ~ real_time_point"))))
  
  # Comparaisons par paires test de Dunn
  pwc <- mat_macro %>% 
    dunn_test(formula(paste0(gene," ~ real_time_point")), p.adjust.method = "bonferroni") 
  print(pwc)
  
  pwc <- pwc %>% add_xy_position(x = "real_time_point")
  pwc %>% 
    mutate(xmin = xmin-1,
           xmax = xmax-1) -> pwc
  print(ggboxplot(mat_macro, 
                  x = "real_time_point", 
                  y = gene, 
                  color = "REPONSE", 
                  palette = c("cornflowerblue","brown3", "chartreuse4")) +
          stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
               caption = get_pwc_label(pwc))) +
    labs(title = "Titre en Ac en fonction des groupes",
         y = "Titre en Anticorps") +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(size=20)) +
    theme(legend.position="right")
  
}


## Test stat R vs RP à VT1/2/3/4------------------------------------------------
### Viral load V1 ###
V1 <- mat_macro[mat_macro$real_time_point == "VT1",]

V1 <- V1%>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

for (gene in commun) {
  # Test stat
  res.wilcox <- V1 %>%
    wilcox_test(formula(paste0(gene," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  print(ggboxplot(V1, 
                  title = paste0("Gene expression of ", gene, " at VT1"),
                  x = "REPONSE", 
                  y = paste0(gene), 
                  ylab = gene,
                  xlab = "Groups", 
                  add = "point") + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}

### Viral load V1 ###
V2 <- mat_macro[mat_macro$real_time_point == "VT2",]

V2 <- V2 %>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

for (gene in commun) {
  # Test stat
  res.wilcox <- V2 %>%
    wilcox_test(formula(paste0(gene," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  print(ggboxplot(V1, 
                  title = paste0("Gene expression of ", gene, " at VT2"),
                  x = "REPONSE", 
                  y = paste0(gene), 
                  ylab = gene,
                  xlab = "Groups", 
                  add = "point") + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}


### Viral load V3 ###
V3 <- mat_macro[mat_macro$real_time_point == "VT3",]
V3 <- V3 %>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

for (gene in commun) {
  # Test stat
  res.wilcox <- V3 %>%
    wilcox_test(formula(paste0(gene," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  print(ggboxplot(V1, 
                  title = paste0("Gene expression of ", gene, " at VT3"),
                  x = "REPONSE", 
                  y = paste0(gene), 
                  ylab = gene,
                  xlab = "Groups", 
                  add = "point") + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}

### Viral load V4 ###
V4 <- mat_macro[mat_macro$real_time_point == "VT4",]
V4 <- V4 %>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

for (gene in commun) {
  # Test stat
  res.wilcox <- V4 %>%
    wilcox_test(formula(paste0(gene," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  print(ggboxplot(V1, 
                  title = paste0("Gene expression of ", gene, " at VT4"),
                  x = "REPONSE", 
                  y = paste0(gene), 
                  ylab = gene,
                  xlab = "Groups", 
                  add = "point") + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}


# Other genes-------------------------------------------------------------------
a <- as.data.frame(gene_nano)
gene <- c("IL6", "IL6R", "IL6ST", "OAS1", "OAS2", "OAS3", "OASL", "RNASEL","PLCG2" )
T_F <- colnames(mat_pat_clean_sans_REA) %in% gene
mat_macro <- mat_pat_clean_sans_REA[, T_F == T]

mat_macro <- cbind(mat_macro, mat_pat_clean_sans_REA[,736 : 743])

T_F <- mat_macro$REPONSE %in% "Other"
mat_macro <- mat_macro[T_F == F, ]

ggplot(mat_macro, aes(x = REPONSE, y = IL6, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = IL6R, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = IL6ST, fill = real_time_point)) +
  geom_violin()

ggplot(mat_macro, aes(x = REPONSE, y = OAS1, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = OAS2, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = OAS3, fill = real_time_point)) +
  geom_violin() +

  ggplot(mat_macro, aes(x = REPONSE, y = OASL, fill = real_time_point)) +
  geom_violin()

ggplot(mat_macro, aes(x = REPONSE, y = RNASEL, fill = real_time_point)) +
  geom_violin()+
  
  ggplot(mat_macro, aes(x = REPONSE, y = PLCG2, fill = real_time_point)) +
  geom_violin()


## Test stat---------------------------------------------------------------------
T_F <- mat_macro$REPONSE %in% "NR"
mat_macro <- mat_macro[T_F == F, ]
summary(mat_macro)
for (gene in gene) {
  print(mat_macro %>% 
          group_by(real_time_point) %>%
          get_summary_stats(gene, type = "common"))    #### problème avec titre col ####
  
  # Test kruskal
  res.kruskal <- mat_macro %>% kruskal_test(formula(paste0( gene," ~ real_time_point")))
  print(res.kruskal)
  
  # Taille de l’effet
  print(mat_macro %>% kruskal_effsize(formula(paste0(gene," ~ real_time_point"))))
  
  # Comparaisons par paires test de Dunn
  pwc <- mat_macro %>% 
    dunn_test(formula(paste0(gene," ~ real_time_point")), p.adjust.method = "bonferroni") 
  print(pwc)
  
  pwc <- pwc %>% add_xy_position(x = "real_time_point")
  pwc %>% 
    mutate(xmin = xmin-1,
           xmax = xmax-1) -> pwc
  print(ggboxplot(mat_macro, 
                  x = "real_time_point", 
                  y = gene, 
                  color = "REPONSE", 
                  palette = c("cornflowerblue","brown3", "chartreuse4")) +
          stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
               caption = get_pwc_label(pwc))) +
    labs(title = "Titre en Ac en fonction des groupes",
         y = "Titre en Anticorps") +
    theme(axis.title.x = element_blank(),
          plot.title = element_text(size=20)) +
    theme(legend.position="right")
  
}

















T_F <- mat_macro$REPONSE %in% "NR"
mat_macro <- mat_macro[T_F == F, ]
summary(mat_macro)
print(mat_macro %>% 
        group_by(real_time_point) %>%
        get_summary_stats(RNASEL, type = "common"))    #### problème avec titre col ####

# Test kruskal
res.kruskal <- mat_macro %>% kruskal_test(formula(paste0( 'RNASEL'," ~ real_time_point")))
print(res.kruskal)

# Taille de l’effet
print(mat_macro %>% kruskal_effsize(formula(paste0('RNASEL'," ~ real_time_point"))))

# Comparaisons par paires test de Dunn
pwc <- mat_macro %>% 
  dunn_test(formula(paste0('RNASEL'," ~ real_time_point")), p.adjust.method = "bonferroni") 
print(pwc)

pwc <- pwc %>% add_xy_position(x = "real_time_point")
pwc %>% 
  mutate(xmin = xmin-1,
         xmax = xmax-1) -> pwc
print(ggboxplot(mat_macro, 
                x = "real_time_point", 
                y = "RNASEL", 
                color = "REPONSE", 
                palette = c("cornflowerblue","brown3", "chartreuse4")) +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
             caption = get_pwc_label(pwc))) +
  labs(title = "RNASEL en fonction des groupes",
       y = "RNASEL") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=20)) +
  theme(legend.position="right")





