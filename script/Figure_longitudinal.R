
library(grDevices)

# 2-ouverture des fichier-------------------------------------------------------
palette <- grDevices::blues9

rm(list = ls())
load("data/1.4_mat_pat_clean_FIGURE.rds") #ouverture de la svg

T_F <- mat_pat_clean$time_point %in% c("REA", "T")
mat_pat_clean <- mat_pat_clean[T_F == F, ]

mat <- mat_pat_clean[,737:743]
mat <- arrange(mat, jours_prelevement)
mat$jours_prelevement <- as.character(mat$jours_prelevement)
mat$jours_prelevement <- as.numeric(mat$jours_prelevement)

ggplot(mat, aes(x = jours_prelevement)) +
  geom_histogram(bins = 29L,color = "gray50", fill="white") +
  theme_classic() +
  geom_vline(xintercept = c(6.5, 13.5, 20.5), linetype="longdash", color = "#d18975", linewidth = 1) +
  coord_cartesian(xlim = c(1.5, 28.5), ylim = c(0.5,16)) +
  labs(title="Longitudinal analysis of SARS-CoV-2 infected individuals",
       x = "days after onset of symptoms",
       y = "Number of patients") +
  annotate(geom="text", 
           x=3, 
           y=15, 
           label= paste0("VT1"),
           color="#08306B",
           size = 6) +
  annotate(geom="text", 
           x=10, 
           y=15, 
           label= paste0("VT2"),
           color="#2171B5",
           size = 6) +
  annotate(geom="text", 
           x=16.5, 
           y=15, 
           label= paste0("VT3"),
           color="#6BAED6",
           size = 6) +
  annotate(geom="text", 
           x=25, 
           y=15, 
           label= paste0("VT4"),
           color="#C6DBEF",
           size = 6)


tab <- count(mat, jours_prelevement)
tab <- rename(tab, number_patient = n)

ggplot(tab, aes(x = jours_prelevement, y = number_patient))+
  geom_line(color = "grey30") +
  geom_point(shape = 21, fill="#5DA5DA", color="grey30", size = 2) +
  theme_classic() +
  geom_vline(xintercept = c(6.5, 13.5, 20.5), linetype="longdash", color = "grey70", linewidth = 1) +
  coord_cartesian(xlim = c(1.5, 28.5), ylim = c(0.5,16)) +
  scale_y_continuous(breaks=seq(0, 21, 2)) +
  scale_x_continuous(breaks=seq(0, 30, 1)) +
  labs(title="Longitudinal analysis of SARS-CoV-2 infected individuals",
       x = "days after onset of symptoms",
       y = "Number of patients") +
  annotate(geom="text", 
           x=3, 
           y=15, 
           label= paste0("VT1"),
           color="#08306B",
           size = 6) +
  annotate(geom="text", 
           x=10, 
           y=15, 
           label= paste0("VT2"),
           color="#08519C",
           size = 6) +
  annotate(geom="text", 
           x=16.5, 
           y=15, 
           label= paste0("VT3"),
           color="#4292C6",
           size = 6) +
  annotate(geom="text", 
           x=25, 
           y=15, 
           label= paste0("VT4"),
           color="#6BAED6",
           size = 6)


"#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"

