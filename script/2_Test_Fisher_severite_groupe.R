# test stat Fisher car une des cellules du tableau de contingence sont inférieures à 5
# Hypothèses
# H0 : les variables sont indépendantes, il n'y a pas de relation entre les deux variables catégorielles. Connaître la valeur d'une variable n'aide pas à prédire la valeur de l'autre variable
# H1 : les variables sont dépendantes, il existe une relation entre les deux variables catégorielles. Connaître la valeur d'une variable permet de prédire la valeur de l'autre variable

rm(list=ls())

dat <- data.frame("R" = c(2,25),
                  "RP" = c(2,7),
                  row.names = c("long_covid", "mild"),
                  stringsAsFactors = F)

mosaicplot(dat, 
           main = "Répartion des patients en fonction de leur groupe et de leur sévérité",
           color = T)

test <- fisher.test(dat)
test
test$p.value

# Visualisation
x <- c()
for (row in rownames(dat)) {
  for (col in colnames(dat)) {
    x <- rbind(x, matrix(rep(c(row, col), dat[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
colnames(df) <- c("severite", "groupe")
df

library(ggstatsplot)
ggbarstats(
  df, severite, groupe,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
)


# Test de proportion mild et long covid ----
## savoir si deux proportions mesurées sont identiques ou non
# H0 les proportion sont significativement différente
# H1 les proportion ne sont pas significativement différente

mat<-matrix(c(2,27,2,10),2)

dimnames(mat) <- list(c("long_covid","total") ,c("R","RP"))
prop.test(mat)

mat<-matrix(c(2,10,2,27),2)

dimnames(mat) <- list(c("long_covid","total") ,c("R","RP"))
prop.test(mat)






