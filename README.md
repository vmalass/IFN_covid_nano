# Sommaire 
- [NanoString_Covid](# NanoString_Covid)
- [Recap projet Covid données NanoString](# Recap projet Covid données NanoString)
- [Definitions des groupes](## Definitions des groupes)
- [On refait l’analyse avec les nouveaux groupes](## On refait l’analyse avec les nouveaux groupes)
- [Analyse des fréquences cellulaires à V1, V4 et 6 mois](## Analyse des fréquences cellulaires à V1, V4 et 6 mois)
- [PCA pour séparation des groupes](### PCA pour séparation des groupes)
- [Récapitulatif nombre de patients par expériences](### Récapitulatif nombre de patients par expériences)
- [Corrélation](### Corrélation)
- [Test stat Krukal-Wallis et Dunn entre R et RP](### Test stat Krukal-Wallis et Dunn entre R et RP)

# NanoString_Covid
Data issu des données NanoString du papier ([Early nasal type I IFN immunity against SARS-CoV-2 is compromised in patients with autoantibodies against type I IFNs](https://doi.org/10.1084/jem.20211211)).  
800 gènes étudiers sur 72 patients avec 12 REA covid+, 47 covid, 7 Hcov et 7 témoins.  

# Recap projet Covid données NanoString
- Echelle de temps VT1 J1 à J6, VT2 J7 à J13, VT3 J14 à J20, VT4 J21 à …
- Doublon n°24 en VT2
- Suppression des patients n° 28, 56, 64 et 66 car 1 prélèvement
- Choix MAD (tous les patients dans au moins 2 sur 3 gènes) :
    - NR : 49 / 50 
    - RP : 4 / 10 / 15 / 21 / 23 / 26 / 29 / 48 / 60 / 61
- PCA à partir de tous les prelevements et gènes avec les times points : la PC1 drive 28% de la variabilité des données, séparation des VT1 des autres avec quelques VT2 dans le nuage.

- DE entre VT1 et T à partir de tous les gènes : signature de 65 gènes. On utilise cette liste pour la suite de l'analyse.

- EN VT1 :
    - Clustering heatmap en VT1 à partir de la signature 
        - NR (49 / 50) se regroupe avec les T mais un R (13) est dans le même cluster. On identifie un cluster qui à la même tendance (2 / 22 / 51 / 62), on les identifie comme NR-.

- En VT2 :
    - Clustering heatmap en VT2 à partir de la signature 
        - RP : 4 / 10 / 15 / 21 / 26 / 29 / 48 / 60 / 61  —> 23 non clusterisé. On identifie également une faible persistance de la signature dans un cluster (19 / 23 / 24 / 52 / 59), on les identifies comme RP- (sauf le 24 qui à deux prélèvements en VT2)

*** Le patient 60 est particulier, avec le clustering en VT1 il est dans les NR- et en VT2 avec les RP. Nous avons de créer un groupe pour l'identifier (A)

- PCA à partir de la signature (DE VT1vsT) avec les times points :
    - PC1 81% de la varabilité
        - On plot la PC1 en fonction du jour après début des symptômes. On observe deux courbes avec les patients R qui décroit plus rapidement au cours du temps que les patients RP.

- Dans le [papier](https://doi.org/10.1084/jem.20211211) dont les data proviennent, ils mettent en évidence la réponse IFN avec l'évolution de la covid-19. Nous utilisons la signature IFN obtenue sur GSEA (fichier : data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.tsv) [lien GSEA](http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_INTERFERON_SIGNALING.html).

- Corrélation entre la PC1 de l’ACP obtenue à partir du DE VT1vsT et la signature IFN. Nous obtenons un corrélation de 0.99 avec une pvalue de 1.31e-117.

## Definitions des groupes :
| Groupes        | numéro patient      | nombre de patients |
| ------|-----|-----|
| NR | 13 / 49 / 50  |	3 |
| NR- |	2 / 22 / 51 / 62 |	4 |
| R	| X |	23 |
| RP-	| 19 / 23 / 52 / 59 |	4 |
| RP	| 4 / 10 / 15 / 21 / 26 / 29 / 48 / 61 |	8 |
| A |	60 |	1 |
| T |	X |	7 |

NR : non répondeur / NR- : non répondeur faible / R : répondeur / RP- : répondeur persistant faible / RP : répondeur persistant / A : atypique / T : témoin

## On refait l’analyse avec les nouveaux groupes
- EN VT1 :
    - Clustering heatmap en VT1 à partir de la signature avec les groupes finaux
        - NR (13 / 49 / 50) se regroupe avec les T dans le même cluster. Le cluster NR- (2 / 22 / 51 / 62) avec un T et le patient A (60). Un continuum de réponse regroupe des patients R, RP- et RP.

- EN VT2 :
    - Clustering heatmap en VT1 à partir de la signature avec les groupes finaux
      - RP (4 / 10 / 15 / 21 / 26 / 29 / 48 / 61) se regroupe le même cluster avec le patient A (60) mais également un T. Le cluster RP- (19 / 23 / 52 / 59) avec un patient R (24 qui a 2 prélèvements à VT2). Un continuum de réponse regroupe des patients R, NR-, NR et T.
      


## Analyse des fréquences cellulaires à V1, V4 et 6 mois
### PCA pour séparation des groupes
L'objectif est de séparer au mieux nos groupes en fonction de leur réponse pour obsrevé si il y a des corrélations. Les échantillons analysés par cytométrie à V1 n'ont pas subi la réorganisation avec les nouveaux times points (VT...). Nous avons fait une PCA à partir de la signature (DE VT1vsT) sur les prélèvements V1 pour ne pas perdre de patient. Nous utilisons la PC1(75% de la variabilité) comme axe des x et nous affichons en axe des y les fréquences cellulaires.

### Récapitulatif nombre de patients par expériences
| Groupes | NanoString | main subset V1 | stim V4 | stim 6 mois | Commun |
| ------|-----|-----|------|-----|-----|
| NR | 3 | 1 | 1 | 1 | 1 |
| NR- | 4 | 2 | 1 | 1 | X |
| R	| 23 | 17 | 13 | 18 | 12 |
| RP-	| 4 | 2 | 2 | 4 | 2 |
| RP	| 8 | 7 | 5 | 7 | 5 |
| A | 1 | X | X | 1 | X |

### Corrélation
Nous avons regardé deux choses, les corrélations entre la PC1 et les fréquences mais également le jour de prélèvement après symptômes et les fréquences cellulaire. L'objectif était de regarder si les corrélations étaient bien dues au groupe de patients et non au temps.

V1 supervisé : 10 populations cellulaires significatives
V1 non supervisé : 7 clusters significatifs

### Test stat Krukal-Wallis et Dunn entre R et RP
V1 supervisé : 9 populations cellulaires significatives
V1 non supervisé : 11 clusters significatifs

V4 non stim supervisé : 0 population cellulaire significative
V4 stim supervisé : 1 population cellulaire significative
V4 non stim non supervisé : 4 clusters significatifs
V4 stim non supervisé : 1 cluster significatif
V4 antigène spé stim non supervisé : 5 clusters significatifs

6 mois supervisé : 3 populations cellulaires significatives

victor.malassigne@inserm.fr