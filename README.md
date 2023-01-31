# NanoString_Covid
Data issu des données NanoString du papier ([Early nasal type I IFN immunity against SARS-CoV-2 is compromised in patients with autoantibodies against type I IFNs](https://doi.org/10.1084/jem.20211211)).  
800 gènes étudiers sur 72 patients avec 12 REA covid+, 47 covid, 7 Hcov et 7 témoins.  

# Recap projet Covid données NanoString
- Echelle de temps VT1 J1 à J6, VT2 J7 à J13, VT3 J14 à J20, VT4 J21 à …
- Doublon n°24 en VT2
- Suppression des patients n° 28, 56, 64 et 66 car 1 prélèvement (analyse de base n°24 car pas de V1 mais avec échelle VT plusieurs non pas de VT1 … Donc gardé)
- Choix MAD (tous les patients dans au moins 2 sur 3 gènes) :
    - NR : 49 / 50 
    - RP : 4 / 10 / 15 / 21 / 23 / 26 / 29 / 48 / 60 / 61
- EN VT1 :
    - Clustering heatmap en VT1 sur HVG 
        - NR : 49 / 50 se regroupe avec les T. Le NR et R mélangé
        - 2 DE effectué entre NRvsR à partir de tous les gènes et des HVG
            - 44 gènes commun (81%) on décide de prendre les gènes uniques entre les deux DE soi 54 gènes
                - A partir des 54 gènes on fait une PCA sur tous les times points (VT). On récupère les valeurs de la PC1 pour regarder PC1 vs jours de prélèvement de tous les times points.
- En VT2 :
    - Clustering heatmap en VT2 sur HVG
        - RP : 4 / 10 / 15 / 21 / 26 / 29 / 48 / 61  —> 23 et 60 non clusterisé
        - 2 DE effectué entre RvsRP à partir de tous les gènes et des HVG avec comme RP (4 / 10 / 15 / 21 / 26 / 29 / 48 / 61)
            - 51 gènes communs (75%) on décide de prendre les gènes uniques entre les deux DE soi 68 gènes
                - A partir des 68 gènes on fait une PCA sur tous les times points (VT). On récupère les valeurs de la PC1 pour regarder PC1 vs jours de prélèvement de tous les times points.
- DE entre TvsNR : 1 gène DE LILRA3 (mono & B high expression + dentritic & NK low expression) 
- Différence génique entre les deux DE NRvsR et RvsRP :
    - 43 gènes commun entre les deux DE sur 79 soit 54% communs

## Definitions des groupes :
| Tables | Are | Cool | |----------|:-------------:|------:| | col 1 is| left-aligned | $1600 |	
Groupes	numéro patient	nombre de patients
NR	49 / 50 	2
NR-	2 / 13 / 62 	3
R	X	25
RP-	19 / 23 / 52 / 59 	4
RP	4 / 10 / 15 / 21 / 26 / 29 / 48 / 61 	8
A	60 	1
T	X	7

## On refait l’analyse avec les nouveaux groupes
- 2 DE entre NRvsR à partir de tous les gènes et des HVG : 
    - 55 gènes communs (76%) on décide de prendre les gènes uniques entre les deux DE soi 72 gènes uniques
    - Visualisation par heatmap : on clusterise les 2 NR avec T et les 3 NR- juste à coté des T. R, RP- et RP mélangé en VT1
- 2 DE entre RvsRP à partir de tous les gènes et des HVG : 
    - 51 gènes communs (70%) on décide de prendre les gènes uniques entre les deux DE soi 72 gènes uniques
    - Visualisation par heatmap : on clusterise les RP ensemble avec un T et le patient 60 (A). Les patients RP- clusterise ensemble avec le 24 (mais vu avec PC1 vs jours de prélèvement que classé dans R)
- Différence génique entre les deux DE NRvsR et RvsRP :
    - 45 gènes en commun entre les deux DE sur 81 soit 55% communs

## Corrélation entre les différentes PC1
Dans l’analyse nous avons générés différentes « PC1 » issu de DE entre NRvsR, RvsRP à différents times points (VT1, 2) mais également avec une signature IFN obtenus sur GSEA (fichier data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.tsv) [lien GSEA](http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_INTERFERON_SIGNALING.html) 
Toutes les corrélations sont proches de 1, donc les signatures sont similaires on choisit pour la suite les gènes DE NRvsR unique entre le DE de tous les gènes et HVG.

## Analyse des fréquences cellulaires à V1, V4 et 6 mois
Pour cela nous sommes partie des prélèvements V1 et fait une PCA à partir des gènes DE NRvsR. 
Nous regardons les fréquences en fonction de la PC1 qui correspond à une réponse IFN et nous permet de séparer les différents groupes. On calcule la corrélation ainsi que la p-value associé avec un cor_test.
10 populations cellulaires ressortent significatives à V1.
Pas de population significative à 6 mois.

victor.malassigne@inserm.fr