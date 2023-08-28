# Abstract 

This study examines the transcriptomic responses of COVID-19 patients at different stages of infection. The study reveals distinct gene expression dynamics that allow patients classification. Two main groups were investigated, based on the duration of the transcriptomic interferons response: responders (R) and persistent responders (RP). Quantifications of circulating interferons unveil an increased production of interferon gamma in persistent responders. Moreover, notable variations in the frequency and roles of CD8+ TEMRA or exhausted T lymphocytes are observed between the R and RP groups, along with higher long-term IFN-γ production in RP. The findings underscore the significance of the interferon response while delineating distinct nuances in patient immune dynamics through these two groups. This contributes to a better understanding of responses to COVID-19 infection.

![RSAD2 Kinetic](result/Kinetic_RSAD2.pdf)

# NanoString_Covid
Data issu des données NanoString du papier ([Early nasal type I IFN immunity against SARS-CoV-2 is compromised in patients with autoantibodies against type I IFNs](https://doi.org/10.1084/jem.20211211)).  
800 gènes étudiés sur 72 patients avec 12 REA covid+, 47 covid, 7 Hcov et 7 témoins.  

# Recap projet Covid données NanoString
- Echelle de temps VT1 J1 à J6, VT2 J7 à J13, VT3 J14 à J20, VT4 J21 à …
- Suppression des patients n° 28, 56, 64 et 66 car 1 prélèvement
- Choix MAD (tous les patients dans au moins 2 sur 3 gènes) :
    - 3 groupes NR / R / RP
- PCA à partir de tous les prelevements et gènes avec les times points : la PC1 drive 28% de la variabilité des données, séparation des VT1 des autres avec quelques VT2 dans le nuage.

- DE entre VT1 et T à partir de tous les gènes : signature de 65 gènes. On utilise cette liste pour la suite de l'analyse. NR trop de peu de patients on utilise que les groupes R et RP

- Heatmaps généré à partir des gènes issus du DE VT1 vs T:
    - VT1
    - VT2

On identifie 23 R et 8 RP que l'on utilisera pour la suite de l'analyse.


- PCA à partir de la signature (DE VT1vsT) avec les times points :
    - PC1 81% de la varabilité
        - On plot la PC1 en fonction du jour après début des symptômes. On observe deux courbes avec les patients R qui décroit plus rapidement au cours du temps que les patients RP. Pas de biais du au temps lors du prélèvement

- Dans le [papier](https://doi.org/10.1084/jem.20211211) dont les data proviennent, ils mettent en évidence la réponse IFN avec l'évolution de la covid-19. Nous utilisons la signature IFN obtenue sur GSEA (fichier : data/REACTOME_INTERFERON_SIGNALING.v2022.1.Hs.tsv) [lien GSEA](http://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_INTERFERON_SIGNALING.html).

- Corrélation entre la PC1 de l’ACP obtenue à partir du DE VT1vsT et la signature IFN. Nous obtenons un corrélation de 0.99 avec une pvalue de 1.31e-117. Les deux listes de gènes donnent la même classification.


# Analyse des fréquences cellulaires à V1, V4 et 6 mois
Principale différence : monocyte, CD8 TEMRA, DC, CD8 TEMRA CD57+, CD8 Antigène spécifique.
On compare différente fréquence de population cellulaire entre le groupe R et RP.

# Analyse complémentaire
Dosage IFN gamma, beta et lambda (soit type I, II et III).
Dosage anticorps IgG anti RBD
Dosage IFN gamma à 6 mois


victor.malassigne@inserm.fr