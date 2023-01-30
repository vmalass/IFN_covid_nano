# NanoString_Covid
Data issu des données NanoString du papier ([Early nasal type I IFN immunity against SARS-CoV-2 is compromised in patients with autoantibodies against type I IFNs](https://doi.org/10.1084/jem.20211211)).  
800 gènes étudiers sur 72 patients avec 12 REA covid+, 47 covid, 7 Hcov et 7 témoins.  

Recap projet Covid données NanoString
- Echelle de temps VT1 J1 à J6, VT2 J7 à J13, VT3 J14 à J20, VT4 J21 à …
- Doublon n°24 en VT2
- Suppression patient n° 28, 56, 64 et 66 car 1 prélèvement (analyse de base n°24 car pas de V1 mais avec échelle VT plusieurs non pas de VT1 … Donc gardé)
- Choix MAD (tous au moins 2 fois sur 3 gènes) :
    - NR : 49 / 50 
    - RP : 4 / 10 / 15 / 21 / 23 / 26 / 29 / 48 / 60 / 61
- Clustering heatmap sur HVG 
    - NR : 49 / 50 
    - RP : 4 / 10 / 15 / 21 / 26 / 29 / 48 / 61  —> 23 et 60 non clusterisé
- DE entre NR (49 et 50) et R : 3 patients proche 

victor.malassigne@inserm.fr