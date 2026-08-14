# Screening en deux étages avec agrégation de réglages

Note de méthode issue de la campagne étendue (`results/grid/`). Elle décrit la
seule construction testée qui améliore le screen par indice de queue sans
changer ce qu'il cible, mesure ce qu'elle apporte, et délimite où elle
n'apporte rien.

## 1. Résumé

Deux modifications indépendantes du screen proposé, chacune motivée par un
diagnostic mesuré, et qui se composent :

1. **Agrégation de rangs sur une grille de réglages** — classer par le rang
   médian du score sur une grille 3×3 en (a, b) au lieu d'un réglage unique.
2. **Présélection par quantile SIS** — restreindre le classement final aux d₁
   coordonnées les mieux classées par quantile SIS.

Ensemble, à n=2000, p=1000, ρ=0.25, elles font passer la récupération exacte
de 0.200 à 0.500 sur M2 et de 0.775 à 0.950 sur M1, et rendent la méthode
égale ou supérieure à quantile SIS sur les deux critères dans les quatre
modèles. À n=5000, p=2000, le gain se concentre sur M2 (0.975 contre 0.300
pour quantile SIS) et l'agrégation de réglages devient inutile.

## 2. La méthode

Pour un échantillon (X, Y) de taille n en dimension p, et un budget de sortie d :

```
Étage 1 — présélection
  pour chaque coordonnée j :
      w_j = utilité quantile SIS à tau = 0.95
            (B-splines cubiques, 3 degrés de liberté, He-Wang-Hong 2013)
  S = les d1 coordonnées de plus grande utilité w_j          (d1 = 25 par défaut)

Étage 2 — classement par indice de queue, agrégé sur les réglages
  pour chaque (a, b) dans {0.30, 0.35, 0.40} x {0.05, 0.10, 0.15} :
      alpha = n^-a ;  h = n^-b / 2
      s_j(a,b) = score de Hill local (score_coordinate_cpp) pour tout j
      r_j(a,b) = rang de s_j(a,b) par ordre CROISSANT
  r_j = médiane des r_j(a,b) sur les 9 réglages
  sortie = les d coordonnées de S de plus petit r_j
```

Le score de l'étage 2 est marginal : restreindre à S ne change aucun score, il
change seulement la compétition. Les deux étages sont donc indépendants et
peuvent s'utiliser séparément.

## 3. Pourquoi chacun des deux marche

Les deux corrigent des causes distinctes, ce qui explique qu'ils se composent
au lieu de se recouvrir.

**L'étage 1 abaisse le seuil de valeur extrême.** Le classement se fait par
score croissant, donc chaque coordonnée active doit battre le *minimum* des
p−4 scores nuls. Pour p = 1000 ce minimum se situe à environ 3.09 écarts-types
sous la médiane des nuls, et 3.29 pour p = 2000. Or la séparation réalisée,
mesurée sur M2 à n=5000 (`code/R/diag_tuning_snr.R`), vaut 3.10 au réglage du
papier et au mieux 3.50 sur toute la grille de réglages : **le screen travaille
sur son seuil de détection**, d'où des échecs tout ou rien. Ramener le champ à
d₁ = 25 abaisse ce seuil à 1.92 écart-type, soit une marge confortable.

**L'étage 2 moyenne le bruit de réglage.** Quand le score à réglage unique
échoue, seuls 33 % à 87 % des neuf réglages échouent avec lui
(`code/R/test_stability.R`) : dans une fraction substantielle des réplications
ratées, un autre réglage réussit. Le rang médian récupère cette information.
Cette part de la bimodalité de R_max est donc du bruit d'estimation, et non une
propriété du jeu de données.

**Ce qui ne marche pas, et pourquoi.** Le bagging sur demi-échantillons dégrade
tout (Sure-4 de 0.775 à 0.550 sur M1) : le score dépend trop fortement de n —
la campagne donne 0.282 de Sure-20 à n=1000 contre 0.654 à n=2000 — si bien que
les demi-échantillons échouent 54 % à 92 % du temps *même quand le score complet
réussit*. Agréger des classements faibles donne un classement faible.

## 4. Résultats

40 réplications par cellule, graines de la campagne, donc comparaisons
appariées. Erreur type d'une probabilité : jusqu'à 0.079.

### n = 2000, p = 1000, rho = 0.25 — Sure-4 / Sure-20

| méthode | M1 | M2 | M3 | M4 |
|:---|---:|---:|---:|---:|
| screen proposé, réglage unique | 0.775 / 0.950 | 0.200 / 0.550 | 0.175 / 0.500 | 0.050 / 0.325 |
| screen proposé + 9 réglages | 0.900 / 1.000 | 0.375 / 0.700 | 0.225 / 0.575 | 0.150 / 0.450 |
| quantile SIS tau=.95 | 0.800 / 1.000 | 0.225 / 0.975 | 0.300 / 0.750 | 0.075 / 0.700 |
| deux étages, d1=25 | 0.875 / 1.000 | 0.275 / 0.925 | 0.350 / 0.825 | 0.150 / 0.725 |
| **deux étages + 9 réglages** | **0.950 / 1.000** | **0.500 / 0.975** | **0.350 / 0.825** | 0.125 / **0.800** |

C'est la seule règle testée qui égale ou dépasse quantile SIS sur les deux
critères à la fois dans les quatre modèles. Sur M2 elle double la récupération
exacte à Sure-20 identique.

### n = 5000, p = 2000, rho = 0.25 — Sure-4 / Sure-20

| méthode | M1 | M2 | M3 | M4 |
|:---|---:|---:|---:|---:|
| screen proposé, réglage unique | 1.000 / 1.000 | 0.900 / 0.975 | 0.675 / 0.850 | 0.400 / 0.700 |
| screen proposé + 9 réglages | 1.000 / 1.000 | 0.950 / 1.000 | 0.650 / 0.800 | 0.400 / 0.700 |
| quantile SIS tau=.95 | 1.000 / 1.000 | 0.300 / 1.000 | **1.000 / 1.000** | **0.900 / 1.000** |
| deux étages, d1=25 | 1.000 / 1.000 | 0.925 / 1.000 | 0.800 / 0.975 | 0.575 / 1.000 |
| **deux étages + 9 réglages** | 1.000 / 1.000 | **0.975 / 1.000** | 0.775 / 1.000 | 0.525 / 0.975 |

**Le gain est régime-dépendant.** À n=5000 l'agrégation de réglages n'apporte
plus rien — le SNR y vaut 4.45 à 9.46 contre 2.34 à 5.02 à n=2000, il n'y a
plus de bruit de réglage à moyenner — et quantile SIS seul redevient meilleur
sur M3 et M4. Le seul modèle où la méthode combinée domine à toutes les tailles
est M2, celui qui porte des variables d'échelle.

## 5. Coût

Par réplication à n=2000, p=1000, monocœur : le score de Hill local coûte
0.41 s, quantile SIS à un tau 3.35 s, la génération et les rangs 0.55 s. La
méthode combinée coûte donc 9 × 0.41 + 3.35 + 0.55 ≈ 7.6 s, soit **1.9 fois
quantile SIS seul** et 7.9 fois le screen proposé seul. L'étage 2 pourrait
n'être calculé que sur les d₁ survivants, ce qui ramènerait son coût à
2 % de sa valeur ici — non implémenté, les scores étant calculés sur les p
coordonnées pour permettre la comparaison avec le classement non restreint.

L'argument « implémentation en O(p n log n) » du screen proposé ne survit pas à
l'étage 1, qui est dominé par p régressions quantiles.

## 6. Limites

- **d₁ est un paramètre de réglage supplémentaire**, dont l'optimum dépend du
  modèle. La méthode ne peut jamais récupérer une coordonnée que l'étage 1 a
  écartée : la probabilité que les quatre actives survivent est un plafond dur,
  mesurée entre 0.800 et 1.000 selon le modèle à d₁ = 25.
- **Mesuré sur deux cellules seulement** (n=2000/p=1000 et n=5000/p=2000, à
  rho=0.25), 40 réplications. Les écarts inférieurs à 0.10 ne sont pas
  interprétables : le 0.125 contre 0.150 de M4 est du bruit.
- **La grille de réglages n'a pas été optimisée.** Neuf points ont été choisis
  autour du réglage du papier ; ni sa taille ni son étendue n'ont été réglées.
- **Le Sure-20 reste en retrait sur M2 à n=2000** pour le screen agrégé seul
  (0.700 contre 0.975) : l'agrégation réduit la queue de R_max sans la
  supprimer.

## 7. Reproduction

```sh
KAPPA=0.20 Rscript code/R/test_stability.R 88 \
  results/grid/comparison_cells/g_M*_n2000_p1000_r025.rds
KAPPA=0.20 Rscript code/R/test_pipeline_agg.R 88 \
  results/grid/comparison_cells/g_M*_n2000_p1000_r025.rds
KAPPA=0.20 Rscript code/R/test_pipeline_agg.R 88 \
  results/grid/comparison_cells/g_M*_n5000_p2000_r025.rds
```

Les réplications sont rejouées depuis les graines stockées dans les cellules de
la campagne ; `code/R/test_pipeline.R` vérifie par assertion que le rejeu
reproduit exactement les rangs enregistrés.

## 8. Constructions écartées

Cinq autres pistes ont été testées sur les mêmes données et n'ont pas été
retenues ; le détail est dans les journaux de `results/`.

| construction | verdict |
|:---|:---|
| contraste de l'enveloppe au lieu du niveau | gagne sur M2 (0.88 contre 0.50 de Sure-20), s'effondre sur M3 (0.00) |
| pente de régression quantile extrême | immunité à l'échelle confirmée, variance rédhibitoire (0.175 de Sure-20 sur M2) |
| régression sur les dépassements | viable (0.925 sur M1) mais jamais meilleure ; même compromis biais-variance reparamétré |
| Hill à biais réduit de second ordre | gagne sur M2 (0.450 contre 0.200), dégrade M3 et M4 ; corrige un biais de second ordre là où l'obstacle est un biais de mélange |
| union / intersection des deux écrans | aucun gain à taille de sortie égale |
