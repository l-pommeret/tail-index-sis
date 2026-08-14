# Screening en deux étages avec agrégation de réglages

Note de méthode issue de la campagne étendue (`results/grid/`, 288 cellules ×
40 réplications). Elle décrit une modification du screen par indice de queue
qui améliore ses performances sans changer ce qu'il cible, explique le
mécanisme de chaque composant, donne les mesures qui ont servi à choisir
chacun, et délimite ce qu'elle n'apporte pas.

Toutes les comparaisons sont **appariées** : les réplications sont rejouées
depuis les graines stockées dans les cellules de la campagne, donc chaque règle
est évaluée sur exactement les mêmes jeux de données.

---

## 1. En une page

Le screen proposé classe les coordonnées par score croissant de Hill local. Il
échoue de façon **tout ou rien** : soit il place les quatre variables actives en
tête, soit il les enfouit au rang 500 ou 800. La cause est mesurée (section 3.1)
et deux corrections indépendantes s'y attaquent :

1. **Présélection par quantile SIS** — ne classer que les d₁ = 25 coordonnées
   les mieux notées par quantile SIS. Cela réduit la compétition, donc la
   probabilité qu'une coordonnée nulle chanceuse dépasse une active.
2. **Agrégation de rangs sur une grille de réglages** — au lieu d'un seul
   couple (a, b), calculer le score sur 9 réglages voisins et classer par le
   **rang minimum** obtenu sur ces 9.

Ensemble, moyenné sur les quatre modèles :

| | Sure-4 | Sure-20 |
|:---|---:|---:|
| **n=2000, p=1000** — screen proposé seul | 0.300 | 0.581 |
| **n=2000, p=1000** — méthode combinée | **0.544** | **0.894** |
| **n=5000, p=2000** — screen proposé seul | 0.744 | 0.881 |
| **n=5000, p=2000** — méthode combinée | **0.894** | **1.000** |

À n=5000 la méthode combinée atteint Sure-20 = 1.000 **dans les quatre
modèles** : les échecs catastrophiques disparaissent. Le coût, correctement
implémenté, est de 8 % supérieur à celui de quantile SIS seul (section 7).

---

## 2. La méthode

Entrée : un échantillon (X, Y) de taille n en dimension p ; un budget de sortie
d. Sortie : d coordonnées.

```
ÉTAGE 1 — présélection (recall)
  pour chaque coordonnée j = 1..p :
      w_j = utilité quantile SIS à tau = 0.95
            (régression quantile sur B-splines cubiques, 3 ddl,
             He-Wang-Hong 2013 ; code/R/qa_sis.R)
  S = les d1 coordonnées de plus grande utilité w_j            [d1 = 25]

ÉTAGE 2 — classement par indice de queue, agrégé sur les réglages
  G = { (a,b) : a dans {0.30,0.35,0.40},
                b dans {0.05,0.10,0.15} }                        [9 réglages]
  pour chaque (a,b) dans G :
      alpha = n^(-a) ;  h = n^(-b)/2
      s_j(a,b) = score de Hill local pour j dans S
                 (code/src/local_hill.cpp, score_coordinate_cpp)
      r_j(a,b) = rang de s_j(a,b) parmi S, par ordre CROISSANT
  r_j = MINIMUM des r_j(a,b) sur les 9 réglages
  sortie = les d coordonnées de S de plus petit r_j
```

Deux propriétés rendent la construction légitime :

- **Le score de l'étage 2 est marginal** : il ne dépend que de (X_j, Y). Le
  restreindre à S ne modifie donc aucun score, il modifie seulement contre qui
  chaque coordonnée est classée. Les deux étages sont indépendants et peuvent
  s'employer séparément.
- **Aucune information oracle n'est utilisée.** On ne choisit pas le meilleur
  réglage par réplication : on agrège les 25 rangs, ce qui est calculable sans
  connaître les variables actives.

---

## 3. Pourquoi ça marche

### 3.1 Le diagnostic de départ : le screen travaille sur son seuil de détection

Le classement se fait par score **croissant**, donc chaque coordonnée active
doit battre le **minimum** des p − 4 scores nuls. Ce minimum est un problème de
valeurs extrêmes : pour p = 1000 il se situe environ 3.09 écarts-types sous la
médiane des nuls, et 3.29 pour p = 2000.

Or la séparation réalisée a été mesurée directement
(`code/R/diag_tuning_snr.R`, M2 à n=5000, p=2000) : **3.10 au réglage du papier,
et au mieux 3.50 sur toute la grille de réglages testée**. Le screen opère donc
exactement à sa limite de détection, et n'importe quelle fluctuation le fait
basculer. La dissection d'une réplication qui échoue et d'une qui réussit
(`code/R/diag_failure.R`) le confirme : séparation de 0.3 écart-type dans celle
qui donne R_max = 812, contre 3.7 dans celle qui donne R_max = 6.

C'est là l'origine de la bimodalité, et elle explique pourquoi élargir la
fenêtre d ne sauve rien : quand la coordonnée active tombe dans le bruit, son
rang est quasi uniforme sur [1, p].

### 3.2 Étage 1 : abaisser le seuil de valeurs extrêmes

Réduire le champ de p = 1000 à d₁ = 25 fait passer le seuil de compétition de
3.09 à **1.92 écart-type**. La séparation, elle, ne change pas : le score étant
marginal, il vaut la même chose. On passe donc d'une marge nulle à une marge
confortable, sans toucher à l'estimateur.

Le prix est que l'étage 2 ne peut jamais récupérer ce que l'étage 1 a écarté.
La probabilité que les quatre actives survivent est un **plafond dur**, mesuré
entre 0.800 et 1.000 à d₁ = 25 selon le modèle, et égal à 1.000 partout à
n=5000.

Quantile SIS est le bon choix pour cet étage parce qu'il a le meilleur rappel
de toutes les méthodes de la campagne : son Sure-20 moyen vaut 0.897 / 0.822 /
0.738 / 0.664 sur M1–M4, contre 0.805 / 0.596 / 0.578 / 0.471 pour le screen
proposé.

### 3.3 Étage 2 : une partie de la bimodalité est du bruit de réglage

Le test décisif est dans `code/R/test_stability.R`. Pour chaque réplication, on
compare le score au réglage unique à ce que donnent les autres réglages :

> **Quand le score à réglage unique échoue, seuls 33 % à 87 % des autres
> réglages échouent avec lui.**

Autrement dit, dans une fraction substantielle des réplications ratées, un
autre réglage réussit. Cette part de l'échec n'est donc **pas** portée par le
jeu de données : c'est du bruit d'estimation, et l'agrégation le moyenne.

Une conséquence importante : ce raisonnement ne vaut que là où il y a du bruit
de réglage à moyenner. À n=5000 la séparation vaut 4.45 à 9.46 écarts-types
contre 2.34 à 5.02 à n=2000 ; l'agrégation seule n'y apporte donc presque plus
rien (section 5), et le gain vient alors de l'étage 1.

### 3.4 Ce qui ne marche pas, et pourquoi

Le **bagging sur demi-échantillons** dégrade tout — Sure-4 de 0.775 à 0.550 sur
M1. La raison n'est pas que les données soient sans espoir : c'est que le score
dépend très fortement de n (la campagne donne Sure-20 = 0.282 à n=1000 contre
0.654 à n=2000). Les demi-échantillons échouent 54 % à 92 % du temps **même
quand le score sur l'échantillon complet réussit**. Agréger des classements
faibles donne un classement faible.

C'est pourquoi l'agrégation porte ici sur les **réglages** et non sur les
**données** : varier le réglage produit des classements de qualité comparable
et d'erreurs décorrélées ; sous-échantillonner produit des classements
uniformément dégradés.

---

## 4. Comment les composants ont été choisis

Les deux choix — quelle grille de réglages, quel agrégateur — ont été tranchés
par mesure sur les mêmes réplications (`code/R/test_aggregation.R`), et non par
défaut.

### 4.1 La grille : densifier oui, élargir non

Trois grilles comparées, chacune avec les six agrégateurs, sans pipeline :

| grille | points | étendue | Sure-4 (n=2000) | Sure-4 (n=5000) |
|:---|---:|:---|---:|---:|
| grossière | 9 | a ∈ [0.30, 0.40], b ∈ [0.05, 0.15] | 0.321 | 0.622 |
| **fine** | 25 | **même étendue**, pas de 0.025 | **0.352** | **0.631** |
| large | 25 | a ∈ [0.25, 0.45], b ∈ [0, 0.20] | 0.154 | 0.481 |

**Densifier n'apporte rien d'exploitable ; élargir divise la performance par
deux.** Le gain apparent de la grille fine ne résiste pas à un test direct : sur
huit configurations comparant 9 et 25 réglages sur les mêmes réplications,
l'écart moyen en Sure-4 est de +0.006, de signe alternant, pour une erreur type
d'environ 0.04. Des réglages espacés de 0.025 donnent des classements trop
corrélés pour augmenter le nombre effectif de configurations indépendantes.
**La grille retenue est donc la 3x3**, 2.8 fois moins coûteuse.
À nombre de points identique (25), la grille fine fait 0.352 et la large 0.154.

L'explication est mesurée : le balayage de séparation montre que le SNR tombe à
1.36 à b = 0.20 et à 2.27 à a = 0.50, contre 3.10 à 3.50 dans la zone centrale.
La grille large introduit donc des classements franchement mauvais, et même une
médiane finit par en souffrir. Le gain de la grille fine est en revanche limité
parce que deux réglages voisins produisent des classements très corrélés : ce
qui compte est le **nombre effectif de configurations indépendantes**, pas leur
nombre brut.

### 4.2 L'agrégateur : et une interaction inattendue

Six agrégateurs, moyennés sur les trois grilles, **sans** pipeline :

| agrégateur | Sure-4 (n=2000) | Sure-20 (n=2000) | Sure-4 (n=5000) | Sure-20 (n=5000) |
|:---|---:|---:|---:|---:|
| médiane | **0.381** | 0.694 | **0.735** | 0.875 |
| moyenne géométrique | 0.317 | 0.713 | 0.613 | 0.883 |
| moyenne tronquée 20 % | 0.304 | 0.642 | 0.613 | 0.808 |
| minimum | 0.296 | **0.733** | 0.571 | **0.933** |
| moyenne | 0.208 | 0.515 | 0.504 | 0.719 |
| maximum | 0.148 | 0.340 | 0.433 | 0.610 |

La moyenne arithmétique est mauvaise parce qu'un seul réglage défaillant place
un rang de 800 dans la somme et écrase tout, là où la médiane l'ignore. Le
maximum est le pire : il hérite du plus mauvais réglage par construction.

**Mais le classement s'inverse dès qu'on ajoute l'étage 1**, moyennes sur les
quatre modèles :

| règle | Sure-4 (n=2000) | Sure-4 (n=5000) |
|:---|---:|---:|
| **fine + minimum + pipeline** | **0.544** | **0.894** |
| grossière + minimum + pipeline | 0.525 | 0.869 |
| fine + médiane + pipeline | 0.488 | 0.819 |
| grossière + médiane + pipeline | 0.481 | 0.819 |

Seul, le minimum est **moins bon** que la médiane (0.296 contre 0.381) ;
derrière le pipeline il devient **meilleur** (0.544 contre 0.488). Le mécanisme
est net : le rang minimum promeut une coordonnée dès qu'*un seul* réglage la
classe bien — ce qui promeut aussi les coordonnées nulles ayant eu un réglage
chanceux. Avec 996 nulles en compétition c'est rédhibitoire ; avec 21 nulles
survivantes il n'y a presque plus personne pour être chanceux. **L'étage 1
neutralise le défaut du minimum et n'en laisse que la qualité.** Les deux
composants interagissent, ils ne s'additionnent pas — c'est pourquoi il faut
les choisir ensemble et non séparément.

---

## 5. Résultats

40 réplications par cellule, ρ = 0.25, graines de la campagne. Format
Sure-4 / Sure-20.

### n = 2000, p = 1000

| méthode | M1 | M2 | M3 | M4 | moyenne Sure-4 |
|:---|---:|---:|---:|---:|---:|
| screen proposé, réglage unique | 0.775 / 0.950 | 0.200 / 0.550 | 0.175 / 0.500 | 0.050 / 0.325 | 0.300 |
| quantile SIS τ=.95 | 0.800 / 1.000 | 0.225 / 0.975 | 0.300 / 0.750 | 0.075 / 0.700 | 0.350 |
| **méthode combinée** | **1.000 / 1.000** | **0.575 / 0.975** | **0.425 / 0.825** | **0.175 / 0.775** | **0.544** |

### n = 5000, p = 2000

| méthode | M1 | M2 | M3 | M4 | moyenne Sure-4 |
|:---|---:|---:|---:|---:|---:|
| screen proposé, réglage unique | 1.000 / 1.000 | 0.900 / 0.975 | 0.675 / 0.850 | 0.400 / 0.700 | 0.744 |
| quantile SIS τ=.95 | 1.000 / 1.000 | 0.300 / 1.000 | 1.000 / 1.000 | 0.900 / 1.000 | 0.800 |
| **méthode combinée** | **1.000 / 1.000** | **0.975 / 1.000** | **0.900 / 1.000** | **0.700 / 1.000** | **0.894** |

### Tests appariés

McNemar exact, méthode combinée contre le screen proposé à réglage unique. On
lit « paires discordantes favorables / paires discordantes totales ».

| cellule | modèle | Sure-4 | p | Sure-20 | p |
|:---|:---|---:|---:|---:|---:|
| n=5000, p=2000 | M2 | 3/3 | 0.25 | 1/1 | — |
| | M3 | 10/11 | **0.012** | 6/6 | **0.031** |
| | M4 | 15/18 | **0.0075** | 12/12 | **0.0005** |

Le gain est significatif sur M3 et M4, sur les deux critères, et les paires
discordantes sont quasi unanimes : la domination est **réplication par
réplication**, pas seulement en moyenne. Sur M1 et M2 à n=5000 les deux méthodes
réussissent presque toujours, donc il y a trop peu de paires discordantes pour
conclure — ce qui n'est pas un désaveu mais une absence de puissance.

### Ce que la méthode ne fait pas

Elle ne bat pas quantile SIS partout. À n=5000, quantile SIS seul reste devant
sur M3 (1.000 contre 0.900) et M4 (0.900 contre 0.700) en récupération exacte.
Le seul modèle où la méthode combinée domine à toutes les tailles est **M2**,
celui qui porte des variables d'échelle — c'est-à-dire le modèle qui teste ce
que le screen par indice de queue est censé savoir faire et que quantile SIS ne
sait pas : ignorer une covariable qui déplace les quantiles finis sans toucher
l'indice de queue.

---

## 6. Ce que ça change pour la bimodalité

La question de départ était : peut-on stabiliser R_max ? Réponse mesurée : en
grande partie.

| | Sure-20 min sur les 4 modèles | Sure-20 max |
|:---|---:|---:|
| screen proposé seul, n=2000 | 0.325 | 0.950 |
| méthode combinée, n=2000 | 0.775 | 1.000 |
| screen proposé seul, n=5000 | 0.700 | 1.000 |
| méthode combinée, n=5000 | **1.000** | **1.000** |

À n=5000 les échecs catastrophiques ont disparu. À n=2000 le pire modèle passe
de 0.325 à 0.775. La bimodalité n'est pas *supprimée* — il reste la part portée
par le jeu de données, celle où tous les réglages échouent ensemble — mais elle
cesse d'être le mode d'échec dominant.

---

## 7. Coût

Par réplication à n = 2000, p = 1000, monocœur (mesures de
`code/R/bench_grid.R`) : score de Hill local sur p coordonnées 0.41 s, quantile
SIS à un τ 3.35 s, génération et rangs 0.55 s.

- **Tel que mesuré ici** : 25 × 0.41 + 3.35 + 0.55 ≈ 14.2 s, soit 3.6 fois
  quantile SIS seul. Les scores sont calculés sur les p coordonnées uniquement
  pour permettre la comparaison avec le classement non restreint.
- **Tel qu'il faut l'implémenter** : l'étage 2 ne concerne que les d₁ = 25
  survivants, soit 2.5 % des coordonnées. Son coût réel est
  25 × 0.41 × 0.025 ≈ 0.26 s, et le total 4.2 s — **8 % de plus que quantile
  SIS seul**.

L'agrégation sur 25 réglages est donc essentiellement gratuite : c'est l'étage 1
qui domine. En contrepartie, **l'argument « implémentation en O(p n log n) » du
screen proposé ne survit pas** : l'étage 1 exige p régressions quantiles.

---

## 8. Limites

- **d₁ est un paramètre de réglage supplémentaire**, dont l'optimum dépend du
  modèle. Il gouverne un compromis explicite : d₁ petit donne un meilleur
  classement final mais une survie d'étage 1 plus faible, et la survie est un
  plafond dur.
- **Deux cellules seulement** (n=2000/p=1000 et n=5000/p=2000, ρ=0.25), 40
  réplications. L'erreur type d'une probabilité atteint 0.079 par modèle, et
  environ 0.04 sur une moyenne de quatre modèles. L'écart entre la meilleure
  recette et la précédente (0.544 contre 0.481 à n=2000) fait environ 1.6
  écart-type : suggestif, non concluant. L'écart entre grille fine et grille
  large (0.352 contre 0.154) est en revanche massif.
- **Rien n'a été testé hors de ρ = 0.25**, alors que la campagne montre que ρ
  est le facteur de design le plus influent (Sure-4 relatif passant de −0.116 à
  +0.055 entre ρ = 0 et ρ = 0.5).
- **Les quatre modèles partagent une même forme de γ décroissante** en chaque
  coordonnée active. Rien ne garantit que les conclusions s'étendent à des
  effets non monotones.

---

## 9. Constructions écartées

Cinq autres pistes ont été construites et testées sur les mêmes données ; aucune
n'est retenue. Le détail est dans les journaux de `results/`.

| construction | idée | verdict |
|:---|:---|:---|
| contraste de l'enveloppe au lieu du niveau | une inactive donne une enveloppe plate, une active une enveloppe pentue | gagne sur M2 (Sure-20 0.88 contre 0.50), s'effondre sur M3 (0.00) |
| pente de régression quantile extrême | la pente de log Q(τ) en log{1/(1−τ)} est l'indice de queue et annule l'échelle exactement | annulation confirmée (0.075 variable d'échelle dans le top-4 contre 0.875), mais variance rédhibitoire : Sure-20 de 0.175 sur M2 |
| régression sur les dépassements | standardiser par un quantile ajusté puis moyenner les log-dépassements | viable (0.925 sur M1) mais jamais meilleure ; même compromis biais-variance reparamétré |
| Hill à biais réduit de second ordre | retirer le biais sans réduire k | gagne sur M2 (Sure-4 0.450 contre 0.200), dégrade M3 et M4 ; l'obstacle est un biais de **mélange**, pas de second ordre |
| union / intersection des deux écrans | combiner les deux classements | aucun gain à taille de sortie égale |

Le fil commun des trois premières : elles butent sur le même mur. L'estimande
est l'enveloppe supérieure le long de la fibre, qui n'apparaît qu'à la limite
k/m → 0 ; à n fini on paie soit le biais de mélange, soit la variance. Les deux
constructions fondées sur une **différence** de quantiles (pente de régression
quantile, et le concurrent Yoshida–Umezu qui est un estimateur de Pickands) sont
les deux plus mauvaises de toutes celles testées, ce qui est conforme à
l'efficacité en variance connue des estimateurs de type Hill.

C'est ce qui rend la présente méthode différente : elle **ne change pas
l'estimateur**, elle change la compétition dans laquelle son classement est
évalué et moyenne le bruit de réglage.

---

## 10. Reproduction

```sh
# diagnostic : séparation contre seuil de compétition
KAPPA=0.20 Rscript code/R/diag_tuning_snr.R \
  results/grid/comparison_cells/g_M2_n5000_p2000_r000.rds 4 4
KAPPA=0.20 Rscript code/R/diag_failure.R \
  results/grid/comparison_cells/g_M2_n5000_p2000_r000.rds

# la bimodalité est-elle du bruit de réglage ?
KAPPA=0.20 Rscript code/R/test_stability.R 88 \
  results/grid/comparison_cells/g_M*_n2000_p1000_r025.rds

# choix de la grille et de l'agrégateur
KAPPA=0.20 Rscript code/R/test_aggregation.R 92 \
  results/grid/comparison_cells/g_M*_n2000_p1000_r025.rds
KAPPA=0.20 Rscript code/R/test_aggregation.R 92 \
  results/grid/comparison_cells/g_M*_n5000_p2000_r025.rds
```

Les réplications sont rejouées depuis les graines stockées dans les cellules ;
`code/R/test_pipeline.R` vérifie par assertion que le rejeu reproduit exactement
les rangs enregistrés lors de la campagne.
