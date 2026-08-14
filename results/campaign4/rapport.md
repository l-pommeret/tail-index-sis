# Campagne de comparaison sur la suite de modèles figée

Généré le 2026-08-14 21:07 depuis `results/campaign4`.

## 1. Protocole

- **Modèles** : M1, M2, M3, M4 (`code/R/generate4.R`). Ensemble actif en indice de queue A_gamma = {1,2,3,4} dans les quatre ; ensemble actif en échelle A_scale = {5,…,24} pour M2 seul, dont les coordonnées déplacent les quantiles conditionnels finis sans toucher l'indice de queue.
- **Design** : n = 2000, p ∈ {500, 1000, 2000}, ρ ∈ {0, 0.25}.
- **Réplications** : 100 par cellule et par modèle, soit 2400 jeux de données. Erreur type d'une probabilité : au plus 0.050.
- **Comparaisons appariées** : les huit règles sont évaluées sur exactement les mêmes jeux de données.

Règles comparées :

| # | règle | description |
|---|:---|:---|
| 1 | Screen proposé | score de Hill local au réglage publié (a\*, b\*) = (0.30, 0.10) |
| 2 | Screen agrégé | rang **minimum** du score sur la grille 3×3, a ∈ {0.30,0.35,0.40} × b ∈ {0.05,0.10,0.15} |
| 3 | Yoshida–Umezu | screening de Pickands conditionnel, réglage publié h = 1, k = ⌊0.072n⌋ |
| 4–7 | Quantile SIS | B-splines cubiques à 3 ddl, τ = 0.90, 0.95, 0.975, 0.99 |
| 8 | Pipeline | quantile SIS (τ=0.95) top-25, puis rang minimum sur la grille 3×3 |

Critère : **Sure-d** = probabilité que les quatre variables actives en indice de queue soient toutes classées parmi les d meilleures. d = 4 est la récupération exacte.

## 2. Résultats à ρ = 0.00

Format Sure-4 / Sure-20.

### M1

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.310 / 0.850 | 0.200 / 0.580 | 0.180 / 0.520 |
| **Screen, 9 réglages, rang min** | 0.420 / 0.920 | 0.220 / 0.750 | 0.110 / 0.740 |
| Yoshida–Umezu | 0.020 / 0.110 | 0.010 / 0.040 | 0.030 / 0.040 |
| Quantile SIS τ=.90 | 0.020 / 0.410 | 0.010 / 0.200 | 0.010 / 0.080 |
| Quantile SIS τ=.95 | 0.500 / 0.950 | 0.270 / 0.850 | 0.190 / 0.810 |
| Quantile SIS τ=.975 | 0.230 / 0.890 | 0.160 / 0.780 | 0.080 / 0.690 |
| Quantile SIS τ=.99 | 0.010 / 0.470 | 0.010 / 0.250 | 0.000 / 0.070 |
| **Pipeline + 9 réglages min** | 0.550 / 0.970 | 0.430 / 0.910 | 0.360 / 0.830 |

### M2

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.020 / 0.210 | 0.000 / 0.100 | 0.000 / 0.060 |
| **Screen, 9 réglages, rang min** | 0.020 / 0.380 | 0.000 / 0.230 | 0.000 / 0.160 |
| Yoshida–Umezu | 0.000 / 0.010 | 0.000 / 0.000 | 0.000 / 0.000 |
| Quantile SIS τ=.90 | 0.000 / 0.000 | 0.000 / 0.010 | 0.000 / 0.000 |
| Quantile SIS τ=.95 | 0.000 / 0.030 | 0.000 / 0.010 | 0.000 / 0.020 |
| Quantile SIS τ=.975 | 0.000 / 0.160 | 0.000 / 0.130 | 0.000 / 0.130 |
| Quantile SIS τ=.99 | 0.000 / 0.180 | 0.000 / 0.030 | 0.000 / 0.010 |
| **Pipeline + 9 réglages min** | 0.060 / 0.430 | 0.060 / 0.410 | 0.050 / 0.280 |

### M3

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.020 / 0.150 | 0.010 / 0.080 | 0.000 / 0.060 |
| **Screen, 9 réglages, rang min** | 0.020 / 0.270 | 0.010 / 0.130 | 0.000 / 0.030 |
| Yoshida–Umezu | 0.000 / 0.050 | 0.000 / 0.040 | 0.000 / 0.010 |
| Quantile SIS τ=.90 | 0.000 / 0.030 | 0.000 / 0.020 | 0.000 / 0.030 |
| Quantile SIS τ=.95 | 0.040 / 0.490 | 0.020 / 0.300 | 0.000 / 0.160 |
| Quantile SIS τ=.975 | 0.000 / 0.320 | 0.010 / 0.100 | 0.000 / 0.020 |
| Quantile SIS τ=.99 | 0.000 / 0.040 | 0.000 / 0.010 | 0.000 / 0.000 |
| **Pipeline + 9 réglages min** | 0.080 / 0.550 | 0.040 / 0.400 | 0.000 / 0.170 |

### M4

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.160 / 0.640 | 0.130 / 0.440 | 0.070 / 0.480 |
| **Screen, 9 réglages, rang min** | 0.180 / 0.800 | 0.160 / 0.730 | 0.050 / 0.620 |
| Yoshida–Umezu | 0.030 / 0.090 | 0.020 / 0.090 | 0.010 / 0.070 |
| Quantile SIS τ=.90 | 0.030 / 0.310 | 0.020 / 0.240 | 0.000 / 0.140 |
| Quantile SIS τ=.95 | 0.330 / 0.860 | 0.250 / 0.830 | 0.130 / 0.760 |
| Quantile SIS τ=.975 | 0.150 / 0.880 | 0.080 / 0.760 | 0.050 / 0.570 |
| Quantile SIS τ=.99 | 0.000 / 0.310 | 0.000 / 0.110 | 0.000 / 0.040 |
| **Pipeline + 9 réglages min** | 0.410 / 0.890 | 0.400 / 0.860 | 0.270 / 0.810 |

## 2. Résultats à ρ = 0.25

Format Sure-4 / Sure-20.

### M1

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.820 / 1.000 | 0.760 / 0.930 | 0.690 / 0.940 |
| **Screen, 9 réglages, rang min** | 0.880 / 0.990 | 0.830 / 0.970 | 0.800 / 0.970 |
| Yoshida–Umezu | 0.220 / 0.430 | 0.130 / 0.350 | 0.170 / 0.330 |
| Quantile SIS τ=.90 | 0.290 / 0.910 | 0.210 / 0.790 | 0.170 / 0.730 |
| Quantile SIS τ=.95 | 0.900 / 1.000 | 0.830 / 0.990 | 0.740 / 1.000 |
| Quantile SIS τ=.975 | 0.700 / 1.000 | 0.660 / 1.000 | 0.500 / 0.960 |
| Quantile SIS τ=.99 | 0.110 / 0.850 | 0.050 / 0.570 | 0.010 / 0.470 |
| **Pipeline + 9 réglages min** | 0.960 / 1.000 | 0.890 / 1.000 | 0.900 / 1.000 |

### M2

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.360 / 0.690 | 0.290 / 0.740 | 0.230 / 0.570 |
| **Screen, 9 réglages, rang min** | 0.480 / 0.970 | 0.480 / 0.890 | 0.330 / 0.780 |
| Yoshida–Umezu | 0.010 / 0.100 | 0.000 / 0.070 | 0.010 / 0.030 |
| Quantile SIS τ=.90 | 0.000 / 0.050 | 0.000 / 0.040 | 0.000 / 0.040 |
| Quantile SIS τ=.95 | 0.060 / 0.640 | 0.010 / 0.570 | 0.040 / 0.630 |
| Quantile SIS τ=.975 | 0.190 / 0.930 | 0.130 / 0.870 | 0.170 / 0.840 |
| Quantile SIS τ=.99 | 0.050 / 0.790 | 0.050 / 0.490 | 0.010 / 0.260 |
| **Pipeline + 9 réglages min** | 0.740 / 0.980 | 0.690 / 0.950 | 0.600 / 0.930 |

### M3

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.210 / 0.560 | 0.170 / 0.500 | 0.150 / 0.440 |
| **Screen, 9 réglages, rang min** | 0.190 / 0.790 | 0.170 / 0.620 | 0.110 / 0.580 |
| Yoshida–Umezu | 0.050 / 0.190 | 0.010 / 0.140 | 0.000 / 0.040 |
| Quantile SIS τ=.90 | 0.150 / 0.650 | 0.040 / 0.410 | 0.040 / 0.350 |
| Quantile SIS τ=.95 | 0.430 / 0.970 | 0.340 / 0.860 | 0.290 / 0.830 |
| Quantile SIS τ=.975 | 0.160 / 0.820 | 0.030 / 0.730 | 0.030 / 0.530 |
| Quantile SIS τ=.99 | 0.010 / 0.220 | 0.000 / 0.140 | 0.000 / 0.060 |
| **Pipeline + 9 réglages min** | 0.370 / 0.980 | 0.390 / 0.870 | 0.260 / 0.850 |

### M4

| règle | p = 500 | p = 1000 | p = 2000 |
|:---|---:|---:|---:|
| Screen proposé, réglage unique | 0.750 / 0.930 | 0.690 / 0.920 | 0.750 / 0.890 |
| **Screen, 9 réglages, rang min** | 0.730 / 0.970 | 0.700 / 0.950 | 0.700 / 0.950 |
| Yoshida–Umezu | 0.310 / 0.570 | 0.210 / 0.490 | 0.160 / 0.450 |
| Quantile SIS τ=.90 | 0.350 / 0.930 | 0.280 / 0.840 | 0.280 / 0.810 |
| Quantile SIS τ=.95 | 0.840 / 0.990 | 0.780 / 0.980 | 0.730 / 1.000 |
| Quantile SIS τ=.975 | 0.550 / 0.990 | 0.470 / 0.990 | 0.390 / 0.920 |
| Quantile SIS τ=.99 | 0.070 / 0.740 | 0.040 / 0.480 | 0.010 / 0.300 |
| **Pipeline + 9 réglages min** | 0.920 / 1.000 | 0.830 / 0.980 | 0.850 / 1.000 |

## 3. Le mécanisme sur M2

M2 est le seul modèle porteur de variables d'échelle. Nombre moyen des 20 coordonnées de A_scale retenues dans le top-24, et nombre de variables actives en indice de queue dans le top-4 :

| ρ | règle | p = 500 | p = 1000 | p = 2000 |
|:---|:---|---:|---:|---:|
| ρ=0.00 | **Screen, 9 réglages, rang min** | 4.4 / 2.02 | 4.2 / 1.84 | 3.1 / 1.71 |
| ρ=0.00 | Quantile SIS τ=.95 | 19.9 / 0.19 | 19.9 / 0.10 | 19.8 / 0.13 |
| ρ=0.00 | Yoshida–Umezu | 9.7 / 0.26 | 8.6 / 0.14 | 7.0 / 0.13 |
| ρ=0.25 | **Screen, 9 réglages, rang min** | 2.9 / 3.37 | 2.0 / 3.39 | 1.4 / 3.18 |
| ρ=0.25 | Quantile SIS τ=.95 | 19.2 / 1.79 | 19.2 / 1.72 | 18.5 / 1.85 |
| ρ=0.25 | Yoshida–Umezu | 9.7 / 0.61 | 7.7 / 0.45 | 6.4 / 0.34 |

Format : |top-24 ∩ A_scale| / |top-4 ∩ A_gamma|.

Proportion moyenne de scores Yoshida–Umezu non définis (ratio de Pickands négatif ou indéterminé) : 0.000, maximum 0.000 sur une cellule.

## 4. Lecture d'ensemble

**Le pipeline est la meilleure règle en Sure-20 dans 19 des 24 cellules** (modèle × dimension × corrélation). Décompte des cellules où chaque règle est la meilleure : pipeline + 9 min 19 ; quantile .95 2 ; quantile .975 2 ; screen seul 1.

**M2 sépare les deux familles.** À ρ = 0.25, moyenné sur les trois dimensions, le screen agrégé atteint 0.430 de récupération exacte contre 0.037 pour quantile SIS à τ=.95, et le pipeline 0.677. La cause est mesurée dans la section 3 : quantile SIS retient environ 19 des 20 variables d'échelle dans son top-24, donc il ne place que 1.79 variable active en indice de queue dans son top-4 contre 3.31 pour le screen agrégé.

**L'agrégation de réglages gagne presque partout.** Le rang minimum sur neuf réglages améliore le Sure-20 du screen de 0.121 en moyenne sur les 24 cellules, avec un gain allant de -0.030 à 0.290. Il le dégrade dans 2 cellules seulement, au pire de 0.030, soit moins que l'erreur type de 0.050.

**Yoshida–Umezu est dernier dans 15 des 24 cellules.** Aucun de ses scores n'est indéfini ici, donc son retard tient à sa variance et non à des ajustements ratés : c'est un estimateur de Pickands, fondé sur des différences de quantiles, là où le score proposé moyenne des espacements logarithmiques.

**Le réglage de τ pèse plus que le choix de la famille.** L'écart moyen de Sure-20 entre τ=.95 et τ=.99 vaut 0.402, à comparer aux écarts entre méthodes. Le τ optimal n'étant pas connu en pratique, une comparaison à quantile SIS réglé à sa meilleure valeur avantage celui-ci d'un montant qui n'existe pas dans une application réelle.

**La corrélation reste le facteur de design dominant.** Passer de ρ = 0 à ρ = 0.25 fait gagner 0.389 de Sure-20 au screen agrégé, en moyenne sur modèles et dimensions. Les quatre coordonnées actives étant adjacentes dans la structure AR(1), conditionner sur l'une déplace aussi les trois autres et creuse la dépression de l'enveloppe.

## 5. Réserves

- 100 réplications donnent une erreur type d'au plus 0.050 ; les écarts inférieurs à 0.10 dans une cellule isolée restent à interpréter avec prudence.
- Une seule taille d'échantillon (n = 2000) et deux corrélations. Les conclusions ne s'étendent pas telles quelles à n plus grand, où la campagne précédente montrait que l'agrégation cesse d'apporter.
- Le pipeline introduit un paramètre supplémentaire, la taille d1 = 25 de l'ensemble présélectionné, dont l'optimum dépend du modèle ; il ne peut jamais récupérer une coordonnée que l'étage 1 a écartée.
- L'étage 1 étant p régressions quantiles, l'argument de coût en O(p n log n) du screen proposé ne s'applique pas au pipeline.

