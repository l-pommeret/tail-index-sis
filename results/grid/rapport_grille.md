# Comparaison des methodes de screening sur grille etendue

Genere le 2026-08-14 15:19 a partir de `results/grid/comparison_cells` (14 cellules sur 288).

> **Rapport partiel** : la campagne est en cours, 14 cellules sur 288 sont calculees. Les vues agregees de la section 3 ne portent que sur les cellules disponibles et bougeront encore. Les tables completes de la section 4 sont definitives cellule par cellule.

**Prochaine echeance : batch de 7 cellules attendu vers 15:32 ; fin de campagne estimee vers 17:18.**

Debit mesure 33.9 core-heures par heure de calcul (soit 34 coeurs solo-equivalents) ; 70.8 core-heures restantes sur 92.2.

## 1. Protocole

Les quatre modeles M1-M4 de la Draft 3 (`code/R/generate3.R`), les memes
estimateurs et le meme reglage que `code/R/run_draft3.R comparison`, mais sur
un espace de design elargi.

- **Modeles** : M1, M2, M3, M4 ; ensemble actif en indice de queue A_gamma = {1,2,3,4} (pour M2, les variables d'echelle A_scale = {5,6,7,8} n'agissent que sur les quantiles finis, avec kappa = 0.20).
- **Tailles** : n = 1000, 2000, 5000.
- **Dimensions** : p = 200, 500, 1000, 2000.
- **Dependance** : X ~ AR(1) gaussien de correlation rho = 0, 0.20, 0.25, 0.30, 0.40, 0.50.
- **Replications Monte Carlo** : 40 par cellule, soit 14 cellules et 560 jeux de donnees simules.

Methodes comparees (identiques a celles de la Draft 3) :

- **Tail-index SIS** (propose) : score de Hill local sur rangs empiriques,
  a* = 0.30, b* = 0.10, soit alpha = n^-a*, h = n^-b*/2, epsilon = 0.05 ; les coordonnees sont classees par score **croissant**.
- **Yoshida--Umezu** : screening de Pickands conditionnel, reglage publie h = 1, k = floor(0.072 n).
- **Quantile SIS** (He, Wang et Hong 2013) : B-splines cubiques a 3 degres de liberte, tau = 0.90, 0.95, 0.975, 0.99.

Criteres :

- **Sure-d** = probabilite que les 4 variables actives soient toutes classees parmi les d meilleures (d = 4, 10, 20, 30, 50, 100 ; d = 4 est la recuperation exacte).
- **E(Rmax)**, **Med(Rmax)** = moyenne et mediane du pire rang actif.

Erreur type Monte Carlo d'une probabilite avec 40 replications : au plus 0.079 (et 0.047 pour une probabilite de 0.9).
Flux de graines : 12000019 + cellule*10007 + r*101, disjoint de tous les flux de la Draft 3.

## 2. Controle : cellule de reference de la Draft 3

La cellule n = 2000, p = 1000, rho = 0.25 n'est pas encore calculee ; ce controle apparaitra des qu'elle sera disponible.

## 3. Vues agregees

Moyennes non ponderees des Sure-d sur les cellules concernees (6 cellules par modele).

### 3.1 Modele M1

Profil Sure-d, moyenne sur toutes les cellules du modele :

| Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 |
|:---|---:|---:|---:|---:|---:|---:|
| Tail-index SIS | 0.979 | 0.996 | 0.996 | 0.996 | 0.996 | 0.996 |
| Yoshida--Umezu | 0.425 | 0.662 | 0.688 | 0.696 | 0.717 | 0.738 |
| Quantile SIS t=.90 | 0.858 | 0.988 | 0.996 | 0.996 | 0.996 | 1.000 |
| Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.975 | 0.996 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.99 | 0.875 | 0.988 | 1.000 | 1.000 | 1.000 | 1.000 |

Sure-20 par taille d'echantillon (moyenne sur p et rho) :

| Methode | n = 5000 |
|:---|---:|
| Tail-index SIS | 0.996 |
| Yoshida--Umezu | 0.688 |
| Quantile SIS t=.90 | 0.996 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

Sure-20 par dimension (moyenne sur n et rho) :

| Methode | p = 2000 |
|:---|---:|
| Tail-index SIS | 0.996 |
| Yoshida--Umezu | 0.688 |
| Quantile SIS t=.90 | 0.996 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

Sure-20 par correlation AR(1) (moyenne sur n et p) :

| Methode | rho = 0 | rho = 0.2 | rho = 0.25 | rho = 0.3 | rho = 0.4 | rho = 0.5 |
|:---|---:|---:|---:|---:|---:|---:|
| Tail-index SIS | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Yoshida--Umezu | 0.100 | 0.525 | 0.750 | 0.775 | 0.975 | 1.000 |
| Quantile SIS t=.90 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.99 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

### 3.2 Modele M2

Profil Sure-d, moyenne sur toutes les cellules du modele :

| Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 |
|:---|---:|---:|---:|---:|---:|---:|
| Tail-index SIS | 0.812 | 0.887 | 0.912 | 0.921 | 0.950 | 0.950 |
| Yoshida--Umezu | 0.054 | 0.150 | 0.183 | 0.225 | 0.296 | 0.358 |
| Quantile SIS t=.90 | 0.000 | 0.996 | 0.996 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.95 | 0.217 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.975 | 0.779 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.99 | 0.638 | 0.979 | 1.000 | 1.000 | 1.000 | 1.000 |

Sure-20 par taille d'echantillon (moyenne sur p et rho) :

| Methode | n = 5000 |
|:---|---:|
| Tail-index SIS | 0.912 |
| Yoshida--Umezu | 0.183 |
| Quantile SIS t=.90 | 0.996 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

Sure-20 par dimension (moyenne sur n et rho) :

| Methode | p = 2000 |
|:---|---:|
| Tail-index SIS | 0.912 |
| Yoshida--Umezu | 0.183 |
| Quantile SIS t=.90 | 0.996 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

Sure-20 par correlation AR(1) (moyenne sur n et p) :

| Methode | rho = 0 | rho = 0.2 | rho = 0.25 | rho = 0.3 | rho = 0.4 | rho = 0.5 |
|:---|---:|---:|---:|---:|---:|---:|
| Tail-index SIS | 0.675 | 0.875 | 0.975 | 0.975 | 1.000 | 0.975 |
| Yoshida--Umezu | 0.025 | 0.150 | 0.275 | 0.175 | 0.200 | 0.275 |
| Quantile SIS t=.90 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.975 |
| Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.99 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

### 3.3 Modele M3

Profil Sure-d, moyenne sur toutes les cellules du modele :

| Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 |
|:---|---:|---:|---:|---:|---:|---:|
| Tail-index SIS | 0.375 | 0.550 | 0.662 | 0.713 | 0.738 | 0.738 |
| Yoshida--Umezu | 0.062 | 0.125 | 0.150 | 0.163 | 0.287 | 0.350 |
| Quantile SIS t=.90 | 0.438 | 0.713 | 0.812 | 0.875 | 0.912 | 0.925 |
| Quantile SIS t=.95 | 0.900 | 0.963 | 0.988 | 0.988 | 0.988 | 1.000 |
| Quantile SIS t=.975 | 0.650 | 0.900 | 0.950 | 0.963 | 0.988 | 1.000 |
| Quantile SIS t=.99 | 0.075 | 0.413 | 0.550 | 0.688 | 0.788 | 0.925 |

Sure-20 par taille d'echantillon (moyenne sur p et rho) :

| Methode | n = 5000 |
|:---|---:|
| Tail-index SIS | 0.662 |
| Yoshida--Umezu | 0.150 |
| Quantile SIS t=.90 | 0.812 |
| Quantile SIS t=.95 | 0.988 |
| Quantile SIS t=.975 | 0.950 |
| Quantile SIS t=.99 | 0.550 |

Sure-20 par dimension (moyenne sur n et rho) :

| Methode | p = 2000 |
|:---|---:|
| Tail-index SIS | 0.662 |
| Yoshida--Umezu | 0.150 |
| Quantile SIS t=.90 | 0.812 |
| Quantile SIS t=.95 | 0.988 |
| Quantile SIS t=.975 | 0.950 |
| Quantile SIS t=.99 | 0.550 |

Sure-20 par correlation AR(1) (moyenne sur n et p) :

| Methode | rho = 0 | rho = 0.2 |
|:---|---:|---:|
| Tail-index SIS | 0.450 | 0.875 |
| Yoshida--Umezu | 0.000 | 0.300 |
| Quantile SIS t=.90 | 0.650 | 0.975 |
| Quantile SIS t=.95 | 0.975 | 1.000 |
| Quantile SIS t=.975 | 0.900 | 1.000 |
| Quantile SIS t=.99 | 0.300 | 0.800 |

### 3.4 Modele M4

Profil Sure-d, moyenne sur toutes les cellules du modele :

| Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 |
|:---|---:|---:|---:|---:|---:|---:|
| Tail-index SIS |  NaN |  NaN |  NaN |  NaN |  NaN |  NaN |
| Yoshida--Umezu |  NaN |  NaN |  NaN |  NaN |  NaN |  NaN |
| Quantile SIS t=.90 |  NaN |  NaN |  NaN |  NaN |  NaN |  NaN |
| Quantile SIS t=.95 |  NaN |  NaN |  NaN |  NaN |  NaN |  NaN |
| Quantile SIS t=.975 |  NaN |  NaN |  NaN |  NaN |  NaN |  NaN |
| Quantile SIS t=.99 |  NaN |  NaN |  NaN |  NaN |  NaN |  NaN |

Sure-20 par taille d'echantillon (moyenne sur p et rho) :

| Methode | n =  |
|:---|
| Tail-index SIS |
| Yoshida--Umezu |
| Quantile SIS t=.90 |
| Quantile SIS t=.95 |
| Quantile SIS t=.975 |
| Quantile SIS t=.99 |

Sure-20 par dimension (moyenne sur n et rho) :

| Methode | p =  |
|:---|
| Tail-index SIS |
| Yoshida--Umezu |
| Quantile SIS t=.90 |
| Quantile SIS t=.95 |
| Quantile SIS t=.975 |
| Quantile SIS t=.99 |

Sure-20 par correlation AR(1) (moyenne sur n et p) :

| Methode | rho =  |
|:---|
| Tail-index SIS |
| Yoshida--Umezu |
| Quantile SIS t=.90 |
| Quantile SIS t=.95 |
| Quantile SIS t=.975 |
| Quantile SIS t=.99 |

### 3.5 Composition du top-4

Nombre moyen, parmi les 4 coordonnees les mieux classees, de variables actives en indice de queue (A_gamma = {1,2,3,4}) et de variables d'echelle (A_scale = {5,6,7,8}), moyenne sur toutes les cellules.

| Modele | Methode | |top4 inter A_gamma| | |top4 inter A_scale| |
|:---|:---|---:|---:|
| M1 | Tail-index SIS | 3.98 | 0.00 |
| M1 | Yoshida--Umezu | 2.64 | 0.39 |
| M1 | Quantile SIS t=.90 | 3.85 | 0.06 |
| M1 | Quantile SIS t=.95 | 4.00 | 0.00 |
| M1 | Quantile SIS t=.975 | 4.00 | 0.00 |
| M1 | Quantile SIS t=.99 | 3.87 | 0.00 |
| M2 | Tail-index SIS | 3.74 | 0.02 |
| M2 | Yoshida--Umezu | 0.92 | 0.60 |
| M2 | Quantile SIS t=.90 | 1.63 | 2.37 |
| M2 | Quantile SIS t=.95 | 3.08 | 0.93 |
| M2 | Quantile SIS t=.975 | 3.78 | 0.21 |
| M2 | Quantile SIS t=.99 | 3.58 | 0.02 |
| M3 | Tail-index SIS | 3.03 | 0.00 |
| M3 | Yoshida--Umezu | 0.98 | 0.04 |
| M3 | Quantile SIS t=.90 | 3.11 | 0.04 |
| M3 | Quantile SIS t=.95 | 3.89 | 0.00 |
| M3 | Quantile SIS t=.975 | 3.58 | 0.00 |
| M3 | Quantile SIS t=.99 | 2.06 | 0.00 |
| M4 | Tail-index SIS | NaN | NaN |
| M4 | Yoshida--Umezu | NaN | NaN |
| M4 | Quantile SIS t=.90 | NaN | NaN |
| M4 | Quantile SIS t=.95 | NaN | NaN |
| M4 | Quantile SIS t=.975 | NaN | NaN |
| M4 | Quantile SIS t=.99 | NaN | NaN |

## 4. Tableaux complets

Une table par (modele, n, p) ; lignes = rho x methode. Sure-d pour d = 4, 10, 20, 30, 50, 100.

### 4.1 Modele M1

#### M1, n = 5000, p = 2000

| rho | Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 | E(Rmax) | Med(Rmax) |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **0** | Tail-index SIS | 0.900 | 0.975 | 0.975 | 0.975 | 0.975 | 0.975 | 9.1 | 4.0 |
|  | Yoshida--Umezu | 0.000 | 0.050 | 0.100 | 0.100 | 0.125 | 0.125 | 734.2 | 805.5 |
|  | Quantile SIS t=.90 | 0.500 | 0.925 | 0.975 | 0.975 | 0.975 | 1.000 | 7.1 | 4.5 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 0.550 | 0.925 | 1.000 | 1.000 | 1.000 | 1.000 | 5.7 | 4.0 |
| **0.2** | Tail-index SIS | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.1 | 4.0 |
|  | Yoshida--Umezu | 0.350 | 0.525 | 0.525 | 0.550 | 0.600 | 0.675 | 256.3 | 9.0 |
|  | Quantile SIS t=.90 | 0.925 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.1 | 4.0 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 0.900 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.2 | 4.0 |
| **0.25** | Tail-index SIS | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Yoshida--Umezu | 0.475 | 0.700 | 0.750 | 0.775 | 0.800 | 0.825 | 135.7 | 5.0 |
|  | Quantile SIS t=.90 | 0.950 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 0.900 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.2 | 4.0 |
| **0.3** | Tail-index SIS | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Yoshida--Umezu | 0.575 | 0.775 | 0.775 | 0.775 | 0.800 | 0.825 | 68.1 | 4.0 |
|  | Quantile SIS t=.90 | 0.900 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.1 | 4.0 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 0.925 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.1 | 4.0 |
| **0.4** | Tail-index SIS | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Yoshida--Umezu | 0.650 | 0.925 | 0.975 | 0.975 | 0.975 | 0.975 | 31.4 | 4.0 |
|  | Quantile SIS t=.90 | 0.925 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.1 | 4.0 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
| **0.5** | Tail-index SIS | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Yoshida--Umezu | 0.500 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.8 | 4.5 |
|  | Quantile SIS t=.90 | 0.950 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |

### 4.2 Modele M2

#### M2, n = 5000, p = 2000

| rho | Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 | E(Rmax) | Med(Rmax) |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **0** | Tail-index SIS | 0.425 | 0.650 | 0.675 | 0.675 | 0.775 | 0.775 | 93.6 | 5.5 |
|  | Yoshida--Umezu | 0.000 | 0.025 | 0.025 | 0.025 | 0.075 | 0.125 | 548.0 | 592.5 |
|  | Quantile SIS t=.90 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 8.0 | 8.0 |
|  | Quantile SIS t=.95 | 0.100 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 5.9 | 6.0 |
|  | Quantile SIS t=.975 | 0.600 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.6 | 4.0 |
|  | Quantile SIS t=.99 | 0.175 | 0.925 | 1.000 | 1.000 | 1.000 | 1.000 | 6.7 | 6.0 |
| **0.2** | Tail-index SIS | 0.725 | 0.850 | 0.875 | 0.900 | 0.950 | 0.950 | 15.0 | 4.0 |
|  | Yoshida--Umezu | 0.050 | 0.125 | 0.150 | 0.200 | 0.275 | 0.350 | 326.3 | 150.0 |
|  | Quantile SIS t=.90 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 7.8 | 8.0 |
|  | Quantile SIS t=.95 | 0.275 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 5.8 | 6.0 |
|  | Quantile SIS t=.975 | 0.800 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.3 | 4.0 |
|  | Quantile SIS t=.99 | 0.675 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 5.0 | 4.0 |
| **0.25** | Tail-index SIS | 0.900 | 0.950 | 0.975 | 1.000 | 1.000 | 1.000 | 5.0 | 4.0 |
|  | Yoshida--Umezu | 0.050 | 0.175 | 0.275 | 0.325 | 0.350 | 0.400 | 394.5 | 340.0 |
|  | Quantile SIS t=.90 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 7.7 | 8.0 |
|  | Quantile SIS t=.95 | 0.300 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 5.5 | 6.0 |
|  | Quantile SIS t=.975 | 0.775 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.2 | 4.0 |
|  | Quantile SIS t=.99 | 0.650 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 4.9 | 4.0 |
| **0.3** | Tail-index SIS | 0.875 | 0.925 | 0.975 | 0.975 | 0.975 | 0.975 | 11.3 | 4.0 |
|  | Yoshida--Umezu | 0.050 | 0.150 | 0.175 | 0.225 | 0.325 | 0.425 | 347.1 | 163.5 |
|  | Quantile SIS t=.90 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 7.7 | 8.0 |
|  | Quantile SIS t=.95 | 0.175 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 5.8 | 6.0 |
|  | Quantile SIS t=.975 | 0.775 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.3 | 4.0 |
|  | Quantile SIS t=.99 | 0.650 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.5 | 4.0 |
| **0.4** | Tail-index SIS | 0.975 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 4.2 | 4.0 |
|  | Yoshida--Umezu | 0.050 | 0.175 | 0.200 | 0.225 | 0.300 | 0.350 | 377.9 | 214.5 |
|  | Quantile SIS t=.90 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 7.8 | 8.0 |
|  | Quantile SIS t=.95 | 0.125 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 6.0 | 6.0 |
|  | Quantile SIS t=.975 | 0.750 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.4 | 4.0 |
|  | Quantile SIS t=.99 | 0.800 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.3 | 4.0 |
| **0.5** | Tail-index SIS | 0.975 | 0.975 | 0.975 | 0.975 | 1.000 | 1.000 | 4.7 | 4.0 |
|  | Yoshida--Umezu | 0.125 | 0.250 | 0.275 | 0.350 | 0.450 | 0.500 | 235.6 | 94.0 |
|  | Quantile SIS t=.90 | 0.000 | 0.975 | 0.975 | 1.000 | 1.000 | 1.000 | 8.4 | 8.0 |
|  | Quantile SIS t=.95 | 0.325 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 5.4 | 5.0 |
|  | Quantile SIS t=.975 | 0.975 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.99 | 0.875 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.2 | 4.0 |

### 4.3 Modele M3

#### M3, n = 5000, p = 2000

| rho | Methode | Sure-4 | Sure-10 | Sure-20 | Sure-30 | Sure-50 | Sure-100 | E(Rmax) | Med(Rmax) |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **0** | Tail-index SIS | 0.125 | 0.300 | 0.450 | 0.500 | 0.550 | 0.550 | 211.4 | 34.0 |
|  | Yoshida--Umezu | 0.000 | 0.000 | 0.000 | 0.000 | 0.100 | 0.150 | 755.8 | 684.5 |
|  | Quantile SIS t=.90 | 0.150 | 0.450 | 0.650 | 0.750 | 0.825 | 0.850 | 73.0 | 12.5 |
|  | Quantile SIS t=.95 | 0.800 | 0.925 | 0.975 | 0.975 | 0.975 | 1.000 | 6.6 | 4.0 |
|  | Quantile SIS t=.975 | 0.350 | 0.800 | 0.900 | 0.925 | 0.975 | 1.000 | 9.4 | 5.0 |
|  | Quantile SIS t=.99 | 0.025 | 0.150 | 0.300 | 0.500 | 0.675 | 0.875 | 55.4 | 31.0 |
| **0.2** | Tail-index SIS | 0.625 | 0.800 | 0.875 | 0.925 | 0.925 | 0.925 | 18.1 | 4.0 |
|  | Yoshida--Umezu | 0.125 | 0.250 | 0.300 | 0.325 | 0.475 | 0.550 | 287.5 | 69.5 |
|  | Quantile SIS t=.90 | 0.725 | 0.975 | 0.975 | 1.000 | 1.000 | 1.000 | 4.9 | 4.0 |
|  | Quantile SIS t=.95 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.0 | 4.0 |
|  | Quantile SIS t=.975 | 0.950 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 4.1 | 4.0 |
|  | Quantile SIS t=.99 | 0.125 | 0.675 | 0.800 | 0.875 | 0.900 | 0.975 | 15.6 | 6.0 |

### 4.4 Modele M4

## 5. Cout de calcul

Secondes par replication et par methode, monocoeur, moyenne sur modeles et rho. Quantile SIS est chronometre par valeur de tau. Derniere colonne : proportion moyenne de scores Yoshida--Umezu non definis.

| n | p | Tail-index SIS | Yoshida--Umezu | Quantile SIS (par tau) | YU non defini |
|---:|---:|---:|---:|---:|---:|
| 5000 | 2000 | 4.53 | 25.24 | 73.50 | 0.000 |

