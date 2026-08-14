# Comparaison des methodes de screening sur grille etendue

Genere le 2026-08-14 15:03 a partir de `results/grid/comparison_cells` (7 cellules sur 288).

> **Rapport partiel** : la campagne est en cours, 7 cellules sur 288 sont calculees. Les vues agregees de la section 3 ne portent que sur les cellules disponibles et bougeront encore. Les tables completes de la section 4 sont definitives cellule par cellule.

## 1. Protocole

Les quatre modeles M1-M4 de la Draft 3 (`code/R/generate3.R`), les memes
estimateurs et le meme reglage que `code/R/run_draft3.R comparison`, mais sur
un espace de design elargi.

- **Modeles** : M1, M2, M3, M4 ; ensemble actif en indice de queue A_gamma = {1,2,3,4} (pour M2, les variables d'echelle A_scale = {5,6,7,8} n'agissent que sur les quantiles finis, avec kappa = 0.20).
- **Tailles** : n = 1000, 2000, 5000.
- **Dimensions** : p = 200, 500, 1000, 2000.
- **Dependance** : X ~ AR(1) gaussien de correlation rho = 0, 0.20, 0.25, 0.30, 0.40, 0.50.
- **Replications Monte Carlo** : 40 par cellule, soit 7 cellules et 280 jeux de donnees simules.

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
| Tail-index SIS | 0.425 | 0.650 | 0.675 | 0.675 | 0.775 | 0.775 |
| Yoshida--Umezu | 0.000 | 0.025 | 0.025 | 0.025 | 0.075 | 0.125 |
| Quantile SIS t=.90 | 0.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.95 | 0.100 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.975 | 0.600 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Quantile SIS t=.99 | 0.175 | 0.925 | 1.000 | 1.000 | 1.000 | 1.000 |

Sure-20 par taille d'echantillon (moyenne sur p et rho) :

| Methode | n = 5000 |
|:---|---:|
| Tail-index SIS | 0.675 |
| Yoshida--Umezu | 0.025 |
| Quantile SIS t=.90 | 1.000 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

Sure-20 par dimension (moyenne sur n et rho) :

| Methode | p = 2000 |
|:---|---:|
| Tail-index SIS | 0.675 |
| Yoshida--Umezu | 0.025 |
| Quantile SIS t=.90 | 1.000 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

Sure-20 par correlation AR(1) (moyenne sur n et p) :

| Methode | rho = 0 |
|:---|---:|
| Tail-index SIS | 0.675 |
| Yoshida--Umezu | 0.025 |
| Quantile SIS t=.90 | 1.000 |
| Quantile SIS t=.95 | 1.000 |
| Quantile SIS t=.975 | 1.000 |
| Quantile SIS t=.99 | 1.000 |

### 3.3 Modele M3

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
| M2 | Tail-index SIS | 3.05 | 0.05 |
| M2 | Yoshida--Umezu | 0.28 | 0.47 |
| M2 | Quantile SIS t=.90 | 1.02 | 2.98 |
| M2 | Quantile SIS t=.95 | 2.80 | 1.20 |
| M2 | Quantile SIS t=.975 | 3.60 | 0.33 |
| M2 | Quantile SIS t=.99 | 2.88 | 0.12 |
| M3 | Tail-index SIS | NaN | NaN |
| M3 | Yoshida--Umezu | NaN | NaN |
| M3 | Quantile SIS t=.90 | NaN | NaN |
| M3 | Quantile SIS t=.95 | NaN | NaN |
| M3 | Quantile SIS t=.975 | NaN | NaN |
| M3 | Quantile SIS t=.99 | NaN | NaN |
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

### 4.3 Modele M3

### 4.4 Modele M4

## 5. Cout de calcul

Secondes par replication et par methode, monocoeur, moyenne sur modeles et rho. Quantile SIS est chronometre par valeur de tau. Derniere colonne : proportion moyenne de scores Yoshida--Umezu non definis.

| n | p | Tail-index SIS | Yoshida--Umezu | Quantile SIS (par tau) | YU non defini |
|---:|---:|---:|---:|---:|---:|
| 5000 | 2000 | 4.57 | 25.49 | 73.12 | 0.000 |

