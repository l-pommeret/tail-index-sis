# Audit du départage d'ex æquo — instructions de run (calculateur)

## Le problème

Dans `campaign5.R`/`campaign6.R`/`real_crime5.R`, le classement agrégé est
calculé ainsi :

- rangs par réglage : `rank(S[,k], ties.method = "first")` (ex æquo →
  indice de colonne) ;
- agrégation : `amin = min` des 9 rangs — le minimum de 9 permutations
  produit des ex æquo **structurels** (jusqu'à 9 coordonnées partagent la
  même valeur) ;
- ordre final : `order(amin, seq_len(p))` — **départage par indice**.

Les vraies actives étant les colonnes 1–4, tout ex æquo actif/inactif est
résolu en faveur de l'actif : Sure-4/Sure-20/E(Rmax) de la règle agrégée
sont potentiellement biaisés vers le haut. Les 5 autres règles ont des
scores continus (ex æquo exacts de mesure nulle — vérifié par le
diagnostic) et toutes les méthodes sont **marginales** (le score de j ne
dépend que de (X_j, Y)), donc équivariantes par permutation de colonnes
au bit près : seule la règle de classement peut faire fuir les labels
dans les résultats.

## La correction (`code/R/rank_rules.R`)

Rangs par réglage `ties.method = "average"` (invariant), puis clé
lexicographique **(min, médiane, moyenne, max)** des 9 rangs — le min
reste la clé primaire, la philosophie min-rank est inchangée — puis
départage résiduel par tirage uniforme **seedé** (jamais l'indice).
Flux de seeds dédiés, jamais utilisés ailleurs :

- tie-break : `973000019 + cell_id*100003 + model_index*10007 + r*307`
- permutation de colonnes : `987000037 + cell_id*100003 + model_index*10007 + r*307`

Les seeds de données sont inchangées (flux campaign5) : la comparaison
avant/après est appariée réplication par réplication.

## À lancer (dans cet ordre)

```sh
# 0. Test d'invariance par permutation (~10 min, 1 cœur).
#    Vérifie bitwise l'équivariance des scores et l'invariance de la
#    nouvelle règle ; signale si l'ancienne règle dépend des labels.
Rscript tests/test_rank_invariance.R

# 1. Diagnostic des ex æquo AVANT correction (25 reps × 12 cellules +
#    contrôle des ex æquo chez les concurrents ; ~30 min sur 90 cœurs).
Rscript code/R/diag_ties6.R results/diag_ties6 25 90 5

# 2. Campagne appariée avant/après, 1000 reps × 12 cellules.
#    9 passes de score par rep, PAS de YU/quantile SIS (inchangés, cf.
#    en-tête de campaign7_rankfix.R) -> rapide (~1-2 h sur 90 cœurs).
#    Checkpoints par blocs : relançable tel quel après interruption.
C7_BLOCK=270 Rscript code/R/campaign7_rankfix.R results/campaign7 A1,A2,A3,B1 1000 90

# 3. Données réelles (une passe, ~1 min).
Rscript code/R/real_crime6.R results/draft5_real_crime_rankfix.rds
```

L'étape 2 imprime en fin de run un **crosscheck** : les métriques de
l'ancienne règle recalculées doivent coïncider à l'identique (max|diff| = 0)
avec les rmax stockés de campaign5 (A1/A2 = M1/M3) et campaign6 (A3/B1).
Si ce n'est pas 0, ne rien conclure et me le signaler.

## Sorties à pousser

- `results/diag_ties6/` (`diag.csv`, `competitor_ties.csv`)
- `results/campaign7/` (`summary.csv`, `before_after.csv`, `tie_stats.csv`
  — les `c7_*.rds` aussi si possible)
- `results/draft5_real_crime_rankfix.rds`
- les logs console (crosscheck + test d'invariance)

## Ce qui est mesuré

- `before_after.csv` : pour chaque modèle × p, Sure-4/Sure-20/E(Rmax)/
  Med(Rmax) ancien vs nouveau vs bras permuté, différences appariées et
  leurs erreurs types, fraction de réplications où Rmax change.
- `tie_stats.csv` : nb moyen de coordonnées ex æquo dans amin, fréquence
  d'un groupe d'ex æquo mixte (actif+inactif) chevauchant les positions
  4 et 20 — les seules configurations où le départage peut agir.
- Bras `prm_*` : nouvelle règle + permutation aléatoire des colonnes à
  chaque réplication (implémentée par permutation des lignes de la
  matrice de scores, bitwise-équivalente au recalcul sur X permuté —
  c'est ce que vérifie le test 0).

## Protocole inchangé

Modèles, n, p, rho, d, epsilon, grille, N9, alpha/h, méthodes
concurrentes, nombre de réplications et seeds de données : identiques.
Seuls ajouts : les deux flux de seeds ci-dessus (tie-break, permutation).
Aucun résultat publié n'est modifié par ces scripts ; ils produisent le
matériel avant/après pour décider des mises à jour du manuscrit.
