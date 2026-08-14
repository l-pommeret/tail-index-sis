# Décisions de design — journal

Décisions prises à partir des mesures de la session du 2026-08-14, avec les
chiffres qui les motivent. Les scripts nommés sont dans `code/R/`.

---

## D1. M4 : porter le taux de décroissance de γ de 0.35 à **0.80**

**Décision retenue.**

```r
# code/R/generate3.R, gamma_model3()
# avant
M4 = 0.55 * exp(-0.35 * s4 - 0.20 * u[,1]*u[,2] - 0.15 * u[,3]*u[,4])
# après
M4 = 0.55 * exp(-0.80 * s4 - 0.20 * u[,1]*u[,2] - 0.15 * u[,3]*u[,4])
```

**Pourquoi.** Pour γ = c·exp(−β·s₄ − interactions), conditionner sur u_j = t et
maximiser sur les trois autres coordonnées actives les place en 0, où les deux
termes d'interaction s'annulent. L'enveloppe supérieure vaut donc exactement
c·exp(−βt), et la séparation de population est

    Δ/c = 1 − {e^(−βε) − e^(−β(1−ε))} / {β(1−2ε)},   ε = 0.05

soit 0.157 pour β = 0.35 contre 0.373 pour M1 et M2 à β = 1. M4 était donc
2.4 fois moins détectable que M1 par construction, ce qui explique qu'il soit le
modèle le plus faible dans toute la campagne.

**Mesuré** (`code/R/test_m4_beta.R`, n=2000, p=1000, ρ=0.25, 40 réplications,
Sure-4 / Sure-20) :

| β | Δ/c | screen 9 régl. min | quantile SIS τ=.95 | pipeline + 9 régl. min |
|---:|---:|:---|:---|:---|
| 0.35 (publié) | 0.157 | 0.100 / 0.500 | 0.075 / 0.650 | 0.225 / 0.700 |
| 0.55 | 0.233 | 0.275 / 0.875 | 0.500 / 0.950 | 0.575 / 0.975 |
| 0.70 | 0.284 | 0.625 / 0.975 | 0.700 / 0.975 | 0.800 / 1.000 |
| **0.80** | **0.315** | **0.700 / 0.975** | **0.675 / 0.975** | **0.875 / 0.975** |
| 0.90 | 0.345 | 0.775 / 0.975 | 0.875 / 1.000 | 0.900 / 1.000 |
| 1.00 | 0.373 | 0.925 / 1.000 | 0.750 / 1.000 | 0.950 / 1.000 |

β = 0.80 place M4 dans le régime où le screen et quantile SIS sont au coude à
coude (0.700 contre 0.675 en Sure-4, 0.975 des deux côtés en Sure-20) et où le
pipeline se détache nettement (0.875). Le modèle reste donc **discriminant**
entre les trois règles, ce que β ≥ 0.95 lui ferait perdre par saturation.

**Vérifié aussi à ρ = 0** (n=2000, p=1000, 40 réplications) — le choix tient aux
deux valeurs de ρ, ce qui est un bon signe de robustesse :

| bras | screen 9 régl. min | quantile SIS τ=.95 | pipeline + 9 régl. min |
|:---|:---|:---|:---|
| M1 (référence) | 0.200 / 0.800 | 0.275 / 0.875 | 0.425 / 0.900 |
| M4 publié (β=0.35) | 0.000 / 0.075 | 0.000 / 0.175 | 0.025 / 0.150 |
| **M4 (β=0.80)** | **0.175 / 0.800** | 0.275 / 0.825 | 0.300 / 0.850 |

β = 0.80 fait passer M4 de 0.075 à 0.800 en Sure-20, c'est-à-dire d'un modèle où
rien ne fonctionne à un modèle exploitable, et le cale exactement sur M1
(0.800 contre 0.800). La dynamique entre valeurs de d y est même meilleure qu'à
ρ = 0.25 : écart Sure-20 − Sure-4 de 0.625 contre 0.275, sans aucune saturation.

**Ce qui ne change pas.** Les coefficients d'interaction restent à 0.20 et 0.15.
Ils vivent dans l'intérieur de la surface, pas sur l'enveloppe, donc ils
n'affectent pas la détectabilité — vérifié en comparant β=1 à interactions
publiées (0.925 / 1.000) et à interactions mises à l'échelle ×2.86
(0.825 / 1.000) : écart dans le bruit. M4 conserve donc sa vocation de modèle
non additif.

**Borne supérieure à ne pas franchir : β ≈ 1.5.** Au-delà, γ_min devient si petit
(3×10⁻⁶ à β=3) qu'une partie de la population cesse d'être à queue lourde :
Y = V^(−γ)·ℓ ≈ ℓ. Le Hill local n'a alors plus de queue à estimer sur cette
portion de fibre et s'effondre — Sure-4 de 0.125 à β=3 — alors que la séparation
de population continue de croître (Δ/c = 0.703). Quantile SIS, qui compare des
quantiles finis, n'est pas affecté (0.775). C'est une limite du **cadre de
modélisation**, pas de l'estimateur, et elle mérite d'être signalée dans le
papier.

**Conséquence à assumer.** M4 cesse d'être le modèle « difficile » de l'étude.
Si ce rôle doit être conservé, il faut soit garder une variante à β = 0.35 ou
0.55 en plus, soit reporter la difficulté sur un autre axe (ρ = 0, n = 1000).

---

## D2. Agrégation de réglages : grille **3×3**, agrégateur **minimum**

**Décision retenue** pour l'étage 2 de la méthode combinée
(voir `results/methode_combinee.md`).

Grille a ∈ {0.30, 0.35, 0.40} × b ∈ {0.05, 0.10, 0.15}, classement par le
**rang minimum** sur les 9 réglages.

**Pourquoi 9 et non 25.** Densifier la grille au pas de 0.025 ne change rien :
sur huit configurations comparées (deux tailles × deux agrégateurs × avec et
sans pipeline), l'écart moyen en Sure-4 est de **+0.006**, de signe alternant,
pour une erreur type d'environ 0.04. Des réglages espacés de 0.025 produisent
des classements trop corrélés pour augmenter le nombre effectif de
configurations indépendantes. La grille à 9 points coûte 2.8 fois moins et fait
aussi bien.

**Pourquoi ne pas élargir.** Une grille 5×5 couvrant a ∈ [0.25, 0.45] et
b ∈ [0, 0.20] fait **chuter** le Sure-4 moyen de 0.352 à 0.154 à nombre de points
égal : elle inclut des réglages dont le rapport séparation/bruit mesuré tombe à
1.36 (b = 0.20) contre 3.10 à b = 0.10, et même une médiane finit par en
souffrir.

**Pourquoi le minimum.** Seul, il est moins bon que la médiane (0.296 contre
0.381 en Sure-4). Derrière la présélection par quantile SIS, il devient meilleur
(0.544 contre 0.488) : le rang minimum promeut une coordonnée dès qu'un seul
réglage la classe bien, ce qui promeut aussi les nulles chanceuses — rédhibitoire
face à 996 concurrentes, sans effet face aux 21 qui survivent à l'étage 1.

---

## D3. Points de design de référence

Mesurés sur la campagne (288 cellules) et les tests ciblés :

- **ρ est le facteur le plus influent**, loin devant n et p : Sure-20 moyen du
  screen de 0.319 à ρ=0 contre 0.822 à ρ=0.5. Mécanisme : les quatre
  coordonnées actives étant adjacentes dans la structure AR(1), conditionner sur
  l'une déplace aussi les trois autres et creuse la dépression de l'enveloppe.
  **ρ = 0.30 à 0.40 est un choix de référence plus favorable que 0.25**, et
  défendable.
- **ρ = 0 est le régime le plus défavorable**, à toutes les tailles testées :
  même à n = 4000 et p = 200, le screen agrégé plafonne à 0.838 de Sure-20 contre
  1.000 pour quantile SIS. Seul M2 y reste discriminant.
- **n = 1000 à ρ = 0 est inexploitable** : meilleur Sure-20 de toute la grille
  0.650, la plupart des cellules sous 0.20, pour les quatre règles.

---

## D4. Variante M2 en amas corrélé — **exploratoire, à ne pas présenter comme pré-enregistrée**

Le facteur d'échelle de M2 en somme indépendante lie l'attractivité de chaque
leurre (∝ κ) au dégât infligé au Hill (∝ κ²s) : augmenter le nombre de
covariables d'échelle tue les deux écrans à κ constant, ou affaiblit les leurres
au point que quantile SIS **s'améliore** (Sure-4 de 0.225 à 0.600) si l'on
apparie la variance. La variante en amas — un facteur latent unique observé par
s proxys de corrélation λ — sépare les deux axes et fait de s un vrai cadran.

À n=2000, p=500, ρ=0.30, λ=0.7, s=20, variance appariée : screen agrégé
**0.550 / 0.975** contre **0.050 / 0.750** pour quantile SIS, qui place 19.1 des
20 variables d'échelle dans son top-24.

**Réserve.** κ = 0.20 provient d'un pilote pré-enregistré avec règle de sélection
fixée avant inspection ; ni le passage à l'amas, ni s = 20, ni λ = 0.7 n'ont ce
statut. À présenter comme une variante exploratoire (un M5), pas comme faisant
partie du protocole initial.

---

## Décisions en attente

- Recalibrer κ par modèle pour que M3 et M4 supportent la perturbation d'échelle
  de M2 : le même κ appliqué à un contraste deux à trois fois plus faible écrase
  le signal (Sure-20 de M3 : 0.750 → 0.325 en substituant ℓ₂ à ℓ₁).
- Pilote pré-enregistré sur λ, sur le modèle de celui qui a fixé κ.
- Généraliser le code : `sv_l2` code en dur les colonnes 5:8,
  `simulate_score_streaming3` suppose que tout tient dans les 8 premières
  colonnes, les drivers codent `%in% 5:8`.
