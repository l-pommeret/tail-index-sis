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

## D4. M2 : facteur d'échelle en **amas corrélé**, s = 20, λ = 0.7

**Décision retenue.**

```r
# code/R/generate3.R
LAMBDA  <- 0.7          # charge du facteur latent sur chaque proxy
S_SCALE <- 20L          # A_scale = {5, ..., 24}
# kappa_c calibre pour que sd(log ell) egale celle du M2 publie au meme (n,p,rho) :
#   kappa_c = sd_cible * sqrt(12)     -- sd(log ell) = kappa_c / sqrt(12) exactement
# mesure : 0.1377 a rho=0.25 -> kappa_c = 0.477
#          0.143  a rho=0.30 -> kappa_c = 0.495
#          0.116  a rho=0    -> kappa_c = 0.402

# dans simulate_dataset3(), pour M2 uniquement :
f <- rnorm(n)
for (j in 5:(4 + S_SCALE))
  z[, j] <- LAMBDA * f + sqrt(1 - LAMBDA^2) * rnorm(n)
# ... puis
y <- V^(-gamma) * exp(-V / 2) * exp(kappa_c * (pnorm(f) - 0.5))
```

**Pourquoi changer.** Le facteur d'échelle publié est une somme de contributions
indépendantes, log ℓ = κ·Σ<sub>j∈A_scale</sub>(u_j − ½), ce qui lie l'attractivité
de chaque leurre (∝ κ) au dégât infligé au Hill local (∝ κ²s) — un seul
coefficient gouverne les deux. Conséquence mesurée
(`code/R/test_m2_scale20.R`) : passer s de 4 à 20 à κ constant multiplie
sd(log ℓ) par 2.4 et **tue les deux écrans** (Sure-20 ≈ 0.10 partout) ; et
apparier la variance par κ ∝ 1/√s affaiblit chaque leurre au point que
**quantile SIS s'améliore**, de 0.225 à 0.600 en Sure-4. Le nombre de
covariables d'échelle ne peut donc pas servir de cadran dans cette
paramétrisation.

**Ce que l'amas corrige.** Un facteur latent unique F porte toute la nuisance et
les s coordonnées d'échelle en sont des proxys de corrélation λ. Var(log ℓ) vaut
κ_c²/12 quels que soient s et λ, donc la difficulté pour le screen est figée,
pendant que chaque proxy garde une corrélation marginale λ avec la nuisance et
reste individuellement attractif pour un écran par quantiles. s et λ deviennent
deux cadrans indépendants de la difficulté.

**Mesuré** (`code/R/test_m2_cluster20.R`, sd(log ℓ) = 0.138 dans tous les bras,
n=2000, p=1000, ρ=0.25, 40 réplications, Sure-4 / Sure-20) :

| bras | screen 9 régl. min | quantile SIS τ=.95 | top-24 de quantile ∩ A_scale |
|:---|:---|:---|---:|
| s=4 publié | 0.325 / 0.775 | 0.050 / 0.975 | 3.58 / 4 |
| s=20, λ=0.5 | 0.400 / 0.825 | 0.175 / 1.000 | 14.45 / 20 |
| **s=20, λ=0.7** | **0.425 / 0.925** | 0.050 / 0.550 | 19.02 / 20 |
| s=20, λ=0.9 | 0.425 / 0.900 | 0.000 / 0.075 | 20.00 / 20 |

Le screen est **plat en λ** (0.825 / 0.925 / 0.900) pendant que quantile SIS
s'effondre (1.000 / 0.550 / 0.075). Le mécanisme se lit au chiffre près à λ=0.9 :
quantile SIS met les 20 variables d'échelle dans son top-24, son Sure-20 vaut
0.075 mais son Sure-30 vaut 1.000 — les actives se rangent exactement derrière
les leurres, à la marche d = 4 + s.

**Pourquoi λ = 0.7 et non 0.9.** L'effet est déjà décisif (0.925 contre 0.550) et
la corrélation intra-amas vaut λ² = 0.49, du même ordre que le ρ = 0.5 déjà
présent dans l'étude, donc défendable comme réaliste. À λ = 0.9 elle monte à
0.81, plus spectaculaire mais plus facile à contester.

**Meilleur point de design mesuré** : n=2000, p=500, ρ=0.30 — screen agrégé
**0.550 / 0.975** contre **0.050 / 0.750** pour quantile SIS, qui place 19.1 des
20 variables d'échelle dans son top-24.

**Réserve de présentation.** κ = 0.20 provient d'un pilote pré-enregistré dont la
règle de sélection était fixée avant inspection ; ni le passage à l'amas, ni
s = 20, ni λ = 0.7 n'ont ce statut — ils ont été choisis en connaissant les
résultats. À présenter comme une variante exploratoire (un M5 ajouté à M2), pas
comme faisant partie du protocole initial. Un pilote pré-enregistré sur λ, sur
le modèle de celui qui a fixé κ, lèverait cette réserve.

**Conséquence sur le code.** A_scale passe de {5,…,8} à {5,…,24} : les drivers
qui codent `%in% 5:8` pour la composition du top-4 doivent être mis à jour, et
`simulate_score_streaming3` suppose actuellement que toutes les coordonnées
pertinentes tiennent dans les 8 premières colonnes.

---

## Décisions en attente

- Recalibrer κ par modèle pour que M3 et M4 supportent la perturbation d'échelle
  de M2 : le même κ appliqué à un contraste deux à trois fois plus faible écrase
  le signal (Sure-20 de M3 : 0.750 → 0.325 en substituant ℓ₂ à ℓ₁).
- Généraliser le code : `sv_l2` code en dur les colonnes 5:8,
  `simulate_score_streaming3` suppose que tout tient dans les 8 premières
  colonnes, les drivers codent `%in% 5:8`.
