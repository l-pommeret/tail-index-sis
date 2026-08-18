# Modèles de langue et valeurs extrêmes — rapport complet

*Tout ce qui a été fait pour confronter le criblage par indice de queue à des
données de modèles de langue : le pont théorique, cinq volets expérimentaux,
un bug qui en a invalidé un, et ce qui subsiste.*

Ce document est le point d'entrée. Deux documents de détail l'accompagnent :
[`AMAS_ECHELLE.md`](AMAS_ECHELLE.md) pour l'amas de nuisance à vérité exacte, et
[`../wild/APPLICATION_PERPLEXITE.md`](../wild/APPLICATION_PERPLEXITE.md) pour
l'application finale au format de l'article. Le rapport
[`RAPPORT_COMPLET.md`](RAPPORT_COMPLET.md) reste valable pour le détail du
premier volet, **à l'exception de sa section 8.4** (voir §7 ci-dessous).

---

## 1. La question

L'article crible les covariables dont dépend l'**indice de queue** `γ(x)` d'une
réponse à variation régulière. Les applications naturelles de cette idée sont
rares : il faut une réponse dont la lourdeur de queue ne soit pas un accident du
jeu de données, et des covariables nombreuses, continues et denses.

Les modèles de langue en offrent une, et elle a un statut que le taux de crimes
violents de Communities and Crime n'a pas.

---

## 2. Le pont théorique

Soit `S_t = −log p_θ(x_t | x_<t)` la **surprisal** du jeton `t`. Si la queue
supérieure de `S` est de type exponentiel, `P(S > x) ~ C e^{−λx}`, alors

```
Y = e^S    vérifie    P(Y > y) ~ C y^{−λ}
```

soit Pareto d'indice `γ = 1/λ`, **exactement**. L'indice de queue de `Y` est donc
le taux de décroissance de la distribution de surprise du modèle. Un `γ(x)`
grand signifie : *ce profil de covariables admet des jetons catastrophiquement
surprenants*. C'est, mot pour mot, le critère de curation de données.

On prend `Y = exp(max_t S_t)` sur une fenêtre de `T` jetons — le maximum d'un
échantillon à variation régulière l'est aussi, avec le même indice.

**Deux précautions de mesure, systématiques dans tous les volets.**

*L'amorce.* Les premiers jetons d'un document ont une surprisal artificiellement
haute, le modèle n'ayant pas de contexte. On saute 64 jetons.

*La longueur.* `E[max_t S_t] = γ log T + O(1)` : la longueur est une **nuisance
d'échelle**, elle déplace tous les quantiles finis sans entrer dans l'indice. On
score exactement `T` jetons par document, ce qui la neutralise par construction.

---

## 3. Volet A — première application (CommonCrawl, n = 20 000)

`Y = exp(max_t S_t)` sous un modèle de la famille Pythia, `T = 512` fixé,
`p = 272` covariables de surface, `2nαh = 232` extrêmes locaux par fenêtre.

### Ce qui a été vérifié avant tout criblage

Pas de plafond numérique sur les scores ; `γ > 0` marginalement et par strates ;
fp32 suffisant ; pas de dégénérescence par une classe de jetons ubiquitaire.

### Le criblage

| | `Y` | `Y_Δ` (écart petit/grand modèle) |
|---|---:|---:|
| seuil de permutation familial 5 % | 1,19056 | 1,19202 |
| score minimal observé | 1,19729 | 1,19711 |
| **sous le seuil familial** | **0 / 272** | **0 / 272** |
| sous le seuil par covariable 5 % | 31 (11,4 %) | 23 (8,5 %) |
| fréquence de sélection, médiane / max | 0,05 / 0,25 | 0,05 / 0,25 |

**Résultat nul.** Aucune covariable ne franchit le seuil familial. L'apport de
ce volet est méthodologique : sans la calibration par permutation, on aurait lu
un classement — car il en existe toujours un, signal ou pas.

### γ varie-t-il ?

Diagnostic indépendant de l'écran : dispersion inter-strates de `γ̂` comparée à
celle de strates aléatoires de mêmes tailles.

- loi nulle sur 400 découpages : dispersion médiane **0,0372** (prédiction
  théorique `γ/√k` = 0,0427 — la calibration est saine) ;
- observée sur 276 covariables : médiane **0,0379**. La covariable typique ne
  disperse **pas du tout** `γ` ;
- mais 33 covariables à `p < 0,05` pour 14 attendues, dont un cas net :

| covariable | dispersion | p | γ min → max |
|---|---:|---:|---|
| `rate_space` | **0,1232** | < 0,0025 | 1,208 → 1,498 |
| `mean_word_len` | 0,0850 | 0,0025 | 1,241 → 1,446 |

Pour `rate_space`, les IC de Q1 et Q3 ne se recouvrent pas, et l'amplitude de
0,29 sur `γ` vaut **1,8 nat** sur la perte du pire jeton à `T = 512`. Le profil
n'est pas monotone, donc « plus d'espaces, queue plus lourde » serait faux.

---

## 4. Volet B — la longueur, nuisance d'échelle à vérité analytique

Le premier cas où l'on sait, indépendamment de l'écran, ce qu'une covariable
fait et ne fait pas.

### La vérification directe

Sur cinq strates de longueur, en mode libre :

| | |
|---|---|
| régression de `med(M)` sur `log T` | pente **1,022** nats/nat, résidu max **0,121** |
| déplacement de `med(M)` d'une strate à l'autre | **2,239 nats** |
| `γ` par strate (Hill) | 1,240 à 1,579, **IC tous recouvrants** |
| `ρ(M, longueur)`, mode libre vs `T = 512` fixe | **+0,404** vs **−0,043** |
| Hill(20 %), mode libre vs fixe | 1,300 vs 1,383 — **même indice** |

La pente de 1,022 sur `log T`, avec un résidu maximal de 0,121, est exactement
la prédiction `E[max_t S_t] = γ log T + O(1)`. La localisation bouge de 2,24
nats, l'indice ne bouge pas.

### La conséquence au criblage

| | mode fixe (T neutralisé) | mode libre (T varie) |
|---|---:|---:|
| proxys de longueur dans le top-5 de quantile SIS | — | **5 / 5** aux trois τ |
| proxys de longueur dans le top-5 de l'écran de queue | — | **0 / 5** |
| rang médian sur 16, écran / quantile 0,95 | — | **13,0** / **3,0** |

La probabilité que les cinq proxys occupent les cinq premières places de 16 par
hasard est de `1/C(16,5) ≈ 2×10⁻⁴`.

### Une borne de validité, mesurée

À `n = 2 000` (`2nαh = 65`), l'écran mettait au contraire la longueur **en
tête**, avec une fréquence de sélection de 1,00. L'effet disparaît à
`n = 20 000` (`2nαh = 232`). Le mécanisme est identifiable : conditionner sur un
proxy de longueur homogénéise l'échelle dans la fenêtre et déplace le **biais de
second ordre** du Hill local, pas l'indice ; ce biais décroît avec le nombre
d'extrêmes par fenêtre.

C'est une limite honnête de la méthode, et elle est quantifiée.

---

## 5. Volet C — l'amas d'échelle greffé, à vérité exacte

La vérité du volet B repose sur une approximation asymptotique
(`P(max > y) ~ θTc·y^{−1/γ}`). Ce volet construit une nuisance dont le statut
est **exact**.

### Le lemme de Breiman

Si `Y` est à variation régulière d'indice `−1/γ` et si `V > 0` est **bornée** et
**indépendante de `Y`**, alors `Y·V` est à variation régulière du **même**
indice, et `P(YV > y) ~ E[V^{1/γ}]·P(Y > y)`. Seule l'échelle change. Ce n'est
pas un développement de second ordre : c'est une identité de queue.

### La greffe

```
F  ~ N(0,1)                            facteur latent, indépendant de tout
V  = exp{ κ (Φ(F) − 1/2) }             borné dans [e^{−κ/2}, e^{κ/2}]
z_j = λF + √(1−λ²) ε_j,  j = 1..20     les 20 proxys observés, λ = 0,7
```

`V` est borné et indépendant de `Y` **par construction**, pas par vérification.
Les 20 proxys déplacent donc tous les quantiles finis de `Y·V` sans entrer dans
son indice — exactement. Le facteur latent est tiré une seule fois et partagé
par les cinq bras : seule la **force** `κ` change.

### Résultat

Proxys classés dans les **20 premières** places, sur 20 greffés parmi 292
covariables, `n = 20 000`. Le hasard en donnerait 1,4.

| κ | queue, sélectionnée | queue, agrégée | quantile 0,90 | quantile 0,95 | quantile 0,99 |
|---:|---:|---:|---:|---:|---:|
| 0,0 | 0 | 0 | 0 | 0 | 0 |
| 0,5 | 2 | 1 | 0 | 1 | 1 |
| 1,0 | 3 | 2 | **10** | **11** | 4 |
| 2,0 | **2** | **2** | **19** | **19** | **16** |
| 4,0 | **4** | **5** | **20** | **20** | **20** |

**Entre κ = 0,5 et κ = 2, quantile SIS passe de 1 proxy à 19 ; la méthode de 2 à
2.** À κ = 2 la nuisance ne vaut pourtant que 28 % de l'écart-type de la réponse,
pour une corrélation de rang de 0,21.

Trois réserves, portées au document de détail : la méthode se dégrade quand même
à κ = 4 (4 proxys sur 20, dont la première place) ; l'agrégation sur réglages
n'aide pas ici et nuit à κ = 4 ; et **aucun proxy ne franchit le seuil familial
5 % dans aucun bras** — la contamination est dans les rangs bruts, jamais dans
la procédure calibrée.

---

## 6. Volet E — l'application perplexité, n = 100 000

*(Le volet D est traité en §7, il est invalidé.)*

C'est le volet qui produit un résultat **positif**. Détail complet dans
[`APPLICATION_PERPLEXITE.md`](../wild/APPLICATION_PERPLEXITE.md) ; l'essentiel
ici.

### Le renversement de design

Les covariables de surface se heurtent au filtre de continuité : sur une fenêtre
de 2 600 caractères, les taux de ponctuation rare et les trigrammes à `df < 0,98`
portent un atome de zéros supérieur à 2 %. **73 colonnes sur 512 survivent.**

Une NLL est un flottant. On prend donc comme covariables le **profil de
surprisal d'un petit modèle** (`pythia-70m`, `160m`) sur le même span de 512
jetons — 49 quantiles, quatre moments, statistiques par bloc et par tranche
lexicale — plus leur écart. **264 sur 264 passent le filtre.**

`Y` ne dépend que de `pythia-410m`, `X` que des deux petits modèles : aucun
quantile de `X` n'est fonction du vecteur de NLL dont `Y` est le maximum.

### La réponse est Pareto, vérifié

| k/n | 0,01 | 0,02 | 0,05 | 0,10 | 0,15 | 0,20 | 0,50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| γ (Hill) | 1,227 | 1,188 | 1,222 | 1,271 | 1,312 | 1,349 | 1,615 |

Plateau entre 1,19 et 1,35 sur un facteur 20 en `k`. QQ exponentiel des excès :
`r = 0,9993 / 0,9996 / 0,9993`, de pente 1,22–1,27 cohérente avec Hill.

### Le régime, et un point de méthode

| réglage | α | 2nαh | nαhΔ² | log(pn) | rapport |
|---|---:|---:|---:|---:|---:|
| `a* = 0,35` | 0,018 | 316 | 5,47 | 17,09 | 0,32 |
| `a = 0,20` | **0,100** | **1 778** | **30,76** | 17,09 | **1,80** |

Le réglage `a* = 0,35` a été choisi **à n = 2000**, où il donne `α = 0,070`.
Comme `α = n^{−a}`, garder `a` fixe fait chuter `α` : 0,018 à `n = 10⁵`. **Ce
qui se transporte d'un `n` à l'autre est `α`, pas `a`.** À `α ≈ 0,1`, la
condition `log(pn) = o(nαhΔ²)` est satisfaite pour la première fois dans le
projet, simulations comprises.

### Le criblage

| | `a = 0,20` | `a = 0,35` |
|---|---:|---:|
| score minimal | **0,93896** | **0,89381** |
| seuil familial 5 % | 1,26160 | 1,15892 |
| sous le seuil familial | **52** | **8** |
| sous 5 % par covariable | 142 (attendu 13) | 24 (attendu 13) |
| sous 1 % | 106 (attendu 3) | 12 (attendu 3) |

Les quatre premières sont identiques aux deux réglages et aux trois règles :
`s160m_max`, `s70m_max`, `s160m_kurt`, `s70m_kurt`, avec des fréquences de
sélection de 1,00 / 1,00 / 0,90 / 0,75 sur 20 sous-échantillons.

Profil de `γ` local le long du rang de `s160m_max` :

```
u        0,10  0,20  0,30  0,40  0,50  0,60  0,70  0,80  0,90
gamma    0,83  0,84  0,84  0,85  0,92  0,95  0,95  1,02  1,28    moyenne 0,942
```

contre un `γ` global de 1,271 et des témoins entre 1,25 et 1,27.

**Le maximum de surprisal d'un petit modèle gouverne l'indice de queue de la
surprisal catastrophique d'un grand modèle sur le même texte.** Le kurtosis et
l'asymétrie suivent : c'est la *forme* de la queue du profil qui porte
l'information, pas son niveau.

### La comparaison des règles

Les quatre porteuses sont les quatre plus fortes variations de `γ` au diagnostic
indépendant. Les deux témoins sont les surprisals **moyennes** : elles déplacent
le quantile 0,95 de `log Y` de 1,4 nat sans bouger `γ` (étendue du profil de
Hill 0,163 contre 0,449 pour le maximum).

| règle | porteuses top-4 | rang `s160m_mean` | rang `s70m_mean` |
|---|:---:|---:|---:|
| **queue, sélectionnée** | **4/4** | **177** | **241** |
| **queue, agrégée** | **4/4** | **170** | **189** |
| Yoshida–Umezu, son réglage | 4/4 | 100 | 125 |
| quantile 0,99 | 4/4 | 43 | 45 |
| quantile 0,95 | 3/4 | 24 | 27 |
| quantile 0,90 | 2/4 | 23 | 26 |
| Yoshida–Umezu, h = 2 | 2/4 | 15 | 21 |

`p = 264`, le hasard donnerait 132. **Seules les règles de queue prennent les
quatre porteuses et relèguent simultanément les témoins au-delà du rang 170.**

---

## 7. Le bug d'ordre mémoire, et le volet D qu'il invalide

### Nature

`numpy.tofile` écrit **toujours** en ordre C (par lignes) ; `matrix()` de R lit
en ordre Fortran (par colonnes). Une matrice de rangs transmise sans
transposition arrive brouillée :

```
lecture telle quelle : colonne 1 -> 63 031 valeurs distinctes sur 100 000
lecture transposée   : colonne 1 -> 100 000 valeurs distinctes, min 1, max n
```

**L'écran lisait du bruit, silencieusement.** Le symptôme était un résultat nul
parfaitement crédible : rangs médians au hasard, aucune significativité,
fréquences de sélection plates.

### Portée

| touché | non touché |
|---|---|
| **Volet D**, criblage des activations internes (7 configurations, toutes à 0 covariable significative) — `RAPPORT_COMPLET.md` §8.4 | Volet A (passe par `screen_llm.R`, qui lit un CSV) |
| Tous les criblages « wild » nuls avant correction | Volets B et C |
| | Tous les diagnostics EVT et de variation de `γ` (entièrement en Python) |

Le volet D affirmait qu'aucune direction du flux résiduel n'annonce un échec
catastrophique à venir, sur les couches 3, 6, 12, 18, 24 séparément, sur trois
couches jointes et en ACP. **Cette conclusion est sans valeur** et doit être
retirée. Vu le résultat du volet E, il y a une réelle chance que ces activations
portent du signal ; le criblage est à relancer.

### Correctif

`np.ascontiguousarray(R.T).tofile(...)` dans `export_ppl_ranks.py`,
`assemble_wild.py` et `export_ranks.py`. Et surtout, un **garde-fou** dans
`screen_wild.R` et `real_wild.R` : chaque colonne est vérifiée comme permutation
de `1..n` avant tout scorage, pour qu'une erreur de ce type ne puisse plus
passer en silence.

---

## 8. Bilan

### Ce qui soutient la méthode

**La spécificité, démontrée trois fois, avec une vérité connue indépendamment de
l'écran à chaque fois.**

| volet | vérité | résultat |
|---|---|---|
| B — longueur | asymptotique (`E[max] = γ log T`) | quantile SIS 5/5, écran de queue 0/5 |
| C — amas greffé | **exacte** (lemme de Breiman) | à κ = 2, quantile SIS 19/20, écran 2/20 |
| E — perplexité | vérifiée par profil de Hill local | témoins de localisation aux rangs 177 et 241, contre 23 et 27 |

**Un résultat positif, calibré et stable** (volet E) : significatif au seuil
familial aux deux réglages, quatre premières identiques aux trois règles,
fréquence de sélection de 1,00.

**Une borne de validité mesurée** (volet B) : à `2nαh = 65` l'écran se fait
prendre par la longueur, à `2nαh = 232` il ne s'y fait plus prendre.

### Ce qui ne la soutient pas

**La puissance supérieure n'est établie nulle part.** Au volet E, les quatre
porteuses sont aussi trouvées par Yoshida–Umezu à son réglage et par quantile
0,99. Ce qui distingue la méthode est ce qu'elle **rejette**, pas ce qu'elle
trouve seule.

**L'accord de rang global est défavorable.** Contre une mesure propre de
*variabilité* de `γ` — l'étendue du profil de Hill local — les règles de queue
obtiennent −0,357 et −0,370, quand quantile 0,90 atteint +0,487. L'explication
tient à la définition du score, qui est la **moyenne** du `γ` local et non sa
variabilité : l'écran retient ce qui abaisse le `γ` moyen, signature du biais de
mélange. Les deux fonctionnelles coïncident au sommet et divergent dans le corps
du classement. À noter que le critère par strates est lui-même contaminé : il
corrèle **+0,504 avec la variation de `γ` et +0,506 avec celle du quantile
0,95**, ce qui avantage mécaniquement les écrans par quantile.

**Le volet A est nul**, et le volet D est à refaire.

### Recommandation de rédaction

Présenter cet ensemble comme une **validation de spécificité sur données réelles
à vérité contrôlée**, la puissance restant établie par les simulations. C'est ce
que les chiffres portent, et c'est plus défendable en évaluation qu'une
affirmation de supériorité que le tableau des sept règles ne soutient pas.

---

## 9. Ce qui reste à faire

1. **Relancer le volet D** — criblage des activations internes, avec l'export
   corrigé. C'est le chantier le plus prometteur.
2. **La marche semi-synthétique** : texte réel (Wikipédia propre) plus une
   fraction de documents corrompus dont la sévérité est pilotée par quatre
   covariables plantées. Donnerait un `Sure-4` calculable sur du vrai texte — le
   chaînon qui manque entre les simulations et une vérité inconnue.
3. **L'échelle des modèles** : vérifier que l'ensemble actif du volet E est
   stable à travers la suite Pythia, comme le bloc `N₉` l'est à travers le
   réglage.
4. **Corriger `RAPPORT_COMPLET.md` §8.4**, qui affirme encore un résultat nul
   sur les activations internes.

---

## 10. Fichiers

| | |
|---|---|
| `code/py/extract_nll.py` | réponse `Y = exp(max NLL)`, `T` fixe, amorce |
| `code/py/extract_ppl_same_span.py` | covariables de perplexité, petit → grand modèle |
| `code/py/doc_covariates.py`, `build_covariates.py` | covariables de surface (volet abandonné) |
| `code/py/encode_docs.py` | coordonnées d'encodeur |
| `code/py/export_ppl_ranks.py`, `export_ranks.py` | pont vers R **(corrigés)** |
| `code/py/evt_gate.py` | plateau de Hill, QQ exponentiel |
| `code/py/gamma_varies_ppl.py` | `γ` varie-t-il, indépendamment de l'écran |
| `code/py/inject_scale_cluster.py` | amas de Breiman, volet C |
| `code/R/screen_wild.R` | criblage + calibration par permutation **(garde-fou)** |
| `code/R/real_wild.R` | les sept règles en concurrence **(garde-fou)** |
| `code/R/screen_llm.R`, `screen_internal.R` | volets A et D |

Graines : tirage 20260817 ; départage d'ex æquo 811000033 (criblage) et
973000019 (comparaison) ; permutations `505000 + b` ; sous-échantillonnage 707.
