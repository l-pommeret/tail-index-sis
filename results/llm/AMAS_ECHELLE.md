# Un amas de nuisance d'échelle à vérité exacte, greffé sur des données réelles

*Criblage par indice de queue contre quantile SIS, sur 20 000 documents de CommonCrawl.*

---

## 1. Ce que l'expérience cherche à établir

La méthode de l'article sélectionne les covariables qui font varier **l'indice de queue** γ(x). Un
écran par quantile conditionnel sélectionne celles qui font varier **un quantile élevé**. Les deux
coïncident souvent, et c'est pour cela que la distinction est difficile à démontrer : sur données
réelles, on ne sait en général pas laquelle des deux choses une covariable donnée déplace.

L'expérience construit donc un cas où la réponse est connue **par construction et sans
approximation** : un bloc de covariables dont on peut prouver qu'elles ne portent aucune
information sur γ, tout en déplaçant tous les quantiles finis de la réponse.

Le résultat, en une ligne : à partir d'une nuisance de force modérée, **quantile SIS remplit ses
vingt premières places de ces covariables ; le criblage par indice de queue en retient deux.**

---

## 2. La construction

### 2.1 Le lemme de Breiman

Si `Y` est à variation régulière d'indice `−1/γ` et si `V > 0` est **bornée** et **indépendante de
`Y`**, alors `Y·V` est à variation régulière du **même** indice, et

```
P(Y·V > y)  ~  E[V^(1/γ)] · P(Y > y),      y → ∞
```

Seule la constante d'échelle change. Ce n'est pas un développement asymptotique de second ordre à
contrôler, ni une condition de variation régulière étendue : c'est une identité de queue, exacte
sous les deux seules hypothèses citées, toutes deux vraies par construction ci-dessous.

### 2.2 Ce qui est greffé

```
F  ~  N(0,1)                                     facteur latent, indépendant de tout
V  =  exp{ κ (Φ(F) − 1/2) }                      borné dans [e^(−κ/2), e^(κ/2)]
Y' =  Y · V            soit    M' = M + log V    avec M = log Y
z_j = λ F + √(1−λ²) ε_j,     j = 1..20           les 20 proxys observés de F
```

avec `λ = 0,7`, donc une corrélation intra-amas de `λ² = 0,49`.

`V` est borné (donc `E[V^(1/γ)] < ∞` sans condition sur γ) et indépendant de `Y` (tiré d'un
générateur qui n'a pas vu les données). Les deux hypothèses de Breiman sont satisfaites **par
construction**, pas par vérification empirique. Les 20 coordonnées `z_j` déplacent donc tous les
quantiles finis de `Y'` et n'entrent pas dans son indice de queue — exactement.

Le facteur latent `F` et les proxys `Z` sont tirés **une seule fois** et partagés par les cinq
bras : d'un bras à l'autre, seule la **force** κ de la nuisance change, jamais sa direction. Les
cinq bras sont donc appariés.

### 2.3 Les données porteuses

`Y = exp(max_t S_t)`, où `S_t = −log p_θ(x_t | x_<t)` est la surprisal par jeton, mesurée en fp32
sur des documents de CommonCrawl. `Y` garde sa vraie queue et sa vraie structure de dépendance :
on ne simule ni la réponse, ni les 292 covariables réelles au milieu desquelles l'amas est greffé.

| | |
|---|---|
| documents | `n` = 20 000 |
| covariables réelles retenues | 272 (sur 298, après filtre de continuité) |
| proxys d'échelle greffés | 20 |
| **total** | **`p` = 292** |
| réglage | `α = n^(−0,30) = 0,0512`, `h = n^(−0,15)/2 = 0,1132` |
| extrêmes locaux effectifs `2nαh` | 232 |

Le filtre de continuité écarte 6 covariables (`len_sents`, `len_paras`, `len_lines`, `rate_digit`,
`rate_nonascii`, `max_word_len`) qui ont moins de 100 valeurs distinctes ou un atome de masse
supérieure à 2 % : les fenêtres de rangs y dégénèrent.

### 2.4 Le balayage

| κ | `sd(log V)` | `sd(M')` | `ρ(z_0, M')` | `ρ` max sur les 20 |
|---:|---:|---:|---:|---:|
| 0,0 | 0,0000 | 1,950 | 0,011 | 0,011 |
| 0,5 | 0,1444 | 1,956 | 0,063 | 0,063 |
| 1,0 | 0,2887 | 1,973 | 0,114 | 0,114 |
| 2,0 | 0,5775 | 2,037 | 0,210 | 0,211 |
| 4,0 | 1,1549 | 2,273 | 0,365 | 0,369 |

κ = 0 est le témoin : aucune nuisance, les 20 proxys sont du bruit pur. À l'autre extrémité,
κ = 4 reste modeste au regard de la dispersion de la réponse — `sd(log V)` y vaut la moitié de
`sd(M')`, et la corrélation de rang d'un proxy avec la réponse plafonne à 0,37.

---

## 3. Résultats

Nombre de proxys d'échelle classés dans les **20 premières** places, sur 20 greffés parmi 292
covariables. Un écran parfaitement spécifique en retient 0 ; le hasard en retiendrait 1,4.

| κ | queue, sélectionné | queue, agrégé | quantile .90 | quantile .95 | quantile .99 |
|---:|---:|---:|---:|---:|---:|
| 0,0 | 0 | 0 | 0 | 0 | 0 |
| 0,5 | 2 | 1 | 0 | 1 | 1 |
| 1,0 | **3** | **2** | **10** | **11** | **4** |
| 2,0 | **2** | **2** | **19** | **19** | **16** |
| 4,0 | **4** | **5** | **20** | **20** | **20** |

Rang **médian** des 20 proxys (le hasard donne 146) :

| κ | queue, sélectionné | queue, agrégé | quantile .90 | quantile .95 | quantile .99 |
|---:|---:|---:|---:|---:|---:|
| 0,0 | 175,0 | 126,2 | 205,0 | 211,0 | 140,0 |
| 0,5 | 152,5 | 125,2 | 48,5 | 64,5 | 84,5 |
| 1,0 | 167,5 | 181,0 | 21,0 | 19,5 | 37,0 |
| 2,0 | 164,0 | 145,8 | 11,5 | 11,5 | 11,5 |
| 4,0 | 143,0 | 66,8 | 10,5 | 10,5 | 10,5 |

### Lecture

**La bascule s'amorce à κ = 1,0 et s'achève à κ = 2,0 ; elle ne touche que quantile SIS.** À κ = 1,
quantile .95 fait passer les proxys d'un rang médian de 64,5 à 19,5 et en amène 11 sur 20 dans les
vingt premières places ; le criblage par indice de queue en est à 3, avec un rang médian de 167,5,
soit au-delà du hasard. À κ = 2, les trois niveaux de quantile placent presque tout l'amas en
tête — rang médian 11,5, c'est-à-dire que la moitié des 20 proxys occupent les 11 premières places
sur 292 — pendant que la méthode reste à 2 proxys et un rang médian de 164.

Le contraste se lit donc sur toute la plage, pas seulement à ses extrémités : **entre κ = 0,5 et
κ = 2, quantile SIS passe de 1 proxy à 19, la méthode de 2 à 2.**

**Le comportement de quantile SIS est le comportement correct pour ce qu'il estime.** Les proxys
déplacent réellement les quantiles élevés de `Y'` ; un écran qui vise les quantiles a raison de les
trouver. Le résultat n'est pas que quantile SIS se trompe, mais que **les deux écrans répondent à
deux questions différentes, et que seul le second répond à celle de l'indice de queue.**

**La dégradation de la méthode est lente mais réelle.** À κ = 4 elle retient 4 proxys sur 20 dans
son top 20, et l'un d'eux prend la première place. La méthode n'est pas immune à une nuisance
d'échelle assez violente — l'estimateur de Hill local reste un estimateur en échantillon fini, et
un déplacement d'échelle de grande amplitude finit par contaminer les fenêtres. Ce qu'établit
l'expérience est un **écart de sensibilité d'un ordre de grandeur**, pas une immunité.

**L'agrégation sur les réglages n'aide pas ici, et nuit à κ = 4.** Elle retient 2 proxys à κ = 2
comme la règle sélectionnée, et 5 contre 4 à κ = 4, avec un rang médian qui tombe de 143 à 66,8.
C'est cohérent avec ce qu'on observe en simulation : l'agrégation stabilise `R_max` quand le signal
est présent mais dispersé entre réglages, et n'a rien à stabiliser quand il n'y a aucun signal de
queue à trouver — elle ne fait alors que mettre en commun le bruit de neuf passes.

### Ce que la calibration par permutation dit

**Aucun proxy n'atteint le seuil familial à 5 % dans aucun des cinq bras, pour aucun écran.** Le
nombre de covariables sous le seuil **par covariable** à 5 % ne croît pas avec κ (33, 30, 36, 24,
14 sur 292, contre 15 attendues au hasard) et les proxys n'y sont jamais surreprésentés (2, 2, 4,
2, 4 — pour 20 proxys sur 292, le hasard en donnerait 1 à 2). Les fréquences de sélection par
sous-échantillonnage sont le seul indicateur calibré qui bouge : médiane de 0,05 pour les proxys
comme pour les autres à κ ≤ 1, puis 0,10 et 0,15 pour les proxys à κ = 2 et 4 contre 0,05 pour les
autres.

Autrement dit, **la contamination visible dans les rangs ne franchit jamais le seuil de
signification.** C'est le comportement souhaitable : la méthode ne déclare pas de découverte là où
il n'y en a pas. Mais cela signifie aussi que le classement brut, pris sans calibration, est plus
contaminé que la procédure complète — un point à retenir pour l'usage.

---

## 4. Portée et limites

**Ce que l'expérience établit.** Sur données réelles, avec une vérité exacte et non asymptotique,
le criblage par indice de queue est nettement moins sensible qu'un écran par quantile à une
nuisance qui déplace l'échelle sans toucher à l'indice. C'est une validation de **spécificité**.

**Ce qu'elle n'établit pas.** Elle ne dit rien de la **puissance** : elle ne montre pas qu'il
existe dans ces données une covariable portant réellement γ que la méthode trouverait et que
quantile SIS manquerait. Cette question ne se tranche pas ici, faute de vérité sur qui porte γ dans
un corpus réel ; elle relève des simulations.

**Le rôle de la greffe.** L'amas est artificiel, et c'est délibéré : c'est ce qui rend la vérité
exacte. La contrepartie est que le mécanisme de nuisance n'est pas issu des données. Une nuisance
d'échelle *naturelle* est documentée séparément — la longueur du document, dont le statut d'échelle
découle de `E[max_t S_t] = γ log T + O(1)` — mais son statut repose sur une approximation
asymptotique, `P(max > y) ~ θ T c y^(−1/γ)`, là où celui de `V` repose sur une identité.

Les deux expériences se complètent : l'une est exacte mais greffée, l'autre naturelle mais
asymptotique, et elles concluent dans le même sens.

---

## 5. Reproduction

```sh
python code/py/inject_scale_cluster.py results/llm/main_joint.parquet \
       results/llm/cluster 20 0.7          # 5 bras, S = 20 proxys, lambda = 0,7
sh code/sh/run_cluster.sh                  # criblage des 5 bras, 88 coeurs, 1000 permutations
```

| | |
|---|---|
| graine d'injection | 424242 |
| graine des permutations | `70000 + b` |
| graine du sous-échantillonnage | 909 |
| permutations par bras | 1000 |
| sous-échantillons | 20 |

Sorties : `results/llm/cluster/k*.parquet` (données), `results/llm/cluster_screen_k*/` (criblages),
`results/llm/cluster_k*.log` (journaux).
