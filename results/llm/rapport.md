# Application LLM : indice de queue de la surprise par token

Exécution du protocole sur FineWeb et Pythia. Le rapport suit l'ordre des
décisions, et signale à chaque étape ce qui a été vérifié plutôt que supposé.

---

## 1. En bref

L'application demandait d'abord de régler un piège : la perplexité moyenne se
concentre et tombe dans le domaine de Gumbel, donc γ = 0 et l'écran classerait
du bruit. La réponse retenue est le **maximum de surprise exponentié**, dont le
mécanisme est exact.

Les prérequis sont **tous vérifiés par la mesure** : pas de plafond numérique,
γ > 0 marginalement et dans les seize strates de covariables, fp32 suffisant.
La longueur se comporte exactement comme la nuisance d'échelle que la théorie
des valeurs extrêmes annonce, et cela se vérifie sur préfixes emboîtés à
composition fixe.

Le criblage, lui, **ne trouve rien** — et trois mesures indépendantes le disent.
Ce n'est pas un échec : c'est le résultat que seule la calibration par
permutation permet d'énoncer, là où un tableau de rangs aurait suggéré une
découverte.

Enfin, un **faux positif mesuré** : quand la longueur varie, l'écran la remonte
en tête avec une stabilité parfaite, alors qu'elle n'affecte prouvablement pas
γ. C'est l'écart entre l'estimande et ce que classe l'estimateur à n fini,
observé sur données réelles avec une vérité analytique.

---

## 2. Pourquoi le maximum de surprise, et pas la perplexité

Avec S_t = −log p_θ(w_t | w_<t), la perplexité moyenne
PPL = exp{T⁻¹ Σ_t S_t} est une moyenne de T surprises : sous mélange, log PPL
se concentre, PPL est approximativement lognormale, donc **γ = 0**. Le piège est
qu'un marginal lourd peut coexister avec des conditionnelles lognormales — le
diagnostic de Hill marginal rassurerait à tort, (C1) serait violée et ξ_j ≡ 0.

La réponse retenue est

    Y = exp( max_{t ≤ T} S_t ) = 1 / min_t p̂_θ(w_t | w_<t).

Le mécanisme est exact : si P(S > s | x) ~ c(x) e^{−s/γ(x)}, alors
P(Y > y | x) ~ c(x) y^{−1/γ(x)}, variation régulière **conditionnelle** avec
γ(x) > 0 par construction, et le facteur sous-exponentiel devient exactement le
ℓ_j(s,u) du modèle. Y > 0 automatiquement, sans filtrage.

Une seconde réponse est extraite dans la même passe :

    Y_Δ = exp( max_t ( S^petit_t − S^grand_t ) ) = max_t p_grand / p_petit,

dont l'indice est le taux de la queue du rapport de vraisemblance maximal.
Les deux modèles appartiennent à la même famille, **tokenizers vérifiés
identiques par assertion**, donc les deux surprises sont alignées token par
token et la différence a un sens.

---

## 3. Dispositif

| | |
|:---|:---|
| corpus | FineWeb, dump `CC-MAIN-2025-18` |
| modèle évalué | `EleutherAI/pythia-410m-deduped` |
| modèle de référence | `EleutherAI/pythia-1.4b-deduped` |
| précision | fp32 de bout en bout, logits castés avant softmax |
| déterminisme | batch = 1 |
| burn-in | 64 tokens, tokens spéciaux exclus |
| analyse primaire | T = 512 tokens notés **exactement** |
| analyse secondaire | T libre, la longueur devient la nuisance |
| n | 20 000 documents retenus sur 50 405 lus |
| covariables | 27 de surface + 256 composantes d'ACP d'un encodeur séparé gelé |

**Corpus.** FineWeb applique des filtres heuristiques et de déduplication, pas
de classifieur de qualité ni de seuil KenLM — contrairement à FineWeb-Edu et à
CCNet, qui tronquent exactement la queue supérieure de Y. Le dump 2025 est très
postérieur à la collecte du Pile (2020), donc la non-appartenance est acquise
par construction plutôt que vérifiée. **Réserve** : une date de crawl récente ne
garantit pas un contenu récent ; une page ancienne peut être recrawlée.

**Covariables.** Calculées depuis le texte brut, jamais depuis la passe forward
qui produit Y. Les embeddings viennent de `roberta-large`, gelé, d'une autre
famille que les Pythia. Filtre de continuité : au moins 100 valeurs distinctes
**et** aucun atome portant plus de 2 % de la masse — il écarte
`rate_digit`, `rate_nonascii`, `len_sents`, `len_paras`, `len_lines` et
`max_word_len`, c'est-à-dire les taux massivement nuls et les comptages entiers
de faible cardinalité. 272 covariables sur 278 sont retenues. Le départage des
ex æquo est aléatoire à graine fixée, car la géométrie des fenêtres de rangs en
dépend.

**Régime.** Au réglage (a\*, b\*) = (0.30, 0.15) et n = 20 000 :
α = 0.0512, h = 0.1132, et **2nαh = 232** extrêmes locaux par fenêtre.

---

## 4. Les prérequis, vérifiés

### 4.1 Pas de plafond numérique

max_t S s'étale de **2,35 à 27,92 nats** (médiane 13,89, q99 19,36). La part de
documents à moins d'un nat du maximum global vaut **0,0001**, soit deux
documents. Les écarts entre les vingt plus grands maxima restent réguliers
(médiane 0,143 nat). Une queue tronquée se verrait comme une accumulation sous
la borne ; il n'y en a pas.

### 4.2 γ > 0, marginalement

| réponse | Hill (k=10 %) | ξ (GPD) | IC95 profilé |
|:---|---:|---:|:---|
| Y, petit modèle | 1,276 | 1,159 | **[1,073 ; 1,251]** |
| Y, grand modèle | 1,291 | — | — |
| Y_Δ | 1,338 | 1,176 | **[0,919 ; 1,501]** (pilote) |

L'IC de vraisemblance profilée exclut 0 très largement. Sur le pilote à
n = 2 000 il valait [0,93 ; 1,54] ; à n = 20 000 il se resserre à
[1,07 ; 1,25], comme attendu.

### 4.3 γ > 0 **par strates** — le test qui porte sur (C1)

Quartiles de `gzip_ratio`, `rate_digit`, `ttr` et `byte_entropy`, seize strates
de 5 000 documents : **ξ ∈ [1,14 ; 1,29]**, bornes basses d'IC de 1,02 à 1,16.
Aucun intervalle ne touche 0. Le scénario redouté — marginal lourd par mélange,
conditionnelles lognormales — est écarté.

Un fait ressort et oriente tout le reste : **γ est remarquablement constant**
d'une strate à l'autre, tous les IC se recouvrant.

### 4.4 fp32 suffisant

Recalcul des 200 plus grandes surprises en fp64, modèle compris :

| contrôle | résultat |
|:---|:---|
| reproduction du stockage fp32 | écart max **0,00e+00** |
| écart fp64 − fp32 | médiane −3,6e−06, max **8,1e−04 nat** |
| corrélation de rang | Kendall **0,9998**, Spearman **0,99999** |
| rangs identiques | 196 / 200 |
| vingt premiers identiques | **20 / 20** |

Cinq ordres de grandeur séparent l'erreur numérique de l'échelle des valeurs.

---

## 5. La longueur est une nuisance d'échelle — vérifié, pas supposé

Pour T variables à variation régulière faiblement dépendantes,
P(max > y | x) ≈ θ T c(x) y^{−1/γ(x)} : **T entre dans l'échelle et jamais dans
l'indice**. En logarithme, E[max_t S | x] = γ(x) log T + O(1).

**Test à composition fixe.** Sur 600 documents, en scorant des préfixes
emboîtés T ∈ {128, 256, 512, 1024} du *même* document :

| T | 128 | 256 | 512 | 1024 |
|---:|---:|---:|---:|---:|
| E[max_t S] | 12,245 | 13,040 | 13,964 | 14,782 |

Pente **1,2315 nat/nat**, résidu maximal 0,041 nat, IC95 intra-document
**[1,149 ; 1,314]** — qui **contient les deux estimations indépendantes de γ**
(Hill 1,24–1,30). La croissance est monotone dans 1,000 des cas, comme l'exige
un maximum sur un ensemble croissant.

**Strates de longueur.** L'indice reste indiscernable d'une constante — Hill de
1,240 à 1,579, les cinq IC se recouvrant — pendant que la localisation se
déplace de **2,24 nats**.

C'est la structure de M2 sur données réelles, avec un statut **analytique** :
la longueur arrive naturellement avec son amas de proxys corrélés (caractères,
octets, mots, espaces, lignes, taille gzip), sans construction artificielle.

---

## 6. Le criblage ne trouve rien, et c'est un résultat

Mode fixe, n = 20 000, p = 272.

| | Y | Y_Δ |
|:---|---:|---:|
| seuil de permutation familial 5 % | 1,19056 | 1,19202 |
| **score minimal observé** | **1,19729** | **1,19711** |
| sous le seuil familial 5 % | **0 / 272** | **0 / 272** |
| sous le seuil par covariable 5 % | 31 (11,4 %) | 23 (8,5 %) |
| fréquence de sélection, médiane | 0,05 | 0,05 |
| fréquence de sélection, maximum | 0,25 | 0,25 |
| covariables sélectionnées systématiquement | **0** | **0** |

Trois mesures indépendantes concordent :

1. **Aucune covariable ne bat la loi nulle de permutation.** Le score le plus
   bas observé est *au-dessus* du 5ᵉ percentile de la loi nulle du minimum.
2. **Aucune sélection stable** sur 20 sous-échantillons disjoints : fréquence
   médiane 0,05, maximum 0,25.
3. **La corrélation de Spearman entre les classements de Y et de Y_Δ vaut
   0,016.** Les deux réponses portent sur les mêmes documents et les mêmes
   covariables ; si l'une avait du signal, un accord partiel serait attendu.
   Zéro est la signature de deux classements de bruit indépendants. Pour
   comparaison, l'accord entre l'écran indice de queue et quantile SIS à
   τ = .95 vaut 0,428.

L'excès de 11,4 % sous le seuil par covariable (contre 5 % attendus) n'est pas
un test propre : les composantes d'embedding sont fortement corrélées entre
elles.

**Ce que la dégénérescence aurait donné, et pourquoi ce n'est pas le cas.** Si
une classe de déchets ubiquitaire — hashes, base64, tableaux de nombres —
atteignait le maximum global à toutes les valeurs de tous les attributs, alors
ξ_j ≡ γ\* et l'écran ne verrait rien, exactement comme ici. Le contrôle écarte
cette explication : sur les 300 documents de plus grand max_t S, on compte
**291 tokens distincts, dont 278 alphabétiques** contre 12 non-ASCII, 5 chiffres
et 5 de ponctuation — des mots ordinaires dans des contextes imprévisibles
(`' hereby'`, `' Olymp'`, `' emphasizes'`). Et ces documents sont répartis
presque uniformément sur les quartiles de toutes les covariables (0,18 à 0,31
contre 0,25 pour l'uniforme).

**Lecture.** Les attributs de surface et d'embedding prédisent le **niveau** de
la surprise, pas son **indice de queue**. C'est cohérent avec la constance de γ
par strate observée en 4.3. C'est un énoncé qui soutient la thèse des deux
cibles distinctes — à condition de ne pas le présenter comme une démonstration
de puissance.

---

## 7. Un faux positif mesuré, avec vérité analytique

Le même écran, en mode libre, où T varie :

| | mode FIXE (T neutralisé) | mode LIBRE (T varie) |
|:---|---:|---:|
| rang médian du bloc longueur | **75,4 %** (dernier quart) | **18,8 %** (tête) |
| fréquence de sélection | 0,05 – 0,15 | **1,00** |
| rang médian, quantile SIS τ = .95 | 26,8 % | 18,8 % |

En mode fixe l'écran relègue les cinq proxys de longueur aux rangs 190 à 215
sur 272 et ne les sélectionne presque jamais — comportement correct, T y étant
constant. En mode libre **il les remonte en tête avec une stabilité parfaite**,
alors que la section 5 établit que la longueur n'affecte pas γ.

Le mécanisme est identifiable : le score n'estime pas l'enveloppe supérieure de
population mais une moyenne de Hill locaux à k fini. Conditionner sur un proxy
de longueur homogénéise l'échelle dans la fenêtre, ce qui déplace le **biais de
second ordre** de l'estimateur sans rien changer à l'indice. **L'écran classe une
variation de biais, pas une variation de γ.**

Deux nuances. D'abord, le test de permutation ne déclare rien de significatif en
mode libre non plus : l'écran range la longueur en tête d'un classement dont
aucun élément ne bat le bruit. Sans calibration on aurait lu ce classement comme
un résultat. Ensuite, quantile SIS est **moins** déplacé que l'écran indice de
queue entre les deux modes (26,8 % → 18,8 % contre 75,4 % → 18,8 %) : la
dissociation attendue — quantile SIS piégé, écran immunisé — ne se produit pas
dans ce sens-là sur ces données.

---

## 8. Un point méthodologique sur l'agrégation

L'agrégation par **rang minimum** sur les neuf réglages, validée dans l'étude de
simulation *derrière une présélection*, produit ici un artefact net :
`emb_000` passe du rang 229 au réglage unique au rang 3 après agrégation. C'est
exactement le mode d'échec identifié en simulation — le rang minimum promeut une
coordonnée dès qu'un seul réglage la classe bien, ce qui est sans danger face à
25 survivantes mais promeut les nulles chanceuses face à 272 coordonnées sans
signal. **Sans présélection, la médiane est la règle appropriée** ; les deux sont
rapportées.

---

## 9. Limites

- **Une seule paire (corpus, modèle).** γ est une propriété de la paire, pas du
  corpus seul. La stabilité du jeu actif entre familles de modèles est une
  question ouverte, et intéressante.
- **La date de crawl n'est pas la date du contenu.** Le dump 2025 réduit le
  risque de contamination sans l'annuler.
- **Vérité partielle par sources non mesurée.** Composer le corpus de K sources
  connues et estimer γ par source sur un échantillon indépendant reste à faire ;
  c'est ce qui donnerait un ancrage positif au résultat négatif.
- **Étage 2 sur hold-out non fait.**
- **Le résultat négatif n'est pas une preuve d'absence.** Il établit qu'aucun
  attribut *de ce jeu* ne module γ de façon détectable *à ce n*, avec la
  puissance que la calibration par permutation quantifie.

---

## 10. Reproduction

```sh
python code/py/extract_surprisals.py results/llm/main_fixe.parquet 20000 fixe
python code/py/extract_surprisals.py results/llm/main_libre.parquet 20000 libre
python code/py/diagnostics.py       results/llm/main_fixe.parquet
python code/py/check_fp64.py        results/llm/main_fixe.parquet 200
python code/py/nested_prefixes.py   results/llm/nested.csv 600
python code/py/length_is_scale.py   results/llm/pilot_libre.parquet results/llm/pilot_fixe.parquet
python code/py/embed_covariates.py  results/llm/main_embed.parquet 20000 256
Rscript  code/R/screen_llm.R        results/llm/main_fixe.csv results/llm/screen 92 1000
```
