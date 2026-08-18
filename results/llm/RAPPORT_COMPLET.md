# Application LLM — rapport détaillé

Ce document décrit intégralement le travail mené sur l'application « indice de
queue de la surprise par token » : ce qui a été mesuré, sur quelles données,
avec quels modèles, et ce que chaque résultat établit ou n'établit pas pour
l'article.

Il remplace `results/llm/rapport.md`, dont les sections 6 et 7 ont été écrites
avant les mesures à grande échelle et disent le contraire de ce qu'on observe
maintenant.

---

## 1. Pourquoi cette application, et ce qu'elle doit à l'article

L'article propose un écran de sélection de variables dont la cible est
**l'indice de queue conditionnel** γ(x), et non un quantile fixe. Sa thèse
centrale, portée par le modèle de simulation M2 et par la Section 4.3, est que
les deux cibles se dissocient : une covariable peut déplacer tous les quantiles
finis sans entrer dans l'exposant limite, auquel cas un écran par quantiles la
retient à tort et l'écran par indice de queue ne la retient pas.

Le défaut de l'application existante (Communities and Crime) est qu'aucune
vérité n'y est connue : on lit des rangs sans pouvoir dire s'ils signifient
quelque chose. L'application LLM a été construite pour fournir ce qui manque —
**une nuisance d'échelle dont le statut est établi analytiquement**, sur données
réelles.

---

## 2. Le point de départ : un piège dans la définition de la réponse

L'idée naturelle serait de prendre la perplexité moyenne d'un document,

    PPL = exp{ T⁻¹ Σ_t S_t },   S_t = −log p_θ(w_t | w_<t).

C'est une moyenne de T surprises : sous mélange, log PPL se concentre, PPL est
approximativement lognormale, donc **γ = 0** et l'on est dans le domaine de
Gumbel, sans variation régulière. Pire, le marginal peut *paraître* lourd par
mélange sur des types de documents hétérogènes pendant que chaque loi
conditionnelle reste lognormale : un diagnostic de Hill marginal rassurerait à
tort, la condition (C1) serait violée et ξ_j ≡ 0 — l'écran classerait du bruit.

**Réponse retenue : le maximum de surprise, exponentié.**

    Y = exp( max_{t ≤ T} S_t ) = 1 / min_t p̂_θ(w_t | w_<t).

Le mécanisme est exact : si P(S > s | x) ~ c(x)·e^{−s/γ(x)}, alors
P(Y > y | x) ~ c(x)·y^{−1/γ(x)} — variation régulière **conditionnelle** avec
γ(x) > 0 par construction, et le facteur sous-exponentiel devient exactement le
ℓ_j(s,u) du cadre de l'article. Y > 0 automatiquement, donc aucun filtrage de
zéros.

Une queue exponentielle pour S est ce qu'on attend d'un vocabulaire zipfien :
p_(r) ∝ r^{−s} donne P(S > u) ∝ e^{−u(s−1)/s}, soit γ = s/(s−1). Le
conditionnement du modèle raidit ce calcul unigramme — on mesure γ ≈ 1,24, ce
qui correspond à s ≈ 4,2.

**Seconde réponse, extraite dans la même passe :**

    Y_Δ = exp( max_t ( S^petit_t − S^grand_t ) ) = max_t p_grand / p_petit,

dont l'indice est le taux de la queue du rapport de vraisemblance maximal.
Interprétation : le pire excès de perte du petit modèle sur la référence croît
comme γ_Δ·log T.

---

## 3. Les données

### 3.1 Corpus

**FineWeb, dump `CC-MAIN-2025-18`**, lu en flux depuis le Hub.

Deux écueils à éviter, et pourquoi ce choix les évite :

**Le filtrage par perplexité est disqualifiant.** CCNet et CC-100 (buckets
KenLM), FineWeb-Edu et les classifieurs de qualité tronquent exactement la queue
supérieure de Y : c'est une sélection sur la réponse, elle invaliderait tout.
FineWeb applique des filtres heuristiques et de la déduplication, **pas** de
seuil de perplexité ni de classifieur de qualité.

**La contamination.** Le Pile, sur lequel Pythia est entraîné, a été collecté
jusqu'en 2020. Un crawl de 2025 est très postérieur, donc la non-appartenance
est acquise par construction plutôt que vérifiée péniblement.
**Réserve honnête** : une date de crawl récente ne garantit pas un *contenu*
récent — une page ancienne peut être recrawlée. Le risque est réduit, pas
annulé.

**Pas de découpage en morceaux.** Un document du corpus donne une ligne ; on ne
découpe jamais un document en plusieurs observations, ce qui violerait
l'indépendance supposée par la théorie.

### 3.2 Volumes

| jeu | n | lus | usage |
|:---|---:|---:|:---|
| pilote T fixé | 2 000 | 5 064 | mise au point, diagnostics |
| pilote T libre | 2 000 | 2 584 | première vue sur la longueur |
| **principal T fixé** | **20 000** | 50 405 | analyse primaire |
| **principal T libre** | **20 000** | 25 853 | démonstration de la nuisance d'échelle |
| préfixes emboîtés | 600 | — | effet pur de T, composition fixe |
| **états internes** | **50 000** | 194 049 | criblage des activations |

Le rendement de lecture (40 % en mode fixe, 26 % pour les internes) est imposé
par le filtre de longueur : il faut 577 tokens pour l'analyse primaire et 833
pour le découpage passé/futur.

---

## 4. Modèles et protocole de mesure

| | |
|:---|:---|
| modèle évalué | `EleutherAI/pythia-410m-deduped` |
| modèle de référence | `EleutherAI/pythia-1.4b-deduped` |
| précision | **fp32 de bout en bout**, logits castés avant softmax |
| déterminisme | **batch = 1** |
| burn-in | 64 tokens écartés en tête |
| tokens spéciaux | exclus |
| contexte | 2 048 |

**Pourquoi Pythia et non les Qwen en cache.** Le protocole demande un modèle à
données d'entraînement publiques, pour *vérifier* l'appartenance plutôt que la
supposer. Les données d'entraînement des Qwen ne sont pas publiques.

**Pourquoi deux modèles de la même famille.** Les tokenizers sont vérifiés
identiques **par assertion dans le code** : sans cela, S^petit et S^grand ne
seraient pas alignées token par token et la différence n'aurait aucun sens.

**Pourquoi fp32.** Les plus grandes surprises sont précisément les observations
que la précision réduite dégrade le plus, et ce sont elles qui pilotent Hill.
La vérification est en section 6.4.

### 4.1 Les deux découpages

**Analyse primaire, T fixé.** Exactement 512 tokens notés par document, après le
burn-in. La longueur est neutralisée par construction : aucun confondant.

**Analyse secondaire, T libre.** Le document entier jusqu'à la limite de
contexte. La longueur devient la nuisance d'échelle à étudier.

---

## 5. Les covariables

**Règle absolue (point 5 du protocole) : jamais issues de la passe forward qui
produit Y.** Si X était une fonction déterministe des mêmes logits que Y, la loi
conditionnelle pourrait dégénérer.

**22 statistiques de surface**, calculées depuis le texte brut :

- *amas longueur* (8) : caractères, octets, mots, phrases, paragraphes, espaces,
  lignes, taille gzip ;
- *typographie* (7) : taux de chiffres, majuscules, ponctuation, non-ASCII,
  espaces, retours à la ligne, alphabétiques ;
- *lexique* (5) : TTR, taux d'hapax, longueur moyenne et maximale des mots,
  longueur moyenne de phrase ;
- *compression* (2) : ratio gzip, entropie par octet.

**256 composantes d'ACP** des embeddings de `roberta-large`, gelé — une famille
distincte des Pythia qui produisent Y.

**Filtre de continuité**, plus strict que dans l'application crime : au moins
100 valeurs distinctes **et** aucun atome portant plus de 2 % de la masse. Il
écarte `rate_digit`, `rate_nonascii`, `len_sents`, `len_paras`, `len_lines` et
`max_word_len` — les taux massivement nuls et les comptages entiers de faible
cardinalité. **272 covariables sur 278** sont retenues.

**Départage des ex æquo** : aléatoire à graine fixée. Le score travaille sur des
fenêtres de rangs, donc un départage par ordre d'apparition créerait une
géométrie artificielle.

### 5.1 Régime atteint

Au réglage (a\*, b\*) = (0,30 ; 0,15) issu de la phase de réglage sur la suite
de modèles figée :

| jeu | n | p | α | h | **2nαh** | n/p |
|:---|---:|---:|---:|---:|---:|---:|
| primaire | 20 000 | 272 | 0,0512 | 0,1132 | **232** | 74 |
| internes | 50 000 | 1 024 | 0,0389 | 0,0987 | **384** | 49 |
| internes, 3 couches | 50 000 | 3 072 | 0,0389 | 0,0987 | **384** | 16 |

Le protocole visait ≈ 230 extrêmes locaux par fenêtre ; c'est atteint et
dépassé.

---

## 6. Les prérequis, vérifiés avant tout criblage

### 6.1 Pas de plafond numérique

max_t S s'étale de **2,35 à 27,92 nats** (médiane 13,89). La part de documents à
moins d'un nat du maximum global vaut **0,0001**, soit deux documents. Les
écarts entre les vingt plus grands maxima restent réguliers (médiane 0,143 nat).
Une queue tronquée se verrait comme une accumulation sous la borne ; il n'y en a
pas.

### 6.2 γ > 0, marginalement

| réponse | Hill (k = 10 %) | ξ (GPD) | IC95 profilé |
|:---|---:|---:|:---|
| Y, pythia-410m | 1,276 | 1,159 | **[1,073 ; 1,251]** |
| Y, pythia-1.4b | 1,291 | 1,206 | **[1,113 ; 1,305]** |
| Y_Δ | 1,338 | 1,176 | [0,919 ; 1,501] |

L'IC de vraisemblance profilée sur ξ exclut 0 très largement. Sur le pilote à
n = 2 000 il valait [0,93 ; 1,54] ; à n = 20 000 il se resserre à [1,07 ; 1,25].

### 6.3 γ > 0 **par strates** — le seul test qui porte sur (C1)

γ marginal > 0 est compatible avec γ(x) = 0 partout, par mélange. Le test
conditionnel est donc indispensable. Sur les quartiles de `gzip_ratio`,
`rate_digit`, `ttr` et `byte_entropy`, soit **seize strates de 5 000 documents** :
ξ ∈ [1,14 ; 1,29], bornes basses d'IC de 1,02 à 1,16. **Aucun intervalle ne
touche 0.** Le scénario redouté est écarté.

### 6.4 fp32 suffisant

Recalcul des 200 plus grandes surprises en **fp64, modèle compris** :

| contrôle | résultat |
|:---|:---|
| reproduction du stockage fp32 | écart max **0,00e+00** |
| écart fp64 − fp32 | médiane −3,6e−06, max **8,1e−04 nat** |
| corrélation de rang | Kendall **0,9998**, Spearman **0,99999** |
| rangs identiques | 196 / 200 |
| vingt premiers identiques | **20 / 20** |

Cinq ordres de grandeur séparent l'erreur numérique de l'échelle des valeurs.

### 6.5 Pas de dégénérescence par une classe ubiquitaire

La Proposition 2.2 dit qu'une coordonnée est invisible si toutes ses fibres
rencontrent l'argmax global. Sur corpus web, une classe de déchets présente à
toutes les valeurs de tous les attributs (hashes, base64) produirait exactement
ce cas, et se manifesterait comme un résultat nul.

Contrôle sur les 300 documents de plus grand max_t S : **291 tokens distincts,
dont 278 alphabétiques** contre 12 non-ASCII, 5 chiffres et 5 de ponctuation —
des mots ordinaires dans des contextes imprévisibles (`' hereby'`, `' Olymp'`,
`' emphasizes'`). Et ces documents sont répartis presque uniformément sur les
quartiles de toutes les covariables (0,18 à 0,31 contre 0,25). **La
dégénérescence n'est pas l'explication.**

---

## 7. La longueur : une nuisance d'échelle à vérité analytique

Pour T variables à variation régulière faiblement dépendantes, d'indice
extrémal θ,

    P(max_{t≤T} > y | x) ≈ θ·T·c(x)·y^{−1/γ(x)},

donc **T entre dans l'échelle et jamais dans l'indice**. En logarithme, avec
M = max_t S_t = log Y :

    E[M | x] = γ(x)·log T + O(1),     Var[M | x] = π²γ(x)²/6.

C'est le point qui rend cette application adaptée à cet article-là : la longueur
arrive avec son amas de proxys corrélés — caractères, octets, mots, espaces,
lignes, taille compressée — **naturellement, sans construction artificielle**, et
son statut est un théorème et non une hypothèse de modélisation.

### 7.1 Vérification à composition strictement fixe

Le test décisif porte sur des **préfixes emboîtés du même document**, ce qui
élimine tout effet de composition — les documents longs ne sont pas des
documents courts allongés.

| T | 128 | 256 | 512 | 1024 |
|---:|---:|---:|---:|---:|
| E[max_t S] | 12,245 | 13,040 | 13,964 | 14,782 |

- pente **1,2315 nat/nat**, résidu maximal 0,041 nat ;
- IC95 intra-document **[1,149 ; 1,314]**, qui **contient les deux estimations
  indépendantes de γ** (Hill 1,24–1,30) ;
- croissance monotone dans **1,000** des 600 documents, comme l'exige un maximum
  sur un ensemble croissant.

La comparaison entre strates de longueur donnait une pente de 1,02 — l'écart
avec γ venait de la composition, et disparaît quand on la tient fixe.

### 7.2 L'indice ne bouge pas, la localisation oui

Sur cinq strates de longueur en mode libre : Hill de 1,240 à 1,579, **les cinq
IC se recouvrant**, pendant que la médiane de M se déplace de **2,24 nats**.

---

## 8. Le criblage

### 8.1 Analyse primaire (T fixé, n = 20 000, p = 272)

| | Y | Y_Δ |
|:---|---:|---:|
| seuil de permutation familial 5 % | 1,19056 | 1,19202 |
| score minimal observé | 1,19729 | 1,19711 |
| **sous le seuil familial** | **0 / 272** | **0 / 272** |
| sous le seuil par covariable 5 % | 31 (11,4 %) | 23 (8,5 %) |
| fréquence de sélection, médiane / max | 0,05 / 0,25 | 0,05 / 0,25 |
| Spearman entre les classements Y et Y_Δ | **0,016** | (les deux réponses) |

Aucune covariable ne franchit le seuil familial. **Sans la calibration par
permutation, on aurait lu un classement qui existe toujours, signal ou pas** —
c'est l'apport principal de cette application par rapport à la Section 5.

### 8.2 La dissociation, en mode libre

| | mode FIXE (T neutralisé) | mode LIBRE (T varie) |
|:---|---:|---:|
| rang médian du bloc longueur, écran indice de queue | 75,4 % | **81,2 %** |
| proxys de longueur dans le top-5 de quantile SIS | — | **5 / 5** aux trois τ |
| proxys de longueur dans le top-5 de l'écran | — | **0 / 5** |
| rang médian, écran / quantile .95 | — | **13,0** / **3,0** sur 16 |

La probabilité que les cinq proxys occupent les cinq premières places de 16 par
hasard est de 1/C(16,5) ≈ 2×10⁻⁴.

**Une borne de validité mesurée.** À n = 2 000 (2nαh = 65), l'écran mettait
au contraire la longueur **en tête**, avec une fréquence de sélection de 1,00.
L'effet disparaît à n = 20 000 (2nαh = 232). Le mécanisme est identifiable :
conditionner sur un proxy de longueur homogénéise l'échelle dans la fenêtre et
déplace le **biais de second ordre** du Hill local, pas l'indice ; ce biais
décroît avec le nombre d'extrêmes par fenêtre.

### 8.3 γ varie-t-il, tout court ?

Trois résultats nuls successifs peuvent venir d'un manque de puissance ou de
l'absence de variation à détecter. La question se tranche **sans l'écran** : on
compare la dispersion inter-strates de γ̂ pour une covariable à celle de strates
**aléatoires** de mêmes tailles, qui ne contiennent que le bruit
d'échantillonnage.

- loi nulle sur 400 découpages aléatoires : dispersion médiane **0,0372**
  (prédiction théorique γ/√k = 0,0427 — la calibration est saine) ;
- observée sur 276 covariables : médiane **0,0379**. La covariable typique ne
  disperse donc **pas du tout** γ ;
- mais 33 covariables à p < 0,05 pour 14 attendues, et **un cas net** :

| covariable | dispersion | p | γ min → max |
|:---|---:|---:|:---|
| **`rate_space`** | **0,1232** | < 0,0025 | 1,208 → 1,498 |
| `mean_word_len` | 0,0850 | 0,0025 | 1,241 → 1,446 |

Pour `rate_space`, les IC de Q1 et Q3 **ne se recouvrent pas** (1,239 contre
1,197). L'amplitude de 0,29 sur γ vaut **1,8 nat** sur la perte du pire token à
T = 512. Le profil n'est **pas monotone** : γ descend de Q1 à Q3 puis remonte,
donc « plus d'espaces, queue plus lourde » serait faux.

**La dissociation joue dans les deux sens :**

| | effet sur le niveau | effet sur l'indice |
|:---|:---|:---|
| bloc longueur | ρ = **+0,42**, capture 5/5 les écrans à quantile | **nul** |
| `rate_space` | ρ = +0,15 | **significatif**, amplitude 0,29 |

### 8.4 Les états internes du modèle

> **⚠ SECTION INVALIDÉE.** Le criblage décrit ci-dessous passait par
> `code/py/export_ranks.py`, qui écrivait la matrice de rangs en ordre C alors
> que `matrix()` de R la lit en ordre Fortran : l'écran lisait une matrice
> brouillée, et donc du bruit. Le résultat nul rapporté ici est sans valeur.
> Le pont est corrigé et un garde-fou vérifie désormais chaque colonne, mais
> **le criblage reste à relancer**. Voir `RAPPORT_LLM_EXTREMES.md` §7.

Question : quelles directions du flux résiduel annoncent que le modèle va
échouer catastrophiquement ?

**Le découpage passé/futur** évite la dégénérescence : X = activations moyennées
sur les tokens 64–320 (le passé), Y = exp(max_t S) sur les tokens 320–832 (le
futur). X est mesurable par rapport au passé, Y par rapport au futur — c'est un
problème de prédiction, pas une tautologie.

Huit couches extraites (3, 6, 9, 12, 15, 18, 21, 24 sur 24), car le flux
résiduel évolue lentement et des couches voisines sont quasi colinéaires — le
même phénomène que sur les grilles de réglage, où densifier n'augmentait pas le
nombre effectif de configurations indépendantes.

| config | p | score min | verdict | sous 5 % (attendu) | freq max |
|:---|---:|---:|:---|---:|---:|
| couche 3 | 1 024 | 1,17505 | non signif. | 48 (51) | 0,15 |
| couche 6 | 1 024 | 1,18010 | non signif. | 49 (51) | 0,15 |
| couche 12 | 1 024 | 1,17537 | non signif. | 62 (51) | 0,15 |
| couche 18 | 1 024 | 1,17399 | non signif. | 65 (51) | 0,15 |
| couche 24 | 1 024 | 1,17945 | non signif. | 52 (51) | 0,20 |
| ACP 1 024 | 1 024 | 1,17876 | non signif. | 58 (51) | 0,15 |
| **3 couches (6/12/18)** | **3 072** | 1,17836 | non signif. | **154 (154)** | 0,15 |

**Aucune profondeur, aucune base.** La configuration a trois couches est
exactement au niveau du hasard — 154 covariables sous le seuil par covariable
pour 154 attendues, 32 sous 1 % pour 31 attendues — alors que les couches
isolees montraient de legers exces (62 pour 51 a la couche 12). Regrouper trois
couches dilue donc l'exces au lieu de le renforcer, ce qui est la signature du
bruit et non d'un signal disperse.

Que l'ACP ne trouve pas plus que le brut
répond à une question précise : il ne s'agit pas d'un signal réel mal aligné sur
la base des neurones — il n'y a pas de direction, privilégiée ou non.

**Limite de conception à assumer** : les activations sont **moyennées** sur 256
tokens. L'état au dernier token du préfixe, ou un maximum plutôt qu'une moyenne,
capteraient autre chose. L'énoncé correct est donc « pas dans la moyenne du flux
résiduel », pas « pas dans le flux résiduel ».

---

## 9. Ce qui soutient la méthode, et ce qui ne la soutient pas

C'est la section qui compte pour l'article.

### 9.1 Ce qui la soutient : la spécificité

Le bloc longueur capture **intégralement** les trois écrans à quantile (5/5 dans
leur top-5) et occupe le dernier quart chez l'écran indice de queue, avec une
vérité **analytique**. Une covariable qui déplace massivement les quantiles
finis sans toucher l'indice ne trompe que l'un des deux. C'est le pendant réel
de M2, avec une nuisance bien plus violente que celle de M2 (2,24 nats de
décalage contre un sd(log ℓ) de 0,138 en simulation).

### 9.2 Ce qui ne la soutient pas : la puissance

On peut tester si le classement d'un écran suit la variation réelle de γ, en
utilisant le critère externe de la section 8.3, calculé **sans aucun écran** :

| écran | ρ(dispersion de γ, rang) | p permutation | dispersion médiane du top-20 |
|:---|---:|---:|---:|
| tail-index (a\*, b\*) | −0,240 | 0,00005 | 0,0534 |
| tail-index agrégé | −0,208 | 0,0007 | — |
| quantile τ=.90 | −0,226 | 0,0001 | — |
| **quantile τ=.95** | **−0,284** | < 0,00001 | **0,0534** |
| quantile τ=.99 | −0,234 | 0,0001 | — |

Les deux familles identifient les covariables qui modulent γ **avec la même
efficacité**, quantile SIS à τ=.95 faisant même légèrement mieux. Il n'y a donc
**aucun résultat montrant que l'écran vise γ mieux que ses concurrents**.

### 9.3 La concordance avec les simulations

C'est exactement ce que dit la campagne à 1 000 réplications : quantile SIS à
τ=.95 devance l'écran sur le Sure-20 dans les quatre modèles, et le seul
avantage de l'écran est la composition du top-4 sur M2 — 0,17 variable
d'échelle contre 2,19. **Même énoncé en simulation et sur données réelles :
l'avantage est la spécificité, pas la puissance.** Cette concordance est en
soi un soutien, à condition que l'article défende cette thèse-là.

**Ce qu'il ne faut pas écrire** : « notre méthode a identifié les attributs qui
gouvernent la queue ». Le test de la section 9.2 dirait le contraire, et c'est
nous qui l'aurions produit.

---

## 10. Un sous-produit qui n'utilise pas l'écran

À signaler parce qu'il est mesuré et qu'il pourrait servir ailleurs, mais il
faut être clair : **il n'utilise pas la méthode de criblage**, seulement
l'estimande et l'estimateur de Hill.

| | pythia-410m | pythia-1.4b | écart |
|:---|---:|---:|---:|
| surprise moyenne | 2,7749 | 2,5324 | **−0,2425 nat (−8,7 %)** |
| … améliorée sur | | | **100,0 % des documents** |
| indice de queue γ | 1,2758 | 1,2910 | **+0,0153** |
| IC95 sur ξ | [1,073 ; 1,251] | [1,113 ; 1,305] | se recouvrent |

Tripler les paramètres améliore le niveau sur les 20 000 documents **sans
exception**, et ne change pas l'indice. Via E[max_t S] = γ·log T + O(1) : la
capacité déplace le O(1), pas le γ. Un seul pas d'échelle, donc à confirmer sur
l'échelle Pythia complète avant d'en faire quoi que ce soit.

---

## 11. Ce qui reste à faire

- **Un second amas de nuisance d'échelle**, indépendant de la longueur et à
  vérité également analytique, pour doubler la seule démonstration qui porte.
- **Vérité partielle par K sources connues** (prose, code, math, multilingue),
  avec γ estimé par source sur un échantillon indépendant.
- **Étage 2 sur hold-out** (Wang–Tsai ou sous-espace de Gardes–Podgorny).
- **Variantes du résumé du préfixe** : dernier token, maximum, pour lever la
  limite de la section 8.4.
- **Variante tokenizer au niveau mot**.
- **n ≈ 100 000** si l'on veut que `rate_space` franchisse le seuil familial :
  Δ/σ croît en n^0,275.
- **Rétroporter la calibration par permutation sur la Section 5** (crime), qui
  n'a aujourd'hui aucune calibration inférentielle de Δ̂_j. Deux issues, toutes
  deux utiles à connaître avant un rapporteur.

---

## 12. Fichiers et reproduction

**Code** — `code/py/` :

| fichier | rôle |
|:---|:---|
| `extract_surprisals.py` | surprises des deux modèles, covariables de surface |
| `extract_internal.py` | activations du préfixe, réponse sur le suffixe |
| `embed_covariates.py` | embeddings d'un encodeur séparé et gelé, ACP |
| `diagnostics.py` | plafond, GPD + IC profilé, strates |
| `check_fp64.py` | recalcul fp64 des 200 plus grandes surprises |
| `nested_prefixes.py` | effet pur de T, composition fixe |
| `length_is_scale.py` | la longueur agit-elle sur l'échelle seule |
| `gamma_varies.py` | γ varie-t-il, calibré par strates aléatoires |
| `export_ranks.py` | pont binaire vers R (les CSV faisaient plusieurs Go) |

**Code** — `code/R/` : `screen_llm.R` (criblage, permutation, sous-échantillons),
`screen_internal.R` (idem pour les activations).

```sh
python code/py/extract_surprisals.py results/llm/main_fixe.parquet  20000 fixe
python code/py/extract_surprisals.py results/llm/main_libre.parquet 20000 libre
python code/py/extract_internal.py   results/llm/internal.parquet   50000
python code/py/embed_covariates.py   results/llm/main_embed.parquet 20000 256
python code/py/diagnostics.py        results/llm/main_fixe.parquet
python code/py/check_fp64.py         results/llm/main_fixe.parquet 200
python code/py/nested_prefixes.py    results/llm/nested.csv 600
python code/py/gamma_varies.py       results/llm/main_joint.parquet max_S_small 400 4
Rscript  code/R/screen_llm.R         results/llm/main_fixe.csv results/llm/screen 92 1000
python code/py/export_ranks.py       results/llm/internal.parquet results/llm/ranks c12
Rscript  code/R/screen_internal.R    results/llm/ranks c12 results/llm/screen_internal 88 1000
```

Les `.parquet` volumineux ne sont pas versionnés ; les résumés, journaux et
tableaux de criblage le sont.
