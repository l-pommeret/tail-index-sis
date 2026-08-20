# Où l'écran par indice de queue domine-t-il ses concurrents ?

*Étude de puissance comparée sur l'application perplexité, dans le plan (n, a).*

Ce document définit les métriques employées, expose les résultats, et énonce
les limites de la démarche. Le protocole de l'application elle-même est décrit
dans [`APPLICATION_PERPLEXITE.md`](APPLICATION_PERPLEXITE.md).

---

## 1. Le problème que ces métriques résolvent

Sur données réelles il n'y a pas de vérité : on ne sait pas quelles covariables
font varier `γ`. Comparer des écrans exige donc un **critère de réussite**, et
tout critère mal choisi rend la comparaison circulaire.

Deux tentatives ont été écartées en chemin, et il faut les mentionner car elles
étaient tentantes.

**Le critère par variation de `γ` — écarté.** On peut mesurer, pour chaque
covariable, la dispersion inter-strates de l'estimateur de Hill, calibrée par
des strates aléatoires. Mais c'est un critère qui désigne *ce qu'un estimateur
naïf de `γ` juge variable*, puis note les écrans sur leur accord avec lui. Il
avantage mécaniquement ceux qui suivent la variation brute de `γ`, et rien n'y
garantit que les covariables ainsi désignées soient actives plutôt que du bruit.
Mesure de sa contamination : il corrèle **+0,504** avec l'étendue du profil de
Hill local et **+0,506** avec celle du quantile 0,95 — il est à moitié un
critère de localisation.

**La stabilité de notre propre écran — écartée aussi.** Définir les covariables
actives comme celles que l'écran par indice de queue classe en tête à tous les
`n` revient à se donner raison d'avance.

---

## 2. Les métriques retenues

### 2.1 La cible

> **cible** = les covariables que **les sept règles** placent dans leur top-10
> à `n = 100 000`, le plus grand échantillon disponible.

Sur `p = 264` covariables, la cible compte **2 éléments** : `s160m_max` et
`s160m_kurt`.

Trois propriétés la rendent utilisable :

- **symétrique** — aucune règle n'y est privilégiée, il faut l'accord des sept ;
- **définie sur le maximum d'information** — au plus grand `n`, là où toutes les
  méthodes sont le plus fiables ;
- **indépendante des échantillons évalués** — on la fixe une fois, puis on teste
  sur des sous-échantillons plus petits.

Elle n'est pas une vérité : c'est un consensus. Sa justification empirique tient
au balayage en `n` (§4.1) — les deux covariables gardent les rangs 1 et 3 sur un
facteur 20 en `n`, ce qu'aucune covariable tirée au hasard ne fait.

### 2.2 Les deux mesures de réussite

> **atteinte** (*hits*) = combien de covariables de la cible la règle place dans
> **son propre** top-10, à la taille `n` considérée. Entre 0 et 2.

> **pire rang** (*worst*) = le plus mauvais rang que la règle attribue à une
> covariable de la cible. **Petit = bien.** Une règle qui classe les deux cibles
> aux rangs 1 et 3 a un pire rang de 3.

Le pire rang est la mesure principale : elle est continue, donc elle discrimine
là où l'atteinte sature à 2 sur 2. C'est l'analogue direct du `Sure-d` des
simulations — la plus mauvaise position d'une active — à ceci près qu'ici les
« actives » sont un consensus et non une vérité.

Chaque cellule est une **moyenne sur des sous-échantillons répétés** : 20 tirages
jusqu'à `n = 5 000`, puis 10, 5 et 3 quand le coût de quantile SIS l'impose.

### 2.3 L'avantage, et ses deux versions

> **avantage** = (pire rang d'un concurrent) − (pire rang de l'écran de queue).
> Positif = l'écran de queue fait mieux.

Deux façons de choisir le concurrent, et **elles ne disent pas la même chose** :

| version | définition | statut |
|---|---|---|
| **oracle** | minimum sur les cinq concurrents, **calculé tirage par tirage** | irréaliste : suppose de savoir d'avance lequel gagnera sur *ce* jeu |
| **fixe** | une règle concurrente nommée, la même partout | ce qu'un praticien fait réellement |

L'écart entre les deux est considérable, et c'est le principal enseignement
méthodologique de cette étude : une comparaison contre un oracle sur les
concurrents est un adversaire que personne ne peut instancier.

---

## 3. Les sept règles comparées

| règle | famille |
|---|---|
| queue, sélectionnée | indice de queue, réglage unique `(a, 0,15)` |
| queue, agrégée | indice de queue, rang minimal sur le bloc `N₉` |
| Yoshida–Umezu, réglage de son article | `h = 1`, `k = ⌊0,072 n⌋` |
| Yoshida–Umezu, `h = 2` | `k = ⌊0,05 n⌋` |
| quantile SIS, `τ = 0,90 / 0,95 / 0,99` | quantile conditionnel |

Les cinq concurrents **ne dépendent pas de `a`**. C'est ce qui rend la grille
abordable : on les calcule une fois par tirage et l'on balaie `a` pour le seul
écran de queue, à sous-échantillon identique — la comparaison est donc appariée.

---

## 4. Résultats

### 4.1 La cible est-elle du signal ou du bruit ?

Sous-échantillons emboîtés, `α = 0,1` tenu constant, 400 permutations par `n`.

| n | 2nαh | seuil familial 5 % | score minimal | marge | nb significatives |
|---:|---:|---:|---:|---:|---:|
| 5 000 | 139 | 1,1715 | 0,9246 | +0,247 | 2 |
| 10 000 | 251 | 1,1576 | 0,9354 | +0,222 | 5 |
| 25 000 | 547 | 1,2167 | 0,9260 | +0,291 | 21 |
| 50 000 | 987 | 1,2344 | 0,9309 | +0,304 | 40 |
| 100 000 | 1 778 | 1,2610 | 0,9390 | +0,323 | 52 |

**Le score minimal est constant à 0,93 sur un facteur 20 en `n`** pendant que le
nombre de significatives passe de 2 à 52 : une estimation qui converge tandis
que le plancher de bruit se resserre. Du bruit dériverait.

Rangs sous l'écran de queue, par `n` :

| | 5 000 | 10 000 | 25 000 | 50 000 | 100 000 |
|---|---:|---:|---:|---:|---:|
| `s160m_max` | 1 | 1 | 1 | 1 | 1 |
| `s70m_max` | 2 | 2 | 2 | 2 | 2 |
| `s160m_kurt` | 3 | 3 | 3 | 3 | 3 |
| `s160m_mean` *(témoin)* | 146 | 175 | 216 | 251 | 259 |

Trois covariables tiennent le top-10 aux cinq tailles ; les témoins de
localisation s'enfoncent de façon monotone. En revanche **le corps du classement
est instable** — Spearman entre `n` consécutifs seulement +0,64 à +0,80, et
`s70m_skew` passe du rang 161 au rang 9. Seul le sommet est fiable.

### 4.2 Avantage contre l'**oracle** sur les concurrents

| a \ n | 1 000 | 2 000 | 5 000 | 10 000 | 25 000 | 50 000 |
|---|---:|---:|---:|---:|---:|---:|
| 0,15 | −34,5 | −16,1 | −8,2 | −0,6 | −0,6 | −0,33 |
| 0,20 | −1,2 | −8,7 | **+1,2** | −0,4 | 0,0 | −0,33 |
| **0,25** | −9,1 | −13,3 | **+1,7** | **+0,4** | 0,0 | −0,33 |
| 0,30 | −0,6 | −11,9 | −0,4 | **+0,2** | 0,0 | −0,33 |
| 0,35 | −13,0 | −8,0 | −7,0 | **+0,5** | 0,0 | −0,33 |
| 0,40 | −38,0 | −20,5 | −23,4 | **+0,5** | 0,0 | −0,33 |

`2nαh` correspondant :

| a \ n | 1 000 | 2 000 | 5 000 | 10 000 | 25 000 | 50 000 |
|---|---:|---:|---:|---:|---:|---:|
| 0,15 | 126 | 205 | 388 | 631 | 1 198 | 1 947 |
| 0,20 | 89 | 140 | 254 | 398 | 722 | 1 133 |
| 0,25 | 63 | 96 | **166** | 251 | 435 | 660 |
| 0,30 | 45 | 65 | 108 | 158 | 262 | 384 |
| 0,35 | 32 | 45 | 71 | 100 | 158 | 224 |
| 0,40 | 22 | 31 | 46 | 63 | 95 | 130 |

**Contre l'oracle, l'écran de queue ne gagne presque jamais**, et quand il gagne
c'est de 0,2 à 1,7 rang. La meilleure configuration est `(a = 0,25 ; n = 5 000)`,
soit `2nαh ≈ 166`.

Deux traits de forme méritent d'être notés. **Les bords sont mauvais** :
`a = 0,15` et `a = 0,40` sont les deux pires lignes partout — trop d'extrêmes
noie le signal local, trop peu ne l'estime plus. Et **au-delà de `n = 25 000` la
mesure sature** : toutes les règles placent les deux cibles dans les trois
premières, l'avantage tombe à zéro par construction. La grille ne dit rien de ce
régime ; il faudrait une cible plus large pour continuer à discriminer.

### 4.3 Avantage contre chaque concurrent **fixe**

Convention `α = 0,1` constante, règle agrégée :

| avantage contre | 1 000 | 1 500 | 2 000 | 3 000 | 5 000 | 10 000 |
|---|---:|---:|---:|---:|---:|---:|
| quantile 0,90 | −1,4 | −0,4 | +5,3 | +5,3 | **+15,5** | **+12,9** |
| quantile 0,95 | −21,3 | +0,8 | +0,1 | +13,8 | **+5,9** | **+5,2** |
| quantile 0,99 | +17,5 | +15,9 | +5,6 | +18,4 | +2,3 | +0,6 |
| Yoshida–Umezu, son réglage | +49,8 | +79,6 | +68,4 | +89,3 | **+101,8** | +65,1 |
| Yoshida–Umezu, `h = 2` | +37,4 | +53,2 | +61,5 | +66,3 | +76,9 | +84,3 |

Convention publiée `a = 0,35`, la même règle :

| avantage contre | 1 000 | 1 500 | 2 000 | 3 000 | 5 000 | 10 000 |
|---|---:|---:|---:|---:|---:|---:|
| quantile 0,90 | +7,8 | −3,4 | +1,9 | −16,6 | +13,5 | +12,4 |
| quantile 0,95 | −12,2 | −2,1 | −3,4 | −8,1 | +3,9 | +4,7 |
| quantile 0,99 | +26,6 | +13,0 | +2,1 | −3,5 | +0,2 | +0,1 |
| Yoshida–Umezu, son réglage | +58,9 | +76,7 | +65,0 | +67,4 | +99,8 | +64,6 |

**Contre un concurrent fixe, l'écran de queue domine dès `n ≈ 2 000–3 000`**, et
l'écart se creuse ensuite.

### 4.4 L'agrégation contre le réglage unique

Pire rang, convention publiée `a = 0,35` :

| | 1 000 | 1 500 | 2 000 | 3 000 | 5 000 | 10 000 |
|---|---:|---:|---:|---:|---:|---:|
| queue, **agrégée** | **38,4** | **33,8** | **24,2** | 35,8 | **5,9** | **4,0** |
| queue, sélectionnée | 49,4 | 49,9 | 32,5 | 45,9 | 17,0 | 5,3 |

**L'agrégation sur le bloc `N₉` bat le réglage unique à cinq tailles sur six**,
et sous les deux conventions. C'est un argument direct et chiffré pour `N₉`,
indépendant de tout le reste.

---

## 5. Réponse à la question posée

**Où l'écran par indice de queue est-il le meilleur ?**

`a ∈ [0,20 ; 0,25]` avec `n ≥ 5 000`, soit **`2nαh ≈ 165–250` extrêmes locaux**.
C'est la seule région où il domine même l'oracle, et c'est aussi là qu'il domine
le plus largement chaque concurrent fixe.

**Existe-t-il un `n` où il détecte mieux à données rares ?** Non. À
`n ≤ 3 000` — soit 32 à 90 extrêmes par fenêtre — quantile SIS fait aussi bien ou
mieux. L'avantage vient *avec* les extrêmes, il ne les précède pas. C'est
cohérent avec la nature de l'estimateur : un Hill local n'a rien à estimer sans
excès, alors qu'une régression quantile utilise tout l'échantillon.

---

## 6. Limites

**La cible n'a que deux éléments.** Sur `p = 264`, c'est peu, et cela sature la
mesure dès `n = 25 000`. Une cible plus large exigerait un critère moins exigeant
que l'accord des sept règles, donc plus discutable.

**Le consensus n'est pas une vérité.** Si les sept règles partagent un biais, la
cible en hérite. La stabilité en `n` (§4.1) réduit ce risque sans l'annuler.

**`a = 0,25` bat `a* = 0,35` sur ce jeu**, à tous les `n`. Le réglage publié a
été calibré sur la suite de simulations, pas ici ; l'écart n'est pas une
contradiction mais il doit être connu.

**Un seul jeu de données, un seul type de réponse.** Rien ici ne dit que
`2nαh ≈ 165–250` soit un seuil transportable. C'est une observation sur
l'application perplexité, pas une règle générale.

---

## 7. Reproduction

```sh
Rscript code/R/wild_nsweep.R      results/wild/ppl100k_ranks results/wild/nsweep 40 400
Rscript code/R/wild_nsweep_all.R  results/wild/ppl100k_ranks results/wild/nsweep_all 40
Rscript code/R/wild_smalln.R      results/wild/ppl100k_ranks \
        results/wild/nsweep_all/nsweep_all.rds results/wild/smalln 40 20
Rscript code/R/wild_smalln.R      results/wild/ppl100k_ranks \
        results/wild/nsweep_all/nsweep_all.rds results/wild/smalln_a035 40 20 fixed_a
Rscript code/R/wild_grid_na.R     results/wild/ppl100k_ranks \
        results/wild/nsweep_all/nsweep_all.rds results/wild/grid_na 40
```

Graines : sous-échantillons emboîtés 31415 ; tirages répétés `90000 + n + r`
(petits `n`) et `120000 + n + r` (grille) ; départage d'ex æquo 811000033 ;
permutations `606000 + b`.

**Avertissement conservé du journal.** `numpy.tofile` écrit en ordre C,
`matrix()` de R lit en ordre Fortran : tout script lisant la matrice de rangs
vérifie que chaque colonne est une permutation de `1..n` avant de scorer.

---

## Annexe — liste exhaustive des covariables

`p = 264` covariables, soit **88 par source × 3 sources**. Toutes sont
calculées sur le **même span de 512 jetons** que la réponse, et toutes passent le
filtre de continuité (≥ 100 valeurs distinctes, aucun atome > 2 %) — c'est
l'intérêt d'une NLL, qui est un flottant.

### Les trois sources

| préfixe | définition |
|---|---|
| `s70m_` | profil de surprisal `S_t = −log p_θ(x_t \| x_<t)` sous **pythia-70m** |
| `s160m_` | idem sous **pythia-160m** |
| `d70_160_` | **écart** `S_t^{70m} − S_t^{160m}`, jeton par jeton |

La réponse `Y = exp(max_t S_t)` est calculée sous **pythia-410m**, qui
n'intervient dans aucune covariable : aucun quantile de `X` n'est fonction du
vecteur de NLL dont `Y` est le maximum.

### Les cinq familles, par source

| famille | nombre | noms | définition |
|---|---:|---|---|
| quantiles | 49 | `q02` … `q98` | quantiles empiriques du profil de NLL, niveaux 0,02 à 0,98 par pas de 0,02 |
| moments | 6 | `mean`, `sd`, `max`, `min`, `skew`, `kurt` | moyenne, écart-type, maximum, minimum, asymétrie, kurtosis du profil |
| moyennes par bloc | 16 | `bm00` … `bm15` | moyenne de la NLL sur le bloc de 32 jetons consécutifs numéro *b* |
| écarts-types par bloc | 16 | `bs00` … `bs15` | écart-type sur le même bloc — dynamique locale de la surprisal |
| tranche lexicale | 1 | `fb00` | NLL moyenne des jetons d'identifiant BPE < 2 515, c'est-à-dire les plus fréquents |
| | **88** | | |

### Liste explicite des quantiles

Les 49 niveaux présents : `q02`, `q04`, `q06`, `q08`, `q10`, `q12`, `q14`, `q16`, `q18`, `q20`, `q22`, `q24`, `q26`, `q28`, `q30`, `q32`, `q34`, `q36`, `q38`, `q40`, `q42`, `q44`, `q46`, `q48`, `q50`, `q52`, `q54`, `q56`, `q57`, `q60`, `q62`, `q64`, `q66`, `q68`, `q70`, `q72`, `q74`, `q76`, `q78`, `q80`, `q82`, `q84`, `q86`, `q88`, `q90`, `q92`, `q94`, `q96`, `q98`.

**Défaut de nommage à connaître.** Le niveau 0,58 est étiqueté `q57` et non
`q58`, parce que `int(0.58*100)` vaut 57 en virgule flottante
(`0.58*100 = 57.99999999999999`). Les 49 quantiles sont bien tous présents ;
seul le nom d'un d'entre eux est décalé. Cela n'affecte aucun résultat, mais
rend la lecture des noms trompeuse entre `q56` et `q60`.

### Une famille qui a échoué

La famille **tranches lexicales** devait compter 20 covariables par source :
la NLL moyenne des jetons dont l'identifiant BPE tombe dans chacune des 20
tranches égales de `[0, 50 304]`, les identifiants GPT-NeoX étant grossièrement
ordonnés par fréquence — c'était un proxy de rareté lexicale.

**19 des 20 ont été écartées**, et pour une raison qui tient à mon
implémentation, pas aux données : la moyenne d'une tranche vaut `NaN` dès qu'un
document y compte moins de 3 jetons, et la colonne entière est supprimée si
**un seul** des 100 000 documents est dans ce cas. Seule `fb00`, la tranche des
jetons les plus fréquents, est toujours peuplée. Une imputation ou un seuil par
colonne aurait conservé la famille.

Elle n'aurait probablement rien changé au résultat — `fb00` est classée 74e, 61e
et 138e selon la source — mais l'inventaire doit dire ce qui manque et pourquoi.

### Les covariables citées dans ce document

| nom | définition | rôle |
|---|---|---|
| `s160m_max` | maximum du profil de surprisal de pythia-160m | **cible**, rang 1 à tous les *n* |
| `s70m_max` | idem pour pythia-70m | rang 2 à tous les *n* |
| `s160m_kurt` | kurtosis du profil de pythia-160m | **cible**, rang 3 à tous les *n* |
| `s70m_kurt` | kurtosis du profil de pythia-70m | 4e à `n = 10⁵`, mais 44e à `n = 5 000` |
| `s160m_mean` | surprisal moyenne sous pythia-160m | **témoin de localisation** |
| `s70m_mean` | surprisal moyenne sous pythia-70m | **témoin de localisation** |
| `s160m_skew` | asymétrie du profil | instable : 77e à `n = 10 000`, 5e à `n = 10⁵` |
| `s70m_skew` | asymétrie du profil | instable : 161e à `n = 5 000`, 9e à `n = 10⁵` |

Les **témoins de localisation** sont ainsi qualifiés sur mesure, non par
hypothèse : le long du rang de `s160m_mean`, le quantile 0,95 de `log Y` se
déplace de 1,4 nat tandis que l'étendue du profil de Hill local n'est que de
0,163 — contre 0,449 pour `s160m_max`.

