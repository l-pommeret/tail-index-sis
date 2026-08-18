# Application : l'indice de queue de la surprisal d'un modèle de langue

*Criblage par indice de queue contre quantile SIS et Yoshida–Umezu, sur 100 000
documents de FineWeb.*

---

## 1. Ce que l'application apporte

L'application Communities and Crime laisse deux choses ouvertes. La réponse y
est à queue lourde par constatation, pas par construction ; et `n = 1993` est
ce qu'il est, donc le régime de taille finie de la section 4.2 ne peut pas être
quitté.

Celle-ci corrige les deux points, et ajoute une dissociation que le jeu
socio-économique ne permet pas de mettre en évidence, faute de savoir quelles
covariables déplacent quoi.

---

## 2. La réponse, et pourquoi son statut est structurel

Soit `S_t = −log p_θ(x_t | x_<t)` la surprisal du jeton `t` sous un modèle de
langue. Si la queue supérieure de `S` est de type exponentiel,
`P(S > x) ~ C e^{−λx}`, alors

```
Y = e^S   vérifie   P(Y > y) ~ C y^{−λ}
```

soit Pareto d'indice `γ = 1/λ` **exactement**. L'indice de queue de `Y` se lit
donc comme le taux de décroissance de la distribution de surprise du modèle, et
`γ(x)` grand signifie « ce profil de covariables admet des jetons
catastrophiquement surprenants » — le critère de curation de données.

On prend `Y = exp(max_t S_t)` sur une fenêtre de `T = 512` jetons : le maximum
d'un échantillon à variation régulière l'est aussi, avec le même indice.

### Vérification

`n = 100 000`, `pythia-410m`, fp32.

| k/n | 0,01 | 0,02 | 0,05 | 0,10 | 0,15 | 0,20 | 0,50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| γ (Hill) | 1,227 | 1,188 | 1,222 | 1,271 | 1,312 | 1,349 | 1,615 |
| écart-type | 0,039 | 0,027 | 0,017 | 0,013 | 0,011 | 0,010 | 0,007 |

Plateau entre 1,19 et 1,35 sur `k/n ∈ [0,01 ; 0,20]`, soit un facteur 20 sur le
nombre d'excès. QQ exponentiel des excès : `r = 0,9993 / 0,9996 / 0,9993` à
`k/n = 0,05 / 0,10 / 0,20`, de pente 1,22–1,27 — cohérente avec Hill, la pente
d'un QQ exponentiel *étant* l'excès moyen.

---

## 3. Le protocole

### Neutralisation de la longueur

`E[max_t S_t] = γ log T + O(1)` : la longueur est une nuisance d'échelle. On
score donc **exactement** `T = 512` jetons par document, après 64 jetons
d'amorce non scorés (les premiers jetons d'un document ont une surprisal
artificiellement haute, le modèle n'ayant pas de contexte). Tous les documents
sont ainsi à `T` identique et la nuisance disparaît par construction.

### Les covariables

`X` = profil de surprisal de `pythia-70m` et `pythia-160m` **sur le même span de
512 jetons**, plus leur écart. Par profil on entend 49 quantiles, les quatre
premiers moments, moyennes et écarts-types sur 16 blocs de 32 jetons, et
moyennes par tranche d'identifiant BPE.

`p = 264` après retrait des colonnes à valeur manquante (88 × 3 blocs).

**Non-dégénérescence.** `Y` ne dépend que de `pythia-410m`, `X` que des deux
petits modèles. Aucun quantile de `X` n'est fonction du vecteur de NLL dont `Y`
est le maximum.

**Pourquoi la perplexité et pas des covariables de surface.** Un premier essai
avec taux de classes de caractères, ponctuation, répétition n-gramme,
compression et trigrammes n'a retenu que 73 colonnes sur 512 au filtre de
continuité : sur une fenêtre de 2 600 caractères, les taux rares et les
trigrammes à `df < 0,98` portent un atome de zéros supérieur à 2 %. Une NLL est
un flottant — **264 sur 264 passent le filtre**. C'est aussi le choix de
features de Wu et al. (arXiv 2509.23488), qui les emploient en moyenne
conditionnelle ; on garde leurs features et on change de régime.

### Le régime

| réglage | α | 2nαh | nαhΔ² | log(pn) | rapport |
|---|---:|---:|---:|---:|---:|
| `a* = 0,35` (manuscrit) | 0,018 | 316 | 5,47 | 17,09 | 0,32 |
| `a = 0,20` | **0,100** | **1 778** | **30,76** | 17,09 | **1,80** |

Point de méthode à noter. Le réglage `a* = 0,35` a été choisi par l'étude de
tuning **à n = 2000**, où il donne `α = 0,070`. Comme `α = n^{−a}`, garder `a`
fixe fait chuter `α` quand `n` grandit : 0,018 à `n = 10⁵`, soit un tout autre
point de fonctionnement. **Ce qui se transporte d'un `n` à l'autre est `α`, pas
`a`.** À `α ≈ 0,1`, la condition `log(pn) = o(nαhΔ²)` est satisfaite pour la
première fois — simulations comprises.

---

## 4. Résultat du criblage

Calibration par 1000 permutations de la réponse ; règles de classement
invariantes par permutation des colonnes (audit du départage d'ex æquo).

| | `a = 0,20` | `a = 0,35` |
|---|---:|---:|
| score minimal observé | **0,93896** | **0,89381** |
| seuil familial 5 % | 1,26160 | 1,15892 |
| sous le seuil familial | 52 | 8 |
| sous 5 % par covariable | 142 (attendu 13) | 24 (attendu 13) |
| sous 1 % par covariable | 106 (attendu 3) | 12 (attendu 3) |

Les quatre premières sont **identiques aux deux réglages et aux trois règles**
(sélectionnée, agrégée, quantile 0,95) :

| rang | covariable | fréquence de sélection | p marginal |
|---:|---|---:|---:|
| 1 | `s160m_max` | 1,00 | 0 |
| 2 | `s70m_max` | 1,00 | 0 |
| 3 | `s160m_kurt` | 0,90 | 0 |
| 4 | `s70m_kurt` | 0,75 | 0 |

La fréquence de sélection est mesurée sur 20 sous-échantillons ; 1,00 signifie
que la covariable entre dans les 20 premières de chacun.

### Lecture

**Le maximum de surprisal d'un petit modèle gouverne l'indice de queue de la
surprisal catastrophique d'un grand modèle sur le même texte.** Profil de γ
local le long du rang de `s160m_max` :

```
u        0,10  0,20  0,30  0,40  0,50  0,60  0,70  0,80  0,90
gamma    0,83  0,84  0,84  0,85  0,92  0,95  0,95  1,02  1,28     moyenne 0,942
```

à comparer au γ global de 1,271 et aux témoins, qui restent entre 1,25 et 1,27
avec une étendue de 0,08 à 0,15.

Le kurtosis et l'asymétrie suivent immédiatement : c'est la **forme** de la
queue du profil de surprisal qui porte l'information, pas son niveau.

---

## 5. Comparaison des règles

Les quatre **porteuses** sont les quatre plus fortes variations de γ selon un
diagnostic indépendant (strates par covariable, Hill par strate, calibration par
strates aléatoires), qui n'utilise aucun des écrans comparés.

Les deux **témoins** sont les surprisals moyennes. Elles déplacent le quantile
0,95 de `log Y` de 1,4 nat sans bouger γ — étendue du profil de Hill 0,163
contre 0,449 pour le maximum. C'est une nuisance de localisation pure, vérifiée
et non postulée.

| règle | porteuses dans le top-4 | rang `s160m_mean` | rang `s70m_mean` |
|---|:---:|---:|---:|
| **queue, sélectionnée** | **4/4** | **177** | **241** |
| **queue, agrégée** | **4/4** | **170** | **189** |
| Yoshida–Umezu, réglage de son article | 4/4 | 100 | 125 |
| quantile 0,99 | 4/4 | 43 | 45 |
| quantile 0,95 | 3/4 | 24 | 27 |
| quantile 0,90 | 2/4 | 23 | 26 |
| Yoshida–Umezu, h = 2 | 2/4 | 15 | 21 |

`p = 264`, le hasard donnerait 132.

**Seules les deux règles de queue prennent les quatre porteuses et relèguent
simultanément les témoins au-delà du rang 170.** Quantile SIS les remonte à
23–27 — ce qui est le comportement *correct* pour ce qu'il estime, puisque les
témoins déplacent réellement les quantiles élevés. Les deux écrans répondent à
deux questions différentes, et seul le second répond à celle de l'indice de
queue.

---

## 6. Limites

**L'accord de rang global est défavorable aux règles de queue.** Contre une
mesure propre de *variabilité* de γ — l'étendue du profil de Hill local, sans
contamination de localisation — les règles de queue obtiennent −0,357 et −0,370,
quand quantile 0,90 atteint +0,487.

L'explication tient à la définition du score, qui est la **moyenne** du γ local
et non sa variabilité : l'écran retient ce qui abaisse le γ moyen, signature du
biais de mélange. Les deux fonctionnelles coïncident au sommet du classement et
divergent dans son corps.

Il faut noter que la contamination joue aussi dans l'autre sens. Le critère par
strates corrèle **+0,504 avec la variation de γ local et +0,506 avec celle du
quantile 0,95** : il est lui-même à moitié un critère de localisation, ce qui
avantage mécaniquement les écrans par quantile.

**La découverte principale n'est pas exclusive.** `s160m_max`, `s160m_kurt` et
`s160m_sd` sont trouvées par les deux familles d'écrans. Ce qui distingue la
méthode ici est ce qu'elle **rejette**, non ce qu'elle trouve seule : c'est un
argument de spécificité, du même ordre que celui de l'amas d'échelle à vérité
exacte, et non de puissance supérieure.

**Le résultat ne dépend pas du franchissement du seuil.** Il est significatif
aussi à `a* = 0,35`, où le rapport n'est que de 0,32 : le Δ réel de ces données
est bien plus grand que la valeur de référence 0,186 utilisée dans le tableau du
régime. Le franchissement est un point de présentation, pas ce qui fait marcher
l'application.

---

## 7. Reproduction

```sh
python code/py/sample_fineweb.py results/wild/main_docs.parquet 120000
python code/py/extract_ppl_same_span.py results/wild/main_docs.parquet \
       results/wild/ppl100k 100000 1          # GPU, ~85 min
python code/py/export_ppl_ranks.py results/wild/ppl100k/ppl_cov.parquet \
       results/wild/ppl100k_ranks
python code/py/gamma_varies_ppl.py results/wild/ppl100k/ppl_cov.parquet 500 4 0.02
WILD_A=0.20 Rscript code/R/screen_wild.R results/wild/ppl100k_ranks \
       results/wild/screen100k_fix 40 1000
Rscript code/R/real_wild.R results/wild/ppl100k_ranks \
       results/wild/ppl100k/ppl_cov_gammavar.csv results/wild/real_wild.rds 40
```

| | |
|---|---|
| corpus | FineWeb `sample-100BT`, fichier `000_00000.parquet` |
| graine de tirage | 20260817 |
| graine du départage d'ex æquo | 811000033 (criblage), 973000019 (comparaison) |
| graine des permutations | `505000 + b`, 1000 permutations |
| graine du sous-échantillonnage | 707, 20 plis |

**Avertissement d'implémentation.** `numpy.tofile` écrit toujours en ordre C ;
`matrix()` de R lit en ordre Fortran. Une matrice de rangs transmise sans
transposition arrive brouillée et le criblage lit du bruit, silencieusement.
`screen_wild.R` et `real_wild.R` vérifient désormais que chaque colonne est une
permutation de `1..n` avant de scorer.
