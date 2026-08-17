# Journal — renforcement mathématique (v2)

Objet : obtenir des preuves plus fortes / plus générales que celles de `manuscript/main-v2.tex`,
avec moins d'hypothèses. Chaque itération est menée par deux sessions GPT en parallèle ;
la vague suivante reçoit ce journal (obstacles + acquis) pour ne pas refaire le même chemin.

Fichier transmis à chaque session : `maths/paper-v2-math.tex`
(main-v2 + introduction, model_population, estimation_theory, appendices, discussion).

## État initial du papier (résumé des points d'appui)

- (C1)–(C2) : variation régulière conditionnelle + régularité de tail uniforme.
- Prop. « Projection » : l'indice projeté est le sup essentiel conditionnel de γ ; (S) sert
  seulement à remplacer l'ess-sup par le max sur toute la fibre.
- Prop. « Detectability » : Δ_j > 0 ⟺ I_ε \ π_j(M) de mesure positive.
- (E1)–(E3↑) : conditions sur les quantiles projetés ; dérivées de conditions primitives en annexe.
- Thm « Score concentration » : max_j |Ψ̂_j − Ψ_j| = O_P(√(log(pn)/(nαh))) + o(Δ_min).
- Thm « Sure screening » + Cor. « Exact recovery » sous log(pn) = o(nαh Δ_min²).

## Cibles connues (a priori, à confirmer/infirmer par les sessions)

- Borne inférieure / optimalité minimax du régime log(pn) ≍ nαh Δ_min² : manquante.
- Distribution asymptotique de Ψ̂_j (pas seulement concentration) ; intervalles / test de Δ_j = 0.
- Affaiblir (E2)–(E3↑) : peut-on remplacer les modules uniformes par une condition de
  second ordre standard (ρ, A(·)) et obtenir un compromis biais-variance explicite ?
- Dépendance : (S) et les fibres — que devient la cible quand le support conditionnel varie ;
  cas des copules à dépendance de queue.
- Indépendance des observations : extension mélangeante / série temporelle.
- Cas non détectable (π_j(M) = I_ε) : impossibilité formelle, ou score alternatif qui le voit.

## Critère d'arrêt de la campagne

On ne s'arrête pas sur « c'est mieux qu'avant ». On s'arrête quand le papier est dans le meilleur
des mondes possibles : la théorie est complète, et ce qui n'y est pas est *démontré* hors de portée.
Rien n'est acquis tant que la preuve n'est pas écrite en entier. Un obstacle démontré infranchissable
compte comme acquis : il devient un théorème d'impossibilité et clôt sa ligne.

**I. La théorie de population est exacte, pas suffisante.**
1. Caractérisation nécessaire *et* suffisante de ce qu'un screen marginal peut voir — la
   détectabilité `Δ_j > 0` comme équivalence, sous des hypothèses minimales (sans (S), avec le
   support conditionnel réel).
2. Le cas non détectable (`π_j(M) = I_ε`, l'exemple diagonal) est réglé : soit un théorème
   d'impossibilité pour toute statistique marginale, soit un score d'ordre supérieur qui le voit,
   avec sa propre théorie.
3. La cible est stable sous dépendance : ce que devient l'enveloppe quand les fibres se déforment,
   énoncé quantitativement (continuité en la copule, pas seulement « ça change »).

**II. L'estimation est optimale, pas seulement suffisante.**
4. Borne inférieure minimax sur une classe contenant les modèles du papier, **appariée** à la borne
   supérieure : le régime `log(pn) ≍ nαh Δ_min²` est le vrai seuil, constantes comprises si possible,
   sinon à facteur log près explicité. Si le seuil actuel n'est pas optimal, on le remplace.
5. Un seuil de transition net (sharp threshold) : en deçà, aucune procédure ne réussit ; au-delà,
   celle du papier réussit.
6. Les hypothèses de nuisance tombent : (E2), (E3↑) remplacées par une condition de second ordre
   standard `(ρ, A(·))` avec compromis biais-variance explicite, et (S) éliminée.
7. **Adaptativité** : choix de `h`, `α`, `d` sans connaître `Δ_min`, `ρ`, ni `s`, avec garantie
   d'oracle. L'agrégation par rangs sur le bloc de neuf tunings reçoit enfin une théorie, au lieu
   d'être une recette.

**III. L'inférence existe.**
8. Loi limite pour `Ψ̂_j`, et pour le maximum sur `p` coordonnées (extrême de champ), donnant un test
   de `Δ_j = 0` valide uniformément sur la classe.
9. Contrôle d'erreur multiple sur le jeu retenu : FDR ou FWER, non asymptotique de préférence.
10. Correction de biais explicite (Hill à biais réduit local), avec le gain de taux correspondant.

**IV. Le domaine de validité est large.**
11. Dépendance temporelle : mélange (β-mixing) ou champ dépendant, mêmes conclusions, taux dégradé
    explicitement quantifié.
12. `γ` non continu / non borné inférieurement, ou support de `Y` avec indice nul : ce qui survit.
13. Robustesse aux contaminations : une fraction `ε_n` d'observations arbitraires, seuil de rupture.

**V. Tout est vérifié.**
14. Chaque théorème retenu a été confié à une session indépendante dont le seul but est de le casser
    (contre-exemple, hypothèse cachée, circularité) et a survécu.
15. Les énoncés sont compilables en LaTeX, numérotés, cohérents entre eux, et cohérents avec les
    simulations et l'application déjà dans le papier ; les constantes annoncées sont vérifiées
    numériquement au moins une fois.
16. **Saturation** : deux vagues consécutives n'apportent plus ni théorème nouveau, ni hypothèse
    affaiblie, ni contre-exemple.

Tant qu'une ligne de I–IV n'est ni démontrée ni démontrée impossible, la campagne continue.

## Phase II — compression (ne s'ouvre qu'après la phase I)

**Porte.** Cette phase ne démarre **pas avant** que les seize lignes du cahier soient acquises ou
démontrées impossibles, et vérifiées. Chercher l'élégance avant d'avoir la vérité produit des preuves
courtes et fausses. Ordre non négociable : d'abord juste et complet, ensuite court et beau.

**Objectif.** Réécrire *toutes* les preuves du papier — annexes comprises, et les acquis des vagues
intégrés — en **moins de dix pages au total**. Aujourd'hui la masse est de l'ordre de quarante :
`appendix_aux` (604 lignes), `appendix_primitive` (482), `appendix_proofs` (114), plus l'addendum de
A (15 p.), le patch de D et le document de F (12 p.). Il ne s'agit pas de couper des résultats mais
de trouver les bons énoncés : la formulation dont tout le reste découle en trois lignes.

**Ce qu'on cherche, concrètement.**
- Le lemme central qui absorbe plusieurs preuves actuelles. Candidats visibles : la représentation
  quantile randomisée + Rényi rend interchangeables plusieurs arguments qui sont écrits séparément ;
  la coloration et le comptage de faux dépassants sont deux formes du même contrôle de dépendance.
- Une seule inégalité de concentration paramétrée, dont les théorèmes de score, de liste et
  d'agrégation soient des instances, au lieu de trois preuves parallèles.
- Le bon niveau d'abstraction pour la géométrie de population : enveloppe = sup essentiel
  conditionnel, tout le reste (détectabilité, invisibilité, ordre d'interaction) en corollaires.
- Éliminer les constantes explicites inutiles : elles gonflent les preuves et ne sont pas certifiées
  (cf. R2). Une constante universelle non nommée vaut mieux qu'un `1/32` invérifiable.
- Supprimer les hypothèses devenues inutiles : chaque hypothèse tombée en phase I raccourcit sa preuve.

**Règle de compression.** Aucune réduction ne doit affaiblir un énoncé. Toute preuve raccourcie
repasse par une session adverse dont le seul travail est de vérifier que le raccourci n'a pas escamoté
un pas. Une preuve courte non vérifiée ne compte pas.

**Règle de conduite.** Ne pas s'arrêter avant que le cahier soit rempli. Pas d'arrêt pour cause de
temps écoulé, de vague décevante, de session qui rend une réponse partielle, ni de « c'est déjà
bien ». Une vague qui échoue produit au minimum une entrée « obstacle » dans ce journal, et la
vague suivante repart de là. Le seul arrêt légitime est le critère lui-même — ou une instruction
explicite de Luc.

## Itérations

### Vague 0 — 2026-08-16
Connexion navigateur testée (Brave, compte Pro). Aucune session lancée pour l'instant.

### Vague 1 — 2026-08-16 (GPT-5.6 Sol, effort Pro)
Deux sessions parallèles, mêmes pièces jointes (`paper-v2-math.tex`), consignes volontairement larges.

- **A — « Strengthening Statistical Theory »**
  <https://chatgpt.com/c/6a80faf7-91c0-83ed-80ae-d268090b4d1c>
  Consigne : pousser les maths plus loin — conclusions plus fortes, hypothèses plus faibles,
  aller où la vraie théorie se trouve, preuves complètes.
- **B — « Proofs in High Dimensional Statistics »**
  <https://chatgpt.com/c/6a80fb18-e408-83eb-b652-2430ad73773e>
  Consigne : lire en referee hostile — quelles hypothèses portent réellement le résultat,
  que reste-t-il démontrable sans elles, puis démontrer le plus fort atteignable.

Durées : 35 min (A), 34 min (B). Les deux ont rendu.

#### Convergences (A et B, indépendamment — signal fort)

1. **(S) est éliminée, remplacée par la géométrie du support réel.** Avec
   `C_j(u) = supp K_j(u,·)` compact et γ continue,
   `ξ_j(u) = ess sup_{v∼K_j(u,·)} γ(u,v) = max_{v∈C_j(u)} γ(u,v)`.
   Preuve (identique des deux côtés) : le max `m` est atteint sur le compact ; tout voisinage d'un
   maximiseur a masse conditionnelle > 0 ; continuité ⟹ ess sup ≥ m−η pour tout η.
   Puis `ξ_j(u) = γ* ⟺ C_j(u) ∩ M_{j,u} ≠ ∅`, et
   `Δ_j > 0 ⟺ λ{u ∈ I_ε : C_j(u) ∩ M_{j,u} = ∅} > 0`, **sans continuité de u ↦ ξ_j(u)**.
   Hypothèse nulle minimale (A : « (NA) », B : « (NS) ») : pour `j ∉ A`, `C_j(u) ∩ M_{j,u} ≠ ∅`
   p.p. (S) l'implique et est strictement plus forte.
   → **critère I.1 atteint** (sous réserve de vérification).

2. **L'atomlessness / la continuité de (E1) tombent complètement.** PIT conditionnel *randomisé*
   (transformée distributionnelle) : avec `Z_ij ∼ U(0,1)` auxiliaires,
   `W_ij = F(Y_i−|U_ij) + Z_ij{F(Y_i|U_ij) − F(Y_i−|U_ij)}`, `V_ij = 1 − W_ij` sont i.i.d.
   uniformes conditionnellement à la colonne, et `Y_i = Q_j(1/V_ij, U_ij)` p.s. Les auxiliaires
   n'entrent jamais dans l'estimateur : pur artefact de preuve. (P1) perd sa clause d'atomlessness.

3. **La cible doit être `D = {j : Δ_j > 0}`, pas `A`.** Les théorèmes se réécrivent sans supposer
   que toute coordonnée active est détectable ; `D = A` devient une condition d'identification
   séparée. Plus honnête et strictement plus fort que la version actuelle.

4. **`log(pn) = o(nh²)` n'est pas nécessaire.** A : l'union DKW ne porte que sur les `p` colonnes,
   d'où `log(2p) = o(nh²)`. B : la condition disparaît — on ancre la fenêtre sur la statistique
   d'ordre `U_(r)j` et on contrôle la largeur par concentration des espacements uniformes
   (`U_(b)−U_(a) ∼ Beta(d, N−d)`, Chernoff), d'où `ρ_n(h,x) = h + 2√(hx/n) + 2x/n = h{1+o(1)}`
   dès que `log(pn) = o(nαh)`. **B est plus fort ici.**

5. **L'équicontinuité de ξ_j tombe.** A : remplacée par variation bornée (`q°_n ≤ C_ε(V_n+Γ)/n`).
   B : supprimée purement — la moyenne des `ξ_j(U_(r)j)` est une intégrale empirique tronquée,
   contrôlée par Hoeffding (Lemme 6), sans aucun module. **B est plus fort ici.**

6. **Le doublement uniforme (P3) tombe**, remplacé par la même quantité des deux côtés : la moyenne
   de déficit *inclinée* `m̄(T) = sup_{j,u} sup_{t≥T} M_{j,u}(t)/B_{j,u}(t) → 0`
   (A : « (P3-L) », B : « (TD) »). (P3) implique `m̄(T) = O(1/T)`, donc rien n'est perdu côté
   simulations, mais le doublement n'est plus présenté comme nécessaire.

#### Apports propres à A

7. **Nouveau score par blocs de rangs disjoints.** B ≍ 1/h blocs disjoints de m ≍ nh rangs, un Hill
   par bloc sur k ≍ αm extrêmes, score = moyenne des B statistiques. Conditionnellement à la colonne,
   les termes de Rényi sont **indépendants** ⟹ MGF exacte ⟹
   `max_j |Ψ̂_j − Ψ_j| = O_P(√(log 2p/(nα)) + log 2p/(nα))`, soit **le facteur h gagné**, avec des
   modules seulement *en moyenne* sur les blocs (plus faible que le sup). Condition de screening :
   `log(2p) = o(nα Δ_D²)` au lieu de `log(pn) = o(nαh Δ_min²)`. Subsiste une condition locale
   `log(p/h) = o(nαh)` (chaque bloc a besoin d'assez d'extrêmes). A signale honnêtement que les
   simulations actuelles utilisent les fenêtres glissantes : le taux `nα` ne peut pas leur être
   attribué sans refaire les campagnes.

8. **Recouvrement exact sans connaître `s`.** `γ̂max = max_j Ψ̂_j`, `Δ̂_j = γ̂max − Ψ̂_j` ; tout seuil
   `2e < λ < Δ_D − 2e` donne `{j : Δ̂_j > λ} = D`. Supprime l'exigence `d = s`.

9. **Impossibilité I — la barrière `1/log n` est informationnelle sous (C1)–(C2) nues.** Construction
   à deux points : `T_n = n^{2γ0}`, `δ_n = a/(2γ0 log n)`, lois identiques sous `T_n`, Pareto
   d'indices différents au-delà. Alors `Δ_1,n = a/(4γ0 log n) > 0` mais
   `‖P_0^{⊗n} − P_1^{⊗n}‖_TV ≤ 1/n`, donc **aucune procédure** (même multivariée, même non
   marginale) ne détecte uniformément les écarts d'ordre `1/log n`. Un seuil d'apparition
   quantifié (type (P2)) est donc *nécessaire*, pas commode.

10. **Impossibilité II — l'invisibilité marginale vaut pour la loi marginale entière.** Avec
    `g(t) = γ0 + a cos 2πt`, comparer `γ0(u) = g(u_2)` et `γ1(u) = g((u_1+u_2) mod 1)` :
    la coordonnée 1 est inactive dans l'un, active dans l'autre, et pourtant la loi de `(U_1, Y)`
    est *identique*. Toute règle marginale a erreur ≥ 1/2 à tout n. Ce n'est donc pas un défaut du
    score moyenné, ni même de l'indice de queue.

11. **Ordre d'interaction exact.** `r*(j) = min{|J| : j ∈ J, Δ_J > 0} ≤ s`, borne **atteinte** :
    pour `γ(u) = g((u_1+⋯+u_s) mod 1)`, `r*(j) = s` pour toute coordonnée active. Donc aucun screen
    de groupe d'ordre fixe ne répare l'invisibilité marginale.
    → **critère I.2 : le cas non détectable est tranché par l'impossibilité** (9–11).

#### Apports propres à B

12. **Impossibilité III — sous dépendance non restreinte, `A` n'est pas identifiable du tout.**
    `U_1 = U_2 = Z` ; `γ^(1) = 1+u_1` et `γ^(2) = 1+(u_1+u_2)/2` coïncident sur la diagonale, donc
    les lois jointes observées sont identiques alors que `A^(1) = {1} ≠ {1,2} = A^(2)`.
    La phrase « unrestricted cross-coordinate dependence » du papier est **incompatible** avec la
    récupération de `A`. Il faut soit (S)/(NS), soit redéfinir la cible (le classement des `Ψ_j`),
    soit définir l'activité modulo le support de `U`.

13. **(C2) est essentielle, mais s'affaiblit en condition de reste *moyennée* (AR).** Preuve de
    projection refaite sous (AR) via la mesure aléatoire inclinée `K_y*` ⟹ `δ_ξ`.
    Contre-exemple montrant que la variation régulière *ponctuelle* ne suffit pas : `Y = B(u) = u^{-1/β}`
    avec prob. 1/2, sinon Pareto(1) ; chaque loi conditionnelle a indice 1 et `L ≡ 1` à terme, mais
    le mélange a indice `1/β > 1`. La cause est l'absence d'un seuil d'apparition commun.

14. **Le biais spatial passe de `h log(pn)` à `h log(1/α)`.** Condition horizontale affine en log
    (HQ) : `|log{Q_j(s,v)/Q_j(s,u)}| ≤ A(ρ) + B(ρ) log s` pour `|u−v| ≤ ρ`, avec `B(ρ) < ξ̲ − b`.
    Plus faible que (E3↑) : aucune borne uniforme en `s → ∞` n'est exigée, la croissance naturelle
    `|ξ_j(v)−ξ_j(u)| log s` est permise. Le gel spatial (Lemme 8) donne
    `|H − Ȟ| ≤ 2A(ρ) + B(ρ){L_k + 2 log(3/α)}`.

15. **La troncature globale `min_{i,j} V_ij ≥ (pn)^{-2}` disparaît**, remplacée par un contrôle
    bilatéral `P{V_(k+1) ∉ [α/3, 3α]} ≤ 2e^{-k/4}` (Chernoff multiplicatif dans les deux sens).

16. **Ce que (P3) cachait : la codimension effective.** Si `H_{j,u}(z) ≍ a z^κ`, alors
    `H(2z)/H(z) → 2^κ`, donc un doublement uniforme **borne κ uniformément**. Or après projection
    sur une coordonnée active, κ ≈ le nombre d'autres coordonnées actives devant approcher
    simultanément leur maximiseur. Transformée de Laplace : `−B'(t)/B(t) ∼ κ/t`, donc le biais de
    variation lente est d'ordre **`κ/log(1/α)`** et la condition honnête est
    `κ_n/log(1/α_n) = o(g_n)`, pas `1/log(1/α_n)`. **C'est ici que la taille de l'ensemble actif
    entre dans la théorie** — le papier la faisait disparaître par une constante de doublement uniforme.

17. **Théorème fini non asymptotique complet** (Théorème 10 de B) avec les quatre termes de
    probabilité `2pe^{−nε²/32} + (6p+2pL)e^{−x} + 2pL e^{−Kt²/4} + 2pL e^{−K/4}` et
    `E_n(x,t) = 2A(ρ) + B(ρ){1+t+2log(3/α)} + Γt + (1+t)b + C_ε Γ{√(x/n) + 1/n}`.

#### Critiques d'audit (B)

18. Les définitions des modèles de simulation (A1–A3, B1) **n'étaient pas dans le fichier transmis**
    — erreur de ma part, `computation_simulation.tex` et `real_data.tex` avaient été omis.
    Corrigé : `maths/paper-v2-full.tex` les inclut désormais. La vérification de
    Corollary~\ref{cor:app-models} n'a donc pas pu être auditée.
19. Le lemme du coin gaussien est une **esquisse, pas une preuve** : enveloppe intégrable uniforme,
    contrôle uniforme de l'inverse gaussien sur la région mobile, erreur intégrée `O(1/r_z)`,
    stabilité sous la frontière non linéaire, extension pondérée, uniformité sur les familles
    conditionnelles — tout cela est affirmé, pas établi.
20. **L'abstract surinterprète** : il annonce que la dimension polynomiale est abordable dès que le
    compte local d'extrêmes domine `log p`, alors que le théorème imprimé exige en plus le
    déplacement des rangs et les biais de quantiles projetés, et que le transfert primitif ajoute
    `h log(pn)`.

#### Tension à trancher (le point le plus intéressant de la vague)

A obtient `√(log 2p / (nα))`, B conserve `√(log(pn) / (nαh))`. Le facteur `h` sépare les deux.
A change l'estimateur (blocs disjoints) et gagne l'indépendance entre blocs ; B garde les fenêtres
glissantes. Question ouverte : le facteur `h` est-il un artefact du chevauchement, ou une vraie
limite ? **Une borne inférieure minimax appariée trancherait**, et c'est exactement le critère II.4.

#### Bilan contre le cahier

- Acquis (à vérifier) : **I.1**, **I.2** (par impossibilité), **II.6** (partiellement : (S), (E1)
  continuité, (E2)→b, (E3↑)→(HQ), (P3)→(TD) — tous affaiblis).
- Progrès partiel : **I.3** (non-identifiabilité sous dépendance libre : la cible dépend du support ;
  reste la continuité quantitative en la copule).
- Non traité : **II.4/II.5** (borne inférieure minimax appariée, seuil net), **II.7** (adaptativité,
  théorie de l'agrégation par rangs), **III** (loi limite, test, FDR, biais réduit),
  **IV** (mélange, contamination, γ irrégulier), **V** (vérification indépendante).

#### Obstacles / dettes

- Aucune borne inférieure ne concerne encore le régime principal `log p` vs `nαh` : celle de A porte
  sur la barrière `1/log n` sous (C1)–(C2) non quantifiées.
- Le score par blocs de A n'a aucune validation numérique.
- La codimension `κ_n` (point 16) doit être calculée pour les modèles du papier ; si `κ_n` croît avec
  `s`, plusieurs énoncés du papier sont à corriger, pas seulement à affaiblir.

### Vague 2 — 2026-08-16 (GPT-5.6 Sol, effort Pro)
Pièce jointe : `maths/paste-wave2.txt` = `paper-v2-full.tex` (manuscrit complet, sections de
simulation incluses cette fois) + le présent journal jusqu'à la fin de la vague 1.

- **C — « Manuscript Review Request »**
  <https://chatgpt.com/c/6a810462-141c-83eb-9f78-0bb47b2be9eb>
  Consigne : casser. Tout énoncé du journal est suspect tant qu'il n'a pas été re-démontré ;
  hypothèses cachées, circularité, constantes qui ne survivent pas, contre-exemples. Ce qui survit,
  dire sous quoi exactement. → sert le critère **V.14**.
- **D — « Statistics Manuscript Review »**
  <https://chatgpt.com/c/6a810503-eaec-83eb-a63a-fe2020f45ec7>
  Consigne : ne pas refaire ce que le journal contient, aller au manquant ; le litige du facteur `h`
  est signalé, une borne inférieure appariée le trancherait. Terrain libre au-delà.
  → vise **II.4/II.5**, et ce que la session choisira dans **III/IV**.

#### Résultat D — le facteur `h` tombe pour l'estimateur **du papier**, avec borne inférieure appariée

C'est le résultat le plus fort de la campagne à ce jour. D confirme le diagnostic de C — le `h` vient
de ce qu'on borne le supremum du profil avant de moyenner — mais va plus loin : **il n'est pas
nécessaire de changer d'estimateur.**

**D1. Concentration directe du score à fenêtres glissantes.** On garde la décomposition *signée*
avant moyennage, au lieu de prendre le sup. On colore les centres de grille par leur résidu modulo
`2d_n+1` : deux fenêtres de même couleur sont **disjointes**, donc leurs variables de Rényi sont
indépendantes conditionnellement à la colonne. Hölder sur les classes de couleur donne
`log E(e^{θS_a} | U_j) ≤ A²v_nθ²/(2{1−A c_n|θ|})` avec
`v_n = (χ_n/L_n²)Σ_r 1/k_r`, `c_n = χ_n/(L_n K_n)`, `χ_n = 2d_n+1`.
Avec la géométrie de rangs déjà démontrée dans le papier (`K_n ≥ nαh/4`, `χ_n ≤ Cnh`, `L_n ≥ c_ε n`) :
`v_n, c_n ≤ C_ε/(nα)`. D'où, à `x ≍ log p` :
`max_j |Ψ̂_j − Ψ_j| = O_P(√(log p/(nα)) + log p/(nα)) + biais`.
Le proxy de variance est transparent : `χ_n/(L_n K_n) ≍ nh/(n·nαh) = 1/(nα)`.

**D2. Borne inférieure appariée pour le même estimateur** (le cœur). Modèle Pareto exact,
colonne de rangs indépendante des réponses, donc `Z_i = log Y_i` i.i.d. exponentielles. Famille à un
paramètre `f_β` qui ne modifie que l'échelle des excès au-dessus de `q_λ = γ log(m/(λk))`, la masse
au-dessus de `q_λ` restant fixée. Identité de Fisher
`Cov(H_m, Σψ(Z_i)) = γ² a_{m,k,λ}` avec `a = E(N∧k)/k ≥ P(N ≥ k)`, Chernoff
`P(N < k) ≤ exp{−(λ−1)²k/(2λ)}`, puis projection de Hoeffding et
`Σ_i N_i² ≥ L_n²m²/n` par Cauchy. Résultat :
`Var(Ψ̂) ≥ γ²/(λnᾱ)·[1 − exp{−(λ−1)²k/(2λ)}]²`, et avec `λ_k = 1 + k^{−1/4}` :
`liminf_n (nα/γ²)·Var(Ψ̂) ≥ 1`.
Borne supérieure finie appariée sur le même modèle (aucun biais spatial ni de variation lente) :
`P{|Ψ̂ − γ| > γ[√(2x/(L_nᾱ)) + x/(L_nᾱ)]} ≤ 2e^{−x}`, union sur `p` sans hypothèse d'indépendance
entre coordonnées. Donc `Ψ̂ − γ ≍_P (nα)^{−1/2}`. **Le litige est clos** : le `h` est réel pour une
statistique locale ou pour le supremum du profil, faux pour le score intégré.

**D3. Le seuil en dimension est `nαΔ² ≍ log(p/d)`, pas `nαhΔ²`** (Fano avec liste). Sous-modèle
plongé dans (C1)–(C2) **avec support de fibre plein** : `B ∼ Bern(q)`, `T = q^{−γ_0}`, et sous `P_j`,
`Y = T exp{(γ_0 − δU_j)E}`. Alors `ξ_j(u) = γ_0 − δu`, `ξ_ℓ ≡ γ_0` pour `ℓ ≠ j`, et `Δ_j = δ/2`.
`D(P_j‖P_0) ≤ qδ²/(3γ_0²)`, information-radius, puis Fano avec liste :
pour tout `Ŝ` de taille `≤ d < p`,
`(1/p)Σ_j P_j{j ∈ Ŝ} ≤ [nqδ²/(3γ_0²) + log 2]/log(p/d)`.
Avec `q = 3α` : région d'impossibilité `nαΔ_j² ≲ log(p/d)`.
**Borne supérieure appariée sur le même sous-modèle** : la statistique linéaire
`T_ℓ = n^{-1}Σ_i B_i(U_{iℓ} − 1/2)log(Y_i/T)`, `Ĵ = argmin_ℓ T_ℓ`, vérifie par Bernstein
`P_j(Ĵ ≠ j) ≤ 2p exp{−nqδ²/(12000γ_0²)}`. Le seuil est donc **atteint à constante absolue près**.
→ **critères II.4 et II.5 atteints** pour la frontière stochastique (restent hors champ : le biais de
variation lente projetée et la barrière d'apparition `1/log(1/α)`, qui ne produisent pas de `h`).

**D4. Vraie garantie d'oracle pour l'agrégation à neuf réglages** — réponse directe à l'objection de
C. Avec `r_j^agg = min_{λ∈Λ} r_j(λ)` et `|Λ| = M` :
`#{j : r_j^agg ≤ R} ≤ MR`. Donc si **pour chaque** `j ∈ D` il existe **un** réglage `λ_j` avec
`r_j(λ_j) ≤ R`, alors les `MR` premières positions de l'agrégat contiennent `D` — et des coordonnées
actives différentes peuvent être sauvées par des réglages différents. Pour `M = 9` : top `9s`.
Version probabiliste : `P{D ⊆ Â^agg_{MR}} ≥ 1 − Σ_{j∈D} P{min_λ r_j(λ) > R}`.
Ce n'est ni un théorème d'intersection ni un oracle complet, mais c'est une **vraie** théorie
partielle de l'agrégation implémentée. → **II.7 partiellement atteint.**

**Conséquences rédactionnelles.** Séparer trois conditions logiquement distinctes :
`log p = o(nαΔ_D²)` (bruit du score) ; `log(pn) = o(nαh)` (chaque fenêtre locale atteint la zone de
travail de Pareto) ; localisation des rangs et biais déterministes `= o(Δ_D)`. La proposition de
profil ponctuel **garde** le taux `nαh` — seul son transfert au score moyenné était dispendieux.
La phrase de l'abstract qui fait de `nαh` la taille d'échantillon effective **du score** doit changer.

#### Résultat C — audit hostile (archivé : `maths/wave2-C-answer.md`)

**Deux affirmations de la vague 1 sont fausses.**

**F1. La suppression de (S) est fausse telle qu'énoncée** (claim 1). Remplacer la fibre entière par
`C_j(u) = supp K_j(u,·)` est correct pour `ξ_j`, mais **continuer à comparer à
`γ* = max_{[0,1]^p} γ` ne l'est pas** : cette quantité dépend de γ *hors du support observé* et
n'est donc pas une fonctionnelle de la loi des données. Contre-exemple : `U_1 = U_2 = Z`,
`γ_1(u) = 1 + (u_1−u_2)²` contre `γ_0 ≡ 1` — lois observées **identiques**, pourtant
`max γ_1 = 2` et `max γ_0 = 1`, d'où `Δ_j = 1` contre `Δ_j = 0`.
*Correction* : prendre `γ*_U = ess sup γ(U) = max_{x ∈ supp U} γ(x)` et
`Δ^U_j = γ*_U − (1/|I_ε|)∫_{I_ε} ξ_j`. Quantité invariante par toute modification de γ hors du
support. L'équivalence `Δ^U_j > 0 ⟺ λ{u : C_j(u) ∩ {v : ι_j(u,v) ∈ M_U} = ∅} > 0` tient alors, et
toujours **sans continuité** de `u ↦ ξ_j(u)`. La cible est `D_U = {j : Δ^U_j > 0}`.

**F2. Le recouvrement sans connaître `s` est faux** (claim 8). `γ̂max = max_j Ψ̂_j` n'estime `γ*_U`
que **s'il existe au moins une coordonnée de gap nul**. Contre-exemple : `p = 1`, `γ(u) = 1 − u/2` —
la seule coordonnée est détectable, mais `Δ̂_1 = 0` identiquement, jamais sélectionnée.
*Correction* : l'énoncé vaut sous `D_U^c ≠ ∅`, et reste non adaptatif (il faut une borne `e`).

**Litige du facteur `h` : tranché.** Le facteur `h` **n'est pas informationnel**. Il provient du
contrôle du *supremum* du profil glissant, ensuite moyenné. Le score par blocs disjoints atteint bien
`O_P(√(log p/(nα)) + log p/(nα))` par Bernstein exponentiel sur `B_n k_n ≍ nα` spacings
indépendants, sous la faisabilité locale `log(p/h) = o(nαh)`. Corollaire : **le régime `nαhΔ²` du
papier n'est pas minimax-optimal tel que formulé**. Une borne inférieure en `nαh` exigerait un autre
problème (signal localisé de largeur `h`), et l'information y vaudrait `nαhδ² ≍ nαΔ²/h`.
Reste ouvert : le score à fenêtres glissantes atteint-il *lui aussi* le taux global `nα` après une
analyse directe de la dépendance ?

**Ce qui survit** (claims 2, 3 corrigé, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16), chacun
re-démontré indépendamment. Notamment `log(pn) = o(nh²)` **supprimée entièrement** (ancrage sur
`U_(r)j`, Chernoff sur `Beta(d, n+1−d)`, et `log(pn) = o(nαh)` implique `log(pn) = o(nh)` puisque
`α ≤ 1`) ; l'équicontinuité tombe sans passer par la variation bornée ; la troncature globale
`(pn)^{-2}` disparaît par le contrôle bilatéral de `V_(k+1)`.

**Nuances ajoutées.**
- *I.3 précisé* : l'enveloppe est stable en **distance de Hausdorff sur les supports**
  (`|max_C γ − max_{C'} γ| ≤ ω_γ(d_H(C,C'))`), et **pas** en variation totale, Wasserstein, ni sous
  convergence faible. Contre-exemple `K_η = (1−η)δ_0 + ηδ_1`, `γ(v) = 1+v` : `d_TV = W_1 = η` mais
  l'ess sup saute de 1 à 2. « Continuité en la copule » sans topologie de support est **faux**.
- *II.7 aggravé* : l'agrégation par rangs sur neuf réglages est un théorème d'**intersection**, pas
  d'oracle. Contre-exemple à quatre coordonnées : le premier réglage est parfait, les rangs minimaux
  valent `(1,2,1,2)`, les départages médian/moyen/max égalisent actif et nul, et la randomisation
  finale peut **exclure une coordonnée active**. Un bon réglage ne suffit pas — contrairement à la
  lecture « robustesse » que le papier en donne.
- *Codimension confirmée* : `κ_{j,n} ≍ s_n − 1` pour les projections actives, `≍ s_n` pour les
  inactives ; condition honnête `κ_n/log(1/α_n) = o(Δ_D)`. Les simulations à quatre variables actives
  restent valides ; toute lecture de l'annexe C comme couvrant `s_n` croissant ne l'est pas.
- *Constantes* : architecture du théorème fini validée, mais les constantes exactes du journal
  (`1/32`, `6p+2pL`, …) **ne sont pas certifiées** ; seule la version à constantes universelles survit.
- *Vérifications numériques* : les exposants `κ` du modèle AR(1) à `ρ = 1/4` et les gaps de population
  à `ε = .05` sont recalculés et **concordent** — `0.1863953328` (A1/B1), `0.1073057160` (A2),
  `0.1575534638` (A3). Le lemme du coin gaussien a la bonne forme limite et la bonne constante, mais
  sa preuve imprimée reste une esquisse (domination, frontière mobile).

**Apport positif (début du critère III).** Pour le score par blocs, sous sous-lissage (`α → 0`,
`nα → ∞`, biais `o((nα)^{-1/2})`), à coordonnée fixée :
`√(nα)(Ψ̂_j − Ψ_j) ⟹ N(0, σ_j²)`, `σ_j² = (1/|I_ε|²)∫_{I_ε} ξ_j(u)² du`,
par Lindeberg conditionnel sur les spacings exponentiels indépendants ; l'aléa des ancres contribue
`n^{-1/2}`, négligeable après multiplication par `√(nα)`. **Manquent** : la loi du maximum sur `p`,
un test uniformément valide de `Δ_j = 0`, un estimateur de `γ*_U`, le contrôle FWER/FDR.

### Vague 3 — 2026-08-16 (GPT-5.6 Sol, effort Pro)
Pièce jointe (par **upload de fichier**, plus par presse-papier) : `maths/paste-wave3.txt`
= manuscrit complet + ce journal + `wave1-A-addendum.tex` + `wave2-D-patch.tex`.

- **E — « Proof Review and Analysis »**
  <https://chatgpt.com/c/6a8111d4-9760-83eb-9563-f89ad171418e>
  Consigne : casser la vague 2, que personne n'a vérifiée — concentration par classes de couleur,
  borne inférieure de variance appariée, seuil de Fano avec liste, garantie d'agrégation, cible de
  population corrigée. → **V.14** pour la vague 2.
- **F — « Statistics Manuscript Analysis »**
  <https://chatgpt.com/c/6a81120c-aab0-83eb-b65b-87b9eb17408f>
  Consigne : aller à ce que le journal déclare non fait — inférence au-delà d'une coordonnée fixée,
  dépendance et contamination, réglage réellement pilotable par les données. → **III**, **IV**, **II.7**.

#### Résultat F — inférence simultanée, réglage honnête, dépendance, contamination
Archivé : `maths/wave3-F-inference.tex` (1 166 lignes, 14 énoncés, 14 preuves) — **compile**,
`maths/wave3-F-inference.pdf`, 12 pages. Travaille sur la cible corrigée `γ*_U` / `D_U`, et **isole
explicitement** ses deux seuls emprunts (la concentration à réglage fixe et la borne de Fano de la
vague 2) au lieu de les re-démontrer en silence.

**III — inférence (le critère est atteint pour l'essentiel).**
- Intervalles de confiance **simultanés et non asymptotiques** pour les `p` gaps corrigés `Δ_j`.
  Point clé : **aucune coordonnée de gap nul n'est requise** — la queue **inconditionnelle** de la
  réponse estime `γ*_U` directement, ce qui contourne la faille F2 de la vague 1 au lieu de la
  contourner par hypothèse.
- FWER non asymptotique, ensemble intérieur familialement valide, et ensemble extérieur de
  sure-screening **honnête**.
- Taille retenue `d` pilotée par les données.
- Dépendance arbitraire entre coordonnées **et** entre réglages candidats : sans effet sur la validité.

**II.7 — réglage (largement atteint).** Réglage **spécifique à chaque coordonnée** par maximisation
de bornes de confiance inférieures, avec garantie d'oracle ; et un théorème de **réglage pilote
arbitraire** par découpage de l'échantillon. Honnêteté explicite : ne vaut que pour des estimateurs
candidats portant un rayon valide à réglage fixe (certificats de biais de queue, spatial, quadrature).
L'adaptation uniforme à un seuil d'apparition non quantifié n'est **pas** revendiquée — le théorème
d'impossibilité à apparition retardée l'interdit sous (C1)–(C2) nues. Cohérent avec la vague 1.

**IV — dépendance temporelle.** Extension par **amincissement à sûreté de dépendance** sous
régularité absolue (β-mélange) : toute la construction se transfère avec l'erreur additive explicite
`(N−1)β(q)`, la taille i.i.d. étant remplacée par `N`. Taille d'échantillon effective explicite.

**IV — contamination (les deux côtés).**
- Le score de Hill glissant implémenté a un **point de rupture d'une seule ligne** (remplacement).
- Plancher informationnel : **aucune méthode** ne bat un plancher de gap en `ε_c/α` sous
  contamination de Huber ou par remplacement fort.
- Une statistique de queue bornée **atteint** cet ordre, conjointement à la frontière propre
  `√(log p/(nα))` : `Δ/γ_0 ≍ ε_c/α + √(log p/(nα))`, avec bornes inférieure et supérieure appariées
  **sur le sous-modèle de queue exact**.
  ⚠️ *Corrigé après audit H* : j'avais écrit `√(ε_c/α)`. Le terme est **linéaire**, sans racine.
  Raison : la variation totale locale entre exponentielles est du **premier ordre** en le changement
  d'échelle relatif — calcul exact `TV{Exp(μ(1−z)), Exp(μ)} = z(1−z)^{(1−z)/z}`, encadré par
  `z/e ≤ TV ≤ z`. Le recouvrement des voisinages de Huber a donc lieu à `qΔ/γ_0 ≍ ε_c`, pas à
  `qΔ²/γ_0² ≍ ε_c`. Erreur de transcription de ma part, sans conséquence sur le document de F.

**Dette laissée ouverte, et signalée par F elle-même** : étendre la borne supérieure robuste du
sous-modèle exact à toute la classe des quantiles projetés demande une construction de seuil local
robuste (seuillage + élagage), non traitée. C'est le prochain objectif naturel de IV.

#### Résultat E — audit de la vague 2 : confirmée, avec trois réparations

**Verdict global** : concentration par classes de couleur **survit** ; borne inférieure de variance
**survit mais la conclusion était fautive** ; Fano avec liste **survit, mais la borne supérieure
appariée imprimée est fausse** ; agrégation **survit exactement, l'interprétation non** ; cible
corrigée **survit après réparation « presque partout »**.

**R1. La cible corrigée exige un « p.p. ».** L'inclusion `C_j(u) ⊆ supp(P_U)` n'est garantie que pour
**presque tout** `u`, pas pour une version ponctuelle arbitraire du noyau conditionnel.
Lemme 1 (preuve par base dénombrable de rectangles ouverts de masse jointe nulle et désintégration)
fournit le `N_j` négligeable. Le théorème d'enveloppe et l'équivalence `Δ^U_j > 0 ⟺ …` se réénoncent
donc p.p. Conséquence conservée : un proxy inactif mais dépendant **peut** appartenir à `D_U`, une
coordonnée structurellement active **peut** avoir un gap marginal nul, et modifier γ hors du support
ne change jamais `D_U`.

**R2. Variance ≠ borne inférieure en probabilité.** L'inégalité finie de D est saine, mais
« aucune borne `o_P((nα)^{-1/2})` n'est possible » **ne découle pas** d'une minoration de variance :
une variable peut être `o_P(a_n)` avec une variance bien plus grande, à cause de valeurs rares.
*Réparation* : la borne supérieure appariée donne un contrôle uniforme du **moment d'ordre 4**
(`E X_n⁴ ≤ C γ⁴/A_n²` par intégration de la queue sous-gamma à deux régimes), puis Paley–Zygmund sur
`X_n²` donne `P(|Ψ̂ − γ| ≥ c_0 γ/√(nα)) ≥ c_1 > 0`. La conclusion tient donc, mais par ce détour.
*Nuance* : l'appariement est **en taux, pas en constante** — la minoration vaut ≈ `γ²/(nα)`, la
majoration ≈ `γ²/(L_n α)` avec `L_n/n → 1−2ε`. Ne pas écrire « optimal à constante près ».
La forme fine avec `Σ_i N_i²` est la meilleure minoration.

**R3. La borne supérieure appariée de Fano était fausse.** `Ĵ = argmin_ℓ T_ℓ` ne renvoie **qu'une**
coordonnée : elle n'apparie `log(p/d)` que si `d = 1`. Pour `d = p^{1/2}` le seuil vaut `(log p)/2`,
pour `d = p/log p` il vaut `log log p` — l'argmin n'atteint ni l'un ni l'autre.
*Réparation* (Théorème 8) : prendre `Ŝ_d` = les `d` plus petites valeurs de `T_ℓ`. Alors
`P_j{j ∉ Ŝ_d} ≤ 2(1 + (p−1)/d)·exp{−nqδ²/(2400γ_0²)}`, par décomposition
`{j ∉ Ŝ_d} ⊆ {T_j ≥ −t} ∪ {#{ℓ ≠ j : T_ℓ ≤ −t} ≥ d}` et Markov sur le compte de faux dépassants —
**sans hypothèse d'indépendance** entre les `T_ℓ`. La frontière `nαΔ² ≍ log(p/d)` est donc bien
atteinte, à constante absolue près. *Détail* : `q = 3α` impose `α < 1/6`, à énoncer.

**R4 (apport nouveau). Le score glissant du papier sait lui aussi exploiter la taille de liste.**
Lemme 9 : si chaque coordonnée vérifie `P{|Ψ̂_j − Ψ_j| ≥ Δ_U/2} ≤ π_n`, alors pour tout `d ≥ s`,
`P{D_U ⊄ Â_d} ≤ [s + (p−s)/(d−s+1)]·π_n`
(compte des faux dépassants + Markov, au lieu d'une séparation par l'erreur maximale). D'où
l'exigence stochastique `nαΔ_U² ≫ log(s + (p−s)/(d−s+1))` : `log(p/d)` pour `s = 1`, et l'on retrouve
`log p` quand `d = s`. Le seuil de Fano avec liste n'est donc pas une curiosité de borne inférieure —
c'est le régime que l'estimateur implémenté atteint dès qu'on cesse de demander la séparation exacte.

**R5. L'agrégation : la proposition est exacte, l'interprétation ne l'est pas.**
`#{j : r_j^agg ≤ R} ≤ MR` est correct, sans indépendance entre réglages. Mais **un réglage parfait ne
garantit pas le top `s`** : contre-exemple à `M = 2`, `p = 4`, cible `{1,2}`, rangs `(1,2,3,4)` et
`(3,4,1,2)` — les rangs minimaux valent `(1,2,1,2)`, les multi-ensembles de rangs sont identiques par
paires donc tous les départages prescrits égalisent, et les deux premières positions sont `{1,3}` :
la coordonnée active 2 est **exclue**. La seule conclusion valide est le top `Ms`.
**Conséquence directe pour le papier** : avec `M = 9` et `s = 4`, le théorème certifie le top **36**,
pas le top 20 utilisé dans les simulations. Certifier le top 20 par cet argument de comptage
exigerait `min_λ r_j(λ) ≤ 2` pour **chaque** coordonnée active. À appeler « couverture par union des
tops avec inflation de liste », jamais « garantie d'oracle à budget constant ».

**R6. Précision sur les couleurs.** `χ_n = 2d_n + 1` est un nombre de couleurs **valide**, pas
nécessairement le nombre chromatique (il peut y avoir moins de centres que de couleurs) ; la preuve
n'exige aucune minimalité. Et le `h` ne disparaît **que** de l'amplitude stochastique du score
intégré : la faisabilité locale `K_n ≍ nαh`, le terme d'échec du seuil (`nαh ≫ log(pn)`), la
localisation des rangs, le biais spatial et la concentration **ponctuelle** du profil
(échelle `(nαh)^{-1/2}`) restent tous en `h`.

- **G — « Manuscript Review Request »** (lancée pendant que F tournait)
  <https://chatgpt.com/c/6a8118aa-2adc-83ed-853f-b44646aafc77>
  Consigne : rendre publiables les deux passages de l'annexe que C a qualifiés d'esquisses —
  l'asymptotique du coin gaussien (domination intégrable uniforme près des faces, frontière non
  linéaire mobile, reste intégré, extension pondérée, uniformité sur les familles gaussiennes
  conditionnelles) et les constantes uniformes du couplage horizontal pour le modèle à échelle
  latente — ou montrer qu'ils sont faux. → **V.15** et dette de rédaction n° 6.

#### Résultat G — les deux esquisses de l'annexe deviennent des preuves
Archivé : `maths/wave3-G-appendix.tex` (1 609 lignes, 5 énoncés, 5 preuves) — **compile**,
`maths/wave3-G-appendix.pdf`, 16 pages. Titre du fichier : « Publication-grade Gaussian-corner
asymptotic and B1 horizontal coupling ». Quatre parties :
1. **Théorème du coin gaussien uniforme** — la version dont C disait qu'elle était affirmée, pas
   établie (domination intégrable près des faces, frontière non linéaire mobile, reste intégré,
   extension pondérée).
2. **Vérification uniforme sur les familles gaussiennes conditionnelles** — l'uniformité que le
   manuscrit invoquait sans la démontrer.
3. **Représentation primitive exacte de la queue pour B1.**
4. **Couplage horizontal uniforme pour B1, avec constantes explicites** — précisément ce que C
   signalait comme « plausible mais non écrit » (dérivées et constantes uniformes en dimension pour
   les moments gaussiens a posteriori).

→ **Dette de rédaction n° 6 levée**, et le volet « énoncés compilables » de **V.15** progresse.
Reste, pour V.15 : la cohérence numérique des constantes de G (non recalculées ici) et l'audit
adverse de G elle-même, qui n'a pas encore eu lieu.

- **H — « Statistical Proof Review »**
  <https://chatgpt.com/c/6a811d17-4800-83eb-b6e4-55661bb48089>
  Consigne : casser le document de F, que personne n'a vérifié — construction de confiance simultanée
  et sa ligne de base inconditionnelle, énoncés familiaux, réglage d'oracle et découpage pilote,
  amincissement β-mélange, plancher de contamination et sa borne supérieure appariée.
  → **V.14** pour la vague 3. Pièce jointe `maths/paste-wave5.txt` (manuscrit + journal + A + D + F).

#### Résultat H — audit de F : l'algèbre tient, l'« honnêteté » non

**Survivent exactement** : l'identification de `γ*_U` par la queue **inconditionnelle** (preuve
complète refaite, par mesure inclinée `K_y* ⇒ δ_{γ*_U}`) ; la décomposition de Rényi du Hill global
et sa borne finie ; l'algèbre max–min des intervalles simultanés (aucune indépendance requise, ni
entre coordonnées, ni entre réglages, ni entre la ligne de base et les scores) ; le **FWER de
l'ensemble intérieur** ; la validité post-sélection à menu fini ; le découpage pilote ;
l'amincissement avec l'erreur exacte `(N−1)β(q)` (lemme d'approximation produit redémontré) ; le
point de rupture d'une ligne ; le plancher de contamination en `ε_c/α`.

**Défauts réels.**

**H1. La bande « honnête » ne l'est pas.** Le rayon annoncé contient `b_0(a)`, qui est une propriété
**inconnue** de la loi inconditionnelle : la variation lente ne la rend ni calculable ni même finie à
un seuil fini donné. Le côté score exige des certificats déterministes, le côté ligne de base
substituait silencieusement la vraie valeur — c'est un rayon **d'oracle**, pas un intervalle de
confiance implémentable. *Réparation* : exiger un certificat connu `b̄_0(a) ≥ b_0(a)`, et énoncer
l'intervalle auto-normalisé `[γ̂_0/(1+t) − b̄_0]_+ ≤ γ*_U ≤ γ̂_0/(1−t) + b̄_0`, ce qui a l'avantage de
supprimer aussi le besoin de connaître `Γ`.
*Et aucune réparation sans certificat n'est possible* : construction à deux points (`T_n = n^{2γ_0}`,
`d_n = c/log n`, lois identiques sous `T_n`, `‖·‖_TV ≤ 1/n`) ⟹ tout intervalle uniformément honnête
sur la classe à variation régulière non quantifiée a une largeur **≥ d_n**, donc pas `o(1/log n)`.
Le découpage d'échantillon n'y change rien.

**H2. L'ensemble extérieur au seuil zéro est généralement vide de contenu.** Même avec des
estimateurs **exacts**, tout rayon strictement positif donne `Δ̄_j > 0` pour une coordonnée nulle,
donc `D̂_out = {1,…,p}`. La formule « taille retenue pilotée par les données » est littéralement vraie
et substantiellement trompeuse. L'énoncé utile est l'emboîtement des ensembles de niveau à seuil
`z > 0`. Sans signal minimal, aucune procédure uniformément valide ne peut exclure des coordonnées
dont le gap peut être arbitrairement proche de zéro.

**H3. Le réglage n'adapte qu'entre candidats certifiés**, pas au biais de queue inconnu. À écrire
explicitement. Idem : `Π_n → 0` doit être énoncé avant de dire que `|D̂_in|` est consistant.

**H4. Décalage d'espace de paramètres dans l'appariement de contamination.** La borne inférieure
compare `P_0` (aucune coordonnée active) à `P_j`, alors que la procédure robuste **renvoie toujours
une coordonnée** : les deux bornes ne portent pas sur le même problème de décision. Deux réparations,
toutes deux fournies : (A) borne inférieure **par paires** `P_j` contre `P_k`, via
`TV(P_j,P_k) ≤ 2qΔ/γ_0`, qui donne le même ordre `ε_c/q` dans la classe à exactement une active ;
(B) ajouter une **sortie nulle** à la procédure, avec le seuil `A_0 = qe^{−2}γ_0Δ_0/6`.

**H5. Portée de la borne robuste, à énoncer.** Sous-modèle de queue exact, plan indépendant, **une
seule** coordonnée active, et surtout **seuil connu** (`q, γ_0, T` entrent dans la statistique) :
c'est un résultat d'oracle. Il ne fournit ni intervalles robustes, ni estimation robuste de `γ*_U`,
ni version robuste du score glissant.

**H6. Conditions manquantes** au corollaire de transition nette : `α < 1/6` (le sous-modèle impose
`q < 1/2`), `log p = o(nα)` (sinon le terme linéaire domine la racine), et un régime de signal
**non vide** — sans quoi le membre de droite peut dépasser le signal maximal admissible `Δ/γ_0 ≤ 1/4`
et l'équivalence n'énonce rien.

**H7. Le mélange doit être uniforme sur la ligne complète** `Z_i = (U_i, Y_i)` en dimension
croissante : une borne coordonnée par coordonnée ne donne pas une borne indépendante de la dimension.
Et découpage pilote + β-mélange **ne se composent pas automatiquement** — il faut un intervalle
temporel entre les deux moitiés, ou un argument d'approximation produit conjoint.

**Apports positifs de H.** Valeur **exacte** de la variation totale entre exponentielles
(`z(1−z)^{(1−z)/z}`, encadrée par `z/e` et `z`), qui remplace Pinsker et fixe la transition de Huber
à un facteur `e` près ; critère géométrique exact d'intersection des voisinages de Huber
(`TV(P,Q) ≤ ε/(1−ε)`) ; et un transfert de remplacement fort **paramétré par `η`** montrant que la
constante `2 − ε_c` imprimée n'est pas optimale — en faisant `η ↓ 0` avec `nε_cη² → ∞`, le seuil de
remplacement fort tend vers le seuil de Huber.

**Restent ouverts, dixit H** : ligne de base robuste inconditionnelle, seuils locaux robustes, et
intervalles de gap simultanés robustes sur **toute** la classe des quantiles projetés. Les deux
premiers sont précisément ce qui a été confié à la session I.

### Vague 4 — 2026-08-16 (lancée sans attendre H)
Pièce jointe : `maths/paste-wave6.txt` = manuscrit + journal + A + D + F + G (288 ko).

- **I — « Mathematical Proof Request »**
  <https://chatgpt.com/c/6a8121d8-5d84-83eb-919f-e4a265840037>
  Consigne : la seule lacune que F déclare elle-même — étendre la borne supérieure robuste du
  sous-modèle de queue exact à toute la classe des quantiles projetés, par une construction de
  seuillage et d'élagage local robuste ; ou prouver qu'elle ne peut exister. → **IV**.
- **J — « Mathematical Proof Review »**
  <https://chatgpt.com/c/6a8121da-fb20-83eb-88aa-9efb42303a35>
  Consigne : casser G — domination intégrable uniforme, frontière non linéaire mobile, reste intégré,
  extension pondérée, uniformité revendiquée sur les familles gaussiennes conditionnelles, et **chaque
  constante explicite**, à recalculer. → **V.14** pour G.

#### Résultat J — audit de G : les deux résultats principaux tiennent, une affirmation est fausse

**Survivent** : l'asymptotique du coin gaussien (exposant, facteur gaussien à variation lente, signe
du terme de moyenne, constante de Dirichlet, ordre de la frontière non linéaire, constante pondérée)
et le fait que le modèle à échelle latente satisfait le couplage horizontal **avec des constantes
indépendantes de `p`**. L'enveloppe intégrable uniforme est valide sous l'hypothèse substantielle
`min_i (Ω1)_i ≥ ν_0 > 0`.

**Quatre défauts.**
1. **Une estimation de face mince manque.** G borne la masse de `g_z` près des faces mais **pas celle
   de la densité limite de Dirichlet** `f_ν`. L'estimation manquante est vraie et immédiate
   (`ν_i ≥ ν_0 > ϑ` ⟹ `∏x_i^{ν_i−1} ≤ C ∏x_i^{ϑ−1}`), mais sans elle la conclusion « uniformément
   sur tout mesurable `B` » n'est pas établie. **Une ligne à insérer.**
2. **Le sandwich de la frontière non linéaire est énoncé globalement** alors que Taylor n'était
   établi que sur un ensemble redimensionné borné. Vrai, mais la preuve demande **deux** arguments de
   bornitude distincts : `‖x‖₁ ≤ c_*^{-1}` sur le domaine réel, `‖x‖₁ ≤ a_-^{-1}` sur le simplexe de
   comparaison, avec `C_0 = (M_D/2)·max{c_*^{-2}, a_-^{-2}}`.
3. **Défaut de version au bord.** La factorisation de queue de B1 n'est prouvée que sur `(0,1)^p`
   (elle passe par `Φ^{-1}`), alors que (P2) est imprimée ponctuellement sur le cube fermé. Le bord
   est de probabilité nulle : ce n'est pas un obstacle probabiliste mais un défaut de version
   conditionnelle régulière. À énoncer « pour `P_U`-presque tout `x` », ou à fixer une version
   arbitraire sur le bord — les constantes sont inchangées.
4. **Une affirmation est fausse** : G explique que le nombre de proxys fixé à `q = 20` intervient
   dans le couplage. **Non.** La dépendance en `q` s'annule exactement : `β{1+(q−1)λ²} = λ` pour tout
   `q ≥ 1`, puisque `β = λ/(σ²+qλ²)` et `σ² = 1−λ²`. Les constantes horizontales restent donc
   indépendantes de la dimension **même si le bloc de proxys grandit**. Phrase à supprimer — et le
   résultat est plus fort que ce que G revendiquait.

**Constantes : toutes valides, plusieurs très lâches.** J recalcule et améliore :
- l'enveloppe de dérivée du reste de B1 passe de `2e⁴e^{1/2−2t}` à **`(32/9)e^{−2t}`** ;
- pour le couplage, la clé est que dans B1 `h` et `m` ne bougent **jamais simultanément** (si `j` est
  un proxy, le bloc actif en est indépendant, donc `h` ne bouge pas ; sinon c'est `m`), donc les deux
  contributions se combinent par un **maximum**, pas par une somme. Les constantes affûtées sont
  environ **2,6 à 3,2 fois plus petites** que celles imprimées.
- Les exposants `κ` des familles conditionnelles et les bornes spectrales par Gershgorin sont
  recalculés et **concordent** avec G.

⚠️ *Prudence de transcription* : les valeurs numériques et l'orientation des fractions
(`3/5` contre `5/3`, etc.) proviennent d'un rendu de page qui **mutile les fractions** — j'ai déjà
commis une erreur de ce type sur le terme de contamination. Avant tout usage dans le manuscrit, les
relire dans le fil <https://chatgpt.com/c/6a8121da-fb20-83eb-88aa-9efb42303a35>, pas dans ce journal.
La capture disque de J a échoué (presse-papier sans focus) ; le doublon produit a été détecté et
supprimé.

**Bilan G après audit** : conserver les deux résultats principaux, insérer les deux lignes de preuve
manquantes, réparer la version au bord, supprimer l'affirmation sur `q = 20`, et remplacer de
préférence les constantes de couplage par les versions affûtées.

#### Résultat I — seuillage robuste : la borne appariée vaut sur **toute** la classe

I démontre d'abord ce qui **ne marche pas** : supprimer un nombre fixe des plus grandes observations
puis appliquer Hill garde une influence non bornée, l'adversaire n'ayant qu'à insérer une observation
extrême de plus qu'on n'en a retirées.

**La construction.** Garder les `k` plus grandes observations locales (le seuil reste **local**),
puis **jeter la valeur numérique du seuil** ; former les **écarts logarithmiques par paires** entre
les `k` retenues ; **winsoriser** chaque écart à un plafond `T` ; moyenner. Chaque statistique locale
vit dans `[0,T]`, et le terme non borné `log(Y/Y_seuil)` a disparu.

**L'identité centrale.** Si `E₁, E₂` sont exponentielles standard indépendantes, alors
`|E₁ − E₂| ∼ Exp(1)`, donc `E[φ_T(ξ|E₁−E₂|)] = m_T(ξ) = ξ(1 − e^{−T/ξ})`. La winsorisation
transforme donc `ξ` par une fonction **connue et strictement croissante** : le score transformé a
**exactement le même ensemble de gaps positifs**, avec `κ_T Δ_j ≤ Δ_{j,T} ≤ Δ_j` et, pour `T = Γ`,
`κ_Γ = 1 − 2/e ≈ 0,264`. On ne perd qu'une constante fixe de force de signal.

**L'arithmétique du point de rupture** — c'est là que la largeur de bande s'annule :
une ligne modifiée touche `O(nh)` fenêtres ; le score moyenne par `1/n` ; et à l'intérieur d'un
top-`k` local elle ne pèse que `O(1/(nαh))` (elle n'entre que dans `k−1` des `C(k,2)` paires).
Produit : `O(1/(nα))`. Donc `nε_c` remplacements de lignes **complètes** donnent `O(ε_c/α)` —
**déterministe**, l'adversaire ayant le droit d'inspecter tout l'échantillon propre avant de choisir.

**Les théorèmes.** Concentration propre à l'échelle `nα` (même Hölder par classes de couleur, avec
différences bornées `2T/k`) ; ligne de base inconditionnelle robuste par la même statistique par
paires appliquée globalement (`≤ 8Tε_c/a`) ; puis **intervalles de gap simultanés robustes sur toute
la classe des quantiles projetés**, avec FWER et recouvrement exact en corollaires. Une inversion
`G_T` de classe C¹ permet même de revenir à `Ψ_j` d'origine, au prix d'un biais de courbure en
`1/k` (Efron–Stein).

**Le seuil.** `Δ_min ≫ ε_c/α + √(log p/(nα)) + biais certifiés (queue, spatial, quadrature)`.
Il **apparie** la borne inférieure du sous-modèle plongé, et les rôles de `h` sont enfin séparés :
`nαh` gouverne la faisabilité locale et le biais spatial, `nα` la fluctuation stochastique intégrée,
`ε_c/α` le plancher de contamination. **Aucun `h` dans la frontière.**

I énonce aussi sa limite : rien de tout cela ne s'étend à (C1)–(C2) nues — l'impossibilité à
apparition retardée de la vague 1 l'interdit. Les trois vagues sont donc cohérentes sur ce point.

→ **Critère IV clos** (mélange, contamination, plancher, borne appariée), sous réserve de l'audit
adverse de I, qui reste à faire.

### Vague 5 — 2026-08-16
Pièce jointe : `maths/paste-wave7.txt` (298 ko) = manuscrit + journal + A + D + F + G.

- **K — « Theory Reconstruction Request »**
  <https://chatgpt.com/c/6a8128e3-ba54-83ed-847a-3a567de1c40e>
  Consigne : écrire la théorie **une fois, correctement** — un énoncé autonome, chaque hypothèse dans
  sa forme survivante la plus faible, chaque résultat dans sa forme corrigée, rien de ce qu'un audit a
  tué ; là où deux vagues démontrent la même chose différemment, garder la meilleure preuve et dire
  pourquoi. C'est le préalable de la **phase II** : on ne comprime pas six documents qui se
  contredisent par endroits, on comprime une théorie unique.
#### Résultat K — théorie unifiée (première version)
K traite les verdicts adverses comme **contraignants** et énonce ses choix de preuve :
- projection par la mesure inclinée conditionnelle, qui **admet les atomes** — donc les lois d'indice
  projeté dégénérées des coordonnées nulles — contrairement à un argument par densité ;
- pour le score glissant, **combinaison des deux réparations** : gel à la statistique d'ordre
  `U_(r)j` (supprime le déplacement DKW *et* toute régularité de `u ↦ ξ_j(u)`) **puis** moyennage de
  l'expansion de Rényi signée par coloration. Strictement préférable à « borner le sup puis
  moyenner », et **sans changer l'estimateur implémenté** ;
- le score par blocs disjoints n'est gardé **que** pour ce qu'il apporte en propre : un TCL simple à
  coordonnée fixée sous hypothèses de biais moyennées ;
- borne inférieure spécifique à l'estimateur **en probabilité**, par variance + moment d'ordre 4 +
  Paley–Zygmund — jamais la variance seule ;
- screening en liste par **comptage de faux dépassements** : garde le gain `log(p/d)` sans exiger la
  séparation uniforme de chaque coordonnée nulle ;
- coin gaussien : la preuve **réparée** (estimation de face mince incluse, deux arguments de
  bornitude pour la frontière) ; B1 énoncé presque partout sur le cube ouvert, constantes
  indépendantes de la taille du bloc de proxys.

**Limite de portée, signalée par K elle-même** : le paquet transmis ne contenait **pas** le résultat
de I (terminé après la construction du bundle). K énonce donc la borne supérieure robuste seulement
sur le sous-modèle d'oracle à queue exacte, « aussi loin que les preuves incluses vont ».
→ **Une seconde passe de reconstruction est nécessaire**, incluant I et L.

⚠️ Document **non capturé sur disque** : le presse-papier a échoué (fenêtre Brave sans focus) et a
reproduit G ; doublon détecté et écarté avant écriture. Le texte reste dans le fil.

- **L — « Mathematical Proof Request »**
  <https://chatgpt.com/c/6a8128e6-4cc4-83eb-8470-76f2abe17dd0>
  Consigne : les deux points que H déclare ouverts sous contamination — un estimateur robuste de la
  ligne de base `γ*_U`, et des intervalles de gap simultanés robustes valides sur **toute** la classe
  des quantiles projetés, pas seulement le sous-modèle plongé. Les construire, ou prouver qu'ils ne
  peuvent exister au taux propre. → fin de **IV**.

#### Résultat L — ligne de base robuste + intervalles robustes sur toute la classe, **et** bornes inférieures appariées

**Construction (annulaire, appariée, élaguée)** — différente de celle de I. Bacs **fixes disjoints**
sur l'échelle des `U` (pas de fenêtres de rangs glissantes) ; dans chaque bac, deux anneaux séparés
`A_k = {k+1,…,2k}` et `B_k = {3k+1,…,4k}` — les `k` plus grandes observations sont **entièrement
jetées** ; on forme les `k²` écarts croisés, on **élague le quart supérieur et le quart inférieur**,
on moyenne, on calibre par une constante `μ_k` exponentielle exacte (encadrée : `log(3/2) ≤ μ_k ≤ log 4`)
et on plafonne à `64Γ`. N'utilise **ni** seuil oracle, **ni** masse de queue oracle, **ni** `γ_0`.

**Résultats.** Ligne de base robuste `γ̂_U^R` avec
`|γ̂_U^R − γ*_U| ≲ b̄_0 + Γ{ε_c/a + √(log(1/δ)/(na)) + log(1/δ)/(na)}` ;
scores robustes simultanés avec
`max_j |Ψ̂_j^R − Ψ_j| ≲ b̄_Ψ + ω̄_Ψ + q̄_B + Γ{ε_c/α + √(log(p/δ)/(nα)) + …}` ;
puis **intervalles de gap simultanés** à couverture familiale, sans aucune indépendance requise —
ni entre coordonnées, ni entre la ligne de base et les scores. Plus : réglage post-sélection sur
menus finis, **par coordonnée**, la validité étant préservée ; et contamination de Huber ramenée au
remplacement fort via la concentration du nombre d'indicateurs (`ε_c^eff`).

**Bornes inférieures appariées, nouvelles.** Par intersection de voisinages de Huber
(`TV(P,Q) ≤ ε/(1−ε)`) et la valeur **exacte** de la variation totale entre exponentielles, plus Fano :
toute famille d'intervalles honnêtes a une largeur `≳ γ_0{ε_c/α + √(log p/(nα))}`, et la ligne de base
`≳ γ_0{ε_c/a + 1/√(na)}`. Le taux atteint est donc **minimax**, aux biais déterministes près.
L reprend aussi la construction à apparition retardée : même à `ε_c = 0`, **aucune** version
uniformément honnête n'existe sans certificat de biais de queue fini. Cohérent avec les vagues 1 et 3.

**⚠️ Tension à trancher entre I et L.** L affirme que les bacs **fixes** sont *nécessaires* sous
remplacement de lignes complètes : l'adversaire contrôlant aussi `U`, un **rang** adverse peut
modifier de **nombreuses** fenêtres glissantes, alors qu'une ligne ne touche que **deux** multi-ensembles
de bacs fixes. Or la construction de I est bâtie sur des **fenêtres de rangs glissantes** et son
arithmétique de rupture suppose qu'une ligne ne touche que `O(nh)` fenêtres — ce qui est vrai pour
une ligne *ordinaire*, mais peut-être pas pour une ligne *adverse* choisie après inspection de
l'échantillon. Les deux résultats concluent au même taux `ε_c/α`, par des chemins incompatibles :
**l'un des deux a une faille, ou les deux modèles d'adversaire diffèrent sans le dire.**
→ à arbitrer en priorité (session N).

#### Résultat N — arbitrage I contre L : chacun a raison sur un point, tort sur l'autre

**L a raison sur la prémisse.** L'affirmation « une ligne modifiée ne touche que `O(nh)` fenêtres »
est **fausse** sous adversaire adaptatif. Adversaire explicite : remplacer la ligne de rang 1 en lui
donnant une covariable **supérieure à toutes les autres**. Son rang devient `n`, et toute la suite
se décale d'un cran : `Z'_q = Z_{q+1}`. Pour **Θ(n)** centres non tronqués, la fenêtre change.
Ce qui reste `O(nh)`, c'est seulement le nombre de fenêtres **contenant** la ligne remplacée.
Donc la preuve de rupture de I, telle qu'énoncée, ne vaut que sous rangs gelés — remplacement de la
**réponse seule**, ou covariables de confiance.

**L a tort sur la conclusion.** Les bacs fixes ne sont **pas nécessaires**. Car ces `Θ(n)`
changements ne sont pas quelconques : ce sont des **translations exactes**, `T(Z'_{I_r}) = T(Z_{I_{r+1}})`.
Sur un bloc maximal de centres où l'identité vaut, la somme **télescope** :
`Σ_{r=u}^{v} {T(Z'_{I_r}) − T(Z_{I_r})} = T(Z_{I_{v+1}}) − T(Z_{I_u})`, borné par l'oscillation `B_0`
**quelle que soit la longueur du bloc**. Seuls `O(d)` centres exceptionnels, près des deux extrémités
du déplacement de rang, exigent la borne grossière `A_0/K`. D'où, pour **une** ligne :
`|S(Z') − S(Z)| ≤ [B_0 + (C_0 d + C_1)A_0/K]/L ≍ B/(nα)` — **aucun facteur de largeur de bande**.
Puis `m` remplacements adaptatifs par inégalité triangulaire **trajectorielle** : `≤ C_ε B ε_c/α`,
sans union bound sur `p`, l'adversaire ayant même le droit de choisir des rangs différents dans
chacune des `p` coordonnées. Seule convention requise : un départage déterministe des ex æquo.

**Le fait structurel, à retenir.** Le télescopage vaut pour le score **intégré** et **échoue** pour un
supremum, une fenêtre isolée, une grille éparse, des poids fortement non uniformes, ou une
fonctionnelle dépendant du centre. C'est **exactement le même phénomène** que l'annulation du `h` par
classes de couleur : le score moyenné est structurellement mieux conduit que le profil ponctuel.
Deux découvertes de la campagne, obtenues indépendamment, sont donc une seule.

**Deux réserves nouvelles.**
- *Sur I* : si le plafond de winsorisation `B_n → ∞`, la borne devient `B_n ε_c/α`. Une constante
  indépendante de `n` exige un plafond **fixe**, une borne globale connue sur l'indice, ou une
  normalisation à dénominateur uniformément minoré.
- *Sur L* : si ses « bacs fixes » sont en réalité des bacs de **quantiles empiriques recalculés sur
  les covariables contaminées**, ils ne sont pas fixes — leurs coupures bougent, et il faut soit des
  coupures de population de confiance, soit un pilote propre, soit un théorème séparé de frontière
  robuste. Et des bacs déterministes sur l'échelle **brute** des covariables réintroduisent les
  conditions de densité marginale que les rangs empiriques servaient précisément à éliminer —
  ce qui coûterait l'un des arguments de vente du papier.

**Verdict** : les deux résultats tiennent, aucune preuve n'est à jeter, mais **celle de I doit être
remplacée** par l'argument de télescopage, et le mot « nécessaires » doit disparaître de L.

- **N — « Statistical Robustness Conflict »**
  <https://chatgpt.com/c/6a813377-07f4-83eb-829c-f6b912eb79a8>
  Consigne : trancher la tension I/L. Sous remplacement de lignes complètes où l'adversaire choisit
  la **covariable** après avoir vu l'échantillon propre, un score robuste à fenêtres de rangs
  glissantes garde-t-il un point de rupture borné, oui ou non ? Preuve, ou adversaire explicite. Et
  dire précisément **quel modèle d'adversaire** chaque construction exige réellement.

#### Résultat M — audit numérique indépendant : **tout est juste**, sauf une phrase

Recalcul intégral et indépendant. **Aucun** gap de population, matrice de précision, exposant de coin,
somme de ligne ou borne spectrale n'est faux.

- **Gaps** : valeurs exactes `0,186395332752104` (A1/B1), `0,107305716040764` (A2),
  `0,157553463789990` (A3). Les trois décimales du manuscrit et les dix décimales du journal sont des
  arrondis corrects. Profils confirmés : `e^{-u}/2`, `e^{-u/2}/2`, `e^{-0,8u}/2`, et `1/2` pour les
  inactives. Détail joli : dans A3 les interactions **disparaissent** de l'enveloppe à une coordonnée,
  parce que leurs partenaires libres maximisent en zéro.
- **Exposants de coin** : `34/15`, `41/15`, `14/5 + λ_j²/(1−λ_j²)`, et `14/5` pour les proxys B1 ;
  encadrement global `34/15 ≤ κ_j ≤ 43/15` confirmé, avec la forme exacte `κ_j = 14/5 + 1/(16^{j−4}−1)`.
- **Orientation des fractions, tranchée** (c'était mon point d'incertitude) : les valeurs propres de
  toutes les matrices de précision conditionnelles vivent dans `[3/5, 5/3]`, et la borne inférieure
  uniforme de somme de ligne est **`ν_0 = 3/5`**, *atteinte* par les lignes médianes — donc exacte, pas
  conservatrice. Les bornes spectrales, elles, ont de la marge (plus petite valeur propre réelle
  ≈ 0,6819, plus grande ≈ 1,5568).
- **Enveloppes de queue A** : coefficient primitif exactement 1, reste et sa dérivée
  **super-exponentiellement** petits, avec enveloppes explicites.
- **B1** : `s_*² = 51/1031`, `β = 70/1031` confirmés ; factorisation de Lambert-W correcte.

**La seule assertion spécifique fausse** reste la phrase sur `q = 20` : la simplification
`β{1+(q−1)λ²} = λ` est exacte pour tout `q ≥ 1`. M note que G **affiche** l'annulation puis affirme
immédiatement que `q = 20` intervient — les deux phrases sont contradictoires. J avait raison.

**Améliorations de constantes, au-delà de ce que J avait trouvé.**
- Facteur exact sur `L_D` : **64/25 = 2,56** (M retrouve indépendamment le « ~2,6 » de J), soit
  `L_D = 1164,64` au lieu de `2981,47`.
- Facteur sur `L_H` : **3,213**, soit `1198,31` au lieu de `3850,31`.
- **Nouveau** : G maximisait `(y−1)e^{−y}` en `y = 3/2` alors que le maximum est en `y = 2`. D'où
  `C_h = 1/4 + 1/(2e²) ≈ 0,3177` au lieu de `0,5289`, et `L_H ≈ 952,29` — un facteur **4,04** sous G.
  Donc même la version affûtée de J n'était pas optimale.
- Le coefficient de couplage des coordonnées PIT passe de `6,826` à `3,345` (la plus grande
  corrélation entre coordonnées distinctes vaut 0,49, dans le bloc de proxys).

→ **V.15 rempli** : les quatre documents compilent, et les nombres sont désormais vérifiés
indépendamment. Aucune constante annoncée n'est fausse ; plusieurs étaient inutilement lâches.

- **M — « Recompute and Verify Models »**
  <https://chatgpt.com/c/6a813167-25cc-83eb-846d-0518218d1bf1>
  Consigne : **tout recalculer de zéro** pour les quatre modèles de simulation — gaps de population
  au niveau de trimming annoncé, matrices de précision conditionnelles et exposants de coin, bornes
  spectrales et de somme de ligne, enveloppes de coefficient de queue et de reste, constantes de
  couplage horizontal — et signaler **chaque** désaccord, si petit soit-il, en disant qui a raison.
  Montrer le calcul, pas seulement le verdict. → **V.15**, second volet (cohérence numérique), le
  seul élément du cahier que personne n'avait pris en charge de façon systématique.

### Vague 6 — 2026-08-16

- **K2 — « Statistics Manuscript Reconstruction »**
  <https://chatgpt.com/c/6a813beb-8aec-83eb-8c21-97a0bb2a986f>
  Seconde passe de reconstruction, avec le journal complet (six vagues). Consigne : écrire la théorie
  entière **une fois, correctement et de façon autonome**, et **re-dériver soi-même** les trois
  résultats dont le LaTeX n'est pas dans le paquet (score winsorisé par paires sur fenêtres
  glissantes ; score annulaire élagué sur bacs fixes avec bornes minimax appariées ; arbitrage
  imposant le remplacement de la preuve de rupture par le télescopage) plutôt que de les citer.
  Re-dériver au lieu de citer est **volontaire** : cela vaut vérification indépendante des trois
  résultats non archivés, et compense l'absence de leur source.

- **O — « Mathematical Proof Validation »**
  <https://chatgpt.com/c/6a813c5a-2188-83eb-a038-7c3e334db827>
  Consigne : **attaquer les trois impossibilités** sur lesquelles la campagne s'appuie le plus —
  barrière `1/log n` sous (C1)–(C2) non quantifiées, invisibilité marginale au niveau de la loi
  entière, non-identifiabilité de `A` sous dépendance libre. Les constructions sont-elles valides,
  les conditions de régularité annoncées sont-elles réellement satisfaites par les contre-exemples,
  l'ordre des quantificateurs est-il correct, et chaque impossibilité couvre-t-elle vraiment la
  classe pour laquelle on l'invoque ? *Motif* : plusieurs lignes du cahier ne sont « acquises » que
  parce que ces trois théorèmes tiennent. Personne ne les a attaqués depuis leur production à la
  vague 1 — c'est le plus gros angle mort de **V.14**.

#### Résultat O — les trois impossibilités tiennent ; **deux de mes glosses sont retirées**

**Aucun des trois contre-exemples ne s'effondre.** Les constructions sont valides, les conditions de
régularité sont réellement satisfaites (le recollement de la barrière `1/log n` est bien sans atome ;
`γ^(1) = γ_0 + a cos{2π(u_1+u_2)}` est **C^∞** sur le carré, donc aucune discontinuité au bord ;
tous les exemples ont `L ≡ c ≡ 1`, soit (C1)–(C2) sous leur forme la plus forte). Mais la campagne —
c'est-à-dire moi — les a utilisées **au-delà de leurs quantificateurs**.

**Correction 1 — la barrière `1/log n` est *uniforme*, pas ponctuelle.** L'énoncé exact est
`∃(P_{0,n},P_{1,n}) ∀(φ_n)`, avec des **tableaux triangulaires** ne dépendant pas de la procédure.
Il exclut `sup_P P^n(erreur) → 0`, **pas** `P^n(erreur) → 0` à `P` fixé : pour un modèle à seuil
**fixe**, un échantillon croissant finit par voir au-delà du seuil et la consistance ponctuelle n'est
pas exclue. Partout où j'ai écrit « ne peut pas être détecté », lire « n'est pas détectable
**uniformément** ». En outre le théorème ne montre **pas** que (P2) est *la* condition nécessaire —
seulement qu'**une** restriction excluant l'apparition arbitrairement retardée l'est ; ni qu'il existe
une barrière universelle `1/log(1/α_n)` : les deux échelles ne coïncident que si
`log(1/α_n) ≍ log n`, ce qui vaut pour `α_n = n^{-a}` et pas en général. L'usage pour la largeur des
intervalles honnêtes, lui, est **validé**.

**Correction 2 — l'invisibilité marginale ne vaut que pour une décision *locale à la coordonnée*.**
Ma phrase « toute règle marginale a erreur ≥ 1/2 » est **trop large** et je la retire. Le théorème
porte sur les rangs mesurables par rapport au seul couple `(U_j, Y)`. Une procédure autorisée à
regarder **une autre** coordonnée distingue les deux modèles, et O donne la statistique explicite :
`T_n = n^{-1} Σ cos(2πU_{i2}) log Y_i`, d'espérance `a/2` sous `P_0` et `0` sous `P_1`.
Un algorithme de screening qui calcule toutes les statistiques marginales **puis** fait dépendre
l'étiquette de la coordonnée 1 de ce qui s'est passé en coordonnée 2 n'est pas couvert.
Second retrait : `Δ_j = 0` **n'implique pas** que la loi marginale entière soit muette.
Contre-exemple : `γ(u) = γ_0 − (u_1−u_2)²` a un indice projeté **plat** (`Δ_1 = 0`), et pourtant
`P(Y > y | U_1 = 1/2) > P(Y > y | U_1 = 0)` — une statistique marginale d'ordre supérieur la voit.
*En compensation*, O fournit un **renforcement** : avec `p = 3`, phases `(u_2+u_3)` contre
`(u_1+u_2+u_3)`, **la liste complète des marginales bivariées de population** ne suffit pas à
identifier `A`. C'est la bonne façon d'énoncer « aucune méthode marginale ».

**Correction 3 — et c'est la plus lourde : (S) telle qu'imprimée n'identifie pas `A`.**
L'impossibilité sous dépendance libre est **confirmée et plus forte** que je l'avais notée : `A`
n'est pas seulement difficile à estimer, ce **n'est pas une fonctionnelle de la loi observable** ;
la définition du papier traite un choix arbitraire de version hors-support comme une vérité
structurelle. O ajoute une **proposition** : dans le sous-modèle Pareto exact, `A` est identifiable
pour tout `γ` admissible **si et seulement si** `supp(P_U) = [0,1]^p` (preuve dans les deux sens).
Or le (S) du manuscrit n'exige le support de fibre plein que pour `u ∈ I_ε` — ce qui **n'implique
pas** le support joint plein. Contre-exemple explicite de O : densité jointe à marginales uniformes,
absolument continue, satisfaisant (S) **exactement comme imprimé**, mais dont le support omet un carré
ouvert `J × J` avec `J ⊂ (0,ε)` ; deux indices y coïncident presque sûrement avec `A = ∅` et
`A = {1,2}`. Et (NS) seule ne répare rien non plus (modèle constant contre `1+(u_1−u_2)²`).
**Donc mon « (S) implique (NA), donc tout va bien » était faux.** Les seules réparations valides :
viser `D_U` (défini par la loi, invariant), ou exiger `supp(U) = [0,1]^p` **plus** `A = D_U` si l'on
tient à `A`. Sur l'exemple diagonal, `D_U = {1,2}` sous les **deux** représentations, avec
`Δ_1 = Δ_2 = 1/2` : la coordonnée 2, proxy observationnel parfait, appartient à la cible identifiable
même si une extension la déclare inactive. C'est exactement l'invariance recherchée.

**Attention supplémentaire** : « activité modulo le support » n'est **pas** bien définie non plus —
sur la diagonale, `1+z` s'écrit aussi bien `g_1(u_1)` que `g_2(u_2)`, donc `{1}` et `{2}` sont deux
représentations minimales. Il faudrait quotienter par les relations déterministes entre coordonnées.
`D_U` évite l'ambiguïté ; c'est un argument de plus pour en faire la cible officielle.

**Dette de rédaction n° 11** : le papier énonce (S) sur `I_ε` et en déduit l'identification de `A`.
C'est insuffisant. Soit viser `D_U`, soit exiger le support joint plein.

- **P — « Statistical Manuscript Review »**
  <https://chatgpt.com/c/6a813c5c-8c98-83ed-a50e-517baec26772>
  Consigne : spécifier les **expériences numériques** que la théorie corrigée exige désormais —
  estimateurs, modèles, tailles, grilles de réglage ; ce que chaque expérience confirmerait ou
  réfuterait ; quels nombres publiés doivent être recalculés ou retirés ; et, là où une expérience
  pourrait contredire un théorème, lequel et comment. *Motif* : la section de simulation du papier
  teste une théorie qui n'est plus la sienne — taille de liste non certifiée, score par blocs jamais
  simulé, scores robustes inédits, taux `nα` et non `nαh`.

#### Résultat P — programme numérique : 8 expériences, et l'inventaire des nombres à retirer

**L'expérience décisive (n° 1).** Modèle `R0` : `log Y = γE` avec `Y` **indépendant** de tous les
`U` — donc **aucun** biais de queue, spatial ou de quadrature. On y isole la variance pure. Le test
qui tranche : comparer `(α,h) ∈ {(0,05 ; 0,20), (0,10 ; 0,10), (0,20 ; 0,05)}`. Les trois ont le
**même `αh`**, donc le même compte local d'extrêmes — mais `nα` varie d'un **facteur 4**. L'ancienne
théorie en `nαh` prédit des variances comparables ; la théorie corrigée prédit un rapport de 4.
Régression descriptive `log Var(Ψ̂) = c − log n − log α + β_h log h` : on attend `β_h → 0` pour les
scores intégrés et `β_h → −1` pour le profil ponctuel. Plus des tests **non asymptotiques** :
comparer le taux de dépassement observé au `2e^{−x}` du théorème pour `x ∈ {3,4,5,6}`.

**Les sept autres** : comparaison propre avec le score par blocs jamais simulé (n° 2, avec balayage
en `n`, absent aujourd'hui) ; agrégation testée **à la taille de liste qu'elle certifie**, avec
assertions dures dans le code (`#{j : r_j^min ≤ R} ≤ 9R` doit tenir à **chaque** réplication — toute
violation est un bug ou une réfutation, sans qualificatif Monte-Carlo) (n° 3) ; transition de liste
`nαΔ² / log(p/d)`, tracée aussi sous l'ancienne normalisation pour montrer laquelle fait collapser
les courbes (n° 4) ; séparation `A` contre `D_U` par deux modèles construits exprès — proxy
structurellement inactif mais marginalement détectable, et actif marginalement invisible (n° 5) ;
codimension cachée avec `s ∈ {2,4,8,16}` **à gap de population constant** `0,186…`, de sorte qu'une
dégradation ne puisse pas être imputée à un signal plus faible (n° 6) ; robustesse avec cinq attaques
dont **R2**, l'attaque adaptative qui teste précisément le télescopage — remonter la ligne de rang 1
au rang `n` et faire exploser sa réponse — en vérifiant que la sensibilité **ne croît pas en `1/h`**
(n° 7) ; couverture simultanée avec la distinction explicite entre diagnostic *oracle* et analyse
*implémentable* (n° 8).

**Nombres publiés à retirer.**
- Le calcul `√(log(pn)/(nαh)) ≈ 0,67` : **caduc**. À `nα ≈ 204`, le terme corrigé vaut `≈ 0,184`.
- Toute phrase faisant de `nαh` la taille effective du score. La valeur `2nαh ≈ 65` reste, mais
  comme **compte local d'extrêmes** seulement.
- La frontière `nαhΔ² ≍ log p` → `nαΔ² ≍ log(p/d)`, plus la condition séparée de faisabilité locale.
- « L'agrégation certifie Sure-20 » : elle certifie le **top 36**. Sure-20 reste publiable comme
  performance empirique, explicitement étiquetée.
- La borne supérieure par argmin pour `d > 1`.
- Tout intervalle dit « honnête » dont le rayon contient le `b_0(a)` inconnu : c'est un oracle.
- L'ensemble extérieur au seuil zéro comme taille pilotée par les données.

**À recalculer** : les heatmaps et le réglage par défaut (aujourd'hui `(a*,b*) = (0,30 ; 0,15)`,
moyenne 0,780), avec sélection du défaut sur un **échantillon pilote indépendant** — pas sur les
200 jeux qui servent à afficher la carte ; les tables d'agrégation, en ajoutant Sure-36 et les
quantiles 90/95 du pire rang ; les temps d'exécution ; et l'application aux données de criminalité
si l'estimateur recommandé change.

**À conserver** : gaps de population, profils projetés, exposants de coin, bornes spectrales, calculs
de queue B1 — tous revérifiés par M. Et le modèle B1 lui-même, qui teste une différence
scientifiquement importante entre screening sur quantile fixe et sur indice de queue.

*Conséquence de la lacune d'archive* : P signale que la constante de calibration `μ_k` du score
annulaire doit être **récupérée dans la preuve de L** avant gel du code ; en attendant, la
précalculer par Monte-Carlo à erreur `< 10^{-4}` et **l'étiqueter comme approximation**.

### Vague 7 — 2026-08-16
Pièce jointe : `maths/paste-wave9.txt` (327 ko), journal incluant O et P.

- **Q — « Structural Set Identification »**
  <https://chatgpt.com/c/6a81437e-3090-83eb-86e8-8e3293bebd31>
  Consigne : faire le travail structurel que la campagne a **contourné**. O a montré que `A` n'est
  identifiable que si le support est le cube entier, que (S) imprimée ne suffit pas, et que
  « l'activité sur le support » ne donne **pas** un ensemble minimal unique (sur la diagonale, `{1}`
  et `{2}` sont deux représentations minimales). La campagne a répondu en abandonnant `A` pour `D_U`.
  Q doit faire l'inverse : définir proprement l'équivalence induite par les **relations
  déterministes** entre coordonnées sur le support, démontrer quelles sont les représentations
  minimales et quand elles sont uniques, caractériser **quel objet quotient** est identifiable, et le
  relier à `D_U`. Dire ce qui est récupérable et ce qui ne l'est pas, avec preuves.
  *Motif* : abandonner `A` est la bonne décision statistique, mais laisse un trou conceptuel — on ne
  sait pas encore *ce qu'on renonce à connaître*.

- **R — « Statistical Manuscript Analysis »**
  <https://chatgpt.com/c/6a814380-7a80-83eb-8375-bca8a52c2092>
  Consigne : l'**application réelle**. Le papier finit sur les taux de criminalité de 1 993 communes
  américaines, `p = 99`, échantillon de queue effectif ≈ 200. Tout l'amont a changé : cible
  support-définie, taux `nα`, liste certifiée plus longue que celle rapportée, intervalles simultanés
  désormais disponibles, point de rupture d'une ligne. Que peut encore affirmer cette section, que ne
  peut-elle plus, que faut-il recalculer — et **l'échantillon de queue effectif est-il seulement assez
  grand** pour que les conditions corrigées tiennent ? Être quantitatif sur le régime.
  *Motif* : c'est la seule partie du papier que la campagne n'a jamais examinée, et c'est celle qui
  parle au lecteur.

#### Résultat Q — le bon objet structurel n'est pas un ensemble de coordonnées

**Thèse.** Abandonner `A` n'était pas nécessaire ; ce qui était faux, c'est de chercher un
**ensemble**. L'objet canonique est une **antichaîne dans un ensemble ordonné quotient** d'information
engendrée par les coordonnées.

**Construction.** Sur le support `S`, préordre déterministe `A ⪯_S B` ⟺ « `B` détermine `A` partout
sur `S` ». Trois formulations équivalentes (Prop. 1) : la condition ponctuelle ; l'existence d'une
unique application **continue** `r_{A←B}` avec `π_A = r_{A←B} ∘ π_B` (`π_B` est une application
quotient d'un compact, d'où la continuité) ; et l'inclusion des tribus engendrées `F_A^S ⊆ F_B^S`.
D'où l'opérateur de clôture `cl_S` et le quotient `Q_S = 2^{[p]}/∼_S`.
Point important : c'est une équivalence sur les **parties**, pas sur les coordonnées — des relations
de groupe comme `U_3 = h(U_1,U_2)` n'induisent aucune équivalence entre singletons.

**Trois minimalités qui se séparent.** Inclusion-minimale ; information-minimale (dans le quotient) ;
générateurs minimaux d'une classe. **Contre-exemple décisif** : `U_1 = X`, `U_2 = V` indépendantes,
`W = (X+V)/2`, `U_3 = H(W)`, `θ = W`. Alors `{3}` et `{1,2}` sont **toutes deux**
inclusion-minimales, mais `{3} ≺_S {1,2}` strictement. Donc quotienter *après coup* la famille des
ensembles inclusion-minimaux ne suffit pas : **la minimisation doit avoir lieu dans le quotient**.

**Unicité, deux critères distincts.** Ensemble brut unique ⟺ l'intersection de toutes les
représentations est elle-même une représentation (Thm 4). Classe quotient unique ⟺ le champ
d'information commun est **engendré par des coordonnées** (Thm 5). Et l'unicité peut échouer même
dans le quotient : avec `U_2 = U_1` ou `1−U_1` à pile ou face et `θ = 2 + a cos(2πu_1)`, l'antichaîne
a **deux** classes incomparables, l'information commune étant `cos(2πU_1)`, qu'aucune partie de
coordonnées n'engendre. Quotienter règle les doublons interchangeables ; cela ne fabrique pas
d'unicité quand deux canaux distincts portent chacun le signal **plus** un supplément différent.

**Pourquoi le support plein est exactement le seuil.** Sur un produit plein, `A ⪯_S B ⟺ A ⊆ B`, et
l'on passe d'un point à l'autre **une coordonnée à la fois** en restant dans le cube — c'est ce
télescopage qui donne l'unique plus petite représentation. Sur un support non produit, les points
intermédiaires peuvent **sortir de `S`** : le raisonnement s'effondre exactement là.

**Lien avec `D_U` (Thm 8), et c'est la vraie surprise.** Avec `M_A = ess sup(Z | F_A)` :
`A` est une représentation suffisante **⟺** `γ*_U − E M_A` atteint son **maximum** `Δ_* = γ*_U − EZ`.
Donc « représentation suffisante » = **gap de groupe maximal**, et **non** « gap positif ».
`D_U` pose une autre question : est-ce que l'information d'une seule coordonnée **exclut** quelque
part les indices globalement maximaux.

**Aucune inclusion dans un sens ni dans l'autre**, avec deux exemples :
- *Structurellement indispensable, gap marginal nul* : `Z = γ_0 − c(U_1−U_2)²` — l'unique
  représentation minimale est `[{1,2}]`, et pourtant `Δ_1 = Δ_2 = 0`, donc `D_U = ∅`.
- *Proxy à gap positif, non suffisant* : `U_2 = H((X+V)/2)`, `Z = 1+X` — la coordonnée 2 n'appartient
  à aucune représentation minimale suffisante, mais `2 ∈ D_U`.

**⚠️~~Défaut du *trimming*~~ — RECTIFIÉ par l'audit T, voir ci-dessous.** Q affirmait que le gap
tronqué n'est pas invariant sur les classes. C'est vrai pour le quotient **presque sûr**, faux pour le
quotient **du support** — donc **faux dans le cadre du papier**. Détail dans le résultat T.

**Identifiable ≠ estimable uniformément (Thm 9).** Mélanger `P_0` (où `U_1 = U_2`) avec `Q`
(indépendantes) à taux `η_n = n^{-2}` donne `‖·‖_TV ≤ 1/n` entre les expériences à `n` points, alors
que les objets structurels **diffèrent**. Aucune procédure ne récupère donc uniformément les
relations déterministes exactes sans hypothèse séparant dépendance déterministe et quasi déterministe.

**Bilan.** Récupérables : `γ*_U`, `Z` modulo `P`, le quotient, l'antichaîne `A_P(Z)`, `D_U` ; et sous
continuité, `S`, `γ|_S`, `Q_S`, `A_S(γ|_S)`. Non récupérables : `γ` hors de `S` ; `A_cube` si
`S ≠ [0,1]^p` ; un représentant **privilégié** dans une classe ; une classe minimale unique quand
l'antichaîne en a plusieurs incomparables ; et les relations déterministes exactes, uniformément.
→ **dette de rédaction n° 12** : remplacer l'ensemble actif structurel par l'antichaîne
`A_S(γ|_S)`, qui **redonne** exactement la définition du papier dans le cas du support produit plein.

#### Résultat R — l'application réelle : les rangs survivent, leur statut change

**Le point qui sauve tout le travail numérique déjà fait** : `Δ_j − Δ_k = Ψ_k − Ψ_j`. Changer la
ligne de base du maximum sur le cube au maximum sur le support **ne change aucun classement**. Tous
les rangs calculés restent des résumés de données valides ; c'est leur **interprétation** qui change.

**Le régime, en chiffres.** `n = 1993`, `p = 99`, réglage retenu `α = 0,102`, `h = 0,160` :
`nα = 204` mais `nαh = 32,6`. Le taux corrigé est donc bien meilleur — le repère à constante unité
passe de **0,611** (ancien calcul du papier) à **0,187**. Mais :
- Les comptes locaux ne sont **pas** uniformes : « environ 65 » est un chiffre d'**intérieur**. Sur les
  neuf réglages, le plus petit `K` local va de **15 à 57**, l'intérieur de 20 à 95. Au réglage
  `(0,40 ; 0,20)`, `log(pn)/(nαh) = 1,168` — franchement non asymptotique.
- **La calibration simultanée échoue pour les neuf réglages** : le théorème exige
  `K+1 ≥ 4 log(8MpL/δ) = 77,4`, or le maximum atteint dans tout le bloc est **58**, et **43** au
  réglage par défaut. Même à un seul réglage il faudrait ≈ 68,7. Cela ne prouve pas que le classement
  est mauvais ; cela prouve qu'**on ne peut pas y attacher la garantie non asymptotique**.
- Le biais est le vrai problème : `1/log(1/α) = 0,439`, et avec une codimension `κ = 3` on est à
  `1,32` avant même les constantes du modèle. `α ≈ 0,10` n'est pas un niveau extrême.
- **Rupture d'une ligne, aggravée par la géométrie** : toutes les largeurs de bande du bloc dépassent
  `ε = 0,05`, donc **chaque** rang de covariable appartient à au moins une fenêtre tronquée — une
  seule réponse aberrante affecte le score de **toutes** les coordonnées. Et agréger neuf classements
  de Hill n'y change rien, puisque les neuf partagent la ligne contaminée.
- Échelle de contamination : `ε_c/α` vaut `0,049` à 0,5 % et `0,098` à 1 % — **du même ordre que
  l'échelle de détection**. Une contamination de 1 % suffit à effacer le signal.
- Inflation de liste à `p = 99` : le top 36 fait **36,4 %** des variables ; si `|D_U| ≥ 11`, la
  garantie top-`9s` devient **vide**.
- Intervalle de ligne de base illustratif à `k = 204` : `[0,249 ; 0,479]` — déjà large **avant**
  d'ajouter les intervalles de score.

**Ce qui survit, avec une hiérarchie que je trouve juste** : (1) *consensus* — pauvreté et structure
familiale biparentale, bien classées par des fonctionnelles de cible différentes ; (2) *candidats
stables de gap de queue* — taux d'hommes jamais mariés et usage des transports publics, dont l'écart
aux screens par quantile est stable sur les neuf réglages ; (3) *candidat sensible au réglage* —
emploi agricole, qui passe du rang 2 au rang 47 ; (4) *promotions d'un seul réglage* — plomberie et
grands logements, dont les rangs au réglage par défaut sont 36 et 46, à ne **pas** présenter comme
des preuves comparables.

**Ce qui ne peut plus être affirmé** : que le screen trouve les variables structurellement actives ;
que les variables de tête « déterminent » l'exposant limite ; que le top 20 est certifié ; que la
stabilité sur neuf réglages démontre une **robustesse** (c'est de la stabilité de réglage, pas de la
résistance à la contamination) ; et que l'application se situe dans un régime où les conditions
suffisantes sont vérifiées. La formule de l'abstract « sépare les prédicteurs de l'exposant limite
des prédicteurs d'un quantile élevé » devient : « met en évidence des désaccords entre deux
classements empiriques dont les cibles de population diffèrent ».

**Conclusion de R, que je fais mienne** : l'application reste scientifiquement utilisable comme
exercice de screening exploratoire transparent et soumis à analyse de sensibilité ; elle n'est pas
soutenable comme découverte certifiée de variables d'indice de queue. → **dette n° 13**.

### Vague 8 — 2026-08-16

- **S — « Statistical Manuscript Analysis »**
  <https://chatgpt.com/c/6a814be9-94a4-83eb-9b83-6b82377f0cd3>
  Consigne : **la loi du maximum sur `p`** — le seul point du cahier qu'aucune vague n'a traité.
  Tous les énoncés simultanés reposent sur des unions bornées : les seuils sont conservateurs et les
  tests ne sont **pas calibrés**. Il faut une approximation gaussienne ou de valeurs extrêmes pour
  `max_j (Ψ̂_j − Ψ_j)`, avec taux d'erreur explicite, ou un schéma multiplicateur démontré ; puis le
  test calibré de gap nul et une procédure FDR. Difficulté réelle : les scores sont **fortement
  dépendants** entre coordonnées (fenêtres de rangs chevauchantes, réponses extrêmes partagées) et
  `p` croît. Si c'est impossible sous les conditions du papier, le démontrer.
  → dernier morceau de **III**.

- **T — « Mathematical Manuscript Review »**
  <https://chatgpt.com/c/6a814bec-61a8-83ed-9470-0921d39015f0>
  Consigne : casser la théorie structurelle de Q, que personne n'a vérifiée — bonne définition des
  quotients, les équivalences en sont-elles vraiment, les contre-exemples font-ils ce qu'on leur
  fait dire, l'argument de continuité est-il correctement employé, et la non-invariance du trimming
  est-elle réelle ou un artefact. → **V.14** pour Q.

#### Résultat T — audit de Q : le cœur tient, **deux quotients avaient été confondus**

**Le diagnostic.** Q a mélangé deux notions : la détermination **exacte partout sur `S`** (quotient
`Q_S`, applications de factorisation **continues**) et la détermination **presque sûre** (quotient
`Q_P`, applications seulement mesurables). Les propriétés ne sont pas les mêmes, et Q passe de l'une
à l'autre sans le dire.

**Survivent** : préordre, clôture, quotient, antichaîne, les deux critères d'unicité, le
contre-exemple `W = (X+V)/2`, le contre-exemple de non-unicité par réflexion (avec l'information
commune identifiée **exactement** comme `σ{cos(2πU_1)}`), et la non-estimabilité uniforme des
relations déterministes.

**Faux : la non-invariance du trimming, dans le cadre du papier.** L'argument est joli.
Si `[{j}]_S = [{k}]_S` avec marginales uniformes, alors `T` est une bijection **continue** de `[0,1]`
— donc strictement monotone — et préserve la mesure. Les deux seules possibilités sont donc
`T(u) = u` ou `T(u) = 1−u`. **Or l'intervalle de trimming `I_ε = [ε, 1−ε]` est symétrique**, donc
préservé par les deux. Conclusion : `[{j}]_S = [{k}]_S ⟹ Δ_j^ε = Δ_k^ε`. Le gap tronqué **est**
invariant sur les classes singleton du quotient du support.
Le contre-exemple de Q exigeait une réarrangement **discontinu** préservant la mesure (demi-rotation),
exclu par la continuité même qui sert à construire `Q_S`. Comme le papier travaille après transformée
intégrale (marginales uniformes) et suppose `γ` continue, **il est dans le cas `Q_S`** : la dette que
j'avais notée **ne s'applique pas**. En revanche la non-invariance est bien réelle sur `Q_P`, avec
demi-rotation et `d` supportée dans `I_ε \ T(I_ε)` — et le gap **non tronqué** reste invariant des
deux côtés.

**Deuxième nuance.** « Gap de groupe maximal ⟺ suffisance » vaut pour la suffisance **presque sûre**,
pas automatiquement pour la factorisation **exacte sur le support**. Contre-exemple : demi-rotation et
`θ(x,y) = 2+y` — `Z` est `σ(U_1)`-mesurable donc gap maximal, mais l'adhérence du graphe contient
`(1/2, 0)` **et** `(1/2, 1)`, donc `θ|_S` ne se factorise pas par `π_{\{1\}}`.

**Troisième nuance.** « Le support plein est le seuil exact » ne vaut que pour l'identifiabilité de
l'ensemble actif de l'**extension ambiante**, avec quantificateur pire-cas sur `γ`. Ce n'est **pas**
un seuil pour que le quotient soit bien posé ni booléen : sur la bande circulaire
`S_a = {d_T(x,y) ≤ a}`, support **propre** du carré à marginales uniformes, on a quand même
`A ⪯ B ⟺ A ⊆ B`.

**Quatrième nuance.** Le contre-exemple `W` de Q ne sépare que **deux** des trois notions de
minimalité ; pour la troisième (générateurs minimaux d'une classe fixée) il faut l'exemple trivial
`U_1 = U_2 = X`, où `{1}` et `{2}` sont deux générateurs minimaux distincts de la même classe.

**Bilan** : la thèse principale de Q survit — l'objet structurel est une antichaîne de classes
d'information, pas un ensemble de coordonnées — mais l'énoncé final doit distinguer explicitement
`Q_S` (continu, identifié sous continuité de `γ`, se réduit à l'ensemble actif du papier sur le cube
plein) de `Q_P` (presque sûr, mesurable). C'est la **troisième fois** qu'un audit corrige ce journal
et non une session : je notais un défaut du papier qui n'en était pas un.

#### Résultat S — loi du maximum sur `p` : impossibilité, puis réparation par multiplicateurs

**1. Aucune loi de valeurs extrêmes universelle.** Sous dépendance inter-coordonnées non restreinte,
il ne peut exister de limite de type Gumbel **sans covariance**. Le maximum doit être approché par le
maximum d'un **vecteur gaussien à covariance arbitraire, estimée sur les données**.

**2. Impossibilité sous la classe actuelle du papier** — et la raison est inattendue. Même une
approximation gaussienne ou par multiplicateurs est **fausse** pour la statistique centrée sur la
population `max_j √(nα)(Ψ̂_j − Ψ_j)`. La cause n'est **ni** le chevauchement des fenêtres, **ni** la
dépendance entre coordonnées : c'est l'absence de restriction **à l'échelle du TCL** sur le biais de
queue à seuil fini. La classe autorise des modèles nuls où la statistique a une **dérive
déterministe tendant vers l'infini**. C'est la même racine que la barrière `1/log n` de la vague 1 :
le certificat de biais manquant, encore lui, mais cette fois il tue la **calibration** et pas
seulement la largeur d'intervalle.

**3. Réparation, sous quatre ingrédients ajoutés** : représentation asymptotiquement linéaire par
ligne ; correction de biais à l'échelle TCL, ou sous-lissage ; variances de coordonnées non
dégénérées ; estimation consistante de la covariance des influences de ligne. Alors un
**multiplicateur gaussien commun aux lignes** est valide — et il préserve **automatiquement** toute
la dépendance engendrée par les réponses extrêmes partagées, les fenêtres de rangs chevauchantes et
la dépendance arbitraire entre coordonnées, puisqu'il rééchantillonne la ligne entière. Erreur
d'approximation générique `C{(log⁷(2pn)/(nν_n))^{1/6} + δ_n^{1/3}}`.

**Contexte bibliographique** (S a consulté la littérature) : bornes gaussiennes sur les rectangles en
grande dimension type Chernozhukov–Chetverikov–Kato ; un travail récent construit un multiplicateur
par lignes pour des estimateurs de Hill en grande dimension sous dépendance transversale non
restreinte, avec la condition caractéristique `k_min/log⁷p → ∞`. **À citer** dans le papier — le
schéma n'est pas sans précédent, ce qui le rend plus crédible mais impose l'attribution.

**FDR** : sous dépendance non restreinte entre coordonnées, la règle step-up valide est
**Benjamini–Yekutieli** avec correction harmonique, pas Benjamini–Hochberg. À corriger partout où le
journal disait « FDR » sans préciser.

→ **III complété** : loi limite à coordonnée fixée (vague 3), intervalles simultanés (vague 3),
et désormais calibration du maximum — conditionnelle, avec l'impossibilité qui explique pourquoi la
version inconditionnelle n'existe pas.

#### Résultat K2 — théorie unifiée produite (94 min de réflexion), **non capturable**

K2 a réfléchi **1 h 34** et produit le document complet : préambule LaTeX, titre « Unified corrected
theory for high-dimensional tail-index screening », résumé annonçant la cible support-définie. Le
texte est dans le fil, mais **impossible à extraire** : le bouton de copie du bloc écrit du **vide**
(bloc virtualisé), le bouton « copier la réponse » n'écrit rien non plus, et `get_page_text` tronque
les blocs de code. Le presse-papier **système** fonctionne (testé) — c'est l'écriture **depuis la
page** qui échoue sans focus de fenêtre.
*Contournement engagé* : demander à K2 de réémettre le document **par tranches de ~6 000 caractères,
en LaTeX brut sans balise de code**, une par message, avec marqueur `PIECE k/N END`, ce qui rend
chaque tranche lisible par `get_page_text`. Coût estimé : ~8 tranches.

- **U — audit de S** (même fil que T)
  <https://chatgpt.com/c/6a814bec-61a8-83ed-9470-0921d39015f0>
  Consigne : attaquer les trois affirmations de S — l'absence de limite EVT sans covariance, la
  construction de **dérive déterministe** censée vivre dans la classe du papier, et la validité du
  multiplicateur par lignes avec son taux en racine sixième. La question centrale : la linéarisation
  par ligne est-elle plausible pour une fonctionnelle de Hill sur **fenêtre de rangs**, où chaque ligne
  entre dans de nombreuses fenêtres ?

#### Résultat U — audit de S : réparations et une condition renforcée

**1. Pas de limite EVT sans covariance : confirmé**, avec un argument net. Coordonnées identiques ⟹
`max_j G_j = Z` (aucune croissance extrême) ; coordonnées indépendantes ⟹ Gumbel après normalisation
`√(2 log p)`. Aucun centrage-échelle déterministe commun ne peut couvrir les deux : les médianes sont
`0` et `√(2 log p)`, les échelles de fluctuation `1` et `1/√(log p)`.
**Et U vérifie que le contre-exemple est bien dans le modèle** — c'est le point que je voulais tester :
sous « dépendance non restreinte », oui ; et même **sous support de fibre plein**, en prenant
`X_ij = Z_i + τ_n ε_ij` avec `τ_n` assez petit, la densité jointe reste partout positive (support
plein) alors que **toutes les permutations de rangs coïncident** avec probabilité tendant vers 1.
Le support plein ne sauve donc pas.

**2. L'obstruction par dérive : correcte, mais la construction de S était invalide.**
Écrire `ℓ_n(s) = s^{d_n}` **n'est pas** une fonction à variation lente (`ℓ_n(ts)/ℓ_n(s) = t^{d_n}`
ne tend pas vers 1) — la construction telle qu'énoncée sortait de la classe. U la **répare par
recollement retardé** : `S_n = n²`, `d_n = 1/log n`, quantile `s^{γ+d_n}` puis `S_n^{d_n} s^γ`.
Vérifications complètes : variation lente à `n` fixé, `b_n(3α) ≤ d_n`, coefficient limite
`S_n^{d_n/γ} = e^{2/γ}` **borné uniformément**. Sur un événement de probabilité `≥ 1 − 1/n`, toutes
les observations sont **avant** le recollement, donc l'échantillon se comporte exactement comme un
Pareto d'indice `γ + d_n` : `√(nα)(Ψ̂ − Ψ) = √(nα)/log n + O_P(1) → ∞`. Cela arrive **déjà à `p = 1`**,
sans aucun rapport avec le chevauchement ni la dépendance.

**3. Objection subtile, et sa réparation.** Le théorème de screening imprimé suppose `Δ_min > 0` :
un null global à indice constant n'est **pas littéralement** dans son domaine — U note que c'est
« déjà un défaut sérieux pour un théorème d'inférence censé calibrer des hypothèses nulles ». Puis il
**enfonce le clou** en plongeant la dérive dans un modèle à **gap positif fixe** (`γ(u) = γ_0 − δu`,
même recollement), avec réglages explicites `α_n = n^{-1/2}`, `h_n = n^{-1/3}` : toutes les conditions
du manuscrit sont vérifiées (`nα_n h_n = n^{1/6} → ∞`, `h_n log n = o(1)`), et la dérive persiste.
Ce n'est donc pas un artefact de null.

**4. Reformulation exigée.** Pas « l'approximation gaussienne est fausse », mais : *une approximation
gaussienne **centrée en zéro** de la statistique **centrée sur la population** est invalide sans
contrôle du biais à l'échelle du TCL*. Elle redevient possible en centrant par `E Ψ̂_j`, en visant un
**pseudo-paramètre à seuil fini**, ou en estimant et corrigeant le biais.

**5. La condition de S était trop faible — correction importante.** Pour un TCL coordonnée par
coordonnée il suffit de `√(nα)·max_j|β_nj| → 0`. Pour approcher un **maximum**, il faut
**`√(nα log p)·max_j |β_nj|/σ_j → 0`** : décaler chaque coordonnée gaussienne de `d` déplace la loi du
maximum de l'ordre de `d√(log p)` par anti-concentration. Le « biais à l'échelle du TCL » de S ne
suffit donc pas.

**6. Linéarisation par ligne : plausible, non établie.** C'était ma question prioritaire. Les
ingrédients existent (PIT conditionnel + Rényi + moyennage signé), et la preuve de borne inférieure
exhibe déjà une **première projection de Hoeffding** `Σ_i (N_i/L_n) h_1(Z_i)`, ce qui atteste
l'existence d'un terme dominant par ligne. Mais — et c'est le point — « une borne inférieure obtenue
de la première projection ne montre pas que le reste dégénéré est négligeable ». Le multiplicateur
reste donc un **théorème de haut niveau non démontré pour le score de Hill sur fenêtres de rangs**.

**7. Taux.** Le terme en racine sixième n'est une borne CCK valide qu'**après** hypothèses de moments
explicites, et ce n'est pas la meilleure disponible aujourd'hui ; le terme en racine cubique n'est
correct **qu'avec son facteur logarithmique**, et la simple consistance de la covariance ne suffit pas.

→ **III** reste donc ouvert sur son dernier point : la calibration exige soit un pseudo-paramètre à
seuil fini, soit une correction de biais explicite, **et** une preuve de linéarisation par ligne qui
n'existe pas encore.

- **V — représentation linéaire par ligne** (fil de T/U)
  <https://chatgpt.com/c/6a814bec-61a8-83ed-9470-0921d39015f0>
  Consigne : combler le trou que U vient d'identifier. Démontrer **ou réfuter** une représentation
  asymptotiquement linéaire du score glissant intégré, avec reste explicite : fonction d'influence
  d'une ligne complète, contrôle des composantes de Hoeffding **dégénérées** d'ordre ≥ 2 sous la
  géométrie de fenêtres chevauchantes et les statistiques d'ordre supérieures partagées, et
  conditions sur `n, α, h, p` rendant le reste négligeable **avec la marge `√(log p)`** qu'exige le
  maximum. Puis, si elle tient, le théorème de multiplicateur avec son **vrai** taux et le
  **pseudo-paramètre à seuil fini** qu'il calibre réellement.

*Note d'archivage — K2.* Les cinq tranches existent dans le fil (`1/30` puis un regroupement
`5/5`, plus un bloc intermédiaire de 34 184 caractères). **La transcription verbatim depuis le texte
de page est abandonnée** : le rendu convertit les mathématiques en Unicode et perd des commandes
(`\esssup{...}`, `\phi_T`, indices), de sorte que le fichier reconstitué ne compilerait pas et
donnerait une **fausse impression d'archive fidèle**. Sont conservées les parties 1 et 5
(`wave6-K2-part1.tex`, contenu de la 5 consigné ci-dessus). L'archive propre attend une capture par
presse-papier, qui exige que la fenêtre du navigateur ait le focus.
*Contenu vérifié de la partie 5* : scores robustes par paires winsorisées (transformée `m_T`, lemme
de sensibilité locale `4T/k_r`, lemme de télescopage avec la remarque sur ses conditions d'échec),
construction annulaire (calibration `log(3/2) ≤ μ_k ≤ log 4`, non-expansivité `ℓ¹` des moyennes
élaguées, stabilité au remplacement), bornes inférieures de contamination (intersection de Huber,
variation totale exacte entre exponentielles, transfert Huber → remplacement fort par troncature
binomiale), **proposition finale « quel modèle d'adversaire exige chaque preuve »** en quatre points,
et chaîne de théorèmes en une page.

- **W — « Audit of Mathematical Journal »**
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  Consigne : auditer **ce journal lui-même**, pas les mathématiques. Chercher les endroits où deux
  entrées se contredisent, où une correction tardive n'a pas été répercutée sur une affirmation
  antérieure, où une conclusion est énoncée plus largement que le résultat qui la porte, où une
  quantité ou une condition est écrite différemment à deux endroits, et où quelque chose noté comme
  acquis dépend d'un résultat affaibli depuis. Classer par risque d'induire en erreur quiconque
  réécrirait le papier **à partir du seul journal**. Citer les passages en conflit.
  *Motif* : trois audits ont déjà corrigé ce journal plutôt qu'une session. Le journal est devenu le
  document de référence de la campagne — il est donc lui-même une source d'erreur, et personne ne
  l'avait encore traité comme telle.

#### Résultat V — la linéarisation par ligne : plausible, **et c'est un théorème neuf à écrire**

**L'annulation de la largeur de bande, une troisième fois — au niveau de la fonction d'influence.**
Pour une log-queue exponentielle, l'influence de premier ordre à fraction `α` est
`χ_α(Z) = α^{-1}[(Z − q_α)_+ − γ·1{Z > q_α}]`, d'espérance nulle et de variance `γ²/α` (par absence
de mémoire). L'influence d'une **ligne** dans le score intégré vaut
`φ_ij = (n/L_n) Σ_{r : i ∈ W_jr} χ_jr/m_r`. Or `#{r : i ∈ W_jr} = O(nh)`, `m_r ≍ nh`, `L_n ≍ n`,
donc `(n/L_n)Σ_r 1/m_r = O(1)` : d'où `Var(φ_ij) ≍ 1/α` et `Var(n^{-1}Σφ) ≍ 1/(nα)`.
**Aucun `h^{-1}` ne subsiste.** C'est la version « fonction d'influence » du même phénomène déjà vu
par classes de couleur (vague 2) puis par télescopage (vague 5). Trois dérivations indépendantes du
même fait : le score **intégré** est structurellement mieux conduit que le profil.

**Six choses restent à démontrer** — et V les liste précisément :
1. Un développement de Bahadur d'ordre intermédiaire **uniforme** sur tous les seuils locaux.
2. **Le développement doit se faire au niveau du score intégré.** Si l'on obtient des influences
   locales d'enveloppe `1/(αh)` et qu'on applique un théorème en grande dimension **avant** de
   moyenner, le paramètre de moment garde le `h` et l'on retombe sur `nαh`. *L'annulation doit
   précéder l'approximation gaussienne* — remarque méthodologique décisive.
3. Les rangs empiriques sont des fonctionnelles **globales** de chaque colonne : il faut un
   développement conditionnel en tableau triangulaire, pas une preuve i.i.d. classique.
4. Le reste doit être contrôlé **à l'échelle du maximum** : `√(log p)·max_j √(nα)|r_nj| = o_P(1)`.
5. Un multiplicateur **implémentable** exige des contributions par ligne **observables** ; les
   variables de Rényi sont des artefacts de preuve, pas des résidus observés. Il faut un
   développement par seuil estimé, du cross-fitting, ou une projection jackknife.
6. Les PIT randomisés : les randomisations auxiliaires propres à chaque coordonnée ne doivent **pas**
   déterminer la covariance de la statistique réelle.

**Verdict de V** : « une linéarisation par ligne au niveau du score est plausible, mais c'est un
**théorème neuf**, pas une conséquence de la preuve de concentration existante. Le chevauchement des
fenêtres n'est pas fatal : une fois la représentation acquise, tout le chevauchement est encodé dans
`φ_ij`. Le difficile est de démontrer cette représentation. »

**Le multiplicateur.** Un multiplicateur **par ligne** (`e_i` commun à toutes les coordonnées) ;
utiliser des `e_ij` indépendants par coordonnée **effacerait la covariance transversale** et serait
invalide.

**Le taux, corrigé sur trois points.**
- `B_n² ≍ 1/α` est justifié : pour l'influence Pareto exacte, `E|X|^{2+k} ≍ α^{-k/2}` et
  `‖X‖_{ψ₁} ≍ α^{-1/2}`. La **forme** du taux en racine sixième est donc correcte.
- Mais une condition de **moment/Orlicz explicite est indispensable** : sans elle, des valeurs rares
  énormes tuent l'approximation **même à `p = 1`**.
- Et le taux **n'est pas optimal** : un résultat plus récent donne `{log⁵(pn)/(nα)}^{1/4}`, et des
  bornes proches de `n^{-1/2}` existent sous conditions plus fortes. « Si le résultat qualifie ce
  taux d'optimal, cette affirmation meurt. »

**Le terme en racine cubique était mal énoncé.** La comparaison gaussienne vaut
`Δ^{1/3}{1 ∨ log(p/Δ)}^{2/3}` — **avec** son facteur logarithmique. Et `Δ_n = o_P(1)` **ne suffit
pas** quand `p → ∞` : il faut `Δ{1 ∨ log(p/Δ)}² = o_P(1)`. Contre-exemple : `Δ = 1/log p` tend vers
zéro sans que la borne suive.

**Dichotomie élégante sur le reste** : si l'on a une borne en haute probabilité `r_n`, la
contribution est **linéaire**, `r_n√(log p)` — pas de racine cubique. La racine cubique n'apparaît
que si l'on ne dispose que d'un moment d'ordre 2, `ρ_n²`, par optimisation de
`ρ_n²/a² + Ca√(log p)`. « Un "reste en racine cubique" générique, sans préciser la norme ni le type
de borne, n'est pas un théorème. »

#### État de l'archive (à jour)

**Sur disque et compilés** : `wave1-A-addendum.tex/.pdf` (15 p.), `wave2-D-patch.tex`,
`wave3-F-inference.tex/.pdf` (12 p.), `wave3-G-appendix.tex/.pdf` (16 p.), plus `wave1-A-answer.md`
et `wave2-C-answer.md`.

**Uniquement dans les fils** (URL au journal, à recapturer) : l'audit H de F, l'audit J de G, le
document de I (seuillage robuste), la théorie unifiée de K.

*Pourquoi* : la capture disque passe par `cmd+c` → `pbpaste`, et le presse-papier système n'est
accessible que si la fenêtre Brave a le **focus**. Sans focus, la copie échoue **en silence** et
`pbpaste` renvoie l'ancien contenu — d'où l'obligation absolue de comparer la capture aux fichiers
existants avant d'écrire (`cmp` contre `maths/*.tex`), garde-fou qui a déjà évité deux doublons
trompeurs. Les contournements testés et **écartés** : `navigator.clipboard.writeText` (refusé hors
focus), téléchargement de liens `sandbox:` (ne part pas), lecture du DOM d'un bloc de code
(virtualisé : 1 704 caractères rendus pour un document de plusieurs dizaines de milliers),
API interne du site (`/backend-api/conversation/<id>` → HTTP 404).
**Seule solution connue** : mettre la fenêtre au premier plan, puis capturer en une passe.

*Méthode d'envoi (corrigée).* Le presse-papier échoue dès que la fenêtre Brave perd le focus système.
L'outil d'**upload de fichier** ne dépend pas du focus : c'est désormais la voie normale.
Vérifier avant envoi : `#pièces jointes == 1` **et** longueur du composeur ≈ celle de la consigne
(un collage résiduel peut arriver en différé et gonfler le composeur à 200 000 caractères).

*Frappe clavier (contournement).* Les événements clavier synthétiques n'atteignent que l'onglet
**actif du navigateur** : dans un onglet d'arrière-plan, `type` ne dépose rien, même après avoir
donné le focus DOM à l'élément. Contournement qui marche partout :
`el.focus(); document.execCommand('insertText', false, texte)` pour écrire, puis
`document.querySelector('#composer-submit-button').click()` pour envoyer — les deux passent par le
DOM et ignorent l'onglet actif. C'est la méthode à utiliser dès qu'on pilote plusieurs sessions en
parallèle. Faire la vérification pièces jointes/longueur **dans le même script** que le clic d'envoi,
avec un `ABORT` si le compte est faux.

#### ⚠️ INSTANTANÉ HISTORIQUE — PÉRIMÉ, NE PAS UTILISER (état arrêté à la vague 2)
> Conservé pour la traçabilité. **Le registre en vigueur est le « Registre au terme de la vague 8 »,
> à la fin de ce fichier.** Plusieurs statuts ci-dessous ont été renversés depuis, et la liste de
> dettes qui suit s'arrête à l'item 10 alors qu'il en existe treize.

#### Bilan du cahier après la vague 2 *(périmé)*

- **Atteints, sous réserve de la vérification adverse restante** : I.1 (version corrigée `Δ^U_j`),
  I.2 (par impossibilité : trois théorèmes négatifs), I.3 (stabilité en distance de Hausdorff sur les
  supports, avec contre-exemple montrant que TV/Wasserstein/faible ne suffisent pas),
  II.4 et II.5 (borne inférieure de variance appariée pour l'estimateur du papier + seuil
  `nαΔ² ≍ log(p/d)` par Fano avec liste, atteint), II.6 ((S), continuité, (E2), (E3↑), (P3), la
  troncature globale et `log(pn) = o(nh²)` — toutes tombées ou strictement affaiblies).
- **Partiels** : II.7 (oracle « union des tops » pour l'agrégation ; manquent le choix pilotable de
  `h, α, d` sans connaître `Δ_D`, et une borne `e` utilisable) ; III (loi limite à coordonnée fixée
  pour le score par blocs ; manquent la loi du maximum sur `p`, un test uniformément valide de
  `Δ_j = 0`, un estimateur de `γ*_U`, FWER/FDR, la correction de biais).
- **Non entamés** : IV (mélange/dépendance temporelle, contamination `ε_n`, `γ` irrégulier ou indice
  nul), V.14 pour les résultats **de la vague 2 eux-mêmes** (personne n'a encore attaqué D).

#### Dettes de rédaction accumulées (à répercuter dans le manuscrit)

1. Cible : `γ*_U` et `D_U`, pas `γ*` sur le cube ni `A`.
2. Abstract : `nαh` n'est pas la taille d'échantillon effective du score ; c'est `nα`.
3. Abstract : la dépendance « non restreinte » est incompatible avec la récupération de `A`.
4. Agrégation : remplacer la justification empirique par la garantie « union des tops » (D4), et
   retirer toute lecture « un bon réglage suffit ».
5. Annexe C : les constantes ne sont pas indépendantes de `s_n` ; expliciter `κ_n`.
6. Lemme du coin gaussien : compléter la preuve (domination uniforme, frontière mobile).
7. Le score par blocs disjoints reste utile (indépendance littérale, hypothèses de biais moyennées
   plus faibles) mais **n'est pas nécessaire** pour retirer le `h` : ne pas l'imposer comme méthode.
8. Agrégation : le théorème certifie le top `9s` = **36**, pas le top 20 des simulations. Soit on
   change la taille de liste rapportée, soit on énonce la condition de sauvetage `min_λ r_j(λ) ≤ 2`,
   soit on présente le top 20 comme un constat empirique explicitement non couvert par le théorème.
9. Ne jamais écrire « optimal à constante près » pour la borne inférieure de variance : l'appariement
   est en **taux**. Et énoncer `α < 1/6` là où le sous-modèle de Fano est invoqué.
10. Screening en liste : ajouter le lemme des faux dépassants (`nαΔ_U² ≫ log(s + (p−s)/(d−s+1))`),
   qui est le bon énoncé quand `d > s` — le papier ne demande la séparation exacte que par habitude.

*Note de méthode.* La transmission par presse-papier ne fonctionne que dans l'onglet au premier plan,
et l'action `type` écrase le presse-papier : coller **d'abord** la pièce jointe, taper la consigne
**ensuite**. Récupération des réponses : `get_page_text` (le bouton « copier » de ChatGPT n'écrit pas
dans le presse-papier système depuis un onglet non focalisé).

*Capture des réponses (résolu).* Les liens de fichiers `sandbox:` de ChatGPT ne se téléchargent pas
depuis ce pilotage, et l'API JS `navigator.clipboard.writeText` est refusée. En revanche **`cmd+a`
puis `cmd+c` dans l'onglet au premier plan fonctionnent** : le texte part dans le presse-papier
système, `pbpaste > fichier` l'archive sans passer par le contexte du modèle. C'est la méthode à
utiliser. Corollaire : **toujours exiger des sessions qu'elles écrivent les preuves dans le corps du
message**, jamais en pièce jointe.
Pour un bloc de code, ne **jamais** mesurer le DOM (`textContent` ne rend que la partie virtualisée
et sous-compte massivement) : cliquer le bouton « copier » du bloc, puis `pbpaste`.

Archivé : `maths/wave1-A-answer.md` (réponse complète de la session A) et
`maths/wave1-A-addendum.tex` — l'addendum complet, 1 168 lignes, 16 énoncés et 16 preuves.
**Compile sans erreur** → `maths/wave1-A-addendum.pdf`, 15 pages. Premier élément du critère V.15
rempli pour ce lot ; reste la vérification adverse (V.14, confiée à la session C) et la cohérence
numérique des constantes.

*Incident.* La session D a reçu le bundle **deux fois** : un collage cru « raté » avait en réalité
créé une pièce jointe invisible à mon sélecteur, et `cmd+a`+`Delete` n'efface que le texte du
composeur, jamais les pièces jointes. Contenu correct mais dupliqué. Avant tout envoi : compter les
pièces jointes réellement attachées au message, ne pas se fier à la longueur du composeur.

---

### Vague 9 — en cours

- **X — IV.12, la ligne jamais traitée** (fil de W)
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  Consigne : que devient la théorie quand `γ` est **irrégulière** — seulement mesurable, non minorée,
  nulle sur un ensemble de mesure positive (réponse **pas** à queue lourde sur une partie de l'espace),
  ou tendant vers zéro au bord ? Que survit-il de l'identité de projection, de la cible
  support-définie, de la concentration du score, du screening ? Contre-exemple là où ça casse, preuve
  et condition la plus faible là où ça tient.
  *Motif* : W a établi que **personne** n'avait traité cette ligne du cahier en neuf vagues — je
  l'avais perdue de vue en marquant IV « clos ».

- **Y — « Statistics Manuscript Audit »**
  <https://chatgpt.com/c/6a8179e2-d0a8-83eb-8801-9b5a89069e28>
  Consigne : livrer les **deux items que l'audit W a démasqués comme substitués** plutôt que faits.
  (i) Le transfert depuis une **vraie classe de second ordre** `(ρ, A(·))` vers tout ce que les
  preuves de score consomment, avec le compromis biais-variance **explicite** pour le réglage —
  demandé par II.6, jamais écrit ; à la place on avait substitué un jeu de conditions directes.
  (ii) L'**estimateur local à biais réduit** et son gain de taux — demandé par III, jamais produit.
  Avec la remarque qui les relie : **c'est le même problème**, une classe de second ordre étant
  exactement ce dont un estimateur à biais réduit a besoin. Plus : ce que la classe de second ordre
  **coûte en généralité** face aux conditions directes.

- **Z — audit de Y** (fil de W/X)
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  Consigne : attaquer l'ensemble du résultat Y — l'asymptotique de Laplace et l'inversion de quantile
  qui donnent `ρ = 0` ; l'affirmation que **toute** correction à rapports fixes est inutile en
  `ρ = 0` et que c'est une obstruction d'identifiabilité ; la condition de log-régularité et
  **est-elle satisfaite par les modèles du manuscrit** ; et le gain de taux annoncé.
  *Motif* : c'est le résultat le plus contre-intuitif de la campagne — il déclare inopérant l'outil
  standard de réduction de biais en valeurs extrêmes. S'il est juste, il change une pratique ; s'il
  est faux, il ferait écrire une bêtise au papier.

#### Résultat Z — audit de Y : **noyau correct, presque toutes les conséquences surévaluées**

**Ce qui survit.** Pour une fibre **diffuse** dont le déficit pondéré est à variation régulière
*authentique*, `H_D(z) ∼ a z^κ L(1/z)`, la transformée de Laplace est à variation régulière d'indice
`−κ` (Karamata), le théorème de densité monotone donne `−tB'(t)/B(t) → κ`, et l'inversion de quantile
donne `ρ = 0`, `A(s) ∼ −κξ/log s`. **Cette partie est juste.**

**Correction 1 — « à peu près `z^κ` » ne suffit pas.** Si l'on suppose seulement `H_D(z) ≍ z^κ`,
**rien ne suit**. Contre-exemple : densité `h(z) = z^{κ−1}{κ + ε cos(log z)}`, dont la transformée
contient `t^{−κ}C cos(log t + φ)`, si bien que `−tB'/B` **ne converge pas**. Il faut la variation
régulière, pas un encadrement en puissance.

**Correction 2 — la condition de reste du manuscrit est trop faible pour que tout cela tienne.**
Le manuscrit suppose `sup_x|r(t,x)| → 0` et `sup_x|∂_t r| → 0`. La projection exige une condition
d'**incrément local** : `sup_{|a|≤M} t·|log R(t+a) − log R(t)| → 0`. Contre-exemple :
`ε(t) = t^{−1/2} sin t` satisfait les conditions imprimées (`ε → 0`, `ε' → 0`) mais son incrément est
d'ordre `t^{−1/2}`, qui **domine** l'incrément de Laplace `κa/t`. Donc
**`H_D ∼ az^κ` ne force PAS le second ordre annoncé dans la classe primitive actuelle du papier.**
*(Pour les quatre modèles de simulation, en revanche, le reste est assez rapide : A1–A3
super-exponentiel, B1 en `O(y^{−2})`.)* Et si le facteur primitif a une pente logarithmique `λ/t` du
même ordre, le coefficient devient `A(s) ∼ −(κ−λ)ξ/log s` : **le coefficient `−κξ` n'est pas garanti.**

**Correction 3 — « la projection force `ρ = 0` » est FAUX tel qu'énoncé.** L'énoncé exige la
**diffusion** (`κ > 0`), alors que la condition de fibre du manuscrit **autorise un atome** en déficit
nul. Deux contre-exemples : fibre dégénérée `D ≡ 0` (alors `B ≡ μ({0})`, aucun facteur à variation
lente, et si le quantile primitif a `ρ < 0` on a `ρ < 0`) ; et `μ = mδ_0 + (1−m)δ_δ`, qui donne
`Pr(Y>y) = my^{−1/ξ} + (1−m)y^{−1/ξ−δ}`, donc **`ρ = −δξ < 0`**.
Énoncé correct : *une masse de déficit quasi maximale **diffuse et à variation régulière** engendre un
terme dominant `ρ = 0`* — et non « la projection force toujours `ρ = 0` ».

**Correction 4 — le no-go n'est valide que pour des poids BORNÉS.** Avec des poids non bornés on
**peut** annuler : `H̃(s) = ((L+ℓ)/ℓ)·H(cs) − (L/ℓ)·H(s)` avec `L = log s`, `ℓ = log c`, a des poids de
somme 1 et **supprime le terme principal**. Prix : poids d'ordre `L`, donc bruit stochastique gonflé
d'autant. Et **connaître `κ`** suffit avec **un seul** seuil : `H(s)/(1 − κ/L) = ξ + o(1/L)`.
Donc l'affirmation « aucune combinaison » est trop large ; il faut dire « aucune combinaison à
rapports fixes **et poids bornés** ».

**Correction 5 — « obstruction d'identifiabilité » est la mauvaise étiquette.** L'identifiabilité
demanderait que deux valeurs de `κ` donnent la même loi observable ; or `κ` **est** identifiable :
`H(cs) − H(s) = a·log c/(log s)² + o(…)`. Le coefficient est présent, simplement **un ordre plus bas**.
C'est de la **faible identification**, pas de la non-identifiabilité — et aucun théorème d'impossibilité
n'en découle sans spécifier l'expérience statistique.

**Correction 6, la plus lourde — le gain annoncé `1/log²` est FAUX, réfuté par le modèle le plus
propre possible.** Prendre `D ∼ Gamma(κ,1)` et `Pr(Y > e^t | D) = e^{−(1/ξ+D)t}` : **aucun** reste
primitif, et la projection est **exacte** : `Pr(Y>e^t) = e^{−t/ξ}(1+t)^{−κ}`. L'inversion soignée donne
`q(L) = ξL − κξ log L + O(1)`, puis
`H(s) = ξ − κξ/log s − κ²ξ·loglog s/(log s)² + …`.
Après séparation en puissance : `H_PS(s) = ξ + (κ²ξ/τ)·loglog s/(log s)² + …`.
Le biais corrigé est donc `≍ loglog s / log² s`, **pas** `O(1/log² s)`. L'origine est l'inversion
elle-même : le terme `−κξ log log s` de `log U(s)` **revient** quand on développe la pente logarithmique.

**Correction 7 — et pour les modèles réels du papier, c'est pire.** Le développement de coin gaussien
`H_D(z) = C z^κ{√(2π)r_z}^{κ−d}e^{−b r_z}{1+O(1/r_z)}` avec `r_z ∼ √(2log(1/z))` introduit des
facteurs en `√(loglog s)`. Le biais après correction devient
`1/(log s·(loglog s)^{3/2})` ou `1/(log s·(loglog s)²)`. Z le confirme par un second contre-exemple
(`q(L) = ξL − a log L − 2c√(log L)`), qui satisfait la log-régularité **faible** et laisse pourtant un
biais `≫ L^{−2}`. Et il vérifie que **les coefficients sont génériquement non nuls** pour les modèles
du papier : `b_j(u) = (4/15)Φ^{−1}(u)` pour `j = 1,4`, `(8/15)Φ^{−1}(u)` pour `j = 2,3` — non nuls sauf
en `u = 1/2` ; et `κ_j − d = −11/15` ou `−4/15`, non nuls aussi.

**Correction 8 — la « log-régularité » n'était jamais définie.** La version faible
`A(s^τ)/A(s) → 1/τ` est **très** insuffisante : elle ne contrôle que l'auxiliaire, pas le reste
`o{A}` que la correction fait remonter. Il faut une condition **directement sur le biais de Hill**
(`H(s) = ξ − a/log s + O(1/log²s)`), et non sur une auxiliaire — car `A` n'est définie qu'à
`A(1+η)` près, de sorte qu'une propriété au **troisième** ordre peut valoir pour une représentation
et pas pour une autre. Sinon, fixer l'auxiliaire **canonique** `A_*(s) = d log U/d log s − ξ`.

**Bilan.** Y garde : le calcul de Laplace, l'inversion sous hypothèses renforcées, le no-go à poids
bornés, et l'idée de séparation en puissance — qui **échappe réellement** au no-go puisque le rapport
`s^{τ−1}` diverge. Y perd : l'universalité de `ρ = 0`, l'étiquette d'impossibilité, le taux
`1/log²`, et la condition de log-régularité telle qu'énoncée. **Le gain reste réel mais plus modeste,
et il faut le recalculer modèle par modèle.**
→ Dans le papier : ne **pas** annoncer de plancher `1/log²`. Annoncer une amélioration
`o(1/log s)` avec le taux exact dépendant de la géométrie du coin.

#### Résultat Y — second ordre et biais réduit : **la projection force `ρ = 0`**

**La découverte.** Proposition 3 : si la masse de déficit près du maximum de fibre se comporte comme
`H(z) ∼ a z^κ L(1/z)`, alors `B(t) ∼ aΓ(κ+1)t^{−κ}L(t)`, et le quantile projeté est de second ordre
avec **`ρ = 0`**, `A(s,u) ∼ −κξ_j(u)/log s`. Autrement dit : **le mécanisme même de projection, sur
les modèles du papier, produit un second ordre logarithmique et non polynomial.**

**Conséquence immédiate et fatale pour la méthode standard.** La réduction de biais classique (deux
seuils en rapport fixe `q`, poids `1/(1−q^{−ρ})`) a un dénominateur qui **s'annule** en `ρ = 0`. Et
Y montre que ce n'est **pas** un accident algébrique : pour **toute** combinaison linéaire finie à
rapports fixes préservant l'indice (`Σ w_ℓ = 1`), le biais résiduel vaut
`A(1/α)·Σ w_ℓ + o(A) = A(1/α)(1+o(1))` — **elle n'enlève rien**. C'est une obstruction
d'identifiabilité : à une seule échelle multiplicative de seuil, `A(t)log x` est indiscernable au
premier ordre d'un petit décalage de `ξ`.
→ **Un jackknife de Hill à rapports fixes ne couvre pas les modèles à coin gaussien du papier.**

**La réparation : seuils séparés en puissance.** Prendre `k_1 ∼ α^s m` au lieu de `qαm`, avec
`c_r = (log t_1/log t_0)^{−λ}`, sous une condition de **log-régularité** de l'auxiliaire
(`A(t^s,u) = s^{−λ}A(t,u) + O(a^{1+δ})`, soit `z ↦ A(e^z,u)` à variation régulière d'indice `−λ` ;
le cas des mélanges de fibres donne `λ = 1`). Gain : le plancher de biais de queue passe de
`1/log(n/Λ)` à **`1/log²(n/Λ)`**. Coût : la fraction plus extrême `α^s`, donc faisabilité locale
`nα^s h ≫ log(pn)`.

**Le transfert demandé (II.6), enfin écrit.** Classe quantitative `Q₂^-(ρ,δ,η;a)` avec forme
canonique `log Q_j = b_j(u) + ξ_j(u)log t + A_j(t,u)/ρ + R_j`, et théorème de transfert vers **tout**
ce que les preuves consomment : le module de queue `b_n(3α) ≤ C{a_α + a_α^{1+δ}}` **cesse d'être une
hypothèse**, la condition horizontale (HQ) suit avec `A(r), B(r) ≤ Cr^η`, et l'on récupère le terme
spatial **plus fin** `h^η log(1/α)` au lieu du `h log(pn)` d'origine.

**Le compromis biais-variance, avec une conséquence nette sur le réglage.** Il n'y a **pas d'optimum
intérieur en `h`** : une fois le faux facteur `h` retiré de la variance du score, **diminuer `h` est
toujours favorable** jusqu'à ce que la faisabilité locale morde. D'où `h* = ℓ_n log(pn)/(nα*)`.
Pour `α` : `α* ≍ (Λ/n)^{1/(2β+1)}` sans correction, `(Λ/n)^{1/(4β+1)}` avec (cas quadratique `δ=1`),
et le taux passe de `(Λ/n)^{β/(2β+1)}` à `(Λ/n)^{2β/(4β+1)}` — **strictement meilleur**.

**Ce que la réduction achète vraiment**, et Y le dit honnêtement : **pas** un changement de la
frontière stochastique irréductible `nαΔ² ≍ Λ_{p,d}` — celle-là était déjà juste. Elle permet de
prendre **`α` bien plus grand à gap donné**, donc de fournir davantage d'extrêmes au problème
stochastique. La frontière de screening passe de `nΔ^{2+1/β} ≫ Λ` à `nΔ^{2+1/(2β)} ≫ Λ`.
Coût en variance : facteur constant `v_ρ(q)` — par exemple **5** pour `ρ = −1`, `q = 1/2`.

**Et une mise en garde que je retiens** (§9) : le second ordre **standard** ne donne que
`reste = o(A(t))`, ce qui n'implique **pas** `O(A^{1+δ})`. Contre-exemple explicite :
ajouter `t^ρ ε(t) sin(log t)` avec `ε ↓ 0` arbitrairement lentement — la loi reste 2RV, mais la
correction à rapports fixes laisse un terme `t^ρ ε(t)` qui décroît plus lentement que toute puissance
de `A`. Donc « une vraie classe de second ordre » et « un gain de taux explicite » sont deux demandes
**presque** identiques et pas tout à fait : *l'auxiliaire de second ordre identifie le terme à
annuler ; c'est un reste de troisième ordre quantitatif qui décide de ce qu'on gagne après annulation.*

**Le coût en généralité (ce que j'avais demandé d'expliciter).** C'est une **spécialisation stricte**,
pas un affaiblissement gratuit. Les conditions directes tolèrent : une auxiliaire changeant de signe,
irrégulière selon `j` et `u`, à vitesses sans rapport entre elles ; des biais logarithmiques ou plus
lents ; des profils contrôlés **en moyenne** seulement ; et des profils d'indice projeté à variation
bornée ou simplement mesurables. La classe de second ordre exige au contraire un seuil d'apparition
commun, une auxiliaire à variation régulière, un `ρ` **borné loin de zéro** dans le régime polynomial,
un module de reste quantitatif, et un contrôle **hölderien horizontal** de l'échelle, de l'indice, de
l'auxiliaire et du reste.
**Recommandation qui en découle** : garder les conditions directes comme cadre principal, et
présenter la classe de second ordre comme un **régime spécialisé** où l'on obtient en plus un réglage
optimal explicite et une réduction de biais — en signalant que le cas `ρ < 0` **ne couvre pas** les
modèles du papier, et que c'est la version `ρ = 0` qui s'applique.

→ **II.6 livré dans son contenu d'origine** ; **III : la moitié « biais » de la calibration est
close** — reste la représentation par ligne et l'influence observable (V).

#### Résultat X — indices irréguliers et nuls : IV.12 enfin traitée

X commence par **séparer quatre questions** que le cahier confondait : la projection tient-elle ?
la cible reste-t-elle définie ? le score vise-t-il encore cette cible ? le biais est-il assez petit ?

**Théorème 1 — projection sous uniformité inclinée.** Le critère exact passe par la loi de mélange
**inclinée par le dépassement**, `Π_y(dx) = F̄_x(y)K(dx)/F̄_K(y)`, et l'identité
`F̄_K(cy)/F̄_K(y) = ∫ {F̄_x(cy)/F̄_x(y)} Π_y(dx)` — qui contient tout le problème, puisque la mesure
`Π_y` **bouge avec `y`**. Deux conditions suffisent : (a) uniformité **sous `Π_y`** de la convergence
des ratios, (b) concentration de `Π_y` sur les indices proches du sup essentiel. **Ni continuité,
ni atteinte du maximum, ni minoration de `γ`.**

**Théorème 2 — Pareto exact avec `γ = 0` autorisé.** Avec `Z = 1/γ − 1/m ∈ [0,∞]` (convention
`1/0 = ∞`), `F̄_K(y) = y^{-1/m}B(log y)` où `B(t) = E e^{−tZ}`. Comme `ess inf Z = 0`, l'inclinée de
`Z` tend vers 0, donc `B(t+a)/B(t) → 1` : `B(log y)` est **à variation lente**. Couvre `γ`
discontinue, `γ = 0` sur un ensemble de mesure **arbitrairement grande**, `γ ↓ 0` au bord, et un
supremum **approché sans être atteint**. Conclusion : **la partie d'indice nul ne contribue en rien
à la queue asymptotique dès que `m > 0`.**

**Deux contre-exemples, dont un très frappant.**
- *Indices conditionnels tous égaux à 1, mélange d'indice `1/β > 1`* : mécanisme d'apparition
  retardée déjà connu.
- **Tous les indices conditionnels sont NULS, et le mélange est Pareto d'indice 1.** Prendre
  `Y | U = u ∼ Exp(u)` : pour tout `u > 0` la queue est exponentielle donc `γ(u) = 0`, mais
  `∫₀¹ e^{−uy}du = (1−e^{−y})/y ∼ 1/y`. Une famille de lois **légères** dont les taux s'approchent de
  la dégénérescence **se mélange en une queue lourde**. Conséquence directe : un estimateur de Hill
  global **n'estime pas nécessairement zéro** quand tous les indices conditionnels sont nuls.
  Réparation minimale : l'instance « indice nul » de la condition inclinée ; suffisant plus simple :
  **variation rapide uniforme** (`F̄_x(y) ≤ y^{−M}` pour tout `M`, uniformément).

**La cible.** Pour `γ` seulement mesurable, `γ*_U = max_{supp}` **n'est plus valide** — seule
`ess sup γ(U)` est invariante par la loi (contre-exemple : deux versions différant en un point de
mesure nulle donnent le même modèle observable et des maxima 0 et 1). Le journal écrivait parfois le
maximum sur le support : **licite seulement sous continuité**. `ξ_j` reste mesurable, via
`{ξ_j > q} = {K_j(u, {γ > q}) > 0}`.

**Théorème 3 — détectabilité sans ensemble de maximiseurs.** L'ensemble `M = {γ = γ*_U}` peut être
**vide** alors que le sup essentiel est positif. La caractérisation par contact avec `M` doit être
remplacée par les **ensembles de niveau quasi maximaux** `M_η = {γ > γ*_U − η}` :
`ξ_j(u) = γ*_U ⟺ K_j(u, M_η) > 0 pour tout η`, et
`Δ_j > 0 ⟺ ∃η > 0 : λ{u ∈ I_ε : K_j(u, M_η) = 0} > 0`.
Contre-exemple montrant la nécessité : `γ(u₁,u₂) = u₂` pour `u₂ < 1` et `0` en `u₂ = 1` — ici
`M = ∅`, et pourtant `Δ₁ = 0`, ce que seule la formulation quasi maximale donne correctement.

**Et une limite honnête de la cible.** Si `γ*_U = 0`, alors `γ(U) = 0` p.s., donc tous les `ξ_j`,
tous les `Ψ_j` et tous les `Δ_j` sont nuls : **`D_U = ∅`**. Ce n'est *pas* un échec d'estimation,
c'est une limite de l'objet : **l'indice de valeurs extrêmes ne distingue pas les vitesses de
décroissance des queues légères**. À écrire tel quel dans le papier.

**Quatre conclusions de X** : aucune minoration de `γ` n'est nécessaire pour la cible ni pour une
borne supérieure additive de concentration ; **la seule mesurabilité ne suffit pas** pour que le
score glissant implémenté vise cette cible ; `γ*_U = 0 ⟹ D_U = ∅` ; et les hypothèses d'indice nul
**ponctuelles** ne survivent pas au mélange sans condition uniforme ou inclinée.

**Suite de X (parties II et III), lue plus tard — et c'est le plus important.**

**Les fenêtres locales prennent des SUPREMA ESSENTIELS, pas des moyennes locales.** Sous les
hypothèses de projection, l'indice de la loi fenêtrée vaut
`ξ_{j,h}(u) = ess sup_{v ∈ W_h(u)} ξ_j(v)`. Le Hill local idéal ne vise donc **pas** une moyenne
mobile de `ξ_j` mais un **supremum essentiel mobile**. Tout le problème de localisation est là.

**Condition spatiale exacte la plus faible (EUR).** Avec l'enveloppe essentielle supérieure
`f^♯(u) = lim_{h↓0} ess sup_{|v−u|≤h} f(v)`, on a toujours `f ≤ f^♯` p.p., et la convergence
`Ψ_{j,h} → Ψ_j` a lieu **si et seulement si** `ξ_j^♯ = ξ_j` p.p. Version quantitative pour le
screening : `Ω_n(h_n) = max_j (1/|I_ε|)∫[ess sup_{|v−u|≤h_n} ξ_j(v) − ξ_j(u)]du = o(Δ_min)`.
**Strictement plus faible que la continuité uniforme, strictement plus fort que la mesurabilité.**

**Discontinuités inoffensives** : profils semi-continus supérieurement ; profils continus par
morceaux à frontières de mesure nulle ; profils hölderiens **même s'annulant** (`ξ(u) = u^β` ne pose
aucun problème structurel) ; et une **région entière d'indice nul** `ξ ≡ 0` sur `[0,a]` ne coûte
qu'une couche limite `O(h)`. **Le fait que `γ` s'annule n'est donc pas le problème.**

**Le problème, c'est l'entrelacement dense — et le contre-exemple est spectaculaire.**
Soit `F` un **ensemble de Cantor gras** (fermé, d'intérieur vide, de mesure `> 0`) et `G` son
complémentaire, ouvert **dense**. Poser `Y|U=u ∼ Exp(1)` si `u ∈ F`, `Pareto(1)` si `u ∈ G`. Alors
`γ = 0` sur `F`, `1` sur `G`, donc `γ*_U = 1` et `Δ_1 = λ(F)/|I_ε| > 0` : **la coordonnée a un gap
de population strictement positif**. Mais tout intervalle contient une part `w > 0` de Pareto, donc
pour **tout** `h > 0` fixé, `ξ_{1,h}(u) = 1` partout, d'où `Ψ_{1,h} = 1` et **gap nul**.
La cible visée par le score glissant est donc **zéro** alors que la vraie vaut `λ(F)/|I_ε|`.
Ce n'est **pas** un problème de reste de queue : les lois conditionnelles sont **exactement**
exponentielles et Pareto. C'est un échec **pur de localisation spatiale**.
Pire, X rend l'échec **adapté à n'importe quel réglage** : par construction en cascade (cellules de
taille `h_{n_k}/8`, fraction retirée `2^{-2k}`), on obtient `w/α_{n_k} → 0`, si bien que les
statistiques d'ordre supérieures sont engendrées par la composante Pareto **uniformément sur les
centres**. Conclusion : **aucun choix déterministe de `(h_n, α_n)` ne rend le score glissant
uniformément valide sur toutes les `γ` mesurables.**

**Zéro n'est pas inconsistant, mais il n'y a aucun taux uniforme.** Le Hill de population tend bien
vers `ξ` même pour `ξ = 0` (Potter + convergence dominée). En revanche, avec
`Q_T(t) = t^a` jusqu'à `T` puis `T^a{1+log(t/T)}`, l'indice **limite est zéro** alors que la cible
de Hill à seuil fini vaut `≈ a`. Version échantillonnale : si `T_n ≫ n`, tout l'échantillon est
**exactement** Pareto d'indice `a` avec probabilité tendante vers 1. Donc
`sup_{P : γ(P)=0} P(|γ̂| > a/2) ↛ 0` : **ni consistance uniforme, ni intervalle honnête, ni théorème
de screening** sur une classe non quantifiée de queues à variation lente. C'est l'**analogue à indice
nul de l'obstruction d'apparition retardée** — la même racine, encore.

**Et une remarque de portée sur l'estimateur** : la théorie de Hill est formulée pour des indices de
Fréchet **positifs**, alors que des procédures de type moment conditionnel couvrent Fréchet, Gumbel
et Weibull. Garder Hill est donc une **restriction substantielle**, à assumer explicitement.

→ **IV.12 traitée**, et bien plus richement que prévu. IV devient : dépendance temporelle et
contamination acquises sous conditions ; indices irréguliers **caractérisés**, avec la condition
spatiale exacte (EUR), la classification des discontinuités inoffensives, et **deux impossibilités
nouvelles** (poussière dense d'indices mélangés, absence de taux uniforme à indice nul).

# REGISTRE AU TERME DE LA VAGUE 8

*Établi après l'audit W du journal lui-même (28 défauts de tenue de registre, dont 6 critiques).
Ce registre **remplace** tout bilan antérieur. En cas de conflit avec une entrée plus ancienne,
c'est celui-ci qui fait foi.*

## Conventions de notation — à respecter partout

L'audit a relevé trois dérives dangereuses. Elles sont closes par les conventions suivantes.

| Symbole | Sens **unique** désormais |
|---|---|
| `γ*_U` | `ess sup γ(U)` — maximum **sur le support observé**. Jamais le maximum sur le cube. |
| `Δ^U_j` | `γ*_U − (1/|I_ε|)∫_{I_ε} ξ_j` — gap **tronqué, support-défini**, d'une coordonnée. |
| `Δ⁰_A` | gap **non tronqué** d'un groupe `A` (invariant sur le quotient d'information). |
| `D_U` | `{j : Δ^U_j > 0}` — cible identifiable du screening marginal. |
| `s_D` | `|D_U|` — taille de la cible de screening. **À utiliser dans « top 9·s »**. |
| `s_A` | taille d'une représentation structurelle minimale. **Distincte de `s_D`** : Q a montré qu'aucune inclusion ne vaut dans un sens ni dans l'autre. |
| `Δ_min` | `min_{j ∈ D_U} Δ^U_j`. Les anciennes notations `Δ_D`, `Δ_U`, `g_n` désignaient cette même quantité — **ne plus les employer**. |

## Statut réel des seize lignes du cahier

**Acquis, sous réserve d'audit adverse déjà effectué**
- **I.1** — détectabilité en équivalence, **version corrigée** (`Δ^U_j`, support-défini). La version
  d'origine, comparant au maximum sur le cube, était **fausse** (C).
- **I.2** — cas non détectable tranché **par l'impossibilité**, avec la portée corrigée par O
  (voir ci-dessous).
- **II.4 / II.5** — borne inférieure de variance appariée pour l'estimateur du papier (D, réparée par
  E via moment d'ordre 4 + Paley–Zygmund), et frontière `nαΔ² ≍ log(p/d)` par Fano avec liste, avec
  **sélecteur top-`d`** (l'argmin de D ne valait que `d = 1`). Appariement **en taux, pas en constante**.

**Partiels — ne pas présenter comme acquis**
- **I.3** — *le critère a été modifié en route*. Demandé : continuité quantitative **en la copule**.
  Obtenu : stabilité en **distance de Hausdorff sur les supports**, avec démonstration que la
  continuité échoue en variation totale, Wasserstein et convergence faible. C'est un théorème de
  remplacement légitime, **pas** ce qui était demandé.
- **II.6** — *contenu modifié en route*. Demandé : transfert vers une condition de second ordre
  standard `(ρ, A)` avec compromis biais-variance explicite. Obtenu : de nombreuses hypothèses
  affaiblies ou supprimées — (S), continuité, (E2), (E3↑), (P3), troncature globale,
  `log(pn) = o(nh²)` — via des conditions directes (HQ), (TD), déficit incliné. **Le transfert
  demandé n'a jamais été écrit.**
- **II.7** — post-sélection sur menu fini et découpage pilote **conditionnels à des certificats de
  biais connus** ; garantie « union des tops » avec inflation de liste. **Pas** d'adaptativité libre
  de nuisances : H a montré que l'adaptation ne vaut qu'entre candidats **déjà certifiés**, et
  l'impossibilité à apparition retardée interdit l'adaptation au biais inconnu.
- **III** — **ouvert**. Loi limite à coordonnée fixée (F) et intervalles simultanés à certificats
  (F, réparés par H) : acquis. Manquent : la calibration du maximum pour le score **réellement
  implémenté** (S a donné un gabarit de haut niveau ; U a montré que sa construction de dérive était
  invalide et l'a réparée ; V a établi que la linéarisation par ligne est **un théorème neuf**, avec
  six lemmes listés), **et** la construction explicite de Hill à biais réduit, que le critère
  exigeait et que personne n'a fournie.
- **IV** — **partiel, pas clos**. Dépendance temporelle (F, avec les deux conditions ajoutées par H :
  mélange uniforme sur la **ligne entière** en dimension croissante, et non-composition automatique
  avec le découpage pilote) et contamination (I/L, avec la preuve de I **remplacée** par le
  télescopage de N) : acquis. **IV.12 — `γ` irrégulier, non minoré, ou indice nul — n'a jamais été
  traité par aucune session.**
- **V.14** — audits rendus pour A, D, F, G, I/L (par arbitrage), Q, S. **Non audités** : la théorie
  unifiée de K2, les résultats de V, et cet audit W lui-même.
- **V.15** — **non rempli**. M a audité les nombres des **quatre modèles de simulation** et n'y a
  trouvé aucune erreur : c'est un résultat solide, mais c'est *son* périmètre. Le critère exige en
  outre la cohérence avec simulations et application — or P impose de retirer sept nombres publiés et
  R de requalifier l'application entière. Les constantes fines des théorèmes non asymptotiques
  (`1/32`, `6p+2pL`) **ne sont pas certifiées** ; seule la version à constantes universelles survit.
- **V.16 — saturation : non atteinte.** Chaque vague depuis le début produit un théorème, une
  réfutation ou une réparation.

## Portée exacte des trois impossibilités (corrigée par O)

1. **Apparition retardée** : impossibilité **uniforme** sur tableaux triangulaires, **pas ponctuelle**
   à loi fixée. Écrire « pas détectable **uniformément** ». Elle montre qu'**une** restriction
   excluant l'apparition arbitrairement retardée est nécessaire — **pas** que (P2) soit *la* condition
   nécessaire. `1/log n` et `1/log(1/α_n)` ne coïncident que si `log(1/α_n) ≍ log n`.
2. **Invisibilité marginale** : ne vaut que pour une décision mesurable par rapport au **seul couple
   `(U_j, Y)`**. Une procédure regardant une autre coordonnée distingue les modèles. Et `Δ_j = 0`
   **n'implique pas** que la loi marginale soit muette. L'énoncé fort correct porte sur la **liste
   complète des marginales bivariées** (`p = 3`, phases décalées).
3. **Non-identifiabilité sous dépendance libre** : confirmée et **plus forte** — `A` n'est pas une
   fonctionnelle de la loi. Réparation : viser `D_U`, ou exiger `supp(U) = [0,1]^p` **plus** `A = D_U`.
   **(S) telle qu'imprimée dans le papier ne suffit pas**, ni (NS).

## Contradictions littérales résolues

- *Litige du facteur `h`* : **clos** par D et confirmé par C, E, N et V (quatre dérivations). La
  phrase de C « reste ouvert : le score glissant atteint-il aussi `nα` » est **périmée**.
- *Preuve de robustesse de I* : le **résultat** survit, **la preuve doit être remplacée** par le
  télescopage de N, avec plafond de winsorisation **fixe** et départage déterministe des ex æquo.
  L'ancienne formule « aucune preuve n'est à jeter » était contradictoire — corrigée ici.
- *Agrégation* : ce n'est **ni** un oracle **ni** une intersection, mais une **couverture par union
  des tops avec inflation de liste** : top `M·s_D`, soit **36** pour `M = 9`, `s_D = 4` — pas 20.
- *Bacs fixes* : **suffisants, non nécessaires** (N). Et si les bacs sont des quantiles empiriques
  recalculés après contamination, ils ne sont pas fixes.
- *Recouvrement sans connaître `s`* : la méthode `max_j Ψ̂_j` de A est **fausse** sans coordonnée de
  gap nul (C). La réparation de F — estimer `γ*_U` par la queue **inconditionnelle** — est une
  construction **différente**, à ne pas présenter comme validant celle de A.

## Divergence numérique à trancher

**RÉSOLU par recalcul direct.** Ce n'était pas une contradiction : les deux sessions parlaient de
**deux designs différents**, et le journal ne le disait pas.

| | `n` | `p` | `nα` | `nαh` | ancien `√(log(pn)/(nαh))` | corrigé |
|---|---|---|---|---|---|---|
| **Simulation** (P) | 2000 | 1000 | 204,5 | 32,70 | **0,666** ≈ 0,67 | `√(log p/(nα))` = **0,184** |
| **Application** (R) | 1993 | 99 | 204,0 | 32,64 | **0,611** | `√(log 2p/(nα)) + lin.` = **0,187** |

Les deux anciens nombres sont **justes**, chacun pour son design ; ils diffèrent parce que
`log(pn)` vaut 14,5 en simulation contre 12,2 pour l'application. Les deux corrigés sont justes
aussi, mais **avec des conventions différentes** : P prend `√(log p/(nα))` seul, R ajoute le terme
linéaire et utilise `2p`. Coïncidence trompeuse : `0,184` et `0,187` se ressemblent alors qu'ils ne
mesurent pas la même chose sur le même design.
**Règle** : écrire systématiquement le design **et** la formule à côté de chaque nombre.

## Statut des fichiers, à ne pas confondre avec le statut mathématique

`wave3-G-appendix.tex` est sur disque **dans sa version d'avant l'audit J** : les deux lignes de
preuve manquantes, la réparation de version au bord et le retrait de l'affirmation sur `q = 20`
**n'ont pas été appliqués au fichier**. Idem pour les constantes affûtées de M. Le fichier compile,
mais il n'est pas à jour. Aucun fichier n'a été patché après audit.


## Dettes de rédaction 11–16 (ajoutées après les vagues 8–9)

11. **(S) sur `I_ε` n'identifie pas `A`** (O). Viser `D_U`, ou exiger `supp(U) = [0,1]^p` **et**
    `A = D_U`. Ne pas présenter (S) comme réglant l'identification.
12. **L'objet structurel est une antichaîne de classes d'information**, pas un ensemble de
    coordonnées (Q, corrigé par T). Distinguer explicitement le quotient du **support** (continu,
    identifié sous continuité de `γ`, redonne la définition du papier sur le cube plein) du quotient
    **presque sûr** (mesurable).
13. **L'application** : classement inchangé (`Δ_j − Δ_k = Ψ_k − Ψ_j`) mais **statut** changé — screening
    exploratoire, pas découverte certifiée. Retirer les sept nombres listés par P ; hiérarchiser
    consensus / stables / sensibles au réglage / promus par un seul réglage.
14. **La condition de reste primitive du papier est trop faible** (Z) pour que le second ordre projeté
    existe : `sup|r| → 0` et `sup|∂_t r| → 0` **n'impliquent pas** la condition d'incrément local
    `sup_{|a|≤M} t|log R(t+a) − log R(t)| → 0`. Contre-exemple `t^{−1/2} sin t`. À renforcer, ou à
    restreindre explicitement aux modèles où le reste est assez rapide (ce qui est le cas de A1–A3 et B1).
15. **Condition spatiale exacte** (X) : remplacer toute hypothèse de continuité du profil par
    `ξ_j^♯ = ξ_j` p.p. (enveloppe essentielle supérieure), avec la version quantitative
    `Ω_n(h_n) = o(Δ_min)`. Signaler que **les fenêtres visent un sup essentiel mobile**, pas une moyenne,
    et donner le contre-exemple de Cantor gras comme borne de validité.
16. **Garder Hill est une restriction** (X) : la théorie couvre les indices de Fréchet positifs ; les
    procédures de type moment conditionnel couvrent Fréchet, Gumbel et Weibull. À assumer explicitement,
    puisque `γ*_U = 0 ⟹ D_U = ∅`.

### Vague 10

- **AA — « Mathematical Assumptions and Proofs »**
  <https://chatgpt.com/c/6a818d98-9998-83ed-9352-e68707a6b12f>
  Consigne : **assembler une classe primitive minimale**. Les audits ont laissé les hypothèses
  primitives éparpillées — condition de reste trop faible (Z, avec contre-exemple `t^{−1/2}sin t`) ;
  condition spatiale exacte par enveloppe essentielle supérieure (X, avec le Cantor gras) ; uniformité
  inclinée pour la projection à indices nuls ou irréguliers (X). Personne n'a réuni tout cela.
  Il faut : l'ensemble d'hypothèses **le plus faible** dont découlent **tous** les résultats survivants,
  la preuve de chaque implication, et pour **chaque** hypothèse un contre-exemple montrant qu'on ne
  peut ni la retirer ni l'affaiblir. Plus : lesquels des quatre modèles de simulation la satisfont.
  *Motif* : c'est le préalable à toute réécriture — le papier ne peut pas énoncer quinze conditions
  directes hétérogènes, et la campagne a maintenant assez de contre-exemples pour savoir lesquelles
  sont **nécessaires**.

*Note de fichiers.* Un en-tête « STATUT D'AUDIT » a été ajouté en tête des quatre `.tex` archivés,
signalant qu'ils sont dans leur version d'origine et que les réparations ne leur ont pas été
appliquées. Sauvegardes en `.bak`. Vérifié : `wave1-A-addendum.tex` compile toujours.

#### Résultat AA — classe primitive minimale : **trois interfaces exactes, dix contre-exemples**

**Le diagnostic de départ, que je fais mien** : ne pas remplacer (P1)–(P4) par une autre liste plate
de conditions de régularité, ce qui remélangerait quatre questions distinctes —
**projection → cible de population → localisation → estimation à seuil fini**. D'où une classe
**stratifiée**, dont les clauses qualitatives sont **nécessaires et suffisantes**, et dont la clause
quantitative est **nécessairement dépendante de la suite** (les constructions à apparition retardée
interdisent toute condition ponctuelle fournissant un taux uniforme).

**Les trois interfaces, chacune exacte.**
- **(TU1)–(TU2), uniformité inclinée** — exacte pour la **projection**. La mesure sous laquelle
  l'uniformité est requise n'est pas `K` mais l'inclinée mobile `Π_{K,y}`. Et AA démontre que
  **(TU2) est *nécessaire* étant donné (TU1)** : conditionnellement à l'approximation des ratios,
  la concentration inclinée **équivaut** à la projection sur le sup essentiel. Ce n'est donc pas une
  hypothèse de commodité.
- **(CTI), condition d'incrément canonique de queue projetée** — AA démontre qu'elle est
  **équivalente** au second ordre projeté (Thm 5), avec `χ_{ρ,ξ}(a)` couvrant `ρ<0` et `ρ=0` d'un
  seul énoncé. C'est le remplaçant exact des deux conditions imprimées `sup|r| → 0`, `sup|∂_t r| → 0`,
  qui ne contrôlent pas la pente **relativement** au terme `1/t` engendré par la projection.
- **(EUR)** — exacte pour la **localisation**, avec l'équivalence prouvée dans les deux sens
  (Thm 3 : `Ψ_{j,h} → Ψ_j ⟺ ξ_j^♯ = ξ_j` p.p.), via l'identité
  `ξ_{j,h}(u) = ess sup_{v ∈ J_h(u)} ξ_j(v)`.

**Et une clause honnête** : pour la synchronisation au seuil de travail (WS), AA note qu'**aucun
audit n'a montré la nécessité** d'une paramétrisation hölderienne ou par couplage particulière. La
bonne hypothèse est donc le **module réel** `q_n(ρ_n, α_n) = o(Δ_min)`, les conditions de couplage
et de Hölder devenant des **lemmes de vérification**. C'est exactement la distinction que la campagne
avait perdue de vue en substituant des conditions les unes aux autres.

**Sharpness, clause par clause (§12), dix contre-exemples** : base de référence par sup essentiel
(diagonale) ; variation régulière ponctuelle insuffisante (seuil d'apparition `u^{-1/β}` non borné) ;
indices nuls non préservés (`Exp(u)` → Pareto) ; condition de reste imprimée insuffisante
(`t^{-1/2} sin t`) ; comparabilité en puissance insuffisante (`κ + ε cos(log z)`) ; diffusion
nécessaire pour `ρ=0` (atome + gap) ; mesurabilité insuffisante (Cantor gras) ; seuil commun
indispensable (recollement retardé) ; faisabilité locale indispensable ; frontière stochastique non
améliorable (Fano avec liste).

**Suppressions explicites (§13)** — à retirer du papier, pas à affaiblir : support de fibre plein ;
atomlessness ; continuité et équicontinuité ; **minoration positive de `γ`** ; doublement uniforme du
déficit — avec un **contre-exemple neuf** : atomes `d_m ↓ 0` de masses telles que
`p_m/Σ_{ℓ>m} p_ℓ → ∞`, où le doublement échoue alors que la moyenne inclinée tend bien vers zéro ;
et `log(pn) = o(nh²)`, artefact de la route DKW abandonnée.

**Vérification des quatre modèles (§14)** : tous satisfont la classe, **branche diffuse `ρ = 0`**.
Explicitement : `Ω_n(h) = O(h)` (par monotonicité des profils `e^{-c_j u}/2`),
`q_n = O(n^{-b} log n)`, exposants de coin `34/15 ≤ κ_j ≤ 43/15`, restes primitifs assez rapides
(A1–A3 super-exponentiels, B1 en `O(e^{-2t})`) pour satisfaire la condition d'incrément renforcée
que la classe imprimée n'avait pas. **Réglages admissibles : `a < 1`, `b > 0`, `a + b < 1`** — la
contrainte `b < 1/2` venait de la route DKW et **n'est plus primitive**. Les neuf paires du bloc
d'agrégation la satisfont.

**Qualification (§15), cohérente avec Z** : les quatre modèles **ne** satisfont **pas** l'énoncé de
troisième ordre uniformément — le facteur de coin gaussien engendre des corrections en `loglog` et
`√loglog`. Donc séparation en puissance ⟹ `o(1/log s)` **dépendant du modèle et de la coordonnée**,
jamais le plancher universel `1/log²`.

**Théorème de transfert primitif minimal** (l'énoncé final) : cinq clauses (i)–(v) plus deux
conditions de vitesse, dont découlent : indices projetés = sups essentiels conditionnels, cible `D_U`,
taux stochastique `1/√(nα)` du score intégré, sure screening top-`d`, et recouvrement exact à
`d = s_D`. **Minimale au sens fort** : TU exacte pour la projection, CTI équivalente au 2RV projeté,
EUR exacte pour la localisation ; tout ce qui est plus fort est superflu, tout ce qui est plus faible
aux trois interfaces critiques est réfuté par un contre-exemple explicite.

→ **II.6 livré**, et l'annexe primitive du papier a désormais un **remplaçant unique** au lieu de
quinze conditions hétérogènes.

- **AB — audit de AA** (fil de W/X/Z)
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  Consigne : attaquer la classe primitive minimale. Les **trois équivalences** en sont-elles vraiment,
  **dans les deux sens** ? L'argument de **nécessité** de la clause de concentration inclinée est-il
  correct ? Les contre-exemples de suppression font-ils ce qu'on leur fait dire ? Et la vérification
  des quatre modèles, avec les exposants admissibles `a<1, b>0, a+b<1`, est-elle juste ?
  *Motif* : c'est le document qui remplacerait l'annexe primitive entière du papier. Une erreur y
  contaminerait tout le reste, puisque tous les théorèmes en dépendraient.

### Verdict de AB (23 min 48 s) — la classe primitive est une *bonne interface*, pas une classe minimale

**Sentence globale.** L'architecture stratifiée est jugée « substantiellement meilleure » que le
bloc primitif plat du manuscrit (qui empaquetait borne inférieure positive + Lipschitz + support de
fibre plein + non-atomicité + doublement uniforme du déficit en un seul lot suffisant). **Mais le
mot « exact » est employé trop librement.** Classification finale, interface par interface :

| Interface | Verdict de AB |
|---|---|
| Projection inclinée (TU1)+(TU2) | implication directe **vraie** ; nécessité de (TU2) **vraie mais seulement conditionnellement à (TU1)** ; la prétention « ⟺ » pour la paire entière est **fausse** |
| Incrément canonique (CTI) | **vraiment équivalent** à 2RV non dégénérée, ρ<0 et ρ=0 ensemble — **mais** seulement après ajout de : RV du premier ordre en couche séparée, non-annulation éventuelle de `a`, convergence **pour tous** les multiplicateurs (pas seulement au pas canonique), et une **branche Pareto exacte** distincte |
| Enveloppe supérieure (EUR) | **exactement équivalent** à la consistance du score de fenêtre **de population à loi fixe** ; **pas** équivalent à la consistance de l'estimateur empirique ; et sa **forme qualitative ne suffit pas en tableau triangulaire** |
| Suppressions | « pour la plupart correctes, mais plusieurs sont des *remplacements* ou des *changements de cible* » |
| Modèles A1–A3, B1 | vérifiés pour la consistance du score au premier ordre et le criblage ; **pas** pour la correction de biais séparée en puissance ni pour l'inférence simultanée |

#### Les trois défauts qui mordent

1. **(EUR) qualitative est insuffisante en tableau triangulaire** — *le défaut le plus sérieux*.
   Contre-exemple des **damiers de plus en plus fins** : on peut avoir `ξ_j^♯ = ξ_j` presque partout
   pour chaque `n`, et pourtant `∫ M_{h_n} f_n → 1` alors que `∫ f_n → 1/2`. L'égalité p.p. ne
   contrôle rien uniformément quand la loi bouge avec `n`. La bonne condition est un **module
   quantitatif** :
   `Ω_n(h_n) = max_{j≤p_n} |I|⁻¹ ∫_I { M_{h_n} ξ_{j,n}(u) − ξ_{j,n}(u) } du → 0`,
   et pour le criblage **`Ω_n(h_n) = o(Δ_min,n)`**.
   *Conséquence pour nous* : partout où le registre dit « (EUR) est la condition spatiale exacte »,
   il faut écrire « (EUR) est exacte à loi fixe ; en tableau triangulaire c'est `Ω_n(h_n)=o(Δ_min)` ».

2. **Consistance de population ≠ consistance d'estimateur.** Même `Ω_n(h_n)→0` ne donne que
   l'approximation de la cible par l'indice asymptotique des lois fenêtrées. Il faut encore, en plus :
   `nα_n h_n → ∞`, `√(log(p_n n)/(nα_n h_n)) → 0`, et `B_n^tail → 0`. (EUR) ne les remplace pas.

3. **La borne inférieure positive sur γ : suppression réelle mais à branche.** Supprimable pour la
   projection de population et les bornes additives supérieures — **mais** les fibres d'indice nul
   exigent une branche séparée, et surtout **la formule du déficit réciproque n'est pas définie quand
   la fibre et l'indice projeté sont tous deux nuls**. Variance et inférence ne survivent pas
   inchangées. AB développe aussi le cas `γ ≥ 0` : si `γ*_U = 0` alors `D_U = ∅`, donc la cible
   elle-même est aveugle aux différences entre queues légères (contre-exemple `e^{-y}` vs `e^{-y²}`,
   les deux d'indice 0, gap nul) — ce n'est pas un échec d'estimation mais **une limite de la cible**.

#### Corrections de détail retenues
- « `a<1` » doit être **`0<a<1`** : l'exposant nul était admis par mégarde. Domaine juste :
  `0<a<1, 0<b<1−a`. La suppression de l'ancien `b<1/2` est en revanche **confirmée**.
- (CTI) exige que la convergence ait lieu **pour tout multiplicateur réel**, localement uniformément,
  pas seulement au pas canonique ; et une **branche Pareto exacte** (`a(t)≡0` éventuellement) sinon
  les queues exactement Pareto sont exclues de la classe.
- La suppression du support de fibre plein est valide **pour `D_U`**, pas pour l'ancien ensemble
  structurel ambiant — c'est un changement de cible, pas une suppression.
- La non-atomicité est supprimable pour la loi de la réponse (quantiles généralisés + PIT randomisé),
  mais **les ex æquo de queue supérieure demandent encore du soin**.
- La continuité n'est pas supprimée : elle est **remplacée** par la régularité essentielle supérieure.

#### Statut
La classe **n'est pas minimale** : deux de ses couches sont des conditions **sur les objets projetés**,
donc pas primitives, et l'une des équivalences annoncées n'est que conditionnelle. Description sûre :
**une classe stratifiée à deux interfaces de population exactes** (l'incrément canonique, et
l'enveloppe à loi fixe), le reste étant suffisant-et-utile mais pas nécessaire.

**Dettes ouvertes par AB** (ajoutées au cahier) :
- (AB-1) remplacer partout (EUR) qualitative par le module `Ω_n(h_n)=o(Δ_min,n)` en régime triangulaire ;
- (AB-2) écrire la branche `γ=0` (déficit réciproque non défini, `D_U=∅`, variance/inférence) ;
- (AB-3) ajouter à (CTI) les qualifications : premier ordre séparé, tous multiplicateurs, non-annulation, branche Pareto exacte ;
- (AB-4) corriger `a<1` en `0<a<1` partout ;
- (AB-5) ne plus annoncer « classe primitive minimale » mais « classe stratifiée à deux interfaces exactes ».

## VAGUE 11 — réparation de la classe, et la dette la plus ancienne

- **AC — réparation de la classe stratifiée** (fil de AA, qui l'a écrite)
  <https://chatgpt.com/c/6a818d98-9998-83ed-9352-e68707a6b12f>
  Les cinq défauts de AB lui sont renvoyés. Priorité explicite : la **couche spatiale en tableau
  triangulaire** sous forme utilisable, avec vérification du module `Ω_n(h_n)=o(Δ_min,n)` pour les
  quatre modèles **et rates explicites** — quelles plages de largeur de bande survivent, lesquelles
  sont perdues. Puis la branche `γ ≥ 0`. Puis : l'une des deux couches sur objets projetés peut-elle
  être remplacée par une condition **vraiment primitive** sur la surface `γ` et les noyaux de fibre,
  ou est-ce impossible ?

- **AD — représentation linéaire par ligne pour la calibration du maximum** (fil de W/X/Z/AB)
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  *C'est la dette la plus ancienne du cahier*, ouverte depuis la vague 3 et jamais close. Sans elle,
  le papier ne peut offrir que Bonferroni ou Benjamini–Yekutieli. Avec elle : lignes indépendantes,
  bootstrap par multiplicateurs gaussiens à ligne commune, borne CCK sur hyper-rectangles, valeur
  critique calibrée sur le **max**. L'obstacle est nommé dans la consigne : le score est bâti sur des
  **rangs empiriques** (une ligne déplace le rang de toutes les autres) et sur une **fenêtre dont
  l'appartenance est elle-même déterminée par les rangs** — l'influence d'une ligne n'est donc pas
  locale. Le télescopage d'édition de rangs déjà établi suggère que la perturbation est petite, mais
  il n'a jamais servi qu'à borner un pire cas, jamais à **extraire un terme linéaire**.
  Consigne : établir la représentation, ou prouver ce qui l'obstrue et quelle calibration plus faible
  reste atteignable.

### Résultat de AC — la couche spatiale réparée, et **trois impossibilités**

AC livre la réparation complète (8 sections) et va au-delà de ce qui était demandé.

#### 1. Le rayon effectif des fenêtres de rangs
Point neuf et important : les fenêtres **de rangs empiriques** n'ont pas le rayon `h_n` mais un rayon
de travail élargi par la fluctuation des rangs,
`r_n = h_n + 2√(h_n x_n / n) + 2x_n/n`, avec `x_n ≍ log(p_n n)`.
C'est à **ce** rayon qu'il faut imposer la condition spatiale, pas à `h_n`.

#### 2. Les plages de largeur de bande — vérifié indépendamment
Condition survivante : `h_n → 0` et `nα_n h_n / log(p_n n) → ∞`, soit **`0 < b < 1−a`**.
La restriction `b < 1/2` du manuscrit est **inutile** (confirmé pour la seconde fois) : dès que
`a < 1/2`, toute la bande `1/2 ≤ b < 1−a` survit. La grille affichée du papier
(`a ∈ {.25,…,.50}`, `b ∈ {0,…,.40}`) a pour somme maximale `0.50+0.40 = 0.90 < 1` : **toute cellule
à `b > 0` est couverte**, y compris les neuf réglages du bloc d'agrégation.

Trois régions sont perdues :
- **`b = 0`** (donc `h = 1/2`) : les fenêtres ne localisent jamais. AC donne les valeurs exactes du
  module, que **j'ai recalculées indépendamment** (avec `ξ_j(u) = 0.5 e^{-cu}`, `c = 1, 0.5, 0.8`,
  et `ε = 0.05` — valeur que je retrouve en reproduisant les gaps `0.186 / 0.107 / 0.158` du papier) :

  | modèle | `Ω(1/2)` (AC) | `Ω(1/2)` (recalcul) | gap | part du gap consommée |
  |---|---|---|---|---|
  | A1/B1 | 0.137713026295563 | **0.137713026295563** | 0.186 | **73.88 %** |
  | A2 | 0.081176584085900 | **0.081176584085900** | 0.107 | **75.65 %** |
  | A3 | 0.117500459573997 | **0.117500459573996** | 0.158 | **74.58 %** |

  Accord aux **15 chiffres** (dernier chiffre de A3 à l'arrondi près). La colonne `b = 0` n'est donc
  pas « techniquement non couverte » : ses fenêtres de population ne conservent qu'**environ un quart
  du gap nominal**. Cela *explique* enfin l'observation empirique du papier — « la largeur de bande
  dégénérée `b = 0` n'est jamais compétitive » — qui n'était jusqu'ici qu'un constat de simulation.
- **`a + b ≥ 1`** : `nα_n h_n` ne diverge plus, le compte d'extrêmes locaux reste borné.
- **gaps rétrécissants** : il faut `r_n = o(Δ_min,n)` ; si `Δ_min,n ≍ n^{-d}` alors **`d < b < 1−a`**,
  et si l'on veut aussi le gel horizontal de queue, `n^{-b} log n = o(Δ_min,n)`.

#### 3. Trois impossibilités (réponse à la question « peut-on rendre les couches primitives ? »)
- **7.1 — la couche spatiale, oui.** En substituant la définition de `ξ_{j,n}`, le module s'écrit
  entièrement en `γ_n` et les noyaux de fibre : `Ω_n^prim(r) = Ω_n(r)` **exactement**. Mais seulement
  en écrivant explicitement la même géométrie d'enveloppe conditionnelle : **aucun substitut ponctuel
  plus faible n'est valable uniformément en tableau triangulaire** — c'est le damier qui l'interdit.
- **7.2 — le second ordre projeté, non.** Contre-exemple explicite : `p = 1`, même loi de `U`, même
  noyau, `γ ≡ ξ > 0`. Modèle 1 : `F̄₁(e^t) = e^{-t/ξ}` (Pareto exact). Modèle 2 :
  `F̄₂(e^t) = exp{−t/ξ + t^{−1/2} sin t}` (dérivée logarithmique négative pour `t` grand, donc
  survie valide). Même indice de premier ordre `ξ` ; mais les incréments valent
  `2t^{−1/2} sin(a/2) cos(t+a/2) + o(t^{−1/2})` et **oscillent** : ceux de `a = π` et `a = π/2` n'ont
  aucun rapport asymptotique stable, donc aucun auxiliaire régulier éventuellement non nul ne peut
  convenir pour tout `a` réel. **`γ` et `K` seuls ne peuvent pas caractériser la 2RV projetée.**
  Ajouter la distribution pondérée du déficit ne suffit pas non plus.
- **7.3 / 7.4** — même verdict pour la régularité horizontale à seuil fini, et — plus surprenant —
  **même la projection d'indice nul n'est pas déterminée par `γ` et `K`**.

*Portée* : ces trois énoncés closent définitivement la question « existe-t-il une classe entièrement
primitive ». **Non**, et la raison est structurelle : le second ordre est une information sur la
**queue conditionnelle**, que la surface `γ` ne porte pas. Le nom correct est donc
**classe stratifiée minimale**, avec une couche primitive (la spatiale) et des couches de queue
irréductiblement non primitives. C'est un résultat négatif, mais c'est le bon type de résultat
négatif : il justifie l'architecture au lieu de la subir.

#### 4. La branche `γ ≥ 0`
Écrite : (TU+)–(EC+) quand `m_K > 0`, (TU0) quand `m_K = 0` ; **interdiction d'utiliser les déficits
réciproques** dans la branche nulle, remplacés par un certificat de Hill à seuil fini ; et si
`γ*_U = 0`, déclarer `D_U = ∅` explicitement plutôt que de laisser la formule diverger.

**Statut des dettes AB** : AB-1 ✅ (et renforcée : rayon `r_n`, pas `h_n`), AB-2 ✅, AB-3 ✅,
AB-4 ✅, AB-5 ✅ (rebaptisée « classe stratifiée minimale », avec preuve que « primitive » est
**impossible**).

- **AE — audit de AC** (fil de la vague 8)
  <https://chatgpt.com/c/6a813beb-8aec-83eb-8c21-97a0bb2a986f>
  Cibles : le **rayon gonflé** `r = h + 2√(hx/n) + 2x/n` est-il du bon ordre, et surtout *suffit-il* —
  ou l'appartenance à une fenêtre de rangs cache-t-elle une dépendance que l'argument de rayon
  escamote ? Le module au rayon **gonflé** (donc plus grand que `h`) reste-t-il `o(Δ_min)` sur toute
  la région annoncée ? Le témoin oscillant `exp{−t/ξ + t^{−1/2} sin t}` est-il une vraie fonction de
  survie, l'indice de premier ordre est-il vraiment inchangé, et l'oscillation exclut-elle **tout**
  auxiliaire ou seulement l'auxiliaire évident ? Enfin — la question qui décide du statut du
  résultat — s'agit-il d'une **impossibilité**, ou seulement de l'échec d'une caractérisation
  particulière ?

### Résultat de AD — la calibration du maximum : **la dette la plus ancienne est close, par la négative**

Ouverte à la vague 3, jamais tranchée. AD la tranche. Verdict : **pour le score glissant brut tel
qu'il est implémenté, le résultat demandé est faux sous la forme énoncée, et ne peut pas valoir
uniformément sur la classe actuelle.** Mais le diagnostic est précis, et il vient avec une issue.

#### Ce qui n'était PAS l'obstacle
Contre l'intuition qui bloquait la campagne depuis huit vagues : **les rangs empiriques ne sont pas
l'obstruction.** AD écrit la fonction d'influence de ligne du fonctionnel de rang de population,
explicitement (51) :
`φ_{j,n}(z) = |I_ε|⁻¹ ∫_{I_ε} (2h_n)⁻¹ [ 1{u−h_n < v_j ≤ u+h_n} χ_{j,u}(x)`
`  + μ_{j,u}(u+h_n){u+h_n − 1(v_j ≤ u+h_n)} − μ_{j,u}(u−h_n){u−h_n − 1(v_j ≤ u−h_n)} ] du`,
soit **l'influence de Hill ordinaire plus deux termes de bord** engendrés par le déplacement des
quantiles marginaux empiriques. Ces deux termes *sont* la correction de rang linéaire. La
perturbation non locale se linéarise donc parfaitement — ce n'est pas elle qui bloque.

#### Les quatre obstacles réels
1. **Impossibilité d'une influence uniformément bornée** (Proposition 1, avec preuve). Même dans le
   modèle nul Pareto exact, la variance de l'influence est d'ordre `1/α_n` et l'influence de Hill brute
   est non bornée, de norme sous-exponentielle d'ordre `1/α_n`. Une représentation à influence
   `O(1)` **contredirait la borne inférieure de fluctuation `1/√(nα_n)` déjà démontrée**. C'est une
   impossibilité interne, cohérente avec nos propres résultats — donc solide.
2. **Le télescopage d'édition de rangs ne suffit pas.** Il contrôle une **différence de premier
   remplacement** ; il ne contrôle **pas les composantes de Hoeffding dégénérées**. C'est exactement
   le point que le journal enregistrait comme non résolu — AD confirme que l'écart était réel et non
   un défaut de rédaction.
3. **La classe actuelle ne fournit pas** de développement de Bahadur uniforme aux quantiles
   intermédiaires, ni de borne à l'échelle du maximum sur le reste non linéaire.
4. **Le centrage de population a sa propre impossibilité** : la classe permet un biais à l'échelle de
   la population qui détruit l'inférence sur le maximum, indépendamment du reste.

#### L'issue constructive (§9) — one-step à échantillon scindé
Geler **toutes** les nuisances sur un échantillon d'entraînement indépendant `I₁`, et évaluer un score
d'influence one-step sur les lignes retenues `I₂` :
`Δ̂_j = Δ^{(1)}_{j,n} + n₂⁻¹ Σ_{i∈I₂} φ̂^{(1)}_{j,n}(Z_i)`.
Conditionnellement à l'entraînement, les `φ̂^{(1)}_{j,n}(Z_i)` sont **des vecteurs de lignes
indépendants par construction** : le reste plugin rang-Hill non analysé **disparaît**. Il ne subsiste
qu'un reste d'estimation de nuisance, contrôlé par deux conditions d'échelle explicites (45)–(46).
Sous elles, le multiplicateur à ligne commune a l'erreur CCK, avec `n₂` au lieu de `n`.
Coût : on change d'estimateur et on perd une fraction constante de l'échantillon. Une version
**bi-croisée** récupérerait l'essentiel de l'efficacité, mais sa preuve doit traiter la dépendance
créée par le recouvrement des ensembles d'entraînement.

#### La calibration exacte mais trop faible (§11)
Sous le nul **fort** `H₀ : Y ⊥ (U₁,…,U_p)`, permuter conjointement les `Y` (même permutation pour
toutes les coordonnées) et recalculer tout le score donne une **distribution de randomisation exacte**
pour n'importe quelle statistique de maximum — elle préserve les colonnes de rangs, toute la
dépendance croisée, les recouvrements de fenêtres, et le fait que les mêmes réponses extrêmes entrent
dans plusieurs coordonnées. **Mais le nul fort est bien plus fort que « `Δ_j = 0` pour tout `j` »** :
des gaps projetés nuls autorisent corps, échelle et second ordre dépendants des covariables, et même
des changements de loi marginale de queue à indice projeté plat. Le test de permutation n'est donc
**pas** un test du bon nul. À mentionner comme diagnostic, jamais comme calibration du résultat.

#### Ce que le papier doit dire
Le score glissant original n'a, sous la théorie présente, que des **garanties simultanées par
concentration** (Bonferroni / Benjamini–Yekutieli). C'est ce qu'il fait déjà — donc **le papier n'est
pas faux**, il est seulement à la limite exacte de ce qui est démontrable pour cet estimateur. Deux
options honnêtes : (i) garder le score brut et assumer BY, en citant l'impossibilité comme
justification ; (ii) ajouter la variante one-step scindée comme **procédure d'inférence** distincte du
criblage, avec le théorème CCK complet.

**Dette la plus ancienne : CLOSE** — non par la preuve espérée, mais par une impossibilité démontrée
plus une construction alternative qui atteint l'objectif. C'est une clôture légitime.

- **AF — audit de AD, avec une objection que j'apporte moi-même**
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  **L'objection centrale est de mon fait, pas du modèle.** L'impossibilité de AD prouve qu'aucune
  représentation n'a d'influence de ligne **bornée en `O(1)`**. Mais *la bornitude n'a jamais été
  requise par le théorème du maximum* : les bornes CCK sur hyper-rectangles ne demandent pas des
  sommants bornés — elles valent pour des enveloppes non bornées sous conditions de moments ou de
  norme sous-exponentielle, au prix d'un taux explicite faisant intervenir cette norme et une
  puissance de `log p`. Or le calcul de AD **donne** cette norme : `1/α_n`, quantité connue et
  contrôlée. L'impossibilité pourrait donc viser une exigence dont personne n'a besoin.
  Question posée : écrire explicitement la condition CCK avec enveloppe d'ordre `1/α_n` — une
  puissance de `log p` sur `n α_n^k` — et dire si elle tend vers zéro **dans la région de réglage
  admissible**. Si oui, **tout le verdict négatif tombe** et la construction à échantillon scindé
  devient inutile.
  Sous-questions : l'écart de Hoeffding dégénéré est-il réel, ou le télescopage peut-il être relevé
  au second ordre par découplage / Efron–Stein ? Et le one-step scindé a-t-il une propriété
  d'**orthogonalité de Neyman** rendant le reste de nuisance du second ordre — sans quoi il n'est pas
  meilleur que le plugin, et il faut le dire.

### Verdict de AE — audit de AC : « en partie solide, surévalué en quatre endroits »

| Claim de AC | Verdict de AE |
|---|---|
| (a) rayon gonflé `r = h + 2√(hx/n) + 2x/n` | **survit**, à condition que la preuve soit *réellement* recentrée sur l'ancre aléatoire `U_(r)j` (et non sur `r/(n+1)`), la condition spatiale posée sur un intervalle tronqué **élargi**, et l'événement négligeable « toutes les ancres restent dedans » ajouté. **L'appartenance aux fenêtres de rangs ne cache aucune dépendance** — objection levée. C'est exactement le changement qui supprime l'ancien déplacement DKW global. |
| (b) `0<b<1−a` | **correct comme région suffisante** pour les modèles A/B à gap fixe, et toute cellule affichée à `b>0` est couverte. **Mais pas « exact »** : ce n'est pas un énoncé de nécessité. |
| (c) valeurs numériques à `b=0` | **numériquement correctes** (accord avec mon propre recalcul). **Correction de terminologie** : les trois nombres sont les **pertes** de gap, pas les gaps retenus ; les fractions **retenues** sont 26.1 %, 24.4 %, 25.4 %. Et ce sont des **valeurs limites**, pas des valeurs exactes à `n` fini. |
| (d) `d<b<1−a` pour gaps rétrécissants | **FAUX pour le théorème complet.** Ces inégalités ne couvrent que la localisation et la faisabilité locale. L'erreur stochastique exige **en plus `2d < 1−a`**, et surtout le biais de queue projeté actuel, en `O(1/log n)`, **exclut tout gap polynomialement rétrécissant `n^{-d}`, `d>0`**, sauf à imposer une condition de second ordre ou de réduction de biais plus forte. |
| (e) couche spatiale « exactement primitive » | **vrai seulement pour la composante d'indice dominante.** La couche spatiale **à seuil fini** dans son entier n'est pas récupérable depuis l'enveloppe seule. Des substituts plus faibles existent néanmoins. |
| (f) témoin oscillant | **valide** : fonction de survie licite, indice de premier ordre bien inchangé, et **aucun auxiliaire régulier standard** ne peut convenir. *Mais* le comparateur Pareto exact n'est membre du second ordre que **par convention dégénérée** ; un comparateur véritablement de second ordre supprime ce défaut. |
| (g) « impossibilité » | **requalifié.** Ce qui est prouvé : *il n'existe pas de caractérisation nécessaire-et-suffisante de la 2RV projetée, ni de la régularité horizontale à seuil fini, qui soit fonction des seuls `(γ, K)`* — deux modèles de mêmes `(γ,K)` donnent des réponses différentes. Ce qui n'est **pas** prouvé : que la 2RV projetée soit impossible ; qu'elle ne puisse être imposée directement sur les quantiles projetés ; qu'elle ne se déduise pas d'hypothèses primitives **plus riches** faisant intervenir `c(x)`, `r(t,x)`, leurs dérivées et des couplages horizontaux ; qu'aucune condition suffisante utile n'existe ; qu'aucune procédure ne puisse la tester sur données complètes. |

> **Rectification à porter au registre** — j'avais écrit après AC « aucune classe entièrement
> primitive n'existe » et « trois impossibilités ». C'est **trop fort**. L'énoncé correct est une
> **non-détermination par descripteurs réduits** : `(γ, K)` ne suffisent pas. Des hypothèses
> primitives plus riches ne sont pas exclues, et restent une piste ouverte.

#### La dette la plus lourde ouverte par AE
**(AE-1) — le biais en `O(1/log n)` tue les gaps rétrécissants.** Tant que le biais de queue projeté
n'est qu'en `1/log n`, aucun `Δ_min ≍ n^{-d}` avec `d>0` n'est atteignable, quelle que soit la
largeur de bande. Le papier ne peut donc énoncer ses théorèmes qu'à **gaps fixes** — ce qu'il fait —
mais toute extension à gaps rétrécissants exige d'abord un **certificat de biais en `o(n^{-d})`**,
donc une correction de biais séparée en puissance. C'est la porte d'entrée obligée vers l'inférence.

**(AE-2)** recentrage explicite sur `U_(r)j` + intervalle élargi + événement d'ancre, à écrire.
**(AE-3)** remplacer le comparateur Pareto exact par un comparateur de second ordre non dégénéré.
**(AE-4)** dire « région suffisante », jamais « exacte », pour `0<b<1−a`.
**(AE-5)** dire « pertes de gap » et non « gaps retenus » pour les trois nombres à `b=0`.

## VAGUE 12 — la porte d'entrée : le biais, et la primitivité

- **AG — certificat de biais polynomial** (fil de AA/AC)
  <https://chatgpt.com/c/6a818d98-9998-83ed-9352-e68707a6b12f>
  Attaque frontale de (AE-1), *la question la plus importante du programme entier*. Le biais de queue
  projeté n'est qu'en `O(1/log n)`, ce qui **exclut tout gap `n^{-d}`, `d>0`, quelle que soit la
  largeur de bande** : tous les théorèmes du papier sont donc confinés aux gaps fixes, et aucun
  réglage n'en sort. Consigne : construire un certificat en `o(n^{-d})` — score à biais réduit,
  seuils séparés en puissance ou mieux ; prouver l'ordre atteint ; dire exactement quelle hypothèse
  de second ordre il exige ; donner la région admissible en `(a,b,d)` ; **chiffrer le coût en
  variance** (la réduction de biais la gonfle toujours) et dire si la région est **non vide**. Si un
  certificat polynomial est impossible dans la classe de second ordre projetée, le prouver et
  identifier **le meilleur taux de gap atteignable**.

- **AH — jusqu'où peut-on être primitif ?** (fil de AE)
  <https://chatgpt.com/c/6a813beb-8aec-83eb-8c21-97a0bb2a986f>
  AE a distingué non-détermination par `(γ,K)` et impossibilité, en laissant explicitement ouvert
  l'enrichissement du descripteur. Consigne : trouver **le plus petit enrichissement** — `c(x)`,
  `r(t,x)`, dérivées, couplages horizontaux — qui **détermine** la 2RV projetée ; ou prouver
  qu'aucun descripteur bâti sur des quantités de queue conditionnelle n'admet de caractérisation
  nécessaire-et-suffisante, et que seules des conditions suffisantes sont disponibles. Exiger : quelle
  des deux est prouvée ; non-vacuité par une classe non triviale ; vérification sur les modèles de
  simulation ; et si le même enrichissement règle la régularité horizontale à seuil fini.
  *Enjeu* : décide si l'annexe primitive du papier peut être **réellement** primitive, ou si elle
  portera à jamais des conditions énoncées sur les objets projetés.

### Verdict de AF — **le verdict négatif de AD est retiré** ; mon objection était fondée

AF ouvre par : *« the previous verdict was too negative. The sharp objection is correct. »*
L'impossibilité de AD n'établissait que ceci : aucune représentation ne peut employer une influence
bornée **uniformément par une constante `O(1)`**. C'est valide contre la formulation littérale à
influence bornée, **mais ce n'est pas une obstruction à une approximation gaussienne ou par
multiplicateurs de type CCK**, qui n'a jamais exigé de sommants bornés.

#### Le calcul, fait explicitement
Après standardisation, l'influence de Hill a une enveloppe **sous-exponentielle d'ordre `α_n^{-1/2}`**
(et non `1/α_n` comme je l'avais supposé — la standardisation en absorbe la moitié). La condition CCK
qui en résulte est une condition sur la **taille d'échantillon de queue effective `nα_n`** :

- taux classique (sixième racine) : `δ = O[{log⁷n / n^{1−a}}^{1/6}] = o(1) ⟺ a < 1`
- taux amélioré (quatrième racine) : `δ = O[{log⁵(p_n n) / (nα_n)}^{1/4}] = o(1) ⟺ a < 1`

**La largeur de bande n'entre pas** dans cette approximation gaussienne à sommants linéaires ; elle ne
revient que dans la faisabilité locale `nα_n h_n ≫ log(p_n n)`, soit `a+b<1`. Donc **dans toute la
région de criblage du manuscrit `0<a<1, b>0, a+b<1`, l'étape CCK n'est pas l'étape limitante** — il y
a « ample room ». Pour une dimension véritablement sous-exponentielle `log p_n ≍ n^κ`, les bornes
informatives sont `κ < (1−a)/7` (ancien théorème) et **`κ < (1−a)/5`** (théorème amélioré).

#### Registre corrigé
| Étape | Statut |
|---|---|
| **Approximation CCK** | **viable** — rétractation explicite de « l'enveloppe croissante obstrue la calibration » |
| **Linéarisation du score brut** | **ouverte et plausible, non réfutée** |
| **Reste de Hoeffding dégénéré** | **vrai trou de preuve**, mais **aucune impossibilité établie** — un argument d'Efron–Stein au second ordre ou de découplage peut réussir. Repère suffisant naturel : **`a+2b<1`**, que *tous* les réglages du manuscrit satisfont. Ce qui manque est la preuve des différences mixtes et de la queue exponentielle, rien d'autre. |
| **One-step à échantillon scindé** | **rétracté comme « rigoureux et atteignable »** : ce n'est une alternative qu'après avoir prouvé une orthogonalité de Neyman cohérente et des taux de nuisance uniformes. **Pas démontrablement préférable au plugin** en l'état. |
| **Honnêteté sur les gaps limites** | **toujours impossible** pour le score non corrigé sur la classe actuelle — c'est (AE-1), le verrou du biais. |

> **Ce que cela change.** La dette de la vague 3 n'est **pas** close par la négative, comme je l'avais
> écrit après AD. Elle est **rouverte et rétrécie** : il ne reste qu'un seul énoncé à prouver — la
> négligeabilité du reste de Hoeffding dégénéré, sous `a+2b<1` — et non plus toute une théorie. Le
> reste de la chaîne (influence explicite, correction de rang, étape CCK) est acquis. C'est le
> meilleur état où cette question ait jamais été.

**Dette (AF-1)** : prouver la négligeabilité des différences mixtes / du reste dégénéré par
Efron–Stein d'ordre deux ou découplage, avec queue exponentielle, sous `a+2b<1`. **C'est désormais le
seul verrou entre le papier et l'inférence simultanée calibrée sur le maximum.**

- **AI — le dernier verrou : le reste de Hoeffding dégénéré** (fil de AD/AF)
  <https://chatgpt.com/c/6a81695f-7884-83eb-bdd0-1064620fd1b0>
  Il ne reste **qu'un seul énoncé non prouvé** entre le papier et l'inférence simultanée calibrée sur
  le maximum. Consigne : montrer que l'erreur de projection de Hájek — la somme de toutes les
  composantes de Hoeffding d'ordre ≥ 2 — est négligeable **uniformément sur les `p` coordonnées**, à
  l'échelle exigée par l'approximation gaussienne, dans la région repère `a+2b<1` que **tous** les
  réglages du manuscrit satisfont.
  *Indication que j'ajoute* : une différence seconde mixte doit être **bien plus petite** qu'une
  différence première, car remplacer une ligne ne déplace le rang de chaque autre ligne que d'**au
  plus une position** — le second remplacement voit donc une fenêtre quasi inchangée. C'est
  précisément la structure qu'un argument de télescopage **jette**, et c'est là qu'il faut creuser.
  Exigences : décomposition complète, borne sur la différence mixte, étape Efron–Stein ou différences
  bornées, **queue exponentielle** (une borne de variance seule ne suffit pas pour un maximum sur
  `p`), union sur `p`, condition finale en `(a,b)`. Si une étape échoue vraiment : isoler le plus
  petit énoncé faux, le contre-exemplifier, et donner la vraie région d'exposants.

### Résultat de AG — le verrou du biais : **impossible dans la classe actuelle, possible dans une classe plus forte**

Réponse en deux temps, et les deux comptent.

#### 1. Impossibilité dans la classe actuelle — démontrée
**Il n'existe aucun certificat de biais polynomial sous la classe de second ordre projetée du
manuscrit.** Raison structurelle : dans le régime de projection diffuse, le paramètre de second ordre
est `ρ = 0` et l'auxiliaire est d'ordre `1/log t` ; la 2RV ordinaire ne dit rien de plus que « le reste
est `o(A(t))` » — elle **ne fournit aucun taux pour ce petit-o**. AG construit un reste explicite,
**à décroissance arbitrairement lente**, à l'intérieur de la même classe 2RV : il converge plus
lentement que **toute** puissance `t^{-β}`. Les corrections logarithmiques finies ne réparent rien, et
**aucun estimateur** ne peut atteindre un taux polynomial uniforme sur la classe non quantifiée.
Ceci confirme (AE-1) et lui donne son statut définitif : ce n'est pas un défaut de preuve, c'est une
**limite de la classe**.

#### 2. La classe renforcée qui, elle, le permet — et la construction optimale
Il faut une hypothèse de **troisième ordre quantitative** (notée `QLH_β`) : tout le biais de
projection non polynomial doit appartenir à une **forme de nuisance unidimensionnelle connue**, avec
un reste polynomial.

Sous cette classe, **la séparation en puissance n'est pas la bonne construction**. La meilleure est
l'**extrapolation à rapport fixe** : deux fractions de queue différant d'un **facteur constant `λ`**,
avec des poids d'ordre `log(1/α)`. Elle donne :
- biais de queue `O{α^β log(1/α)}` ;
- **inflation de variance seulement `O{log²(1/α)}`**, soit une inflation d'écart-type d'ordre
  `log(1/α)` — le prix, incompressible, d'identifier faiblement le coefficient d'un biais en `1/log t`
  à partir de deux seuils voisins.

**Constante optimale, que j'ai vérifiée indépendamment.** L'inflation vaut
`v ∼ [(λ−1)/(log λ)²] · log²(1/α)`, minimisée en :

| quantité | AG | recalcul |
|---|---|---|
| `log λ` | 1.59362 | **1.59362** |
| `λ` | 4.92155 | **4.92155** |
| valeur minimale | 1.54414 | **1.54414** |

Accord exact. (Le choix de `λ` ne change que des constantes, jamais la région d'exposants.)

#### 3. La région admissible `(a,b,d)` — et elle est non vide
Avec `α_n = n^{-a}`, `h_n = n^{-b}/2`, `p_n ≤ n^C`, `Δ_min,n ≍ n^{-d}` :

| contrainte | origine |
|---|---|
| `a + b < 1` | faisabilité locale |
| `2d < 1 − a` | erreur stochastique |
| `d < aβ` | certificat de biais de queue |
| `d < b` | biais spatial et horizontal |

Un `a` admissible existe ssi `d/β < a < 1−2d`, donc **ssi**
> **`d < β / (2β+1)`**

et l'on peut ensuite choisir `d < b < 1−a`. **La région est genuinement non vide** : pour `β = 1`,
tout `d < 1/3` est atteignable.

#### 4. Le prix à payer, et il est réel
**Les quatre modèles de simulation A1–A3 et B1 sont exclus de la classe renforcée**, tels qu'ils sont
analysés aujourd'hui. Leur masse de déficit en coin gaussien engendre, après transformée de Laplace
et inversion de quantile, des facteurs en `√(log log t)`, `log log t / log t` et puissances inverses de
`log log t` ; la correction par séparation en puissance y laisse des biais en
`1/{log t (log log t)^{3/2}}` ou `1/{log t (log log t)²}`, **à coefficients génériquement non nuls**.
Or `n^d / {log n (log log n)^ν} → ∞` pour tout `d>0` : ces résidus ne sont **jamais** `o(n^{-d})`.
**Conclusion honnête : A1–A3 et B1 ne supportent aujourd'hui aucun théorème à gap polynomialement
rétrécissant.** Leur condition de gap reste inverse-logarithmique.

> **Ce que le papier doit en faire.** Le théorème à gaps rétrécissants existe, il est propre, sa
> région est non vide et sa constante optimale est calculée — mais il vit dans une classe que les
> modèles de simulation actuels **ne satisfont pas**. Deux options : (i) énoncer le théorème sous
> `QLH_β` et dire franchement qu'aucun des quatre modèles ne l'illustre ; (ii) **ajouter un cinquième
> modèle de simulation**, construit pour satisfaire `QLH_β`, afin que le théorème ait un témoin
> numérique. La seconde est nettement préférable, et c'est un travail borné.

**Dette (AG-1)** : construire un modèle de simulation satisfaisant `QLH_β` avec `β` explicite, et
l'ajouter à la campagne numérique pour illustrer le théorème à gaps rétrécissants.

### Résultat de AH — **l'annexe primitive peut être réellement primitive**

*« The broad impossibility alternative is false. »* AH tranche dans le sens positif, et c'est la bonne
nouvelle qui manquait depuis AC.

#### Le quotient minimal complet
Une fois l'indice projeté connu, **toute la famille primitive de queue conditionnelle n'entre dans la
queue projetée que par une seule intégrale scalaire de fibre**. Factorisation exacte :
`F̄_{j,u}(e^t) = e^{-t/ξ_j(u)} · A_{j,u}(t)`, avec
`A_{j,u}(t) = ∫ e^{-tD(v)} c{ι_j(u,v)} {1+r(t, ι_j(u,v))} K_j(u,dv)`, où `D(v) = 1/γ(x) − 1/ξ`.
Cette intégrale scalaire **est le quotient complet minimal** de la famille primitive, à recodage
bijectif près. C'est elle qui explique la non-détermination de AC : `(γ, K)` ne déterminent pas
`A_{j,u}`, parce qu'ils ignorent `c` et `r`.

#### Le plus petit enrichissement en amont — et il est vérifiable
Deux raffinements de **quantités déjà présentes** dans la preuve de transfert du manuscrit :
`M_{j,u}(t)/B_{j,u}(t) = κ_{j,u}/t + o(t^{-1})` et `sup_x |∂_t r(t,x)| = o(t^{-1})`.
Version **géométrique directement vérifiable** : la **masse de déficit réciproque pondérée par `c`
doit être à variation régulière en zéro, d'exposant `κ_{j,u} > 0`**, soit `H_{j,u} ∈ RV_{κ_{j,u}}(0)`.

Conséquence — et c'est un point à retenir : le paramètre de second ordre projeté qui en résulte est
**`ρ = 0`**, *pas* `−κ_{j,u}`, avec auxiliaire de quantile de queue
`A^Q_{j,u}(s) ∼ −ξ_j(u) κ_{j,u} / log s`.
Cohérent avec AG : on est bien dans la branche diffuse `ρ=0` à auxiliaire en `1/log`, donc **cet
enrichissement ne donne pas à lui seul le biais polynomial** — il détermine la 2RV, rien de plus.
Les deux résultats se complètent sans se contredire.

#### Les trois niveaux, distingués proprement
| niveau | statut |
|---|---|
| condition additive de second ordre sur `A_{j,u}(t)` | **exactement nécessaire et suffisante** — mais c'est la condition *projetée*, en forme d'intégrale de fibre |
| paire primitive `tM/B → κ`, `tR₁(t) → 0` | **suffisante**, et dans la branche à reste rapide et maximiseur non atomique, c'est la **caractérisation de von Mises** de l'expansion `ρ=0`. `H_{j,u} ∈ RV_κ(0)` y est **équivalent** au comportement de Laplace requis. |
| sur la classe primitive non restreinte | **non nécessaire** : les lois Pareto exactes forment une branche dégénérée, les maximiseurs atomiques ont d'autres corrections dominantes, et un reste conditionnel peut dominer le déficit |

#### Ce que le papier doit faire — division éditoriale nette
Remplacement proposé pour l'annexe : **(P2-SO)** `R₀(t)→0`, `tR₁(t)→0` ; **(P3-RV)**
`H_{j,u}(z) = a_{j,u} z^{κ} L_{j,u}(1/z){1+o(1)}`, `κ>0` ; **(P4)** couplage horizontal de
`1/γ`, `log c`, `log(1+r)`. L'annexe peut alors **prouver** au lieu de supposer
`Q_j(·,u) ∈ 2RV_{ξ_j(u), 0}(−ξ_j(u)κ_{j,u}/log(·))`, uniformément sur les fibres utilisées par le score.

**Et la vérification est complète pour A1–A3 et B1** : le théorème du coin gaussien fournit (P3-RV),
les dérivées de leur reste satisfont (P2-SO), leurs couplages existants donnent (P4). **La classe
contient les quatre modèles de simulation** — contrairement à la classe `QLH_β` de AG.

**Division retenue** : *théorème principal* sur les conditions projetées directes les plus faibles ;
*annexe primitive* comme **théorème de transfert rigoureux** pour une sous-classe naturelle large
contenant les modèles. Ce n'est pas une équivalence universelle, et il ne faut pas le prétendre.

#### La branche atomique et l'horizontal
- §7 : le cas d'un **atome au maximum** est *génuinement différent* — corrections dominantes autres.
  À traiter comme branche séparée, jamais à absorber.
- §8 : **la régularité horizontale à seuil fini ne découle pas** de l'enrichissement vertical. La
  condition de couplage horizontal (P4) du manuscrit **reste nécessaire comme dispositif**. AC avait
  raison sur ce point ; l'enrichissement ne le règle pas.

## VAGUE 13 — confronter AG et AH l'un à l'autre

- **AJ — audit de AG, et le régime intermédiaire** (fil de AG)
  Objection que j'apporte : le témoin « à décroissance arbitrairement lente » de AG est construit **au
  niveau du quantile projeté**. Or AH vient d'établir que la queue projetée ne dépend de la famille
  primitive **que par une intégrale scalaire de fibre**, et que si la masse de déficit réciproque
  pondérée par `c` est à variation régulière en zéro d'exposant `κ>0`, alors la 2RV projetée suit avec
  auxiliaire `∼ −ξκ/log s` — classe qui **contient les quatre modèles**. Question : le témoin de AG
  est-il **atteignable par projection** depuis une famille primitive légitime satisfaisant cette
  condition, ou n'existe-t-il qu'en postulant la pathologie directement au niveau projeté ? Si non
  atteignable, **l'impossibilité ne mord pas sur la classe que le papier utiliserait**.
  Puis la piste que cela ouvre, et que personne n'a explorée : dans la classe primitive, le biais
  dominant est d'ordre **exactement** `1/log t`, à coefficient **identifié**. Peut-on l'estimer et le
  soustraire pour obtenir `1/log²`, et **itérer** jusqu'à `1/log^J` ? Si oui : estimateur, ordre du
  biais après `J` étapes, inflation de variance en fonction de `J`, et condition de gap résultante.
  **Ce régime intermédiaire — entre gaps fixes et gaps polynomiaux — est peut-être la vraie réponse
  pour ce papier.**

- **AK — audit de AH** (fil de AH)
  Cibles : l'étape de Laplace/Tauber donne-t-elle bien `ρ = 0` et non `−κ` (le signe décide de toute
  la théorie de second ordre en aval), et où exactement `κ` survit-il — seulement dans la constante
  auxiliaire ? La caractérisation de von Mises vaut-elle **dans les deux sens** dans la branche
  annoncée, et le facteur à variation lente ne la casse-t-il pas ?
  **Et une question de cohérence entre les deux sessions** : AG trouve pour les mêmes quatre modèles
  des facteurs en `√(log log t)`, `log log t / log t`, puissances inverses de `log log t`, et en conclut
  que leur biais n'est jamais polynomial ; AH conclut que ces mêmes modèles satisfont sa condition de
  variation régulière avec `κ>0`. Les deux ne peuvent tenir que si le facteur à variation lente
  **absorbe toute la structure en `log log`**. Confirmer ou réfuter explicitement, en calculant la
  masse de déficit de A1 et en exhibant son `κ` et son facteur lent.
  Enfin : correction dominante explicite dans la branche à **atome au maximum**, et le mélange
  atome + coin continu est-il couvert par **aucune** des deux branches ?

### Verdict de AK — le transfert primitif survit ; **AG et AH sont compatibles**, et la raison est nette

#### 1. `ρ = 0` confirmé, et l'endroit exact où `κ` survit
Si la masse de déficit réciproque pondérée par `c` est à variation régulière en zéro d'exposant
`κ>0`, alors la survie projetée **et** le quantile supérieur projeté ont bien le paramètre de second
ordre **`ρ = 0`**, pas `−κ`. Le nombre `−κ` est l'**indice de variation régulière du facteur de
Laplace `B(T)` en la variable log-seuil `T = log y`** : `B(T) ∈ RV_{−κ}`. Ce n'est pas le paramètre de
second ordre pour les changements **multiplicatifs** `y ↦ λy`. La codimension effective `κ` ne
survit que dans l'**auxiliaire dominante** : `A^S(y) ∼ −κ/log y`, `A^Q(s) ∼ −ξκ/log s`.

#### 2. La réconciliation AG / AH — trois variables différentes
C'était la question de cohérence. Réponse : **aucune contradiction**, les deux résultats parlent de
variables différentes.

| objet | variable | comportement dominant |
|---|---|---|
| `H(δ)` | `δ ↓ 0` | `δ^κ L(1/δ)` — **polynomial** |
| `B(T)` | `T → ∞` | `T^{−κ} L(T)` |
| `F̄(y)` | `y → ∞` | `y^{−1/ξ} (log y)^{−κ} L(log y)` — **logarithmique** |

`κ > 0` signifie un comportement polynomial **dans la variable de déficit `δ`**, ce qui se traduit en
comportement **logarithmique** dans la variable de réponse `y`. Les facteurs en `√(log log y)` et
puissances de `log log y` relevés par AG **sont réels** — ce sont exactement le facteur à variation
lente évalué en `T = log y` — et ils sont **subordonnés** à l'auxiliaire dominante en `1/log s`.
Les deux sessions ont raison. **AG-1 et la vérification de AH tiennent simultanément.**

#### 3. Calcul explicite pour A1 — vérifié indépendamment
Pour A1, `γ(u) = ½exp(−Σ₁⁴u_a)`, `U_a = Φ(Z_a)`, bloc actif AR(1) gaussien à `ρ = 1/4`. En
conditionnant sur `U_j = u`, le déficit réciproque exact est `D = 2e^u(e^S − 1)`, linéaire au coin,
donc théorème du coin gaussien en dimension libre `d = 3` avec `κ_j = 1ᵀΩ_j1`, où **`Ω_j` est
l'inverse de la covariance _conditionnelle_** de `Z_{−j}` sachant `Z_j` (et non de la covariance
marginale — c'est le point subtil).

**J'ai recalculé `κ_j` moi-même** :

| coordonnée | AK | recalcul |
|---|---|---|
| `j = 1, 4` (extrémités) | 34/15 | **34/15 = 2.266667** |
| `j = 2, 3` (intérieures) | 41/15 | **41/15 = 2.733333** |

Accord exact. (Avec la covariance *marginale* on obtiendrait 11/5 = 2.2 — donc le conditionnement est
bien traité.) Survies projetées résultantes :
`F̄_{1,u}(y) ≍ y^{−1/ξ₁(u)} (log y)^{−34/15} (log log y)^{−11/30} exp{−(4/15) z_u √(2 log log y)}`
`F̄_{2,u}(y) ≍ y^{−1/ξ₂(u)} (log y)^{−41/15} (log log y)^{−2/15} exp{−(8/15) z_u √(2 log log y)}`

#### 4. Les rétrécissements imposés
- **La réciproque de von Mises était trop large.** Elle est bien à double sens dans la branche
  Pareto conditionnelle exacte à coin non atomique. **Avec un reste conditionnel non nul**, le sens
  direct tient sous la condition de dérivée annoncée, mais **la réciproque exige une condition
  explicite de non-annulation**. À corriger dans l'énoncé.
- Le facteur à variation lente **ne peut pas changer `ρ`** (§5) — vérifié.

#### 5. La branche atomique — trois cas, pas deux
| cas | `H` près de zéro | correction projetée dominante |
|---|---|---|
| coin continu | `H(0)=0`, `H(z) ∼ a z^κ L(1/z)` | `A^Q(s) ∼ −ξκ/log s` |
| **atome + coin continu** | `H(z) = h₀ + a z^β L(1/z)` | `A^Q(s) ≍ −(log s)^{−β−1} L(log s)` |
| **atome + trou spectral** | `H(z) = h₀` sur `(0,δ₀)` | `A^Q(s) ≍ s^{−ξδ₀}`, soit **`ρ = −ξδ₀ < 0`** |

Le mélange atome + coin continu **n'est pas** couvert par la condition sans atome, mais **l'est** en
appliquant la variation régulière à `H(z) − H(0)`. Un mélange atomique **non structuré** — aucune
forme imposée à `H(z)−H(0)` — n'a **aucune correction de second ordre identifiée** et relève du seul
premier ordre. L'ancienne condition de doublement du manuscrit autorisait les atomes mais
**n'identifiait pas le taux résiduel**, donc ne pouvait pas identifier l'auxiliaire — c'est
précisément ce que la nouvelle formulation répare.

> **Remarque de portée** : le troisième cas est le seul de toute la campagne qui produise un
> **`ρ < 0` strict**. C'est donc là — et seulement là — que la théorie de second ordre classique à
> `ρ` négatif s'applique, et donc potentiellement là que le biais polynomial de AG devient
> atteignable. **À croiser avec AJ.**

- **AL — la branche à trou spectral : la frontière entre biais logarithmique et polynomial**
  (fil de AH/AK) — *piste ouverte par mon croisement de AK et AG*
  Toute la campagne n'a produit que des branches à `ρ = 0` et auxiliaire en `1/log` — et c'est
  exactement là que AG a prouvé qu'**aucun** certificat de biais polynomial n'existe, ce qui confine
  tous les théorèmes du papier aux gaps fixes. **Mais l'impossibilité de AG porte sur la branche
  `ρ = 0`.** Or la troisième ligne du tableau de AK — atome au maximum **avec trou spectral** — donne
  `ρ = −ξδ₀ < 0` **strict**, avec une auxiliaire qui est une vraie puissance. La théorie classique de
  second ordre s'y applique donc, et le certificat polynomial devrait y être disponible.
  Consigne : dériver le développement projeté exact dans cette branche (`ρ`, auxiliaire, constante) ;
  donner le biais à seuil fini du score de Hill local ; déterminer si l'extrapolation à rapport fixe
  y atteint un certificat **polynomial** ; donner la région admissible en `(a,b,d)` ; dire **quelle
  taille minimale d'atome** est requise ; et surtout — **un trou `δ₀ → 0` avec `n` interpole-t-il
  entre les deux régimes ?** Cela donnerait **la frontière exacte entre biais logarithmique et biais
  polynomial**, donc la condition structurelle précise sous laquelle le criblage à gaps rétrécissants
  est possible. C'est aujourd'hui la principale question ouverte du programme.

## ★★★ RÉSULTAT DE AI — **LE DERNIER VERROU EST OUVERT : LA PREUVE MARCHE** ★★★

*La dette de la vague 3 — la représentation linéaire par ligne pour la calibration du maximum — est
**close par la preuve**, treize vagues après son ouverture.*

#### D'abord : mon indication était fausse
Je suggérais que toute différence seconde mixte gagne un facteur `(nh)^{-1}` sur une différence
première, parce qu'un remplacement de ligne ne déplace chaque rang que d'une position. AI identifie
cet énoncé comme **le plus petit énoncé faux** : il est faux pour Hill, parce que **deux lignes
peuvent déplacer conjointement le seuil local**. L'intuition était « directionnellement juste mais
fausse trajectoire par trajectoire ». *À retenir : aucune borne déterministe de différence mixte ne
peut marcher ici.*

#### La preuve qui marche est **probabiliste**, en cinq temps
1. **Linéariser conditionnellement** chaque statistique de Hill locale à son quantile intermédiaire
   **déterministe** (et non au quantile empirique) ;
2. prouver que son **reste de Hájek local exact** a une échelle de queue **`k^{-3/4}`**, où
   `k ≍ nαh` (Lemme 1 : reste de Hill–Bahadur local ; Lemme 2 : la correction de projection est
   `O(k^{-1})` ; Lemme 3 : queue **exponentielle** du reste dégénéré local exact) ;
3. **moyenner** ces restes dégénérés locaux en réutilisant **le coloriage de fenêtres disjointes**
   déjà établi dans la campagne — la vieille construction par classes de résidus resert ici ;
4. linéariser la **statistique de rang linéaire** restante et borner son reste d'assignation de rang ;
5. **combiner les deux bornes exponentielles**.

#### La représentation obtenue
`φ_{j,n}(Z_i) = n · b_{j,n}(U_{ij}) · h_{m,k}[ −log{1 − F_{Y|j}(Y_i | U_{ij})} ]`
et **la somme de toutes les composantes de Hoeffding d'ordre ≥ 2 est uniformément négligeable à
l'échelle du maximum gaussien** (Théorème 5).

#### La région d'exposants — meilleure que le repère
La preuve délivre
> **`0 < a < 1`, `0 < b < 1/2`, `a + b < 1`**

Le repère `a+2b<1` que AF avait proposé est **strictement plus fort**, donc suffisant. **Les neuf
réglages du manuscrit** (`a ∈ {0.30,0.35,0.40}`, `b ∈ {0.10,0.15,0.20}`) donnent `max(a+2b) = 0.80 < 1`
— **tous couverts**, et ce sont exactement les neuf réglages du criblage agrégé.

> ⚠️ **Note pour le registre** : la contrainte **`b < 1/2` réapparaît ici**, alors que les vagues 8–11
> l'avaient éliminée comme artefact de la route DKW abandonnée. Ce n'est **pas** une contradiction :
> c'est un `b<1/2` **d'origine différente**, issu du contrôle du reste de Hájek, pas de DKW. Le
> registre doit distinguer les deux. Pour le criblage seul, `0<b<1−a` reste valide ; pour
> **l'inférence calibrée sur le maximum**, il faut en plus `b<1/2`.

#### La conséquence, enfin atteinte
Avec `φ_{j,n} = ϕ_{0,n} − ϕ_{j,n}` l'influence de gap :
`Δ̂_j − Δ_{j,n} = n⁻¹ Σᵢ φ_{j,n}(Z_i) + o_P(1/√(nα log p))` **uniformément en `j`**.
La statistique à **multiplicateur gaussien de ligne commune** — un seul `e_i` partagé par toutes les
coordonnées, ce qui **préserve intégralement la dépendance croisée** — a pour erreur
`O[{log⁵(pn)/(nα)}^{1/4}]` plus l'erreur d'estimation de covariance.

> **Et surtout : « No split-sample replacement is required to eliminate the raw score's degenerate
> Hoeffding component. »** Le score **brut, tel qu'implémenté dans le papier**, est calibrable sur le
> maximum. Ni scission d'échantillon, ni changement d'estimateur, ni perte de sample.

#### Bilan de la chaîne d'inférence simultanée
| maillon | statut |
|---|---|
| fonction d'influence explicite (avec correction de rang) | ✅ AD |
| enveloppe sous-exponentielle `α^{-1/2}`, étape CCK viable | ✅ AF |
| **reste de Hoeffding dégénéré négligeable** | ✅ **AI — le dernier verrou** |
| multiplicateur à ligne commune, dépendance croisée préservée | ✅ AI |
| région d'exposants couvrant les neuf réglages du papier | ✅ AI |

**La chaîne est complète.** Le papier peut passer de Benjamini–Yekutieli à une **valeur critique
calibrée sur le maximum**, sans changer d'estimateur. C'est le résultat le plus important de la
campagne.

**Reste à faire** : audit adversarial de AI (échelle `k^{-3/4}`, réutilisation du coloriage, queue
exponentielle du reste d'assignation de rang, conditions d'estimation d'influence et de cohérence de
covariance).

- **AM — audit hostile de AI** (fil de AD/AF/AI)
  Le résultat le plus lourd de la campagne doit subir l'audit le plus dur. Cibles :
  **(1)** d'où vient l'exposant `k^{-3/4}` — vrai ordre ou artefact d'un découpage commode ? Un
  exposant plus faible suffirait-il après moyennage et union sur `p` ? Exiger le calcul.
  **(2) La réutilisation du coloriage** — il avait été introduit pour montrer que le score intégré a
  l'échelle `nα` et non `nαh` ; moyenner des restes **dégénérés** est une autre tâche. *Objection que
  j'apporte* : dans une classe de couleur les fenêtres sont disjointes **en positions de rang**, mais
  **les lignes sous-jacentes sont partagées entre coordonnées et les valeurs de réponse sont les
  mêmes** — les restes peuvent donc être dépendants même à fenêtres disjointes. L'indépendance
  utilisée est-elle réelle, ou seulement **intra-coordonnée** ?
  **(3)** la linéarisation conditionnelle se fait au quantile intermédiaire **déterministe**, mais la
  statistique implémentée utilise l'**empirique** : chiffrer l'erreur de transfert, vérifier qu'elle
  n'est pas absorbée par hypothèse.
  **(4)** borne exponentielle explicite du reste d'assignation de rang, **et le comportement au bord
  de troncature**, où une fenêtre est rognée et le nombre de rangs contributifs n'est pas le nominal.
  **(5) le retour de `b<1/2`** — confirmer que son origine diffère bien de la route DKW éliminée, ou
  découvrir qu'il tombe aussi.
  **(6)** les conditions de cohérence de covariance et d'estimation d'influence, auxquelles la preuve
  renvoie, sont-elles **elles-mêmes atteignables** pour cet estimateur, ou cachent-elles la
  difficulté restante ?

## LIVRABLE — `dossier-preuves.pdf` (31 pages, compile sans erreur)

Document pédagogique rassemblant les résultats de la campagne. Structure par résultat :
**énoncé → « L'idée » (en langage ordinaire) → démonstration → encadré « Ce que cela apporte à
l'article »**. Étiquettes de statut honnêtes (A / A⁻ / V / S / I) pour ne pas laisser passer un
résultat non audité dans l'article.

Chapitres : (1) la cible de population — projection, `ess sup`, `D_U`, trois impossibilités et leurs
portées corrigées ; (2) **l'annulation de `h`** — lemme mgf, coloriage, concentration, corollaire,
avec la remarque sur les trois endroits où `nαh` survit ; (3) la couche spatiale — damiers, module
`Ω_n`, rayon gonflé, calcul explicite à `b=0` ; (4) l'annexe primitive — factorisation de fibre,
témoin oscillant, transfert, `ρ=0` vs `−κ`, `κ` pour A1, trois branches de coin ; (5) la barrière du
biais — impossibilité, `QLH_β`, `λ*`, région `d<β/(2β+1)`, exclusion des modèles ; (6) **l'inférence
simultanée** — influence, CCK, le faux pas des différences mixtes, la preuve en cinq étapes ;
(7) robustesse et bornes inférieures ; (8) **ce qu'il faut changer dans l'article** (10 corrections
obligatoires + 5 ajouts) ; (9) ce qui reste ouvert.

Les trois vérifications numériques indépendantes y sont signalées comme telles (Ω à `b=0`, `λ*`,
`κ = 34/15` et `41/15`).
