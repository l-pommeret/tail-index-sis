# Q003 — Removing the bandwidth from the effective sample size

**Target model:** GPT-5.6 or equivalent, long-form mathematical work.
**Expected output:** LaTeX-ready English, paper notation, complete proofs.
**Companion documents:** the manuscript's Section 3 and Appendices A–B
(restated below in full, so this note is self-contained).

---

## 0. The task in one sentence

The score bound currently reads
$$\max_{j\le p}|\widehat\Psi_j-\Psi_j|
=O_{\mathbb P}\!\left(\sqrt{\frac{\log(pn)}{n\alpha h}}\right)+o(\Delta_{\min}),$$
and we want the stochastic term improved to
$$O_{\mathbb P}\!\left(\sqrt{\frac{\log(pn)}{n\alpha}}\right),$$
i.e. the rank bandwidth $h$ removed from the effective sample size, with
the consequent weakening of the sure-screening condition from
$\log(pn)=o(n\alpha h\Delta_{\min}^2)$ toward
$\log(pn)=o(n\alpha\Delta_{\min}^2)$.

Prove it, or prove that it cannot hold and identify the true rate. Do
not produce a plausible-looking derivation you have not checked: the
whole point is which of the two rates is correct.

---

## 1. The estimator (verbatim)

Data: i.i.d. pairs $(X_i,Y_i)$, $Y_i>0$, $X_i\in\mathbb R^p$, continuous
marginals; latent PIT $U_{ij}=F_j(X_{ij})$. Trimming
$I_\varepsilon=[\varepsilon,1-\varepsilon]$, $0<\varepsilon<1/2$.

Ranks $R_{ij}$ of $X_{ij}$ among $X_{1j},\dots,X_{nj}$,
$\widehat U_{ij}=R_{ij}/(n+1)$. For $u\in I_\varepsilon$,
$$\mathcal W_j(u)=\{i:|\widehat U_{ij}-u|\le h\},\quad
M_j(u)=|\mathcal W_j(u)|,\quad k_j(u)=\lfloor\alpha M_j(u)\rfloor,$$
and, with $Y_{j,u,(1)}\le\cdots\le Y_{j,u,(M_j(u))}$ the local order
statistics,
$$\widehat\xi_j(u)=\frac1{k_j(u)}\sum_{r=1}^{k_j(u)}
\log\frac{Y_{j,u,(M_j(u)-r+1)}}{Y_{j,u,(M_j(u)-k_j(u))}}.$$
The score averages the profile over the deterministic grid
$G=\{r/(n+1):r=1,\dots,n\}\cap I_\varepsilon$, of size $L_n=|G|$:
$$\widehat\Psi_j=\frac1{L_n}\sum_{u\in G}\widehat\xi_j(u),
\qquad
\Psi_j=\frac1{1-2\varepsilon}\int_{I_\varepsilon}\xi_j(u)\,du .$$
The screen ranks coordinates by increasing $\widehat\Psi_j$ and keeps the
$d$ smallest. $\Delta_j=\gamma^\star-\Psi_j$, $\gamma^\star=\max_u\gamma(u)$,
$\Delta_{\min}=\min_{j\in A}\Delta_j$; inactive coordinates have
$\Delta_j=0$.

Window counts are **deterministic**: $M_j(u)$ depends on the ranks only
through the grid position, so it is a nonrandom function of $u$, the same
for every $j$. Write $m_n^-=\lfloor\{h+\min(h,\varepsilon)\}(n+1)\rfloor-1$,
$K_n=\lfloor\alpha m_n^-\rfloor$, $S_n=4n+1$ (the number of distinct
window states per coordinate), $\delta_n=(n+1)^{-1}$.

## 2. The current proof, and exactly where it is lossy

The chain is:

1. **DKW event.** $D_n(\eta)=\{\max_{j,i}|U_{ij}-\widehat U_{ij}|
   \le\eta+\delta_n\}$ has $\mathbb P(D_n(\eta)^c)\le2pe^{-2n\eta^2}$.
   This is what forces $\log(pn)=o(nh^2)$.
2. **Rényi representation.** Conditionally on the column
   $\bm U_j=(U_{1j},\dots,U_{nj})$, the conditional PIT values
   $V_{ij}=1-F_{Y|j}(Y_i\mid U_{ij})$ are i.i.d.\ uniform, and for any
   index set $\mathcal W$ of size $m$ that is deterministic given
   $\bm U_j$, with $k=\lfloor\alpha m\rfloor$ and $V_{(1)}<\cdots<V_{(m)}$
   the ordered PIT values of its members,
   $$L_k=\frac1k\sum_{r=1}^k\log\frac{V_{(k+1)}}{V_{(r)}}$$
   is distributed as the mean of $k$ i.i.d.\ standard exponentials,
   independent of $V_{(k+1)}$.
3. **Frozen-location expansion.** On $D_n(\eta)\cap\{V_{ij}\ge\tau\}
   \cap\{V_{(k+1)}\le3\alpha\}$, for a window realized at $u$,
   $$|H-\xi_j(u)|\le
   2\omega_n^{\uparrow}(h+\eta+\delta_n,\tau,\alpha)
   +\xi_j(u)|L_k-1|+b_n(3\alpha)L_k,$$
   where $b_n(3\alpha)$ is the (E2) modulus and
   $\omega_n^{\uparrow}(r,\tau,\alpha)=\max_j\sup_{u}\sup_{|v-u|\le r}
   \sup_{w\in[\tau,3\alpha]}|\log\{Q_j(w^{-1},v)/Q_j(w^{-1},u)\}|$ is the
   tail-only (E3$^\uparrow$) modulus.
4. **Uniform profile bound.** A union bound over the $pS_n$ admissible
   window states gives, for $\eta>0$, $0<t\le1$, $0<\tau\le3\alpha$,
   $$\mathbb P\Bigl\{\max_{j\le p}\sup_{u\in I_\varepsilon}
   |\widehat\xi_j(u)-\xi_j(u)|
   >2\omega_n^{\uparrow}+\Gamma t+(1+t)b_n(3\alpha)\Bigr\}$$
   $$\le\underbrace{2pe^{-2n\eta^2}}_{\text{ranks}}
   +\underbrace{pn\tau}_{\text{PIT truncation}}
   +\underbrace{2pS_ne^{-K_nt^2/4}}_{\text{spacings}}
   +\underbrace{pS_ne^{-(K_n+1)/4}}_{\text{threshold}} .$$
5. **Profiles to scores.** Pathwise,
   $$|\widehat\Psi_j-\Psi_j|\le
   \sup_{u\in I_\varepsilon}|\widehat\xi_j(u)-\xi_j(u)|
   +\Bigl|\frac1{L_n}\sum_{u\in G}\xi_j(u)-\Psi_j\Bigr|,$$
   the second term being a quadrature error bounded by
   $\omega_\xi(\delta_n)+O(\delta_n)$.
6. Choosing $t_n=\sqrt{(4/K_n)\log\{pS_n(pn)\}}$ makes the spacings term
   vanish and yields $\Gamma t_n=O(\sqrt{\log(pn)/K_n})
   =O(\sqrt{\log(pn)/(n\alpha h)})$ by $K_n\asymp n\alpha h$.

**The lossy step is 5 combined with 4.** The score is an average of
$L_n\asymp n$ profile values, and the proof bounds that average by its
supremum. The supremum genuinely costs $\sqrt{\log(pn)/(n\alpha h)}$;
the average should not.

## 3. The mechanism we believe recovers the factor

Windows have width $2h$ in the rank scale, so **two windows whose
centres are more than $2h$ apart use disjoint index sets**. Conditionally
on the column $\bm U_j$, the pairs $(V_{ij})_i$ are i.i.d.; hence
statistics built on disjoint index sets are **conditionally
independent**. The grid therefore carries about
$$N_h\;\asymp\;\frac{1}{2h}$$
mutually disjoint, hence conditionally independent, windows.

Write $S\asymp2nh$ for the number of grid points inside one window
width, and decompose the grid into $S$ interleaved "separated families"
$G_\tau$, $\tau=1,\dots,S$, each consisting of grid points spaced more
than $2h$ apart, so that all windows within one family are pairwise
disjoint. Then
$$\widehat\Psi_j=\frac1S\sum_{\tau=1}^SA_{j,\tau},
\qquad
A_{j,\tau}=\frac1{|G_\tau|}\sum_{u\in G_\tau}\widehat\xi_j(u),$$
and, by convexity,
$|\widehat\Psi_j-\Psi_j|\le\max_\tau|A_{j,\tau}-\Psi_{j,\tau}|
+\text{quadrature}$, with $\Psi_{j,\tau}$ the corresponding population
average.

Each $A_{j,\tau}$ is an average of $N_h\asymp1/(2h)$ **conditionally
independent** terms, each a mean of $k\asymp2n\alpha h$ exponential
spacings, hence of individual conditional variance $\asymp1/k$. So
$$\operatorname{Var}(A_{j,\tau}\mid\bm U_j)
\;\asymp\;\frac{1}{N_hk}
\;\asymp\;\frac{2h}{2n\alpha h}
\;=\;\frac1{n\alpha},$$
and a Bernstein bound with a union over $pS\asymp2pnh$ pairs
$(j,\tau)$ — costing $\log(pnh)=O(\log(pn))$ — should give
$$\max_{j\le p}\max_{\tau\le S}|A_{j,\tau}-\Psi_{j,\tau}|
=O_{\mathbb P}\!\left(\sqrt{\frac{\log(pn)}{n\alpha}}\right).$$

That is the target. **Verify this, make it rigorous, or refute it.**

## 4. Obstructions you must confront, not skip

1. **The threshold and spacings events are still uniform over all
   windows.** The frozen-location expansion is valid only on
   $\{V_{(k+1)}\le3\alpha\}$, and controlling that for *every* admissible
   window costs $pS_ne^{-(K_n+1)/4}$, which requires
   $\log(pn)=O(K_n)\asymp n\alpha h$. So even if the main term improves,
   a **side condition $\log(pn)=o(n\alpha h)$ may survive**, independent
   of $\Delta_{\min}$. Determine honestly whether it does. If it does,
   the improvement is
   $$\log(pn)=o(n\alpha\Delta_{\min}^2)\ \text{ and }\
   \log(pn)=o(n\alpha h),$$
   which is still a strict gain whenever $\Delta_{\min}^2<h$ — the
   weak-signal regime — and should be stated that way. If the side
   condition can be removed (e.g. by bounding the contribution of the
   bad windows to the *average* rather than excluding them), do so and
   say how.
2. **Sub-exponential tails.** $L_k$ is a mean of $k$ exponentials, so
   $A_{j,\tau}$ is an average of sub-exponential variables. Use a
   Bernstein inequality with the correct sub-exponential (not
   sub-Gaussian) regime, and check which of the two Bernstein branches
   is active at the target deviation. State the regime in which
   $\sqrt{\log(pn)/(n\alpha)}$ is the sub-Gaussian branch.
3. **Bias is untouched.** The terms $2\omega_n^{\uparrow}$ and
   $b_n(3\alpha)L_k$ are biases, not fluctuations; averaging over the
   grid does not reduce them. The $o(\Delta_{\min})$ part of the bound
   must stay exactly as it is. Do not claim any improvement there.
4. **Conditional independence needs care at the boundary.** The
   separated families must be constructed from the *rank* windows, whose
   membership is deterministic given $\bm U_j$ (Lemma A.8 of the
   manuscript). Verify that disjointness of index sets follows from
   separation of centres by more than $2h$ in the rank scale, including
   at the two boundaries of $I_\varepsilon$ where windows are clipped.
5. **Dependence across $j$ is irrelevant but must be handled by the
   union bound only.** Coordinates are dependent; the proof must use no
   independence across $j$, only a union bound, as in the current
   Proposition A.13.
6. **$\Psi_{j,\tau}$ versus $\Psi_j$.** Each separated family averages
   $\xi_j$ over a sparse subgrid, not over $I_\varepsilon$. Quantify
   $\max_\tau|\Psi_{j,\tau}-\Psi_j|$ through the modulus $\omega_\xi$ and
   check it is $o(\Delta_{\min})$ under (E1); a family of spacing $2h$
   has quadrature error $O(\omega_\xi(2h))$, **not** $O(\omega_\xi(1/n))$,
   which is a genuine new requirement — state it explicitly, and say
   whether it forces $\omega_\xi(h)=o(\Delta_{\min})$, i.e. $L_\gamma
   h=o(\Delta_{\min})$ for Lipschitz $\gamma$.

Point 6 is the one most likely to bite. Do not paper over it.

## 5. Deliverables

1. A revised **Proposition (finite-sample score bound)** replacing
   Proposition A.15, with the explicit four-term probability bound and
   the new deviation term.
2. Revised **Theorem 3.2 (score concentration)** and **Theorem 3.3
   (sure screening)** with their exact hypotheses, including every side
   condition that survives.
3. Complete proofs, with the blocking construction stated precisely
   (how the families are indexed, why the windows are disjoint, where
   conditional independence is used).
4. A short subsection **"What the bandwidth still costs"** collecting
   the conditions in which $h$ still appears
   ($\log(pn)=o(nh^2)$ from DKW, $K_n\ge2$, and whichever of the
   threshold/spacings/quadrature conditions survive), with the admissible
   range of $(a,b)$ for $\alpha=n^{-a}$, $h=n^{-b}/2$ compared against
   the current range.
5. A statement of whether the new rate is **optimal**: is
   $\sqrt{\log(pn)/(n\alpha)}$ a lower bound for this score functional,
   or could further structure improve it? A heuristic is acceptable here
   if labelled as such.
6. An explicit numerical comparison at the manuscript's design
   ($n=2000$, $p\in\{500,1000,2000\}$, $\alpha=n^{-0.30}$,
   $h=n^{-0.15}/2$, so $n\alpha h\approx32.7$, $n\alpha\approx204.5$,
   $\log(pn)\approx14.5$ at $p=1000$, gaps $0.186$, $0.107$, $0.158$):
   what do the old and new bounds give, and does the new one become
   informative at this $n$?

## 6. Constraints

- Uniformity over $j\le p$ with $p=p_n\to\infty$ and over the triangular
  array in $n$ is mandatory. No independence across coordinates.
- Do not weaken the target: separation of active from inactive scores
  needs uniform accuracy below $\Delta_{\min}/2$.
- Do not introduce covariate density conditions; the empirical-rank
  construction exists precisely to avoid them.
- Keep the screen unchanged. This is a sharper analysis of the same
  estimator, not a new estimator. If the improvement requires modifying
  the score (for instance averaging over a coarser grid of spacing $h$
  rather than $1/n$), say so explicitly and quantify what changes — that
  would be a different, and reportable, conclusion.
- State clearly at every step what is assumed and what is proved.
