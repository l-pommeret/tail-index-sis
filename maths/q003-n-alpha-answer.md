# Q003 — Removing the bandwidth from the effective sample size

The requested question is whether the stochastic score error can be sharpened from (O_{\mathbb P}{\sqrt{\log(pn)/(n\alpha h)}}) to (O_{\mathbb P}{\sqrt{\log(pn)/(n\alpha)}}), without changing the estimator, while retaining all side conditions honestly.  The answer is **yes**. The improvement is real, but the proof needs two corrections to the proposed argument:

1. One must average the **signed** Rényi-spacing fluctuation, not the absolute profilewise bound.
2. One must center each color class at its own sparse-grid population average. Those sparse-grid averages recombine exactly into the original fine-grid average, so no new assumption (\omega_\xi(h)=o(\Delta_{\min})) is needed.

The condition
[
\log(pn)=o(n\alpha h)
]
does survive, because the present tail assumptions control the frozen quantile expansion only when every local threshold lies inside the tail region (V_{(k+1)}\le 3\alpha). Thus the corrected screening conditions are
[
\boxed{
\log(pn)=o(n\alpha\Delta_{\min}^{2}),\qquad
\log(pn)=o(n\alpha h),\qquad
\log(pn)=o(nh^{2}),
}
]
together with the unchanged bias assumptions.

The full finite-sample, concentration, screening, tuning, optimality, and numerical conclusions requested in Q003 follow. 

---

## 1. Rank-grid notation and the window coloring

Write
[
\delta_n=\frac1{n+1},
\qquad
\mathcal R_n
============

\left{
r\in{1,\ldots,n}:
\varepsilon\le r\delta_n\le 1-\varepsilon
\right},
]
so that
[
G={u_r=r\delta_n:r\in\mathcal R_n},
\qquad
L_n=|\mathcal R_n|.
]
The set (\mathcal R_n) is an interval of consecutive integers.

For (r\in\mathcal R_n), write
[
\mathcal W_{jr}:=\mathcal W_j(u_r),\qquad
m_r:=M_j(u_r),\qquad
k_r:=\lfloor\alpha m_r\rfloor.
]
Because the windows are defined from ranks, (m_r) does not depend on (j) or on the observed values; in particular,
[
m_r\ge m_n^-,
\qquad
k_r\ge K_n:=\lfloor\alpha m_n^-\rfloor.
]

Define the rank-separation integer
[
q_n:=\lfloor 2h(n+1)\rfloor+1,
]
and the number of colors
[
\chi_n:=\min(q_n,L_n).
]
When (q_n\le L_n), color (r\in\mathcal R_n) by its residue modulo (q_n):
[
\mathcal R_{n,c}
================

{r\in\mathcal R_n:r\equiv c\pmod{q_n}},
\qquad c=0,\ldots,q_n-1.
]
When (q_n>L_n), use singleton color classes. In either case there are (\chi_n) nonempty color classes.

The effective sample size appearing below is
[
\boxed{
\mathfrak n_n:=\frac{L_nK_n}{\chi_n}.
}
]

### Lemma 1 — Disjointness within a color class

For every coordinate (j), two distinct windows whose centers belong to the same color class are disjoint:
[
r,s\in\mathcal R_{n,c},\quad r\ne s
\quad\Longrightarrow\quad
\mathcal W_{jr}\cap\mathcal W_{js}=\varnothing.
]

#### Proof

If (q_n\le L_n), then (|r-s|\ge q_n), and therefore
[
|u_r-u_s|
=========

|r-s|\delta_n
\ge q_n\delta_n

> 2h.
> ]
> Suppose that a row (i) belonged to both windows. Then
> [
> |u_r-u_s|
> \le
> |u_r-\widehat U_{ij}|+
> |\widehat U_{ij}-u_s|
> \le2h,
> ]
> a contradiction.

Equivalently, in rank coordinates,
[
|R_{ij}-r|\le h(n+1),
\qquad
|R_{ij}-s|\le h(n+1)
]
would imply
[
|r-s|\le2h(n+1)<q_n.
]
Clipping a window at the lower or upper rank boundary only removes row indices and cannot create an intersection. Thus the assertion remains valid at both boundaries of (I_\varepsilon).

If (q_n>L_n), every class is a singleton and the result is immediate. (\square)

### Lemma 2 — Conditional independence within a color class

Let
[
\mathscr U_j=\sigma(U_{1j},\ldots,U_{nj}).
]
Conditionally on (\mathscr U_j), the variables
[
V_{ij}:=1-F_{Y\mid j}(Y_i\mid U_{ij}),
\qquad i=1,\ldots,n,
]
are independent (\operatorname{Unif}(0,1)). Hence, for fixed (j), all statistics built from the (V_{ij})'s in windows belonging to a common color class are conditionally independent.

No claim of independence across coordinates (j) is made or used.

#### Proof

The rows ((U_{ij},Y_i)) are independent over (i). Conditional on (U_{ij}), continuity of (F_{Y\mid j}(\cdot\mid U_{ij})) gives
[
V_{ij}\mid U_{ij}\sim\operatorname{Unif}(0,1).
]
Consequently, conditional on (\mathscr U_j), the (V_{ij})'s are independent uniforms. By Lemma 1, windows in one color class use disjoint row-index sets, so the corresponding window statistics are conditionally independent. Dependence between different coordinates, including the fact that they share the same (Y_i)'s, is handled solely through a union bound over (j). (\square)

---

## 2. The signed frozen-location expansion

For a grid point (u_r), let
[
V_{jr,(1)}<\cdots<V_{jr,(m_r)}
]
be the ordered conditional PIT values in (\mathcal W_{jr}), and define
[
L_{jr}
======

\frac1{k_r}\sum_{\ell=1}^{k_r}
\log\frac{V_{jr,(k_r+1)}}{V_{jr,(\ell)}}.
]
Conditionally on (\mathscr U_j),
[
L_{jr}\stackrel d=
\frac1{k_r}\sum_{\ell=1}^{k_r}E_{jr\ell},
\qquad
E_{jr\ell}\stackrel{\mathrm{iid}}{\sim}\operatorname{Exp}(1).
]

Set
[
\Omega_n(\eta,\tau)
:=
\omega_n^\uparrow
(h+\eta+\delta_n,\tau,\alpha),
\qquad
\beta_n:=b_n(3\alpha).
]

Let
[
\mathcal D_n(\eta)
==================

\left{
\max_{j\le p}\max_{i\le n}
|U_{ij}-\widehat U_{ij}|
\le\eta+\delta_n
\right},
]
[
\mathcal V_n(\tau)
==================

\left{
\min_{j\le p}\min_{i\le n}V_{ij}\ge\tau
\right},
]
and
[
\mathcal T_n
============

\bigcap_{j\le p}
\bigcap_{r\in\mathcal R_n}
\left{
V_{jr,(k_r+1)}\le3\alpha
\right}.
]

### Lemma 3 — Signed local decomposition

On
[
\mathcal D_n(\eta)\cap\mathcal V_n(\tau)\cap\mathcal T_n,
]
simultaneously for all (j\le p) and (r\in\mathcal R_n),
[
\boxed{
\widehat\xi_j(u_r)-\xi_j(u_r)
=============================

\xi_j(u_r){L_{jr}-1}+R_{jr},
}
]
where
[
\boxed{
|R_{jr}|
\le
2\Omega_n(\eta,\tau)+\beta_n L_{jr}.
}
]

#### Proof

Fix (j,r) and write (u=u_r), (m=m_r), and (k=k_r). For each row in (\mathcal W_{jr}), define the frozen-location observation
[
Z_i=Q_j(V_{ij}^{-1},u).
]
On (\mathcal D_n(\eta)),
[
|U_{ij}-u|
\le
|U_{ij}-\widehat U_{ij}|
+
|\widehat U_{ij}-u|
\le h+\eta+\delta_n.
]
On (\mathcal V_n(\tau)\cap\mathcal T_n), all PIT values that determine the largest (k+1) frozen observations lie in ([\tau,3\alpha]). The local quantile modulus therefore gives the order-statistic comparison
[
\left|
\widehat\xi_j(u)-H^0_{jr}
\right|
\le2\Omega_n(\eta,\tau),
]
where
[
H^0_{jr}
========

\frac1k\sum_{\ell=1}^k
\log
\frac{Q_j(V_{jr,(\ell)}^{-1},u)}
{Q_j(V_{jr,(k+1)}^{-1},u)}
]
is the Hill statistic computed from the frozen observations.

For completeness, body observations with (V_{ij}>3\alpha) cannot invalidate the order-statistic comparison. By monotonicity of the quantile in its tail argument,
[
Q_j(V_{ij}^{-1},U_{ij})
\le
Q_j((3\alpha)^{-1},U_{ij})
\le
e^{\Omega_n}
Q_j((3\alpha)^{-1},u),
]
while
[
Q_j(V_{jr,(k+1)}^{-1},u)
\ge
Q_j((3\alpha)^{-1},u).
]
Thus no body observation can exceed the relevant frozen threshold by more than the multiplicative factor (e^{\Omega_n}). The largest (k+1) order statistics are consequently bracketed by their frozen counterparts multiplied by (e^{\pm\Omega_n}), and the Hill functional changes by at most (2\Omega_n).

Now write
[
Q_j(s,u)=s^{\xi_j(u)}\ell_j(s,u).
]
Then
[
H^0_{jr}
========

\xi_j(u)L_{jr}+\rho_{jr},
]
where
[
\rho_{jr}
=========

\frac1k\sum_{\ell=1}^{k}
\log
\frac{\ell_j(V_{jr,(\ell)}^{-1},u)}
{\ell_j(V_{jr,(k+1)}^{-1},u)}.
]
Because (V_{jr,(k+1)}\le3\alpha), set
[
s=V_{jr,(k+1)}^{-1}\ge(3\alpha)^{-1},
\qquad
t_\ell=
\frac{V_{jr,(k+1)}}{V_{jr,(\ell)}}\ge1.
]
The definition of (b_n(3\alpha)) gives
[
\left|
\log
\frac{\ell_j(t_\ell s,u)}{\ell_j(s,u)}
\right|
\le
\beta_n\log t_\ell.
]
Averaging yields
[
|\rho_{jr}|\le\beta_nL_{jr}.
]
Combining the frozen-location and slow-variation errors gives
[
\widehat\xi_j(u)-\xi_j(u)
=========================

\xi_j(u)(L_{jr}-1)+R_{jr},
\qquad
|R_{jr}|
\le2\Omega_n+\beta_nL_{jr}.
]
(\square)

This signed identity is the step that is missing if one merely averages the old absolute inequality. Averaging
[
|\widehat\xi_j(u)-\xi_j(u)|
]
would retain an expected term of order (k^{-1/2}), and would not recover the factor (h^{-1/2}).

---

## 3. A colored sub-exponential inequality

The following lemma converts disjointness within colors into concentration for the full fine-grid average. It is slightly sharper than taking the maximum over color-class averages: there is no union bound over the colors.

### Lemma 4 — Colored Rényi-spacing Bernstein inequality

Fix (j), condition on (\mathscr U_j), and let (a_r\in[-1,1]) be arbitrary deterministic coefficients. Then, for every (t>0),
[
\boxed{
\mathbb P\left[
\left.
\frac1{L_n}\sum_{r\in\mathcal R_n}
a_r(L_{jr}-1)>t
,\right|,\mathscr U_j
\right]
\le
\exp\left{
-\frac{\mathfrak n_n}{4}(t^2\wedge t)
\right}.
}
]
The same bound holds for the lower tail, and hence
[
\mathbb P\left[
\left.
\left|
\frac1{L_n}\sum_{r\in\mathcal R_n}
a_r(L_{jr}-1)
\right|>t
,\right|,\mathscr U_j
\right]
\le
2\exp\left{
-\frac{\mathfrak n_n}{4}(t^2\wedge t)
\right}.
]

#### Proof

For (E\sim\operatorname{Exp}(1)), define
[
\psi(z)
=======

# \log\mathbb E e^{z(E-1)}

-z-\log(1-z),
\qquad z<1.
]
For (|z|\le1/2),
[
\psi(z)\le z^2.
]

Let the color classes be denoted by (\mathcal R_{n,c}), (c=1,\ldots,\chi_n), and set
[
S_j
===

\frac1{L_n}
\sum_{r\in\mathcal R_n}
a_r(L_{jr}-1).
]
For (\lambda>0), generalized Hölder gives
[
\begin{aligned}
\mathbb E!\left[
e^{\lambda S_j}\mid\mathscr U_j
\right]
&=
\mathbb E!\left[
\left.
\prod_{c=1}^{\chi_n}
\exp\left{
\frac{\lambda}{L_n}
\sum_{r\in\mathcal R_{n,c}}
a_r(L_{jr}-1)
\right}
\right|\mathscr U_j
\right]
\
&\le
\prod_{c=1}^{\chi_n}
\left[
\mathbb E!\left[
\left.
\exp\left{
\frac{\chi_n\lambda}{L_n}
\sum_{r\in\mathcal R_{n,c}}
a_r(L_{jr}-1)
\right}
\right|\mathscr U_j
\right]
\right]^{1/\chi_n}.
\end{aligned}
]
Within each color class the windows are disjoint, so the corresponding (L_{jr})'s are conditionally independent. Thus
[
\begin{aligned}
\log
\mathbb E!\left[
e^{\lambda S_j}\mid\mathscr U_j
\right]
&\le
\frac1{\chi_n}
\sum_{r\in\mathcal R_n}
k_r
\psi\left(
\frac{\chi_n\lambda a_r}{L_nk_r}
\right).
\end{aligned}
]
If
[
0\le\lambda\le\frac{L_nK_n}{2\chi_n}
=\frac{\mathfrak n_n}{2},
]
then
[
\left|
\frac{\chi_n\lambda a_r}{L_nk_r}
\right|
\le\frac12,
]
and hence
[
\begin{aligned}
\log
\mathbb E!\left[
e^{\lambda S_j}\mid\mathscr U_j
\right]
&\le
\frac1{\chi_n}
\sum_{r\in\mathcal R_n}
\frac{\chi_n^2\lambda^2a_r^2}{L_n^2k_r}
\
&\le
\frac{\chi_n\lambda^2}{L_nK_n}
==============================

\frac{\lambda^2}{\mathfrak n_n}.
\end{aligned}
]

If (0<t\le1), take (\lambda=\mathfrak n_nt/2). Chernoff's inequality gives
[
\mathbb P(S_j>t\mid\mathscr U_j)
\le
\exp\left(
-\frac{\mathfrak n_nt^2}{4}
\right).
]
If (t>1), take (\lambda=\mathfrak n_n/2). Then
[
-\lambda t+\frac{\lambda^2}{\mathfrak n_n}
==========================================

-\frac{\mathfrak n_nt}{2}
+
\frac{\mathfrak n_n}{4}
\le
-\frac{\mathfrak n_nt}{4}.
]
This proves the upper-tail claim. Applying it to (-a_r) gives the lower-tail bound. (\square)

The inequality has the correct sub-exponential transition:
[
\exp(-c\mathfrak n_nt^2),\qquad 0<t\le1,
]
and
[
\exp(-c\mathfrak n_nt),\qquad t>1.
]

---

## 4. Revised Proposition A.15

Define the fine-grid quadrature error
[
\mathcal Q_{\xi,n}
:=
\max_{j\le p}
\left|
\frac1{L_n}\sum_{r\in\mathcal R_n}\xi_j(u_r)-\Psi_j
\right|.
]
A standard Riemann-sum estimate gives, for all sufficiently large (n),
[
\mathcal Q_{\xi,n}
\le
2\omega_\xi(\delta_n)
+
\frac{4\Gamma\delta_n}{1-2\varepsilon}.
]

### Proposition A.15 — Finite-sample integrated-score bound

Assume (0<\alpha\le1/3), (K_n\ge2), (0<\tau\le3\alpha), and
[
0<\xi_j(u)\le\Gamma
]
uniformly over (j\le p) and (u\in I_\varepsilon). Then, for every (\eta>0), (\tau\in(0,3\alpha]), and (t>0),
[
\begin{aligned}
&\mathbb P\Bigg[
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|

>

\mathcal Q_{\xi,n}
+
2\Omega_n(\eta,\tau)
+
\beta_n
+
(\Gamma+\beta_n)t
\Bigg]
\
&\quad\le
\underbrace{2pe^{-2n\eta^2}}*{\mathrm{rank\ localization}}
+
\underbrace{pn\tau}*{\mathrm{PIT\ truncation}}
+
\underbrace{
2p\exp\left{
-\frac{\mathfrak n_n}{4}(t^2\wedge t)
\right}
}*{\mathrm{integrated\ spacings}}
+
\underbrace{
pL_n\exp\left{
-\frac23\alpha m_n^-
\right}
}*{\mathrm{local\ threshold}}.
\tag{A.15}
\end{aligned}
]

In particular, since (L_n\le S_n), (\alpha m_n^-\ge K_n), and (K_n\ge1),
[
pL_n e^{-2\alpha m_n^-/3}
\le
pS_n e^{-(K_n+1)/4}.
]
Thus the final term may be replaced by the original manuscript-compatible bound
[
pS_n e^{-(K_n+1)/4}.
]

#### Proof

On the good event
[
\mathcal E_n
============

\mathcal D_n(\eta)\cap
\mathcal V_n(\tau)\cap
\mathcal T_n,
]
Lemma 3 gives
[
D_{jr}
:=
\widehat\xi_j(u_r)-\xi_j(u_r)
=============================

\xi_j(u_r)(L_{jr}-1)+R_{jr},
]
with
[
|R_{jr}|
\le
2\Omega_n+\beta_nL_{jr}.
]
Write (Z_{jr}=L_{jr}-1). Since (L_{jr}=1+Z_{jr}),
[
D_{jr}
\le
2\Omega_n+\beta_n+
{\xi_j(u_r)+\beta_n}Z_{jr},
]
and
[
D_{jr}
\ge
-2\Omega_n-\beta_n+
{\xi_j(u_r)-\beta_n}Z_{jr}.
]
Therefore
[
\frac1{L_n}\sum_rD_{jr}

>

2\Omega_n+\beta_n+(\Gamma+\beta_n)t
]
implies
[
\frac1{L_n}\sum_r
\frac{\xi_j(u_r)+\beta_n}{\Gamma+\beta_n}
Z_{jr}

> t.
> ]
> Similarly,
> [
> \frac1{L_n}\sum_rD_{jr}
> <
> -2\Omega_n-\beta_n-(\Gamma+\beta_n)t
> ]
> implies
> [
> \frac1{L_n}\sum_r
> \frac{\xi_j(u_r)-\beta_n}{\Gamma+\beta_n}
> Z_{jr}
> <-t.
> ]
> Both coefficient sequences lie in ([-1,1]). Lemma 4 consequently yields, conditionally on (\mathscr U_j),
> [
> \mathbb P\left[
> \left.
> \left|
> \frac1{L_n}\sum_rD_{jr}
> \right|

2\Omega_n+\beta_n+(\Gamma+\beta_n)t,
\ \mathcal E_n
\right|\mathscr U_j
\right]
\le
2e^{-\mathfrak n_n(t^2\wedge t)/4}.
]
This bound is uniform in (\mathscr U_j). Integrating and taking a union bound over (j\le p) gives the third probability term. No independence across (j) is used.

For the remaining terms, the coordinatewise DKW inequality gives
[
\mathbb P{\mathcal D_n(\eta)^c}
\le2pe^{-2n\eta^2}.
]
Because every (V_{ij}) is marginally uniform,
[
\mathbb P{\mathcal V_n(\tau)^c}
\le pn\tau.
]

Finally, conditionally on (\mathscr U_j), for a window of size (m_r),
[
V_{jr,(k_r+1)}>3\alpha
]
is equivalent to
[
B_{jr}:=\sum_{i\in\mathcal W_{jr}}
\mathbf 1{V_{ij}\le3\alpha}
\le k_r,
]
where
[
B_{jr}\sim\operatorname{Bin}(m_r,3\alpha).
]
Its mean is (3\alpha m_r), while (k_r\le\alpha m_r). Hence
[
\mathbb P(B_{jr}\le k_r\mid\mathscr U_j)
\le
\mathbb P\left(
B_{jr}\le\frac13\mathbb EB_{jr}
\mid\mathscr U_j
\right)
\le
\exp\left(-\frac23\alpha m_r\right)
\le
\exp\left(-\frac23\alpha m_n^-\right).
]
A union bound over (pL_n) grid windows gives the fourth term.

Since
[
\widehat\Psi_j-\Psi_j
=====================

\frac1{L_n}\sum_rD_{jr}
+
\left{
\frac1{L_n}\sum_r\xi_j(u_r)-\Psi_j
\right},
]
adding (\mathcal Q_{\xi,n}) completes the proof. (\square)

---

## 5. Why there is no new (O{\omega_\xi(h)}) bias

This is the point at which the initially proposed max-over-colors proof can appear to require an additional spatial assumption.

For a color class (c), define
[
A_{j,c}
=======

\frac1{N_c}
\sum_{r\in\mathcal R_{n,c}}
\widehat\xi_j(u_r),
\qquad
B_{j,c}
=======

\frac1{N_c}
\sum_{r\in\mathcal R_{n,c}}
\xi_j(u_r),
]
where (N_c=|\mathcal R_{n,c}|), and let (w_c=N_c/L_n). Then
[
\widehat\Psi_j
==============

\sum_cw_cA_{j,c}
]
and, crucially,
[
\frac1{L_n}\sum_{r\in\mathcal R_n}\xi_j(u_r)
============================================

\sum_cw_cB_{j,c}.
]
Therefore
[
\left|
\widehat\Psi_j-
\frac1{L_n}\sum_r\xi_j(u_r)
\right|
\le
\max_c|A_{j,c}-B_{j,c}|.
]
Only after the color classes have been recombined does one compare the full fine-grid average with (\Psi_j). The resulting quadrature error is (\mathcal Q_{\xi,n}=O{\omega_\xi(\delta_n)+\delta_n}), exactly as before.

It is true that an individual sparse-grid average can satisfy only
[
\max_c|B_{j,c}-\Psi_j|
\le
C_\varepsilon
\left[
\omega_\xi(q_n\delta_n)+
\Gamma q_n\delta_n
\right]
=======

O{\omega_\xi(2h+\delta_n)+h}.
]
This quantity is not necessarily (o(\Delta_{\min})) under the original version of (E1), which controls only (\omega_\xi(\delta_n)). For Lipschitz profiles it would impose
[
L_\xi h=o(\Delta_{\min}).
]

But that individual sparse-grid comparison is **not part of the correct proof**. It would arise only from the unnecessarily strong and generally false centering
[
A_{j,c}\approx\Psi_j
]
for every color separately. The correct centering is (A_{j,c}\approx B_{j,c}), followed by the exact weighted recombination identity.

Thus:

[
\boxed{
\text{The coloring introduces no new population requirement }
\omega_\xi(h)=o(\Delta_{\min}).
}
]

The existing local-quantile bias condition at radius (h+\eta+\delta_n) remains, of course, unchanged.

---

## 6. Revised Theorem 3.2

Let
[
x_n=\log(pn),
\qquad
\eta_n=\sqrt{\frac{x_n}{n}},
\qquad
\tau_n=(pn)^{-2},
]
and define
[
B_n
===

2\omega_n^\uparrow
(h+\eta_n+\delta_n,\tau_n,\alpha)
+
b_n(3\alpha)
+
\mathcal Q_{\xi,n}.
]

### Theorem 3.2 — Uniform integrated-score concentration

Assume the model and quantile regularity conditions underlying the frozen-location expansion, and suppose:

[
0<\alpha_n\le\frac13,\qquad 0<h_n\le1-\varepsilon;
]

[
\log(pn)=o(nh_n^2);
\tag{C1}
]

[
\log(pn)=o(n\alpha_nh_n);
\tag{C2}
]

and (K_n\ge2) eventually. Then
[
\mathfrak n_n\asymp n\alpha_n,
]
with constants depending only on (\varepsilon), and
[
\boxed{
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|
=======================

O_{\mathbb P}\left(
\sqrt{\frac{\log(pn)}{n\alpha_n}}
\right)
+
B_n.
}
\tag{3.2}
]

In particular, if (B_n=o(1)), then
[
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|
=======================

O_{\mathbb P}\left(
\sqrt{\frac{\log(pn)}{n\alpha_n}}
\right)+o(1).
]

#### Proof

First,
[
L_n=(1-2\varepsilon)n+O(1).
]
Moreover,
[
m_n^-
=====

{h_n+\min(h_n,\varepsilon)}(n+1)+O(1)
\asymp nh_n.
]
Condition (C2) implies (n\alpha_nh_n\to\infty), and therefore
[
K_n
===

\alpha_nm_n^-{1+o(1)}
\asymp n\alpha_nh_n.
]

If (q_n\le L_n), then
[
\chi_n=q_n\asymp nh_n.
]
If (q_n>L_n), then (h_n) is bounded below by a positive constant depending on (\varepsilon), (\chi_n=L_n\asymp n), and (K_n\asymp n\alpha_n). In both cases,
[
\mathfrak n_n
=============

\frac{L_nK_n}{\chi_n}
\asymp n\alpha_n.
]

Apply Proposition A.15 with (\eta=\eta_n) and (\tau=\tau_n). The rank term satisfies
[
2pe^{-2n\eta_n^2}
=================

2pe^{-2\log(pn)}
\le\frac{2}{pn^2}
\longrightarrow0,
]
and
[
pn\tau_n=(pn)^{-1}\longrightarrow0.
]
The threshold term satisfies
[
pL_ne^{-2\alpha_nm_n^-/3}
\longrightarrow0
]
by (C2).

Choose
[
t_n=M\sqrt{\frac{x_n}{\mathfrak n_n}}
]
with a fixed sufficiently large (M). Since (C2) implies
[
x_n=o(n\alpha_nh_n)=o(n\alpha_n),
]
we have (t_n=o(1)), so the quadratic Bernstein branch applies. The integrated-spacing probability is bounded by
[
2p\exp\left(-\frac{M^2x_n}{4}\right)=o(1)
]
for (M>2). Hence
[
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|
\le
B_n+
(\Gamma+b_n(3\alpha_n))
O_{\mathbb P}\left(
\sqrt{\frac{x_n}{\mathfrak n_n}}
\right).
]
Because (b_n(3\alpha_n)) is bounded under the stated tail assumptions and (\mathfrak n_n\asymp n\alpha_n), the result follows. (\square)

### Which Bernstein branch is active?

The dimensionless deviation parameter is
[
t_n\asymp
\sqrt{\frac{\log(pn)}{n\alpha_n}}.
]
The sub-Gaussian branch requires (t_n=o(1)), equivalently
[
\log(pn)=o(n\alpha_n).
]
This is automatically implied by the surviving threshold condition
[
\log(pn)=o(n\alpha_nh_n),
]
because (h_n\le1).

If (\log(pn)) were comparable with or larger than (n\alpha_n), the linear sub-exponential branch would instead give a deviation of order
[
\frac{\log(pn)}{n\alpha_n},
]
and the square-root rate would no longer be valid.

---

## 7. Revised Theorem 3.3

Let (A) be the active set, (s=|A|), and suppose
[
\Psi_j\le\gamma^\star-\Delta_{\min}
\quad(j\in A),
\qquad
\Psi_j=\gamma^\star
\quad(j\notin A).
]
Let (\widehat A_d) be the set of indices of the (d) smallest estimated scores.

### Theorem 3.3 — Sure screening at the (n\alpha) scale

Assume the hypotheses of Theorem 3.2 and suppose additionally that

[
\boxed{
\log(pn)
========

o\left(
n\alpha_n\Delta_{\min,n}^{,2}
\right),
}
\tag{S1}
]
and
[
\boxed{
\omega_n^\uparrow
(h_n+\eta_n+\delta_n,\tau_n,\alpha_n)
+
b_n(3\alpha_n)
+
\omega_\xi(\delta_n)
+
\delta_n
========

o(\Delta_{\min,n}).
}
\tag{S2}
]

Then
[
\mathbb P\left(
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|
<
\frac{\Delta_{\min,n}}2
\right)
\longrightarrow1.
]
Consequently, for every (d\ge s),
[
\boxed{
\mathbb P(A\subseteq\widehat A_d)\longrightarrow1.
}
]
If (d=s), then
[
\boxed{
\mathbb P(\widehat A_s=A)\longrightarrow1.
}
]

#### Proof

Theorem 3.2 and (S1) give
[
\sqrt{\frac{\log(pn)}{n\alpha_n}}
=================================

o(\Delta_{\min,n}).
]
Condition (S2) gives (B_n=o(\Delta_{\min,n})). Hence
[
e_n:=
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|
=======================

o_{\mathbb P}(\Delta_{\min,n}).
]
In particular,
[
\mathbb P(2e_n<\Delta_{\min,n})\to1.
]

On the event (2e_n<\Delta_{\min,n}), for every active (j),
[
\widehat\Psi_j
\le
\Psi_j+e_n
\le
\gamma^\star-\Delta_{\min,n}+e_n,
]
while for every inactive (\ell),
[
\widehat\Psi_\ell
\ge
\Psi_\ell-e_n
=============

\gamma^\star-e_n.
]
Since
[
\gamma^\star-\Delta_{\min,n}+e_n
<
\gamma^\star-e_n,
]
every active score is strictly smaller than every inactive score. Thus the first (s) ranked coordinates are exactly (A), and any list of size (d\ge s) contains (A). (\square)

The corrected signal condition is therefore
[
\boxed{
\log(pn)=o(n\alpha\Delta_{\min}^{2}),
}
]
but it must be accompanied by
[
\boxed{
\log(pn)=o(n\alpha h).
}
]

---

## 8. Why the threshold side condition survives

The event
[
V_{jr,(k_r+1)}\le3\alpha
]
ensures that the denominator of the frozen Hill statistic lies in the region where the second-order tail modulus (b_n(3\alpha)) is available. When this event fails, the local threshold lies in the uncontrolled body of the conditional distribution.

Under the present assumptions, the conditional quantile below the (3\alpha) tail level can be arbitrarily small relative to the tail quantiles. Consequently, on a bad-threshold window, the Hill log-ratios can be arbitrarily large. There is no deterministic envelope with which to multiply the fraction of bad windows.

Thus one cannot replace uniform exclusion of bad windows by an average-number-of-bad-windows argument without adding at least one of:

[
\sup_{j,u,w>3\alpha}
\left|
\log Q_j(w^{-1},u)
\right|
<\infty
]
at a quantified rate;

a global body-to-tail quantile envelope;

or a clipped/robustified local statistic.

All of these would strengthen the model or change the estimator. Under the stated class and the unchanged score, the side condition
[
\log(pn)=o(n\alpha h)
]
therefore remains.

The sharper binomial calculation above improves the numerical constant from
[
e^{-(K_n+1)/4}
]
to approximately
[
e^{-(2/3)\alpha m_n^-},
]
but it does not change the asymptotic requirement.

---

## 9. What the bandwidth still costs

The factor (h) disappears from the principal stochastic accuracy, but not from every condition.

### 9.1 Rank localization

Taking
[
\eta_n\asymp\sqrt{\frac{\log(pn)}n}
]
and invoking the existing tail modulus at a radius proportional to (h_n) requires
[
\eta_n=o(h_n),
]
namely
[
\boxed{
\log(pn)=o(nh_n^2).
}
]

If one were willing to formulate the spatial modulus directly at
[
h_n+\sqrt{\frac{\log(pn)}n},
]
the DKW inequality itself would not require (\eta_n=o(h_n)). The stated condition is what permits retention of the manuscript's (O(h_n))-radius formulation.

### 9.2 Existence of local upper order statistics

One needs
[
K_n=\lfloor\alpha_nm_n^-\rfloor\ge2.
]
Since (m_n^-\asymp nh_n), this requires at least
[
n\alpha_nh_n\gtrsim1.
]
The stronger threshold side condition makes (K_n\to\infty).

### 9.3 Uniform threshold control

The unchanged tail-only assumptions require
[
\boxed{
\log(pn)=o(n\alpha_nh_n).
}
]

### 9.4 Spatial bias

The term
[
\omega_n^\uparrow
\left(
h_n+\sqrt{\frac{\log(pn)}n}+\delta_n,
\tau_n,\alpha_n
\right)
]
is untouched. In concrete Lipschitz models it can impose a condition such as
[
L_n^{\mathrm{tail}}h_n=o(\Delta_{\min,n}).
]
This is a property of the original local smoothing bias, not a cost introduced by the coloring.

### 9.5 Population quadrature

Only
[
\omega_\xi(\delta_n)+O(\delta_n)
]
appears. There is no additional (\omega_\xi(h_n)) term.

### 9.6 Main stochastic term

The integrated spacing fluctuation is
[
O_{\mathbb P}\left(
\sqrt{\frac{\log(pn)}{n\alpha_n}}
\right),
]
with no bandwidth factor.

---

## 10. Admissible exponent ranges

Let
[
\alpha_n=n^{-a},
\qquad
h_n=\frac12n^{-b},
\qquad
\Delta_{\min,n}\asymp n^{-d}.
]

### Polynomial (p), or (\log(pn)=O(\log n))

The revised conditions become
[
b<\frac12
]
from DKW,
[
a+b<1
]
from uniform threshold control, and
[
a+2d<1
]
from the new stochastic screening condition.

Thus the revised stochastic admissible region is
[
\boxed{
0<b<\frac12,\qquad
a+b<1,\qquad
a+2d<1.
}
]
If the existing spatial bias is Lipschitz and must satisfy (h_n=o(\Delta_{\min,n})), add
[
\boxed{b>d.}
]

By contrast, the old score analysis required
[
\log(pn)=o(n\alpha_nh_n\Delta_{\min,n}^2),
]
which becomes
[
\boxed{
a+b+2d<1.
}
]
Therefore the old coupled restriction
[
a+b+2d<1
]
is replaced by the two separate restrictions
[
a+b<1,
\qquad
a+2d<1.
]

For fixed gaps, (d=0), the exponent region does not enlarge because the threshold side condition (a+b<1) remains. For shrinking gaps, (d>0), the enlargement is strict.

### Exponentially growing dimension

If
[
\log p\asymp n^\kappa,
]
the revised conditions are
[
\boxed{
\kappa<1-2b,\qquad
\kappa<1-a-b,\qquad
\kappa<1-a-2d.
}
]
The old stochastic condition was
[
\boxed{
\kappa<1-a-b-2d.
}
]

---

## 11. Optimality

There are two distinct questions.

### 11.1 One-coordinate functional lower bound

The scale (1/\sqrt{n\alpha}) is minimax-natural over a class in which the body distribution is unrestricted and only an (\alpha)-fraction of each sample is guaranteed to carry stable tail-index information.

To see this, consider two distributions that:

1. have the same body and the same probability (\alpha) of exceeding a high threshold;
2. conditionally on exceeding that threshold, have exponential log-excess distributions with means (\gamma) and (\gamma+\delta).

The one-observation Kullback–Leibler divergence is
[
\operatorname{KL}(P_{\gamma+\delta},P_\gamma)
=============================================

\alpha,
\operatorname{KL}
{\operatorname{Exp}(\gamma+\delta),
\operatorname{Exp}(\gamma)}
===========================

O(\alpha\delta^2).
]
Thus
[
\operatorname{KL}
(P_{\gamma+\delta}^{\otimes n},P_\gamma^{\otimes n})
====================================================

O(n\alpha\delta^2).
]
Taking
[
\delta\asymp(n\alpha)^{-1/2}
]
keeps the two experiments mutually contiguous. Le Cam's two-point argument therefore prevents uniform estimation of the tail-index score at a smaller order over such a tail-local class.

This argument can be embedded with a profile constant in (u), so that the functional difference in (\Psi_j) is exactly (\delta).

### 11.2 Simultaneous maximum

A (p)-alternative packing suggests the lower scale
[
\sqrt{\frac{\log p}{n\alpha}}.
]
A fully formal packing requires specifying how the coordinate-indexed alternatives are embedded in the joint conditional-tail family, especially under arbitrary cross-coordinate dependence. Thus the (\sqrt{\log p}) assertion is best regarded here as a lower-bound heuristic.

The extra (\log n) inside (\log(pn)) is not intrinsic to the maximum fluctuation. The colored-spacing term itself involves only a union over (p):
[
2p\exp{-c\mathfrak n_nt^2}.
]
The additional (n) enters through uniform rank, threshold, and truncation events. Hence the intrinsic simultaneous scale is plausibly
[
\sqrt{\frac{\log p}{n\alpha}},
]
while the present complete theorem is conveniently stated with (\log(pn)).

Additional parametric structure can improve the rate. For example, in an exact Pareto family valid throughout the distribution, all observations may contain information about (\gamma), and a different estimator can attain an (n^{-1/2}) rate. The (n\alpha) lower bound pertains to the broad tail-local class and to the unchanged tail-fraction score, not to every possible structured submodel.

---

## 12. Numerical comparison at the manuscript design

Take
[
n=2000,
\qquad
\alpha=n^{-0.30}=0.1022565,
\qquad
h=\frac12n^{-0.15}=0.1598879.
]
Then
[
n\alpha=204.513,
\qquad
n\alpha h=32.699.
]

The rate proxies are
[
r_{\mathrm{old}}(p)
===================

\sqrt{\frac{\log(pn)}{n\alpha h}},
\qquad
r_{\mathrm{new}}(p)
===================

\sqrt{\frac{\log(pn)}{n\alpha}}.
]

[
\begin{array}{c|c|c|c}
p & \log(pn) & r_{\mathrm{old}}(p) & r_{\mathrm{new}}(p)\
\hline
500  & 13.8155 & 0.6500 & 0.2599\
1000 & 14.5087 & 0.6661 & 0.2664\
2000 & 15.2018 & 0.6818 & 0.2726
\end{array}
]

Thus the new analysis improves the nominal stochastic scale by the factor
[
\sqrt h=0.3999,
]
or equivalently by approximately (2.50).

At (p=1000), the three gaps and their required half-gaps are
[
\begin{array}{c|c}
\Delta & \Delta/2\
\hline
0.186 & 0.0930\
0.107 & 0.0535\
0.158 & 0.0790
\end{array}
]
whereas
[
r_{\mathrm{new}}(1000)=0.2664.
]
Thus even the constant-free new proxy exceeds every full gap, and exceeds every required half-gap by a large factor.

The corresponding screening-regime ratios are
[
R_{\mathrm{old}}
================

\frac{n\alpha h\Delta^2}{\log(pn)},
\qquad
R_{\mathrm{new}}
================

\frac{n\alpha\Delta^2}{\log(pn)}.
]
At (p=1000),
[
\begin{array}{c|c|c}
\Delta & R_{\mathrm{old}} & R_{\mathrm{new}}\
\hline
0.186 & 0.0780 & 0.4877\
0.107 & 0.0258 & 0.1614\
0.158 & 0.0563 & 0.3519
\end{array}
]
The asymptotic theorem requires the relevant ratio to diverge. The new ratios are six-and-a-quarter times larger, but all remain below (1).

There is also a finite-(n) loss hidden by the asymptotic notation. With the manuscript's (\varepsilon=0.05),
[
L_n=1800,\qquad
q_n=640,\qquad
m_n^-=418,\qquad
K_n=42,
]
and hence
[
\mathfrak n_n
=============

# \frac{1800\cdot42}{640}

118.125.
]
The exact conservative colored rate at (p=1000) is therefore
[
\sqrt{\frac{\log(2,000,000)}{118.125}}
======================================

0.3505.
]
In the simulation models (\Gamma=\gamma^\star=0.5). The canonical choice
[
t=2\sqrt{\frac{\log(pn)}{\mathfrak n_n}}
]
then contributes
[
\Gamma t=0.3505
]
before adding either tail bias or spatial bias.

For a (95%) allocation to the colored-spacing term alone, solving
[
2p\exp(-\mathfrak n_nt^2/4)\le0.05
]
at (p=1000) gives
[
t
=

# 2\sqrt{\frac{\log(40000)}{118.125}}

0.5990,
]
and hence
[
\Gamma t=0.2995.
]
This is still well above every half-gap.

The sharper threshold probability is numerically satisfactory:
[
pL_ne^{-2\alpha m_n^-/3}
\approx
7.6\times10^{-7}
\qquad(p=1000).
]
By contrast, the original much looser expression
[
pS_ne^{-(K_n+1)/4}
]
exceeds (1) and is numerically useless at this design. The sharper Chernoff constant therefore matters for finite-sample bookkeeping, even though it does not change the asymptotic side condition.

Finally,
[
nh^2=51.13,
\qquad
\frac{nh^2}{\log(pn)}=3.52,
\qquad
\frac{n\alpha h}{\log(pn)}=2.25.
]
Taking
[
\eta_n=\sqrt{\frac{\log(pn)}n}=0.0852
]
gives
[
\frac{\eta_n}{h}=0.533.
]
Thus the DKW inflation of the localization radius is not small at (n=2000); the spatial bias must actually be controlled at a radius near
[
h+\eta_n\approx0.245,
]
not merely at (h\approx0.160).

Hence the numerical verdict is:

[
\boxed{
\begin{minipage}{0.88\linewidth}
The bandwidth-free analysis is a substantial and mathematically genuine improvement. It moves the fluctuation into the correct sub-Gaussian branch and improves the nominal screening ratio by a factor (1/h\approx6.25). Nevertheless, at (n=2000) it does not provide an informative uniform sure-screening certificate for any of the reported gaps. Even before adding the unchanged bias terms, the rigorous leading deviation remains larger than every required half-gap.
\end{minipage}
}
]

---

## Final corrected conclusion

For the unchanged integrated sliding-rank Hill score,
[
\boxed{
\max_{j\le p}
|\widehat\Psi_j-\Psi_j|
=======================

O_{\mathbb P}\left(
\sqrt{\frac{\log(pn)}{n\alpha}}
\right)
+
o(\Delta_{\min})
}
]
is valid under the original bias assumptions, together with
[
\boxed{
\log(pn)=o(nh^2),
\qquad
\log(pn)=o(n\alpha h).
}
]
Sure screening follows under the additional signal condition
[
\boxed{
\log(pn)=o(n\alpha\Delta_{\min}^2).
}
]

The factor (h) in the old principal stochastic rate was an artifact of replacing the integrated score by the supremum of the local profile errors. It is not the effective-sample-size rate of the integrated score. The bandwidth remains only in rank localization, local tail occupancy, and the unchanged spatial bias—not in the leading integrated spacing fluctuation.
