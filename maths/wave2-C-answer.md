
Hostile audit: what survives, what fails, and what must be rewritten
Overall verdict

The manuscript’s central population fact—the projected tail index is a conditional essential supremum—and its printed finite-sample concentration theorem are basically sound under the direct projected-quantile assumptions actually used in the proof.

The campaign journal also contains several genuine improvements:

atomlessness can be removed by randomized conditional PITs;
equicontinuity of the projected profiles is unnecessary;
the separate rank-displacement condition log(pn)=o(nh
2
) is unnecessary;
the global PIT truncation (pn)
−2
 can be removed;
the spatial bias can be reduced from hlog(pn) to hlog(1/α);
the uniform doubling condition can be replaced by the directly relevant tilted-deficit condition;
a disjoint-block version of the estimator can attain stochastic order (nα)
−1/2
, rather than (nαh)
−1/2
.

But three campaign conclusions are false as stated.

First, full support cannot simply be removed while retaining

γ
⋆
=
[0,1]
p
max
	​

γ.

That quantity depends on values of γ outside the observed support and is not identifiable. The correct replacement is the essential maximum on the support of U.

Second, the structural active set A is not identifiable under unrestricted covariate dependence. Even after correcting the global maximum, the intrinsically identifiable target is the set of coordinates with a nonzero marginal envelope gap, not the coordinate set on which an arbitrary off-support extension of γ happens to depend.

Third, the proposed estimator

γ
	​

max
=
j
max
	​

Ψ
j
	​

,
Δ
j
	​

=
γ
	​

max
−
Ψ
j
	​


does not recover all detectable coordinates unless at least one zero-gap coordinate exists. When every coordinate is detectable, the largest population score is strictly below the global maximum, and at least one detectable coordinate is assigned estimated gap zero.

The claimed minimax threshold is not established. In fact, the block construction shows that the factor h is not intrinsic to estimation of the integrated score. A minimax theorem involving nαh would require a substantially more restrictive model class, usually one in which the signal itself is localized at spatial scale h.

1. The first fatal defect: the off-support maximum is not identifiable

The journal proposes replacing full fibre support by

C
j
	​

(u)=suppK
j
	​

(u,⋅)

and writing

ξ
j
	​

(u)=
v∈C
j
	​

(u)
max
	​

γ(u,v).

That identity is correct when γ is continuous. What is not correct is to continue comparing ξ
j
	​

(u) with

γ
⋆
=
x∈[0,1]
p
max
	​

γ(x).

Consider

U
1
	​

=U
2
	​

=Z,Z∼Unif(0,1),

and suppose that, conditionally on U, Y is exactly Pareto with tail index

γ
1
	​

(u
1
	​

,u
2
	​

)=1+(u
1
	​

−u
2
	​

)
2
.

On the support {u
1
	​

=u
2
	​

}, the tail index is identically one. Hence the observed law is the same as under

γ
0
	​

(u
1
	​

,u
2
	​

)≡1.

Nevertheless,

[0,1]
2
max
	​

γ
1
	​

=2,
[0,1]
2
max
	​

γ
0
	​

=1.

For either coordinate,

C
1
	​

(u)=C
2
	​

(u)={u},ξ
1
	​

(u)=ξ
2
	​

(u)=1.

Thus the printed definition gives

Δ
j,1
	​

=2−1=1,Δ
j,0
	​

=1−1=0,

even though the two observed data laws are identical.

So the proposed no-(S) target is not merely difficult to estimate. It is not a functional of the data-generating law.

Correct support-intrinsic target

Define

γ
U
⋆
	​

=esssupγ(U).

When γ is continuous on the compact support
S=supp(U),

γ
U
⋆
	​

=
x∈S
max
	​

γ(x).

Then set

Δ
j
U
	​

=γ
U
⋆
	​

−
∣I
ε
	​

∣
1
	​

∫
I
ε
	​

	​

ξ
j
	​

(u)du.

This quantity is invariant under every modification of γ outside the support of U.

Let

M
U
	​

={x∈S:γ(x)=γ
U
⋆
	​

}.

For almost every u,

ξ
j
	​

(u)=γ
U
⋆
	​

⟺{ι
j
	​

(u,v):v∈C
j
	​

(u)}∩M
U
	​


=∅.

Consequently,

Δ
j
U
	​

>0⟺λ[{u∈I
ε
	​

:C
j
	​

(u)∩{v:ι
j
	​

(u,v)∈M
U
	​

}=∅}]>0.
	​


No continuity of u↦ξ
j
	​

(u) is needed for this equivalence: the integrand
γ
U
⋆
	​

−ξ
j
	​

(u) is measurable and nonnegative, so its integral is positive exactly when it is positive on a set of positive measure.

This is the correct version of campaign claim 1.

2. Structural activity is not identifiable under unrestricted dependence

The journal’s second counterexample is correct.

Let

U
1
	​

=U
2
	​

=Z∼Unif(0,1)

and consider

γ
(1)
(u
1
	​

,u
2
	​

)=1+u
1
	​

,γ
(2)
(u
1
	​

,u
2
	​

)=1+
2
u
1
	​

+u
2
	​

	​

.

They agree on the support u
1
	​

=u
2
	​

:

γ
(1)
(Z,Z)=γ
(2)
(Z,Z)=1+Z.

Using the same conditional Pareto law on that support gives identical joint laws of (U,Y), while

A
(1)
={1},A
(2)
={1,2}.

Thus A, as defined by varying coordinates over the full cube, is not identifiable without a support condition.

There is an even deeper issue. Under deterministic relations such as U
1
	​

=U
2
	​

, the same observable tail index 1+Z can be represented using coordinate 1 or coordinate 2. There need not be a unique smallest coordinate representation on the support.

The identifiable screening target is therefore

D
U
	​

={j:Δ
j
U
	​

>0},

not A.

Under full fibre support,

j∈
/
A⟹Δ
j
U
	​

=0,

so D
U
	​

⊆A. If every active coordinate satisfies the argmax-projection condition, then D
U
	​

=A. Without such an identification condition, recovery of A is impossible.

The main screening theorem should therefore be written first for D
U
	​

, followed by a separate corollary saying that it recovers A when D
U
	​

=A.

3. The projection theorem survives

The manuscript’s projection theorem is correct under its uniform tail-remainder assumption.

Let

T=γ(U),ξ=esssup(T∣G),

and write

N
y
	​

=E[c(U)y
−1/T
∣G].

Under (C2),

P(Y>y∣G)=N
y
	​

{1+o(1)}

almost surely, with a deterministic relative error bound.

Introduce the tilted conditional law

K
y
⋆
	​

(A)=
N
y
	​

E[c(U)y
−1/T
1
{T∈A}
	​

∣G]
	​

.

For each δ>0,

K
y
⋆
	​

{T≤ξ−δ}⟶0a.s.

Indeed, on {ξ≥
γ
	​

+δ}, put
r=ξ−δ and s=ξ−δ/2. By the defining property of the conditional essential supremum,

π
δ
	​

=P(T>s∣G)>0

almost surely on that event. Therefore

K
y
⋆
	​

{T≤r}≤
c
−
	​

π
δ
	​

y
−1/s
c
+
	​

y
−1/r
	​

⟶0.

Hence K
y
⋆
	​

⇒δ
ξ
	​

, and

N
y
	​

N
ty
	​

	​

=∫t
−1/g
K
y
⋆
	​

(dg)⟶t
−1/ξ
.

The uniform relative remainder transfers the limit to the projected survival function.

Exact weakening of (C2)

Uniform convergence in (C2) is sufficient but stronger than necessary. It is enough that

E[c(U)y
−1/T
∣G]
E[c(U)y
−1/T
∣r(y,U)∣∣G]
	​

⟶0,

where

F
ˉ
(y∣U)=c(U)y
−1/T
{1+r(y,U)}.

This is the genuinely relevant averaged remainder condition. The tilting proof is unchanged.

Why some uniform or averaged condition is essential

Let U∼Unif(0,1), choose 0<β<1, and conditionally on U=u let

Y={
u
−1/β
,
P,
	​

with probability 1/2,
with probability 1/2,
	​


where P is standard Pareto:
P(P>y)=y
−1
, y≥1.

For every fixed u, once y>u
−1/β
,

P(Y>y∣U=u)=
2
1
	​

y
−1
.

Every conditional distribution therefore has tail index one, with an eventually constant slowly varying factor.

But marginally,

P(Y>y)=
2
1
	​

y
−β
+
2
1
	​

y
−1
∼
2
1
	​

y
−β
.

The projected tail index is 1/β>1, not one.

Thus pointwise conditional regular variation alone does not permit mixture projection. A common tail onset, or the corresponding averaged remainder condition, is indispensable.

4. Support stability: what continuity in the copula can and cannot mean

The envelope is stable under perturbations of the supports, not under ordinary weak, total-variation, or Wasserstein perturbations of the conditional probabilities.

Let C,C
′
 be compact subsets and let γ have modulus of continuity ω
γ
	​

. Then

	​

v∈C
max
	​

γ(v)−
v∈C
′
max
	​

γ(v)
	​

≤ω
γ
	​

{d
H
	​

(C,C
′
)},

where d
H
	​

 is Hausdorff distance.

Proof: if v
⋆
 maximizes γ on C, choose v
′
∈C
′
 within d
H
	​

(C,C
′
). Then

C
max
	​

γ=γ(v
⋆
)≤γ(v
′
)+ω
γ
	​

(d
H
	​

)≤
C
′
max
	​

γ+ω
γ
	​

(d
H
	​

),

and reverse the roles.

Consequently, for two conditional-support families,

∣Ψ
j
	​

−Ψ
j
′
	​

∣≤
∣I
ε
	​

∣
1
	​

∫
I
ε
	​

	​

ω
γ
	​

(d
H
	​

{C
j
	​

(u),C
j
′
	​

(u)})du.

There is no analogous bound in total variation or Wasserstein distance alone. Take

γ(v)=1+v,K
0
	​

=δ
0
	​

,K
η
	​

=(1−η)δ
0
	​

+ηδ
1
	​

.

Then

d
TV
	​

(K
η
	​

,K
0
	​

)=η,W
1
	​

(K
η
	​

,K
0
	​

)=η,

but

K
0
	​

esssup
	​

γ=1,
K
η
	​

esssup
	​

γ=2

for every η>0.

Hence the envelope functional is discontinuous even under weak convergence:
K
η
	​

⇒K
0
	​

, but the essential suprema do not converge.

Campaign target I.3 can therefore be completed only after specifying a support topology or a quantitative near-maximal-mass condition. “Continuity in the copula” without such a qualification is false.

5. Marginal invisibility is genuinely an impossibility

Let U
1
	​

,U
2
	​

 be independent uniforms and let g be a positive, nonconstant, one-periodic continuous function, for example

g(t)=γ
0
	​

+acos(2πt),γ
0
	​

>∣a∣.

Compare

γ
(0)
(u
1
	​

,u
2
	​

)=g(u
2
	​

)

with

γ
(1)
(u
1
	​

,u
2
	​

)=g(u
1
	​

+u
2
	​

),

where periodicity makes the latter well defined without a discontinuity at the boundary. Let Y∣U be exact Pareto with the indicated index.

Under the first model, coordinate 1 is inactive. Under the second it is active. Yet, conditional on U
1
	​

=u, the variable
(u+U
2
	​

)mod1 is uniform, so

L(Y∣U
1
	​

=u)

is identical under the two models. Therefore the entire sample

{(U
i1
	​

,Y
i
	​

):1≤i≤n}

has exactly the same law under both models.

Every test or score based solely on the coordinate-1 marginal experiment has maximal error at least 1/2, for every n.

This proves more than invisibility of the average-envelope score. The entire coordinatewise conditional law is uninformative.

Exact interaction order

For

γ(u
1
	​

,…,u
s
	​

)=g(u
1
	​

+⋯+u
s
	​

),

with independent uniform coordinates and periodic g, conditioning on any proper subset J⊊{1,…,s} leaves at least one independent uniform phase. Hence the entire conditional law of Y∣U
J
	​

 is constant in U
J
	​

.

Only after conditioning on all s coordinates does the tail index vary. Therefore

r
⋆
(j)=min{∣J∣:j∈J, Δ
J
	​

>0}=s.

The upper bound r
⋆
(j)≤s is sharp. No grouped screen of a fixed order smaller than s can repair this form of marginal invisibility.

Campaign claims 10 and 11 survive.

6. The delayed-tail lower bound survives, but only for the naked model

The 1/logn lower bound in the journal is valid under (C1)–(C2) with no quantitative common tail onset.

Fix γ
0
	​

>0 and let

T
n
	​

=n
2γ
0
	​

,q
n
	​

=T
n
−1/γ
0
	​

	​

=n
−2
,δ
n
	​

=
2γ
0
	​

logn
a
	​

.

Let U∼Unif(0,1).

Under the null, the conditional tail index is

γ
0
(n)
	​

(u)=γ
0
	​

+δ
n
	​

.

Under the alternative,

γ
1
(n)
	​

(u)=γ
0
	​

+δ
n
	​

u.

Below T
n
	​

, choose the same conditional distribution under both models, with

P(Y>T
n
	​

∣U=u)=q
n
	​

.

Conditionally on Y>T
n
	​

, let Y/T
n
	​

 be Pareto with the corresponding conditional index.

Both models satisfy conditional regular variation. For y>T
n
	​

,

F
ˉ
(y∣u)=q
n
	​

(
T
n
	​

y
	​

)
−1/γ
r
(n)
	​

(u)
=c
r,n
	​

(u)y
−1/γ
r
(n)
	​

(u)
,

and

c
r,n
	​

(u)=q
n
	​

T
n
1/γ
r
(n)
	​

(u)
	​


is bounded above and away from zero uniformly in n,u, because
δ
n
	​

logT
n
	​

=O(1).

The null gap is zero. Under the alternative, on any symmetric trimmed interval,

Δ
1,n
	​

=δ
n
	​

(1−
∣I
ε
	​

∣
1
	​

∫
I
ε
	​

	​

udu)=
2
δ
n
	​

	​

=
4γ
0
	​

logn
a
	​

.

Couple the two experiments so that U, the exceedance indicator, and the body observation are identical. They differ only if at least one observation exceeds T
n
	​

. Therefore

∥P
0,n
⊗n
	​

−P
1,n
⊗n
	​

∥
TV
	​

≤nq
n
	​

=
n
1
	​

.

Le Cam’s inequality gives, for every test ϕ
n
	​

,

P
0,n
	​

(ϕ
n
	​

=1)+P
1,n
	​

(ϕ
n
	​

=0)≥1−
n
1
	​

.

Thus gaps of order 1/logn cannot be detected uniformly over the bare (C1)–(C2) class.

This does not prove a lower bound under the primitive conditions (P1)–(P4), because those conditions impose a common quantitative tail onset and derivative control. It also does not establish the proposed nαhΔ
2
 threshold.

7. Randomized PITs remove atomlessness completely

Campaign claim 2 is correct.

For a conditional CDF F
j
	​

(⋅∣u), introduce independent
Z
ij
	​

∼Unif(0,1) on an auxiliary probability space and put

W
ij
	​

=F
j
	​

(Y
i
	​

−∣U
ij
	​

)+Z
ij
	​

{F
j
	​

(Y
i
	​

∣U
ij
	​

)−F
j
	​

(Y
i
	​

−∣U
ij
	​

)},V
ij
	​

=1−W
ij
	​

.

Then, conditionally on the column
(U
1j
	​

,…,U
nj
	​

), the V
ij
	​

 are i.i.d. uniforms and

Y
i
	​

=Q
j
	​

(V
ij
−1
	​

,U
ij
	​

)a.s.

This is the standard distributional transform. The auxiliary uniforms are used only in the proof; the Hill statistic computed from the observed Y
i
	​

’s is unchanged.

Consequences:

Continuity of F
Y∣j
	​

 is unnecessary.
The atomlessness clause in (P1) is unnecessary.
The Rényi representation still applies to the randomized uniforms.
Ties in Y create no difficulty because the quantile representation remains valid.

The direct tail-quantile conditions must of course be stated for generalized quantiles.

8. Equicontinuity of ξ
j
	​

 and the nh
2
 condition are unnecessary

The printed proof freezes a rank window at its deterministic rank-grid centre r/(n+1). This forces a DKW displacement of order

logp/n
	​

, which must be o(h).

The stronger argument freezes the same window at the latent central order statistic

A
jr
	​

=U
(r)j
	​

.

If a window contains ranks within h(n+1) of r, then every selected observation lies within a uniform order-statistic spacing of A
jr
	​

. For a spacing involving d adjacent ranks,

U
(r+d)j
	​

−U
(r)j
	​

∼Beta(d,n+1−d).

A beta Chernoff bound yields numerical constants C,c>0 such that, for

ρ
n
	​

(x)=h+C{
n
hx
	​

	​

+
n
x
	​

+
n
1
	​

},
P[
j≤p
max
	​

r∈R
n
	​

max
	​

i∈W
jr
	​

max
	​

∣U
ij
	​

−A
jr
	​

∣>ρ
n
	​

(x)]≤2p∣R
n
	​

∣e
−cx
.

Hence

ρ
n
	​

(x)=h{1+o(1)}

as soon as x=o(nh). The principal local-extreme condition

log(pn)=o(nαh)

already implies log(pn)=o(nh), since α≤1.

There is no need for

log(pn)=o(nh
2
).
The score target needs no modulus of continuity

The frozen target becomes

L
n
	​

1
	​

r∈R
n
	​

∑
	​

ξ
j
	​

(U
(r)j
	​

),

not the deterministic Riemann sum
∑
r
	​

ξ
j
	​

(r/(n+1)).

For any measurable 0≤f≤Γ, the sum of f over the central order statistics differs from

i=1
∑
n
	​

f(U
i
	​

)1
{U
i
	​

∈I
ε
	​

}
	​


only through observations near the two trimming boundaries. Hoeffding bounds for the latter sum and binomial bounds for the two boundary counts give

P[
j≤p
max
	​

	​

L
n
	​

1
	​

r∈R
n
	​

∑
	​

ξ
j
	​

(U
(r)j
	​

)−Ψ
j
	​

	​

>C
ε
	​

Γ{
n
x
	​

	​

+
n
1
	​

}]≤6pe
−cx
.

No continuity, equicontinuity, or bounded variation of the ξ
j
	​

’s is used.

Campaign claims 4 and 5 therefore survive in their stronger B form.

9. Two-sided threshold control removes the global PIT truncation

For m independent uniforms and k=⌊αm⌋≥2,

P{V
(k+1)
	​

∈
/
[α/3,3α]}≤2e
−k/4
.

The upper-tail half is the lemma already printed. For the lower half,

{V
(k+1)
	​

<α/3}={Bin(m,α/3)≥k+1}.

The binomial mean is μ=αm/3, and k+1>3μ. Chernoff’s bound gives

P{Bin(m,α/3)≥3μ}≤exp{−μ(3log3−2)}≤e
−k/4
.

Thus there is no need to control

i,j
min
	​

V
ij
	​


or introduce τ
n
	​

=(pn)
−2
, provided the spatial quantile condition is formulated in the affine-log form below.

Campaign claim 15 survives.

10. The correct spatial condition gives hlog(1/α), not hlog(pn)

Assume that, whenever ∣v−u∣≤r and s≥(3α)
−1
,

	​

log
Q
j
	​

(s,u)
Q
j
	​

(s,v)
	​

	​

≤A
n
	​

(r)+B
n
	​

(r)logs.
(HQ)

Also assume the tail increment bound

	​

log
ℓ
j
	​

(s,u)
ℓ
j
	​

(ts,u)
	​

	​

≤b
n
	​

logt

and

B
n
	​

(r)<
ξ
	​

−b
n
	​

.

The last inequality implies that both

s⟼Q
j
	​

(s,u)s
−B
n
	​

(r)
ands⟼Q
j
	​

(s,u)s
B
n
	​

(r)

are nondecreasing in the working tail.

Let H be the actual Hill statistic and 
H
ˇ
 the statistic after freezing all conditioning values at u. Sorting the upper and lower envelopes gives

∣H−
H
ˇ
∣≤2A
n
	​

(r)+B
n
	​

(r){L
k
	​

+2log
V
(k+1)
	​

1
	​

}.

On

V
(k+1)
	​

∈[α/3,3α],

this becomes

∣H−
H
ˇ
∣≤2A
n
	​

(r)+B
n
	​

(r){L
k
	​

+2log
α
3
	​

}.

Since

H
ˇ
=ξ
j
	​

(u)L
k
	​

+R
k
	​

,∣R
k
	​

∣≤b
n
	​

L
k
	​

,

the event ∣L
k
	​

−1∣≤t gives

∣H−ξ
j
	​

(u)∣≤2A
n
	​

(r)+B
n
	​

(r){1+t+2log
α
3
	​

}+Γt+(1+t)b
n
	​

.
	​


This is the campaign’s claimed spatial improvement, and the derivation is valid. The cancellation in the Hill log-ratios is exactly what removes the need to control the smallest PIT value.

11. A corrected stronger finite-sample theorem for the printed estimator

Let L
n
	​

 be the number of score-grid points and K
n
	​

 the minimum local extreme count. Under randomized PIT quantile representation, bounded measurable projected indices

0<
ξ
	​

≤ξ
j
	​

(u)≤Γ,

the tail condition b
n
	​

, and (HQ), there exist numerical constants C,c>0 such that, for x≥1 and 0<t≤1,

	​

P[
j≤p
max
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣>2A
n
	​

{ρ
n
	​

(x)}+B
n
	​

{ρ
n
	​

(x)}{1+t+2log
α
3
	​

}
+Γt+(1+t)b
n
	​

+C
ε
	​

Γ{
n
x
	​

	​

+
n
1
	​

}]
≤2pe
−c
ε
	​

n
+Cp(L
n
	​

+1)e
−cx
+2pL
n
	​

e
−K
n
	​

t
2
/4
+2pL
n
	​

e
−K
n
	​

/4
.
	​

(*)

The terms are, respectively:

central anchors leaving the enlarged trimmed interval;
order-spacing and trimmed empirical-average errors;
exponential-spacing fluctuations;
the two-sided V
(k+1)
	​

 event.

Choosing

x=C
1
	​

log(pn),t=C
2
	​

K
n
	​

log(pn)
	​

	​

,

and using K
n
	​

≍nαh, gives

j≤p
max
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣=O
P
	​

(
nαh
log(pn)
	​

	​

)+O[A
n
	​

(Ch)+B
n
	​

(Ch)log
α
1
	​

+b
n
	​

].
	​


Only

log(pn)=o(nαh)

is needed for the stochastic part.

Thus the following printed assumptions are not doing real work:

atomlessness/continuity of the conditional CDF;
equicontinuity of ξ
j
	​

;
ω
ξ
	​

(1/n)=o(Δ
min
	​

);
log(pn)=o(nh
2
);
the global PIT truncation (pn)
−2
.

The exact constants 1/32, 6p+2pL, and so forth reported in the journal cannot be certified from the journal summary alone. The universal-constant inequality (∗) is what survives a fresh proof.

12. The score theorem should recover D
U
	​

, not assume D
U
	​

=A

Let

D
U
	​

={j:Δ
j
U
	​

>0},d
U
	​

=∣D
U
	​

∣,Δ
D
	​

=
j∈D
U
	​

min
	​

Δ
j
U
	​

.

By definition, every j∈
/
D
U
	​

 has

Ψ
j
	​

=γ
U
⋆
	​

.

Therefore, on

j
max
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣<Δ
D
	​

/2,

every coordinate in D
U
	​

 precedes every coordinate outside D
U
	​

.

So the correct theorem is

P(D
U
	​

⊆
D
d
	​

)→1,d≥d
U
	​

,

with exact recovery when d=d
U
	​

.

A separate identification corollary may state:

D
U
	​

=A⟹P(
D
s
	​

=A)→1.

This removes the unnecessary assumption that every active coordinate is detectable from the statement of the estimation theorem itself.

13. The proposed automatic gap estimator is false without a baseline coordinate

The journal claims that

γ
	​

max
=
j
max
	​

Ψ
j
	​

,
Δ
j
	​

=
γ
	​

max
−
Ψ
j
	​


allows recovery without knowing s.

This requires

j
max
	​

Ψ
j
	​

=γ
U
⋆
	​

,

equivalently, at least one coordinate with zero gap.

Take p=1 and

γ(u)=1−
2
u
	​

.

The sole coordinate is detectable:

Δ
1
U
	​

>0.

But

γ
	​

max
=
Ψ
1
	​

⟹
Δ
1
	​

=0

identically. It can never be selected by a positive threshold.

The corrected statement is:

If D
U
c
	​


=∅, so that at least one zero-gap coordinate exists, and
max
j
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣≤e, then

∣
Δ
j
	​

−Δ
j
U
	​

∣≤2e.

Any threshold satisfying

2e<λ<Δ
D
	​

−2e

recovers D
U
	​

.

This still requires a usable deterministic or data-driven upper bound e. It is not yet a fully adaptive procedure.

14. Minimum-rank aggregation is consistent only when all tunings are good

The manuscript’s fixed-nine-tuning corollary is correct under the event that every one of the nine tunings strictly separates all detectable and zero-gap coordinates. On that intersection,

r
j
	​

(a,b)≤d
U
	​

(j∈D
U
	​

),r
j
	​

(a,b)≥d
U
	​

+1(j∈
/
D
U
	​

)

for every tuning, so taking minimum ranks preserves separation.

What is false is the stronger adaptive interpretation that one good tuning suffices.

For two active coordinates a
1
	​

,a
2
	​

, two nulls n
1
	​

,n
2
	​

, and two tunings, consider

tuning 1
tuning 2
	​

a
1
	​

1
3
	​

a
2
	​

2
4
	​

n
1
	​

3
1
	​

n
2
	​

4
2
	​

	​


The first tuning is perfect. Minimum ranks are

(1,2,1,2).

The prescribed median/mean/maximum tie-breakers also tie a
1
	​

 with n
1
	​

 and a
2
	​

 with n
2
	​

; final randomization can exclude an active coordinate.

Thus the aggregation theorem presently proved is an intersection theorem, not an oracle theorem. It requires all nine tunings to be asymptotically valid. It gives no guarantee when only one or a subset is well tuned.

Campaign target II.7 remains open.

15. A disjoint-block estimator really does remove the factor h

The journal’s most important constructive claim is essentially correct, after making the estimator and its regularity assumptions explicit.

Partition the central rank range into

B
n
	​

≍h
−1

disjoint consecutive blocks, each containing

m
n
	​

≍nh

observations. In each block calculate one Hill statistic using

k
n
	​

=⌊αm
n
	​

⌋≍nαh

upper observations, and average the B
n
	​

 block statistics:

Ψ
j
	​

=
B
n
	​

1
	​

b=1
∑
B
n
	​

	​

H
jb
	​

.

Conditionally on the j-th covariate column, the PIT variables in different blocks are disjoint and independent. After freezing each block at its central latent order statistic,

H
jb
	​

=ξ
jb
	​

L
jb
	​

+bias
jb
	​

,

where the L
jb
	​

’s are independent means of k
n
	​

 exponentials.

The stochastic part is

S
j
	​

=
B
n
	​

1
	​

b=1
∑
B
n
	​

	​

ξ
jb
	​

(L
jb
	​

−1)=
B
n
	​

k
n
	​

1
	​

b=1
∑
B
n
	​

	​

r=1
∑
k
n
	​

	​

ξ
jb
	​

(E
jbr
	​

−1).

A direct exponential Bernstein bound gives

P[∣S
j
	​

∣>CΓ{
B
n
	​

k
n
	​

x
	​

	​

+
B
n
	​

k
n
	​

x
	​

}
	​

U
j
	​

]≤2e
−x
.

Since

B
n
	​

k
n
	​

≍nα,

a union bound over coordinates yields

j≤p
max
	​

∣S
j
	​

∣=O
P
	​

(
nα
logp
	​

	​

+
nα
logp
	​

).

Under the same affine-log spatial condition as above and, for example, uniformly Lipschitz ξ
j
	​

, the full bound is

j≤p
max
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣=
	​

O
P
	​

(
nα
logp
	​

	​

+
nα
logp
	​

+
n
logp
	​

	​

)
+O[A
n
	​

(Ch)+B
n
	​

(Ch)log
α
1
	​

+b
n
	​

+L
ξ
	​

h].
	​

	​


Local block validity still requires

log(p/h)=o(nαh),

because every individual block must contain enough usable extremes and the two-sided threshold event is union-bounded over pB
n
	​

 blocks.

The main separation condition becomes

logp=o(nαΔ
D
2
	​

),

subject to the local-validity and bias conditions.

This proves that the factor h is not information-theoretic for estimation of the integrated score. It arises from controlling the supremum of the sliding local profile and then averaging that supremum bound.

What remains unresolved is whether the original overlapping-window estimator itself also has the global nα stochastic rate after a direct dependence analysis.

The block estimator is different from the simulation procedure. No simulation or application conclusion can be transferred to it without rerunning the experiments.

16. The printed nαhΔ
2
 regime is not minimax-optimal as presently formulated

A minimax theorem must first specify a statistical model class. Here α and h are principally tuning parameters, not intrinsic parameters of the data-generating law.

Under the broad (C1)–(C2) class, the delayed-tail construction shows that the difficulty can be controlled by an arbitrary tail-onset sequence. There is no universal threshold involving nαh.

Under the primitive model class, the block estimator gives a stronger upper rate nα for the same integrated target. Therefore an nαh lower bound cannot hold on any class broad enough to include the primitive Lipschitz models while also allowing the block estimator’s assumptions.

A lower bound of order nαh could arise for a different problem—for example, estimating a local profile value, or detecting a signal localized to an unknown interval of width h. But then the relationship between the local amplitude and the integrated gap must be specified. If a depression of height δ occupies width h, its integrated gap is of order hδ, and the information scales as

nαhδ
2
≍
h
nαΔ
2
	​

,

not nαhΔ
2
.

Campaign targets II.4 and II.5 have not been reached. The present evidence points away from nαhΔ
2
 as the sharp threshold for the average score.

17. Uniform doubling can be replaced by the tilted-deficit condition

Let

B
j,u
	​

(t)=∫e
−tD
j,u
	​

(v)
c(u,v)K
j
	​

(u,dv),

and

M
j,u
	​

(t)=∫D
j,u
	​

(v)e
−tD
j,u
	​

(v)
c(u,v)K
j
	​

(u,dv).

The quantity actually entering the derivative of the projected slowly varying factor is

B
j,u
	​

(t)
M
j,u
	​

(t)
	​

.

Thus the relevant primitive condition is

m
n
	​

(T)=
j,u
sup
	​

t≥T
sup
	​

B
j,u
	​

(t)
M
j,u
	​

(t)
	​

⟶0.
(TD)

Repeating the inversion argument in Appendix C gives

b
n
	​

(3α)≤C[m
n
	​

{c
γ
	​

log(1/α)}+
R
1
	​

{c
γ
	​

log(1/α)}].

This is strictly weaker than uniform doubling. The doubling assumption is one convenient way to obtain

m
n
	​

(T)=O(T
−1
).

Campaign claim 6 survives, provided the tail-derivative remainder from (P2) is retained. The tilted-deficit condition does not by itself control oscillations in the primitive survival remainder.

18. The effective codimension is real and was hidden in the constants

Suppose

H
j,u
	​

(z)∼a
j,u
	​

z
κ
j,n
	​

L
j,u
	​

(1/z)(z↓0).

Standard Laplace regular variation gives

B
j,u
	​

(t)∼a
j,u
	​

Γ(κ
j,n
	​

+1)t
−κ
j,n
	​

L
j,u
	​

(t)

and

B
j,u
	​

(t)
M
j,u
	​

(t)
	​

∼
t
κ
j,n
	​

	​

.

Therefore the projected tail bias is naturally

log(1/α)
κ
j,n
	​

	​

.

Moreover,

H
j,u
	​

(z)
H
j,u
	​

(2z)
	​

⟶2
κ
j,n
	​

.

A uniform doubling constant therefore forces

n,j,u
sup
	​

κ
j,n
	​

<∞.

If s
n
	​

 active coordinates must approach a corner simultaneously, then typically

κ
j,n
	​

≍s
n
	​

−1

for active projections and κ
j,n
	​

≍s
n
	​

 for inactive ones.

The honest bias condition is consequently

log(1/α
n
	​

)
κ
n
	​

	​

=o(Δ
D
	​

),

not merely

log(1/α
n
	​

)
1
	​

=o(Δ
D
	​

).

This does not invalidate the four-active-variable simulations, where κ is bounded. It invalidates any interpretation of Appendix C as automatically covering growing s
n
	​

 with constants independent of s
n
	​

.

Campaign claim 16 survives.

19. The primitive spatial transfer can also be sharpened

The coupling assumption (P4) implies a quantile comparison of the form

∣logQ
j
	​

(s,v)−logQ
j
	​

(s,u)∣≤A(r)+B(r)logs

with

A(r)=2Γ(L
H
	​

+L
D
	​

ΓC
Q
	​

)r+2Γζ
n
	​

,

and

B(r)=2Γ
2
L
D
	​

r.

Substituting this into the Hill-specific freezing bound gives primitive spatial bias

O(hlog
α
1
	​

+ζ
n
	​

),

not

O{hlog(pn)+ζ
n
	​

}.

The hlog(pn) term in the manuscript comes from controlling every relevant quantile by the smallest PIT among all pn transformed observations. That global minimum is not part of the Hill statistic after the two-sided V
(k+1)
	​

 argument.

Campaign claim 14 survives.

20. Audit of the Gaussian-corner lemma and simulation verification

The Gaussian-corner asymptotic has the correct leading form. Writing

r
z
	​

=−Φ
−1
(z),Ω=Σ
−1
,ν=Ω1,κ=1
⊤
Ω1,b=1
⊤
Ωμ,

the transformed Gaussian density satisfies, after u=zx,

z
κ
(
2π
	​

r
z
	​

)
κ−d
e
−br
z
	​

f
U
	​

(zx)z
d
	​

⟶∣Σ∣
−1/2
e
−μ
⊤
Ωμ/2
i=1
∏
d
	​

x
i
ν
i
	​

−1
	​

.

Integrating over the limiting simplex
∑
i
	​

a
i
	​

x
i
	​

≤1 gives

Γ(κ+1)
∏
i
	​

Γ(ν
i
	​

)a
i
−ν
i
	​

	​

	​

,

so the displayed constant is correct.

The assumptions

i
min
	​

ν
i
	​

>0,

uniform positive definiteness, a uniform C
2
 expansion of D, and a uniformly Lipschitz positive weight are enough to justify the result. One splits the rescaled simplex into the region
min
i
	​

x
i
	​

<e
−
r
z
	​

	​

, whose normalized mass is exponentially negligible, and its complement, on which the normal-quantile expansion is uniform and dominated by an integrable Dirichlet density.

The manuscript’s proof states these domination and moving-boundary steps rather than proving them. The lemma’s statement appears correct, but the printed proof is not yet publication-grade.

For the AR(1) model with ρ=1/4, the claimed exponents check:

κ
1
	​

=κ
4
	​

=
15
34
	​

,κ
2
	​

=κ
3
	​

=
15
41
	​

,

and for inactive coordinates

κ
j
	​

=
5
14
	​

+
1−λ
j
2
	​

λ
j
2
	​

	​

≤
15
43
	​

.

The population gaps also check numerically at ε=.05:

model
A1/B1
A2
A3
	​

Δ
j
	​

0.1863953328
0.1073057160
0.1575534638
	​

	​


For B1, the asserted conditional-tail remainder O(y
−2
) is consistent with the generator because γ≤1/2: the first inversion correction is O(y
−1/γ
), hence uniformly O(y
−2
).

What is still missing is a full proof of the uniform (P4) constants for B1’s posterior Gaussian moments. The assertion is plausible—the latent scale is bounded and the relevant Gaussian posterior means shift linearly—but the derivatives and their dimension-uniform constants are not written down.

The Monte Carlo tables, figures, source code, and included table files are absent from the supplied material. None of the reported simulation frequencies, timings, or real-data ranks can be independently checked from this text.

21. Fixed-coordinate inference becomes available for the block estimator

Although the journal did not supply an inference theorem, the block construction immediately gives one under undersmoothing.

Suppose α→0, nα→∞, the block and tail biases are

o{(nα)
−1/2
},

and the block target approximation has the same order. Then, for fixed j,

nα
	​

(
Ψ
j
	​

−Ψ
j
	​

)⟹N(0,σ
j
2
	​

),

where

σ
j
2
	​

=
∣I
ε
	​

∣
2
1
	​

∫
I
ε
	​

	​

ξ
j
	​

(u)
2
du.

Conditionally on the covariate column, this is a Lindeberg theorem for the independent weighted exponential spacings. Randomness of the empirical block anchors contributes order n
−1/2
, which becomes negligible after multiplication by 
nα
	​

 because α→0.

This is only fixed-coordinate inference. It does not provide:

a limit law for the maximum over p;
a uniformly valid test of Δ
j
	​

=0;
an estimator of the unknown global support maximum;
FWER or FDR control.

The delayed-tail lower bound also shows that no uniformly powerful test against 1/logn alternatives can exist under the bare (C1)–(C2) class.

22. Final campaign ledger
Journal claim	Verdict after audit
1. Replace (S) by actual conditional supports	False as stated. Correct only after replacing the full-cube maximum by the essential maximum on supp(U).
2. Remove atomlessness using randomized PIT	Valid.
3. Target D={j:Δ
j
	​

>0}, not automatically A	Valid after making Δ support-intrinsic.
4. Weaken/remove log(pn)=o(nh
2
)	Valid; it can be removed entirely.
5. Remove equicontinuity of ξ
j
	​

	Valid via latent-order-statistic anchors.
6. Replace uniform doubling by tilted deficit	Valid, together with primitive tail-derivative control.
7. Disjoint blocks give 
logp/(nα)
	​

	Valid for an explicitly defined block estimator under local bias and target-approximation conditions.
8. Estimate γ
⋆
 by max
j
	​

Ψ
j
	​

 and recover without s	False without at least one zero-gap coordinate; not fully adaptive even then.
9. 1/logn impossibility under naked (C1)–(C2)	Valid. Does not cover the stronger primitive class.
10. Entire coordinatewise marginal law can be invisible	Valid.
11. Interaction order can equal s	Valid.
12. A nonidentifiable under unrestricted dependence	Valid.
13. Uniform (C2) essential, averaged replacement possible	Valid.
14. Spatial bias improves to hlog(1/α)	Valid under (HQ) and B<
ξ
	​

−b.
15. Global PIT truncation disappears	Valid using two-sided V
(k+1)
	​

 control.
16. Uniform doubling hides effective codimension	Valid.
17. Exact finite theorem with reported numerical constants	Architecture valid; the exact constants in the journal are not established.
18. Simulation verification omitted from first sessions	Correct objection. Numerical results remain unaudited here.
19. Gaussian-corner proof is only a sketch	Correct objection. Statement appears true under the printed uniformity assumptions, but the proof needs completion.
20. Abstract overstates the dimensional regime	Correct objection. Bias and identification conditions remain essential.
Bottom line

The theory that can presently be defended is:

Projection target=esssup{γ(U)∣U
j
	​

},
	​


with the global reference taken on the observed support, and the identifiable screened set

D
U
	​

={j:γ
U
⋆
	​

−
∣I
ε
	​

∣
1
	​

∫
I
ε
	​

	​

ξ
j
	​

(u)du>0}.
	​


For the printed sliding-window estimator, a stronger proof gives

j
max
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣=O
P
	​

(
nαh
log(pn)
	​

	​

)+O[A(Ch)+B(Ch)log
α
1
	​

+b
n
	​

],
	​


without atomlessness, equicontinuity, log(pn)=o(nh
2
), or global PIT truncation.

For the disjoint-block estimator,

j
max
	​

∣
Ψ
j
	​

−Ψ
j
	​

∣=O
P
	​

(
nα
logp
	​

	​

+
nα
logp
	​

)+bias,
	​


subject to the local feasibility condition
log(p/h)=o(nαh).

Therefore the manuscript’s current stochastic rate is a valid upper bound, but not the strongest available one, and it has not been shown minimax. The two remaining foundational tasks are a properly specified minimax experiment and a genuinely adaptive tuning/selection rule; neither campaign session completed them.