Texte collé.txt
Document
Here is the complete LaTeX source of a statistics manuscript (marginal screening for the conditional tail index in high dimension). I want the mathematics pushed further: stronger conclusions, under weaker assumptions. Go wherever you think the real theory is; give actual proofs, not a plan. Take the time you need.

Mathematical verdict

The manuscript’s next real advance is not another refinement of the existing uniform-profile argument. The current proof controls the maximum error of every local Hill profile over all coordinates and all window states, then bounds the score average by that maximum. This produces the stochastic term

nαh
log(pn)
	​

	​

,

and the main screening theorem additionally assumes that every structurally active coordinate has a positive marginal gap.

I developed a replacement theory with six concrete improvements:

full fibre support is replaced by an exact, substantially weaker argmax-accessibility condition;
continuity and atomlessness of the projected conditional laws are removed completely;
the theorem recovers the identifiable positive-gap set, without assuming every active variable is marginally visible;
the implemented overlapping-window estimator gets a sharper finite-sample theorem and a weaker rank-displacement condition;

a new disjoint-rank-block score attains stochastic score error of order

nα
logp
	​

	​

,

while requiring only averaged spatial and tail regularity;

there are two impossibility theorems, including a full-data lower bound showing that 1/logn signals are not uniformly detectable under unquantified (C1)–(C2).

The complete theorem statements and proofs are in the files at the end.

1. Full fibre support is far stronger than identification requires

The manuscript currently assumes that every conditional fibre has support equal to the entire cube, so that the projected essential supremum becomes a maximum over the full fibre. That can be replaced exactly.

Let

C
j
	​

(u)=suppK
j
	​

(u,⋅)

be the actual conditional support of U
−j
	​

∣U
j
	​

=u. Under (C1)–(C2) and continuity of γ,

ξ
j
	​

(u)=
v∼K
j
	​

(u,⋅)
esssup
	​

γ{ι
j
	​

(u,v)}=
v∈C
j
	​

(u)
max
	​

γ{ι
j
	​

(u,v)}.

The proof is short but decisive. The support is compact, so the maximum m is attained at some v
⋆
	​

∈C
j
	​

(u). Every neighbourhood of v
⋆
	​

 has positive conditional probability, and continuity gives

γ{ι
j
	​

(u,v)}>m−η

throughout a sufficiently small neighbourhood. Therefore the conditional essential supremum is at least m−η for every η>0, hence equals m.

Define

H
j
	​

={u∈I
ε
	​

:({u}×C
j
	​

(u))∩M

=∅},M=argmaxγ.

Then

ξ
j
	​

(u)=γ
⋆
⟺u∈H
j
	​

,

and, without requiring continuity of u↦ξ
j
	​

(u),

Δ
j
	​

>0⟺λ(I
ε
	​

∖H
j
	​

)>0.

This gives the minimal null-identification assumption:

λ(I
ε
	​

∖H
j
	​

)=0for every j∈
/
A.
(NA)

Under (NA), no inactive coordinate has a positive gap. Full support implies (NA), but is much stronger.

There is also an important negative interpretation. If (NA) fails, an inactive variable may restrict the conditional support of the active coordinates so that high-index configurations disappear for some of its values. Its projected tail index then genuinely varies. This is not a false positive caused by estimation: the marginal law of (U
j
	​

,Y) itself contains a tail-index signal.

The exact population target is therefore

D={j:Δ
j
	​

>0}.

Under (NA), D⊆A, and

D=A

if and only if every active coordinate misses the accessible argmax on a set of positive Lebesgue measure.

2. Atomlessness is unnecessary

The continuity clause in (E1) is currently used to manufacture iid conditional uniforms and a quantile representation. The manuscript consequently treats atoms as an obstruction that must be repaired by (P1). This obstruction is only a proof artefact.

Introduce independent auxiliary variables Z
ij
	​

∼U(0,1), used only on the proof probability space, and define the randomized conditional distributional transform

T
ij
	​

=F
j
	​

(Y
i
−
	​

∣U
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
−
	​

∣U
ij
	​

)},V
ij
	​

=1−T
ij
	​

.

Then, for each fixed coordinate j, conditionally on the complete column

U
j
	​

=(U
1j
	​

,…,U
nj
	​

),

the V
ij
	​

 are iid standard uniforms and

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

)almost surely.

For completeness, fix u, write F=F
j
	​

(⋅∣u), and let q=F
−1
(t). Since F(q
−
)≤t≤F(q),

P(T
ij
	​

≤t∣U
ij
	​

=u)
	​

=F(q
−
)+{F(q)−F(q
−
)}
F(q)−F(q
−
)
t−F(q
−
)
	​

=t,
	​


with the zero-jump case interpreted in the evident way. The generalized inverse of the randomized transform equals Y
i
	​

 almost surely.

Thus all the Rènyi-spacing and order-statistic arguments remain valid with arbitrary projected conditional distributions. The estimator is unchanged; only the proof is enlarged by auxiliary uniforms.

The atomlessness clause can therefore be deleted from (P1), and the continuity clause can be deleted from the structural part of (E1).

3. The implemented score admits a sharper theorem

The existing estimator can be handled without assuming that every active coordinate is detectable and without imposing equicontinuity merely to control numerical integration.

Define the actual quadrature error

q
n
∘
	​

=
j≤p
max
	​

	​

∣G∣
1
	​

u∈G
∑
	​

ξ
j
	​

(u)−Ψ
j
	​

	​

.

This is the quantity the proof needs. Requiring

q
n
∘
	​

=o(Δ
D
	​

)

is strictly weaker than requiring a uniform modulus of continuity for every profile.

For

e
n
	​

(η,t,τ)=
R
n
	​

(η,t,τ)=
	​

2ω
n
↑
	​

(h+η+δ
n
	​

,τ,α)+Γt+(1+t)b
n
	​

(3α)+q
n
∘
	​

,
2pe
−2nη
2
+pnτ+2pS
n
	​

e
−K
n
	​

t
2
/4
+pS
n
	​

e
−(K
n
	​

+1)/4
,
	​


the exact finite-sample bound is

P(
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

∣>e
n
	​

(η,t,τ))≤R
n
	​

(η,t,τ).
(1)

The randomized PIT is the only new ingredient needed in the existing profile proof. Everything else follows from the manuscript’s sorting, tail-freezing and exponential-spacing arguments.

On the event in which the score error is at most e
n
	​

,

Ψ
j
	​

−
Ψ
k
	​

≤−(Δ
j
	​

−Δ
k
	​

)+2e
n
	​

.

Consequently,

Δ
j
	​

−Δ
k
	​

>2e
n
	​

⟹
Ψ
j
	​

<
Ψ
k
	​

.
(2)

In particular, when e
n
	​

<Δ
D
	​

/2, every coordinate in D precedes every coordinate in D
c
. Therefore

P(D⊆
A
d
	​

)≥1−R
n
	​

,d≥∣D∣,

and

P(
A
∣D∣
	​

=D)≥1−R
n
	​

.

This is stronger than the present theorem because invisible active coordinates are no longer assumed away. Under (NA), D⊆A; under the additional active detectability condition, D=A.

Sharper rank condition

The DKW term only involves a union over the p coordinate columns, not over the n window states. Hence the first rate condition can be weakened from

log(pn)=o(nh
2
)

to

log(2p)=o(nh
2
).
(3)

The logn factor remains in the exponential-spacing term because the overlapping-window proof still unions over O(n) states. Thus the asymptotic conditions for the implemented estimator become

log(2p)=o(nh
2
),log(pn)=o(nαh),

and

log(pn)=o(nαhΔ
D
2
	​

),

together with the direct bias and quadrature conditions.

Discontinuous profiles are allowed

Equicontinuity can be replaced by bounded variation. If

j≤p
max
	​

TV
I
ε
	​

	​

(ξ
j
	​

)≤V
n
	​

,

then

q
n
∘
	​

≤C
ε
	​

n
V
n
	​

+Γ
	​

.

Thus profiles with jumps are covered. The proof is a bounded-variation midpoint/Voronoi quadrature bound.

4. Exact recovery without knowing s

Suppose at least one zero-gap coordinate exists. Define

γ
	​

max
	​

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
	​

−
Ψ
j
	​

.

On the event

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

∣≤e,

one has

Δ
j
	​

≤2e,j∈D
c
,

whereas

Δ
j
	​

≥Δ
j
	​

−2e,j∈D.

Therefore every threshold satisfying

2e<λ<Δ
D
	​

−2e

gives

{j:
Δ
j
	​

>λ}=D.

This removes the requirement d=s from exact recovery. A deterministic high-probability upper bound from (1) can be used for e.

5. A new score gains the averaging factor

The main positive improvement comes from changing the score estimator, not merely sharpening constants.

Partition the empirical ranks in I
ε
	​

 into B consecutive, disjoint blocks of m ranks. For each coordinate and block, compute one Hill statistic H
jb
	​

 using

k=⌊αm⌋

upper observations, and set

Ψ
j
	​

=
B
1
	​

b=1
∑
B
	​

H
jb
	​

.

Typically,

m≍nh,B≍h
−1
,k≍nαh,

so

Bk≍nα.

Conditional on the j-th covariate column, the blocks contain disjoint sets of observations. Their Rènyi terms

L
jb
	​

=
k
1
	​

r=1
∑
k
	​

log
V
jb,(r)
	​

V
jb,(k+1)
	​

	​


are therefore independent means of k standard exponentials.

Writing

Z
j
	​

=
B
1
	​

b=1
∑
B
	​

ξ
j
	​

(u
b
	​

)(L
jb
	​

−1),

the exact moment-generating function calculation gives, for

V=
k
1
	​

b
∑
	​

a
b
2
	​

,c=
k
1
	​

b
max
	​

a
b
	​

,
logEe
θ∑
b
	​

a
b
	​

(L
b
	​

−1)
≤
2(1−cθ)
θ
2
V
	​

.

With a
b
	​

=ξ
j
	​

(u
b
	​

)/B,

V≤
Bk
Γ
2
	​

,c≤
Bk
Γ
	​

,

and hence

P(∣Z
j
	​

∣>Γ{
Bk
2x
	​

	​

+
Bk
x
	​

})≤2e
−x
.
(4)

Define blockwise moduli b
j,b
	​

 and ω
j,b
	​

, and only require their averages

b
ˉ
n
	​

=
j
max
	​

B
1
	​

b
∑
	​

b
j,b
	​

,
ω
ˉ
n
	​

=
j
max
	​

B
1
	​

b
∑
	​

ω
j,b
	​

.

This is weaker than a supremum over every u.

The resulting finite-sample bound is

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

∣>2
ω
ˉ
n
	​

+2
b
ˉ
n
	​

+q
n,B
	​

+Γ{
Bk
2x
	​

	​

+
Bk
x
	​

}]
≤2pe
−2nη
2
+pnτ+3pBe
−k/4
+2pe
−x
.
	​

(5)

Consequently,

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
nα
log(2p)
	​

	​

+
nα
log(2p)
	​

)+2
ω
ˉ
n
	​

+2
b
ˉ
n
	​

+q
n,B
	​

.
(6)

The weak-signal screening condition becomes

log(2p)=o(nαΔ
D
2
	​

),
(7)

rather than

log(pn)=o(nαhΔ
min
2
	​

).

There remains a separate local feasibility condition

log(p/h)=o(nαh),

because every individual block still needs enough extremes for tail freezing and Hill control. Thus the new theorem improves score variance and weak-gap detectability; it does not pretend that a local Hill statistic can be formed without a diverging local tail count.

This block estimator is new. The manuscript’s current simulations implement overlapping windows, so the nα rate must not be attributed to those numerical results without rerunning them.

6. The deficit-mass doubling condition is unnecessary

The present primitive theory assumes uniform doubling of the deficit mass. That condition is useful for obtaining an explicit O(1/t) rate, but the transfer proof only needs concentration of the exponentially tilted deficit law.

Let

B
j,u
	​

(t)
M
j,u
	​

(t)
	​

=∫e
−tD
j,u
	​

(v)
c{ι
j
	​

(u,v)}K
j
	​

(u,dv),
=∫D
j,u
	​

(v)e
−tD
j,u
	​

(v)
c{ι
j
	​

(u,v)}K
j
	​

(u,dv),
	​


and define

χ
n
	​

(T)=
j≤p
n
	​

max
	​

u∈I
ε/2
	​

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

.

Replace (P3) by

T→∞
lim
	​

n
sup
	​

χ
n
	​

(T)=0.
(P3-L)

For the projected survival factor

A
j,u
	​

(t)=∫e
−tD
j,u
	​

(v)
c{ι
j
	​

(u,v)}{1+r(t,ι
j
	​

(u,v))}K
j
	​

(u,dv),

one obtains directly

	​

A
j,u
	​

(t)
A
j,u
′
	​

(t)
	​

	​

≤3χ
n
	​

(t)+2R
1
	​

(t).
(8)

Moreover,

−(logB
j,u
	​

)
′
(t)=
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

Condition (P3-L) therefore implies

∣logA
j,u
	​

(t)∣=o(t)

uniformly. Quantile inversion gives

T
j,u
	​

(x)=logQ
j
	​

(e
x
,u)≥
2
γ
	​

	​

x

eventually, and differentiation yields

dx
d
	​

logℓ
j
	​

(e
x
,u)=
1−ξ
j
	​

(u)a
j,u
	​

{T
j,u
	​

(x)}
ξ
j
	​

(u)
2
a
j,u
	​

{T
j,u
	​

(x)}
	​

,a
j,u
	​

=
A
j,u
	​

A
j,u
′
	​

	​

.

Combining this with (8),

b
n
	​

(3α
n
	​

)≤6Γ
2
χ
n
	​

(
2
γ
	​

	​

log
3α
n
	​

1
	​

)+4Γ
2
R
1
	​

(
2
γ
	​

	​

log
3α
n
	​

1
	​

).
(9)

The current doubling condition implies χ
n
	​

(T)=O(T
−1
), so none of the simulation verification is lost. But doubling is no longer presented as necessary.

A still weaker sufficient condition follows from

χ
n
	​

(T)≤z+
inf
j,u
	​

H
j,u
	​

(z/2)
D
max
	​

c
+
	​

e
−Tz/2
	​

.

Thus it is enough that every fixed neighbourhood of zero deficit carry uniformly positive fibre mass. No polynomial margin or local doubling shape is required.

7. The 1/logn barrier is information-theoretic under bare (C1)–(C2)

The manuscript correctly identifies the status of order-1/logn gaps as a central unresolved issue. Over the unquantified (C1)–(C2) class, there is a sharp negative answer.

Take p=2, independent uniform U
1
	​

,U
2
	​

, and

T
n
	​

=n
2γ
0
	​

,δ
n
	​

=
logT
n
	​

a
	​

=
2γ
0
	​

logn
a
	​

.

Under the null,

γ
0,n
	​

(u)=γ
0
	​

,
F
ˉ
0,n
	​

(y∣u)=y
−1/γ
0
	​

.

Under the alternative,

γ
1,n
	​

(u)=γ
0
	​

+δ
n
	​

u
1
	​


and

F
ˉ
1,n
	​

(y∣u)={
y
−1/γ
0
	​

,
T
n
−1/γ
0
	​

	​

(y/T
n
	​

)
−1/γ
1,n
	​

(u)
,
	​

1≤y≤T
n
	​

,
y>T
n
	​

.
	​

(10)

Beyond T
n
	​

,

F
ˉ
1,n
	​

(y∣u)=c
n
	​

(u)y
−1/γ
1,n
	​

(u)
,c
n
	​

(u)=T
n
1/γ
1,n
	​

(u)−1/γ
0
	​

	​

,

and

e
−a/γ
0
2
	​

≤c
n
	​

(u)≤1.

Thus (C1)–(C2) hold at every n, with uniformly bounded positive scale constants, but with no common onset rate.

The active coordinate has gap

Δ
1,n
	​

=
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

Yet the null and alternative laws agree exactly below T
n
	​

, and

P(Y>T
n
	​

∣U)=T
n
−1/γ
0
	​

	​

=n
−2
.

Therefore

	​

P
0,n
⊗n
	​

−P
1,n
⊗n
	​

	​

TV
	​

≤n
−1
.

For every screening procedure based on the complete data,

P
0,n
	​

(1∈
A
)+P
1,n
	​

(1∈
/
A
)≥1−n
−1
.
(11)

Hence no estimator—not Hill, not a marginal estimator, and not a full multivariate procedure—can uniformly recover order-1/logn gaps over bare (C1)–(C2). A quantified common onset such as (P2), or a comparable second-order restriction, is information-theoretically necessary.

This lower bound does not rule out order-1/logn detection under (P2) with an explicit bias correction. It settles the unquantified primitive class.

8. Marginal invisibility is also exact at the level of the whole marginal law

There is a second, orthogonal impossibility theorem.

Let

g(t)=γ
0
	​

+acos(2πt),0<a<γ
0
	​

.

Compare

γ
0
	​

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
1
	​

(u
1
	​

,u
2
	​

)=g{(u
1
	​

+u
2
	​

)mod1},

using independent uniform covariates and exact conditional Pareto tails.

Coordinate 1 is inactive in the first model and active in the second. Nevertheless, for every u
1
	​

,

P
1
	​

(Y>y∣U
1
	​

=u
1
	​

)
	​

=∫
0
1
	​

y
−1/g{(u
1
	​

+v)mod1}
dv
=∫
0
1
	​

y
−1/g(v)
dv
=P
0
	​

(Y>y∣U
1
	​

=u
1
	​

).
	​


Thus the complete law of (U
1
	​

,Y) is identical in the two models.

Every decision rule using only the j-th coordinate and Y consequently has worst-case error at least 1/2, for every sample size. This proves that the argmax-projection limitation is not specific to the average score—or even to the tail index. It can persist in the entire coordinatewise marginal experiment.

9. Exact group order

For a coordinate set J, define

ξ
J
	​

(u
J
	​

)=
v
max
	​

γ(u
J
	​

,v),Δ
J
	​

=γ
⋆
−∫ξ
J
	​

(u
J
	​

)du
J
	​

.

Under full support,

Δ
J
	​

>0⟺λ
∣J∣
	​

([0,1]
∣J∣
∖π
J
	​

(M))>0.

The minimal interaction order required to see active coordinate j is

r
⋆
	​

(j)=min{∣J∣:j∈J, Δ
J
	​

>0}.

If s=∣A∣, then

r
⋆
	​

(j)≤s.

This bound is sharp. For

γ(u)=g{(u
1
	​

+⋯+u
s
	​

)mod1},

every projection of the argmax onto fewer than s active coordinates fills the entire projected cube, so

r
⋆
	​

(j)=s

for every active coordinate.

Thus no universally valid fixed-order pairwise or low-order group screen can remove marginal invisibility.

Deliverables

The insertable patch contains all theorem statements and complete proofs:

Insertable LaTeX theory patch
Standalone compilable LaTeX addendum
Compiled 15-page PDF

The cleanest manuscript revision is to make D, rather than A, the unconditional population target; recover A only under null accessibility plus active detectability; delete atomlessness through the randomized PIT; retain the present overlapping estimator as the main implemented method; and present the disjoint-block score as the sharper theoretical variant whose nα score rate requires separate numerical evaluation.

sharper_theory_addendum.pdf
PDF