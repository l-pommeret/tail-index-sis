ttached: a statistics manuscript (LaTeX) followed by the campaign journal of a first round of work on it. The journal records results two independent sessions reached; none of them has been checked by anyone. Your job is to break them. Take every claim as suspect until you have re-proved it yourself: hidden assumptions, circular steps, constants that do not survive, counterexamples. Where a claim survives, say what exactly it survives under. Where it dies, kill it precisely.
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsthm,mathtools,bm,booktabs,graphicx}
\usepackage[numbers,sort&compress]{natbib}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\usepackage[section]{placeins}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{lemma}[theorem]{Lemma}
\theoremstyle{remark}\newtheorem{remark}[theorem]{Remark}
\DeclareMathOperator*{\esssup}{ess\,sup}
\DeclareMathOperator*{\argmax}{arg\,max}
\newcommand{\A}{A}\newcommand{\I}{I_\varepsilon}
\newcommand{\M}{\mathcal{M}}\newcommand{\Pp}{\mathbb P}
\newcommand{\E}{\mathbb E}\newcommand{\R}{\mathbb R}
\newcommand{\gmax}{\gamma^\star}
\title{Tail-Index Sure Independence Screening for\\High-Dimensional Heavy-Tailed Regression}
\author{Anonymous manuscript}\date{}
\begin{document}\maketitle
\begin{abstract}
We propose a marginal screen for the covariates that determine a positive
conditional tail index, designed for problems where the dimension $p$ is large
relative to the number of observations informative about the tail. The screen
rests on one geometric fact: after projection onto a single coordinate, the
conditional distribution is a mixture whose tail is governed by its heaviest
component, so the projected tail index is the upper envelope of the tail-index
surface along the corresponding fibre. Inactive coordinates produce flat
envelopes at the global maximum $\gmax$, while an active coordinate whose
fibres miss the global argmax set depresses the envelope; ranking coordinates
by the average projected index therefore separates active variables from
noise, and we characterize exactly which active coordinates a marginal screen
can see, rather than assuming a generic signal-strength condition. We estimate
the envelope with empirical-rank local Hill profiles. Ranks give every
coordinate the same deterministic window geometry, which removes the covariate
density conditions that fixed-dimension projected-tail theory requires and
permits exact finite-sample union bounds over a diverging number of
coordinates: the maximal score error over $p$ coordinates is
$O_{\Pp}(\sqrt{\log(pn)/(n\alpha h)})$ plus a bias term of order
$o(\Delta_{\min})$, where $\alpha$ is the local tail fraction, $h$ the rank
bandwidth, $n\alpha h$ the effective local extreme sample size, and
$\Delta_{\min}$ the weakest population gap. Sure screening and exact recovery follow when
$\log(pn)=o(n\alpha h\Delta_{\min}^2)$, so a polynomially growing dimension is
affordable whenever the local extreme count dominates the logarithm of the
dimension. The algorithm costs $O(pn\log n)$ time and $O(n)$ working memory,
under a second per screen at $(n,p)=(2000,1000)$. A rank-aggregated version,
taking the minimum rank over a block of nine neighbouring tunings, stabilizes
the ranking against the choice of the smoothing parameters. Simulations across
$p\in\{500,1000,2000\}$, including a model whose scale is driven by a
correlated cluster of twenty inactive covariates, show when screening on
the tail index differs from screening on a fixed high quantile: the
cluster fills the leading positions of moderate-level quantile screens,
whose contamination clears at extreme levels only at a severe variance
cost, while the leading positions of the aggregated tail screen remain
on the active set. An application to violent-crime rates in 1{,}993
U.S. communities, with $p=99$ continuous predictors against an
effective tail sample of about $200$, illustrates the method in the
regime it is designed for, and separates predictors of the limiting
exponent from predictors of a high quantile on real data. The uniform
projected-quantile conditions used by the theorems are derived from
primitive assumptions on the tail-index function, the conditional tails
and the fibre kernels, and verified for the simulation models.
\end{abstract}
\input{sections-v2/introduction}
\input{sections-v2/model_population}
\input{sections-v2/estimation_theory}
\input{sections-v2/computation_simulation}
\FloatBarrier
\input{sections-v2/real_data}
\FloatBarrier
\input{sections-v2/discussion}
\appendix\input{sections-v2/appendix_aux}\input{sections-v2/appendix_proofs}
\input{sections-v2/appendix_primitive}
\bibliographystyle{plainnat}
\bibliography{references}
\end{document}
\section{Introduction}\label{sec:intro}

Let $(X,Y)$ be a random pair with $Y>0$ and $X\in\R^p$, and suppose that,
conditionally on $X=x$,
\begin{equation}\label{eq:intro-rv}
 \frac{\Pp(Y>ty\mid X=x)}{\Pp(Y>y\mid X=x)}
 \longrightarrow t^{-1/\gamma(x)},\qquad y\to\infty,\quad t>0,
\end{equation}
for a positive function $\gamma$, the conditional tail index. Estimating
$\gamma$ is the first step toward extreme conditional quantiles and tail
probabilities \citep{deHaanFerreira2006,hill1975}, and it suffers from a
particularly severe curse of dimensionality: local estimation must find
observations near $x$ and then keep only an upper fraction of their
responses. A one-dimensional window of width $h$ combined with a tail
fraction $\alpha$ leaves an effective sample of order $n\alpha h$, not $n$
\citep{daouiaEtAl2013,gardesStupfler2014}. When $p$ is large, even moderate
sample sizes are effectively small.

This paper asks the screening question: among $p$ candidate covariates,
possibly with $p$ growing with $n$, which ones does $\gamma$ actually depend
on? Sure independence screening \citep{fanLv2008} answers such questions by
ranking coordinates through inexpensive marginal statistics and keeping a
moderate number for subsequent joint analysis. Model-free variants exist for
general dependence \citep{liZhongZhu2012} and for fixed conditional quantiles
\citep{heWangHong2013}, but the limiting exponent in \eqref{eq:intro-rv} is a
different target: a covariate can move the $0.95$ quantile of $Y$ without
affecting $\gamma$, and conversely. Our simulations make this distinction
concrete.

The screen we propose rests on a single geometric mechanism. Transform each
continuous covariate to $U_j=F_j(X_j)$ and let $\xi_j(u)$ denote the tail
index of $Y$ given $U_j=u$ alone. Conditioning on one coordinate leaves the
others free to range over the fibre $\{u\}\times[0,1]^{p-1}$, so the
projected distribution is a mixture, and a mixture of heavy tails is governed
by its heaviest component. Under a full-support condition,
\begin{equation}\label{eq:intro-envelope}
 \xi_j(u)=\max_{v\in[0,1]^{p-1}}\gamma(u,v):
\end{equation}
the projected index is the upper envelope of the tail-index surface along the
fibre. If coordinate $j$ is inactive, every fibre reaches the global maximum
$\gmax=\max_x\gamma(x)$ and the envelope is flat at $\gmax$. If coordinate
$j$ is active and some of its fibres miss the global argmax set, the envelope
dips below $\gmax$ there. Averaging the envelope over a trimmed interval
gives a score $\Psi_j$ whose gap $\Delta_j=\gmax-\Psi_j$ is zero for every
inactive coordinate and positive for every detectable active one. Ranking
coordinates by increasing $\Psi_j$ is the screen.

The same geometry states the method's limitation. An active coordinate is
invisible when every one of its fibres still meets the global argmax set;
this happens, for instance, for $\gamma(u_1,u_2)=\gamma_0-(u_1-u_2)^2$,
whose argmax is the diagonal. Detectability is a property of the argmax
projections, and we characterize it exactly rather than assuming a generic
signal-strength condition.

Estimation follows the envelope. For each coordinate we slide a window over
the empirical marginal ranks, apply a Hill statistic to the upper
$\alpha$-fraction of the responses in the window, and average the resulting
profile. Empirical ranks give every coordinate the same deterministic window
geometry and remove the marginal scales. The price is a rank-displacement
term, which we control simultaneously over $p$ coordinates by a
Dvoretzky--Kiefer--Wolfowitz bound; it is negligible when $\log(pn)=o(nh^2)$.
The main result bounds the maximal score error by
\[
 O_{\Pp}\!\left(\sqrt{\frac{\log(pn)}{n\alpha h}}\right)+o(\Delta_{\min}),
\]
where $\Delta_{\min}$ is the weakest population gap and the second term
collects the tail and spatial biases, and yields sure screening when
$\log(pn)=o(n\alpha h\Delta_{\min}^2)$: high dimension is affordable
exactly when the local extreme count dominates the logarithm of the
dimension, so a polynomially growing $p$ is admissible at a fixed gap.
A blockwise implementation costs $O(pn\log n)$ time and $O(n)$ working
memory beyond the scores.

Our construction adapts the projected-tail framework of
\citet{gardesPodgorny2025}, who estimate a low-dimensional central tail-index
subspace by minimizing an average of projected local Hill estimates over
projection matrices, in fixed dimension. We keep their projection semantics
and local Hill estimator, and change what is needed for screening: the
projections are the $p$ coordinates themselves, $p$ may diverge, windows are
built from empirical ranks with deterministic local counts, and uniformity
over a compact matrix class is replaced by union bounds over coordinates and
window states. Their density assumptions on the covariate are not needed
here; the empirical-rank construction replaces them. Among screening
proposals, \citet{yoshidaUmezu2026} screen on covariate-dependent
extreme-value indices through kernel conditional Pickands contrasts in a
large-bandwidth regime, and parametric tail-index regression \citep{wangTsai2009}, extended to
high dimensions with regularization by \citet{sasakiTaoWang2026},
exploits a parametric link when one is available. Tail dimension reduction methods
\citep{gardes2018,aghbalouEtAl2024,bousebataEtAl2023,arbelEtAl2024} target
related but distinct conditional-tail objects.

The theoretical scope deserves a clear statement. The primitive model
\eqref{eq:intro-rv} identifies the projected index
\eqref{eq:intro-envelope}, while the uniform regularity needed for
growing-dimensional concentration concerns the projected quantiles
themselves, because projection creates logarithmic slow variation even
from exact Pareto conditionals. We state the required conditions where they
act, on the projected quantiles, and then derive them in
Appendix~\ref{app:primitive} from primitive assumptions on $\gamma$, on
the conditional tails and on the fibre kernels, after weakening the
spatial condition to the range of quantiles the proofs actually use. The
same appendix verifies those assumptions for the simulation models.

Section~\ref{sec:population} develops the population theory.
Section~\ref{sec:estimation} defines the estimator and states the screening
results. Section~\ref{sec:simulation} presents the simulation models, the tuning
study together with a rank-aggregated version of the screen built on a
block of nine neighbouring tunings, and a comparison with competing
screens at $p\in\{500,1000,2000\}$ that covers, in the same tables, the
sensitivity to the dimension and to the quantile level of fixed-quantile
screening. Section~\ref{sec:realdata} applies the aggregated screen to
violent-crime rates in U.S. communities.
Section~\ref{sec:discussion} discusses scope and practice. Proofs are in
Appendix~\ref{app:proofs}, supported by the auxiliary results of
Appendix~\ref{app:aux}.
\section{Population theory}\label{sec:population}

\subsection{Model and assumptions}

All quantities below may depend on $n$: the dimension $p=p_n$, the
tail-index function $\gamma$, the active set $\A$, its size $s=|\A|$, and
later the bandwidth $h$, the tail fraction $\alpha$ and the retained size
$d$. The index $n$ is suppressed whenever no ambiguity arises.

\emph{Notation.} A positive measurable function $f$ is slowly varying
at infinity if $f(ty)/f(y)\to1$ as $y\to\infty$ for every $t>0$, and a
survival function $\bar F$ is regularly varying with index $-1/\xi<0$
if $\bar F(ty)/\bar F(y)\to t^{-1/\xi}$ for every $t>0$. Generalized
inverses are $F^{-1}(t)=\inf\{y:F(y)\ge t\}$. For $u\in[0,1]$ and
$v\in[0,1]^{p-1}$, $\iota_j(u,v)\in[0,1]^p$ denotes the vector with
$u$ in coordinate $j$ and $v$ in the remaining coordinates, in order;
we abbreviate $\gamma(u,v)=\gamma\{\iota_j(u,v)\}$ when $j$ is clear
from context. Norms on $[0,1]^{\A}$ are Euclidean. All asymptotic
statements concern a triangular array indexed by $n$, and
$Z_n=O_{\Pp}(a_n)$ means $\lim_{C\to\infty}\limsup_n
\Pp(|Z_n|>Ca_n)=0$.

The observed data are i.i.d.\ pairs $(X_i,Y_i)$ with $Y_i>0$ and
$X_i\in\R^p$ with continuous marginal distribution functions
$F_1,\ldots,F_p$. The theory is formulated on the latent
probability-integral transforms $U_{ij}=F_j(X_{ij})$, so that
$U_i\in[0,1]^p$ has uniform one-dimensional marginals and unrestricted
cross-coordinate dependence; the transforms are never needed by the
algorithm, because the ranks of $X_{1j},\ldots,X_{nj}$ coincide with the
ranks of $U_{1j},\ldots,U_{nj}$ almost surely. With slight abuse of
notation we write $\gamma$ for the tail-index function on the cube,
$\gamma(u)=\gamma_X\{F_1^{-1}(u_1),\ldots,F_p^{-1}(u_p)\}$, where
$\gamma_X$ is the function appearing in \eqref{eq:intro-rv}.

\begin{itemize}
\item[(C1)] \emph{Conditional regular variation.} For every $u$, the
conditional survival function factorizes as
\[
 \bar F(y\mid u)=\Pp(Y>y\mid U=u)=y^{-1/\gamma(u)}L(y,u),
\]
where $L(\cdot,u)$ is slowly varying at infinity,
$\gamma:[0,1]^{p}\to(0,\infty)$ is continuous, and
$0<\underline\gamma\le\gamma(u)\le\Gamma<\infty$ uniformly in $n$ and
$u$.
\item[(C2)] \emph{Uniform tail regularity.} There is a function $c$ with
$0<c_-\le c(u)\le c_+<\infty$ such that
$\sup_u|L(y,u)/c(u)-1|\to0$ as $y\to\infty$.
\end{itemize}
(C1) states that each conditional law is heavy-tailed with index
$\gamma(u)$, as in \eqref{eq:intro-rv}; (C2) strengthens the slow
variation of $L$ to uniform convergence toward a bounded positive limit,
the form needed to control mixtures of the conditional laws.

The active set $\A$ collects the coordinates on which $\gamma$
genuinely depends: $j\in\A$ if and only if there exist $u,v$ differing
only in coordinate $j$ with $\gamma(u)\ne\gamma(v)$. Iterating over the
complementary coordinates shows that $\gamma(u)=g(u_{\A})$ for some
function $g$, and that $\A$ is the smallest coordinate set with this
property. Membership in $\A$ concerns the limiting exponent
only, so a coordinate may shift the body of the conditional distribution, or
its second-order behaviour, without being active. Conditions (C1)--(C2)
parallel conditions (C.1)--(C.2) of \citet{gardesPodgorny2025}, except that
we do not assume the law of $\gamma(U)$ to be absolutely continuous; the
projection argument below handles atoms, which arise as soon as one
coordinate is inactive.

Fix a trimming parameter $0<\varepsilon<1/2$ and write
$\I=[\varepsilon,1-\varepsilon]$. For each coordinate $j$, fix a
pointwise conditional distribution $K_j(u,dv)$ of $U_{-j}$ given
$U_j=u$, and define the conditional law of $Y$ given $U_j=u$ to be
compatible with it through the disintegration identity
\[
 \Pp(Y>y\mid U_j=u)=\int\Pp\{Y>y\mid U=(u,v)\}\,K_j(u,dv);
\]
without such pointwise versions the statements below hold for almost
every $u$. Define
\[
 \xi_j(u)=\esssup\nolimits_{v\sim K_j(u,\cdot)}\gamma(u,v),
\]
the essential supremum of $v\mapsto\gamma(u,v)$ under $K_j(u,\cdot)$;
Proposition~\ref{prop:projection} shows that the conditional law of $Y$
given $U_j=u$ is regularly varying with index $-1/\xi_j(u)$, so
$\xi_j(u)$ is its tail index. One further condition is used only where
stated:

\begin{itemize}
\item[(S)] \emph{Full fibre support.} For every $j$ and $u\in\I$, the
support of $K_j(u,\cdot)$ is all of $[0,1]^{p-1}$. A continuous joint
density, positive on the interior cube, is sufficient.
\end{itemize}

\subsection{Projection turns the surface into an envelope}

Conditioning on $U_j=u$ mixes the conditional laws of $Y$ over the
fibre $\{u\}\times[0,1]^{p-1}$. Heuristically, the slowest-decaying
component of a mixture of heavy tails dictates the tail of the mixture,
so the projected index should equal the largest tail index available on
the fibre, and a localized region of heavy tails should survive
projection intact. The following proposition makes this precise.

For a sub-$\sigma$-field $\mathcal G$, the conditional essential
supremum $\esssup\{\gamma(U)\mid\mathcal G\}$ is the smallest
$\mathcal G$-measurable random variable dominating $\gamma(U)$ almost
surely.

\begin{proposition}[Projection]\label{prop:projection}
Assume (C1)--(C2) and let $\mathcal G\subseteq\sigma(U)$. Then
$\Pp(Y>y\mid\mathcal G)$ is almost surely regularly varying with index
$-1/\xi$, where $\xi=\esssup\{\gamma(U)\mid\mathcal G\}$. In particular
$\xi_j(u)$ is the conditional essential supremum of $\gamma(U)$ given
$U_j=u$. Under (S),
\begin{equation}\label{eq:fibre}
 \xi_j(u)=\max_{v\in[0,1]^{p-1}}\gamma(u,v),\qquad u\in\I.
\end{equation}
\end{proposition}

The essential-supremum identity requires neither full support nor absolute
continuity; (S) enters only to replace the essential supremum by the
maximum over the whole fibre. If dependence removes part of a fibre, the
maximum runs over the conditional support actually available at $u$, and
the population target changes accordingly.

\subsection{Score, gap, and detectability}

Write $\gmax=\max_u\gamma(u)$ for the global maximum and
\begin{equation}\label{eq:score}
 \Psi_j=\frac{1}{1-2\varepsilon}\int_{\I}\xi_j(u)\,du,
 \qquad
 \Delta_j=\gmax-\Psi_j,
 \qquad
 \Delta_{\min}=\min_{j\in\A}\Delta_j.
\end{equation}
The screen ranks coordinates by increasing $\Psi_j$ and keeps the $d$
smallest. Small scores indicate activity: all inactive coordinates share the
flat baseline $\gmax$, while fixing a detectable active coordinate lowers
the attainable envelope. The algorithm uses only the ordering of the
$\Psi_j$ and never estimates $\gmax$.

Let $\M=\argmax_u\gamma(u)$, a nonempty compact set, and let $\pi_j(\M)$ be
its projection onto coordinate $j$.

\begin{proposition}[Detectability]\label{prop:geometry}
Assume (C1)--(C2) and (S).
\begin{enumerate}
\item If $j\notin\A$, then $\xi_j\equiv\gmax$ on $\I$ and $\Delta_j=0$.
\item For every $j$ and $u\in\I$: $\xi_j(u)=\gmax$ if and only if
$u\in\pi_j(\M)$.
\item $\Delta_j>0$ if and only if $\I\setminus\pi_j(\M)$ has positive
Lebesgue measure.
\end{enumerate}
If the restriction of $\gamma$ to the active coordinates has a finite set
of maximizers, every active coordinate is detectable.
\end{proposition}

Claim 3 separates two failure modes that finite-sample experiments tend to
blur. When $\pi_j(\M)$ covers $\I$, no sample size makes coordinate $j$
identifiable from this score; the diagonal example
$\gamma(u_1,u_2)=\gamma_0-(u_1-u_2)^2$, with $\M=\{(t,t):0\le t\le1\}$, is
of this kind, and both of its active coordinates are invisible. When
$\pi_j(\M)$ misses a set of positive measure, detection is a question of
estimation accuracy, quantified by
Section~\ref{sec:estimation}. Figure~\ref{fig:population-profiles}
contrasts the two situations.

\begin{figure}[t]
\centering
\includegraphics[width=.88\linewidth]{figures/population_profiles.pdf}
\caption{Projected envelopes. Left: a monotone active coordinate depresses
the envelope below the inactive baseline $\gmax$. Right: in the diagonal
example every fibre meets the argmax set, so the envelope of an active
coordinate stays flat.}
\label{fig:population-profiles}
\end{figure}

\begin{proposition}[Quantitative margin]\label{prop:margin}
Assume (C1)--(C2), (S), and that for some $u^\star_{\A}\in[0,1]^{\A}$ and $\nu>0$,
\[
 \gamma(u)\le\gmax-c_1\,\|u_{\A}-u^\star_{\A}\|^\nu
 \qquad\text{for all }u,
\]
for some constant $c_1>0$. Then, for $j\in\A$,
\[
 \gmax-\xi_j(u)\ge c_1\,|u-u_j^\star|^\nu,
 \qquad
 \Delta_j\ge c_1\,\E\bigl(|W-u_j^\star|^\nu\bigr),\quad
 W\sim\mathrm{Unif}(\I).
\]
For $\nu=2$ the right-hand side is at least
$c_1\{(1-2\varepsilon)^2/12+(1/2-u_j^\star)^2\}$.
\end{proposition}

Two consequences of \eqref{eq:fibre} help with examples. If
$\gamma(u)=G(\sum_{\ell\in\A}g_\ell(u_\ell))$ with $G$ strictly increasing,
then $\xi_j(u)=G(g_j(u)+\sum_{\ell\in\A,\,\ell\ne j}\max g_\ell)$, and
$\Delta_j>0$ whenever $g_j$ lies below its maximum on a subset of $\I$ of
positive measure. If $\partial_j\gamma\ge c_j>0$ throughout the cube,
arbitrary interactions are allowed and $\Delta_j\ge c_j/2$.

For the benchmark $\gamma(u)=\tfrac12\exp(-u_1-u_2-u_3-u_4)$ used in
Section~\ref{sec:simulation}, $\xi_j(u)=e^{-u}/2$ for $j\le4$ and
$\xi_j\equiv1/2$ otherwise, so $\Delta_j\approx0.186$ at
$\varepsilon=0.05$. The gap is comfortable, yet estimating $\xi_j(u)$ near
its maximum requires sample points whose remaining three active
coordinates are simultaneously small; the simulations keep population
detectability and finite-sample accessibility clearly separated.
\section{Estimation and screening theory}\label{sec:estimation}

\subsection{Empirical-rank local Hill scores}

The estimator follows the envelope construction coordinate by coordinate.
Let $R_{ij}$ be the rank of $X_{ij}$ among $X_{1j},\ldots,X_{nj}$ and
$\widehat U_{ij}=R_{ij}/(n+1)$. For $u\in\I$, form the rank window
\[
 \mathcal W_j(u)=\{i:|\widehat U_{ij}-u|\le h\},\qquad
 M_j(u)=|\mathcal W_j(u)|,\qquad
 k_j(u)=\lfloor\alpha M_j(u)\rfloor,
\]
and, writing $Y_{j,u,(1)}\le\cdots\le Y_{j,u,(M_j(u))}$ for the local order
statistics, the local Hill estimate
\begin{equation}\label{eq:hill}
 \widehat\xi_j(u)=\frac1{k_j(u)}\sum_{r=1}^{k_j(u)}
 \log\frac{Y_{j,u,(M_j(u)-r+1)}}{Y_{j,u,(M_j(u)-k_j(u))}},
\end{equation}
defined when $\alpha M_j(u)>1$. This is the local-fraction Hill statistic of
\citet[Definition~5]{gardesPodgorny2025} evaluated on one-coordinate rank
windows. The score averages the profile over the deterministic grid
$G=\{r/(n+1):r=1,\ldots,n\}\cap\I$,
\begin{equation}\label{eq:emp-score}
 \widehat\Psi_j=\frac1{|G_j^{\rm val}|}\sum_{u\in G_j^{\rm val}}\widehat\xi_j(u),
 \qquad
 G_j^{\rm val}=\{u\in G:\alpha M_j(u)>1\},
\end{equation}
Because the window counts $M_j(u)$ depend on the ranks only through
the grid position, $G_j^{\rm val}$ is the same set for every $j$, and
under the conditions of Theorem~\ref{thm:score} every grid point is
usable, so $G_j^{\rm val}=G$; the convention $\widehat\Psi_j=+\infty$
when no grid point is usable is an implementation safeguard that is
provably never triggered when $K_n\ge2$. The screened set
$\widehat\A_d$ collects the $d$ smallest scores, ties broken by a
seeded auxiliary randomization.

The two localization parameters play different roles: $h$ fixes how many
observations resemble the target rank position, roughly $2nh$ of them, and
$\alpha$ fixes how many of their largest responses enter the Hill statistic,
roughly $2n\alpha h$. Rank windows make these counts deterministic up to
boundary effects, which is what permits exact finite-sample union bounds
over a diverging number of coordinates.

\subsection{Assumptions}

Ranks displace the unobserved transforms $U_{ij}$, windows displace the
conditioning point, and the Hill statistic sees only the upper
$\alpha$-fraction of a projected distribution. The assumptions bound these
three effects uniformly over the coordinates. Throughout this section
every active coordinate is assumed detectable, $\Delta_{\min}>0$ for
every $n$, so that the conditions below are well posed. Let
$Q_j(\cdot,u)$ denote the generalized upper conditional quantile of $Y$
given $U_j=u$, $Q_j(s,u)=F_{Y|j}^{-1}(1-s^{-1}\mid u)$ with
$F^{-1}(t)=\inf\{y:F(y)\ge t\}$, and let $\tau_n=(pn)^{-2}$.

\begin{itemize}
\item[(E1)] \emph{Quantile representation and regularity of the projected
index.} For every $j\le p$ and every $u\in[0,1]$ the conditional
distribution function $F_{Y|j}(\cdot\mid u)$ is continuous, so that
$V_{ij}=1-F_{Y|j}(Y_i\mid U_{ij})$ is conditionally standard uniform and
$Y_i=Q_j(1/V_{ij},U_{ij})$ almost surely. Moreover the projected indices
are equicontinuous on $\I$, uniformly over coordinates and over $n$,
\[
 \omega_\xi(r)=\max_{j\le p}\;
 \sup_{u,v\in\I,\,|u-v|\le r}|\xi_j(u)-\xi_j(v)|
 \longrightarrow0\qquad(r\to0),
\]
and satisfy the grid-scale regularity condition
$\omega_\xi(1/n)=o(\Delta_{\min})$.
\item[(E2)] \emph{Uniform tail regularity.} Writing
$Q_j(s,u)=s^{\xi_j(u)}\ell_j(s,u)$ with $\ell_j(\cdot,u)$ slowly varying,
\[
 \max_{j\le p}\,\sup_{u\in\I}\,
 \sup_{s\ge(3\alpha)^{-1},\,t>1}
 \frac{|\log\{\ell_j(ts,u)/\ell_j(s,u)\}|}{\log t}
 \;=\;o(\Delta_{\min}).
\]
\item[(E3)] \emph{Uniform local quantile regularity.} There exists a
constant $C_0>1$ such that
\[
 \max_{j\le p}\,
 \sup_{u\in\I}\;
 \sup_{v\in[0,1],\,|u-v|\le C_0h}\;
 \sup_{w\in[\tau_n,1-\tau_n]}
 \left|\log\frac{Q_j(w^{-1},v)}{Q_j(w^{-1},u)}\right|
 \;=\;o(\Delta_{\min}).
\]
The theorems use only this one $C_0$; under \eqref{eq:rate-cond} any
$C_0>1$ suffices.
\end{itemize}

Scale nuisances enter through (E3) and not through (E2), and the
condition accommodates them. Writing
$\log Q_j(s,u)=\xi_j(u)\log s+\log\ell_j(s,u)$, the displacement over
$|v-u|\le C_0h$ splits into an index part, bounded by
$2L_\xi C_0h\log(pn)$ when $\xi_j$ is Lipschitz with constant $L_\xi$,
and a scale part which is $O(L_cC_0h)$ whenever the conditional
log-scale profile is Lipschitz in $u$ with constant $L_c$ uniform in
$n$. A covariate that shifts the conditional scale of $Y$ without
touching the index is therefore harmless in the limit, at any bounded
amplitude: both parts are $o(\Delta_{\min})$ as soon as
$h\log(pn)=o(\Delta_{\min})$. Model B1 of
Section~\ref{sec:simulation} is built to make the finite-sample
version of this statement measurable.

(E1) asks for continuity only; strict monotonicity of the conditional
distribution function is not needed, because flat pieces are visited with
probability zero. (E2) is the analogue of condition (H.4) in
\citet{gardesPodgorny2025}, normalized by $\log t$ and required only beyond
the working tail level $(3\alpha)^{-1}$; it accommodates slowly varying
factors such as $\ell(s)\asymp(\log s)^{-2}$, for which the left side is of
order $1/\log(1/\alpha)\asymp1/\log n$ when $\alpha=n^{-a}$, so gaps
$\Delta_{\min}\gg1/\log n$ are covered and smaller gaps are not. (E3)
adapts their condition (H.5): the truncation $\tau_n=(pn)^{-2}$ makes the
probability that any of the $pn$ conditional PIT variables falls below
$\tau_n$ vanish, which is the high-dimensional replacement for their
fixed-dimension truncation. Their density conditions (H.1)--(H.2)
have no counterpart here: empirical-rank windows have deterministic counts,
so no lower bound on a covariate density is required.

These conditions are stated on the projected quantiles rather than on
the primitive model, and deliberately so: projection generates slow
variation of its own, which a condition imposed upstream would have to
anticipate. If $W\sim U(0,1)$ and $\Pp(Y>y\mid W)=y^{-1/(1+W)}$
exactly, then $\Pp(Y>y)\sim4y^{-1/2}/\log y$, so even exact Pareto
conditionals produce a logarithmic factor after projection, and the
left side of (E2) is then of exact order $1/\log(1/\alpha)$. Stating
(E2) downstream measures that factor where it acts, and the resulting
requirement is mild: every gap with $\Delta_{\min}\log(1/\alpha)
\to\infty$ is covered, which for $\alpha=n^{-a}$ admits all gaps
$\Delta_{\min}\gg1/\log n$, and in particular every fixed gap.

Appendix~\ref{app:primitive} derives these conditions from primitive
assumptions. Two steps make this possible. First, the range
$w\in[\tau_n,1-\tau_n]$ in (E3) is wider than the proofs need: only the
$k+1$ largest responses of a window enter a Hill statistic, and
Lemma~\ref{lem:app-tailfreeze} shows that on the working event these are
generated by $w\le3\alpha$. Replacing (E3) by the tail-only condition
\[
 \omega_n^{\uparrow}(C_0h,\tau_n,\alpha)
 =\max_{j\le p}\sup_{u\in\I}
 \sup_{|v-u|\le C_0h}\sup_{w\in[\tau_n,3\alpha]}
 \left|\log\frac{Q_j(w^{-1},v)}{Q_j(w^{-1},u)}\right|
 =o(\Delta_{\min})
 \tag{E3$^{\uparrow}$}
\]
leaves Theorems~\ref{thm:score} and~\ref{thm:sis} unchanged, and removes
any requirement on the body of the conditional distribution. Second,
Theorem~\ref{thm:app-transfer} derives (E1), (E2) and (E3$^{\uparrow}$)
from four primitive conditions --- Lipschitz continuity and
atomlessness, a differentiable uniform tail approximation, uniform
doubling of the fibre-deficit mass for
$D_{j,u}=1/\gamma-1/\xi_j(u)$, and a horizontal coupling of the
primitive tails only --- with
\[
 \omega_\xi(1/n)=O(n^{-1}),
 \qquad
 b_n(3\alpha)=O\{1/\log(1/\alpha)\},
 \qquad
 \omega_n^{\uparrow}=O\{h\log(pn)\}.
\]
Corollary~\ref{cor:app-models} verifies the primitive conditions for the
Gaussian-copula models of Section~\ref{sec:simulation}, whose corner
margins are regularly varying with a slowly varying Gaussian factor, so
that a pure two-sided polynomial margin is neither assumed nor needed.

\subsection{Main results}

\begin{theorem}[Score concentration]\label{thm:score}
Assume (C1)--(C2), (S), (E1)--(E2), (E3$^{\uparrow}$), $\alpha\le1/3$,
$h\le1-\varepsilon$, and
\begin{equation}\label{eq:rate-cond}
 \log(pn)=o(nh^2),
 \qquad
 \log(pn)=o(n\alpha h).
\end{equation}
Then there are deterministic sequences
$a_n=O(\sqrt{\log(pn)/(n\alpha h)})$ and $\beta_n=o(\Delta_{\min})$ with
\[
 \Pp\Bigl(\max_{j\le p}\,|\widehat\Psi_j-\Psi_j|\le a_n+\beta_n\Bigr)
 \longrightarrow1;
\]
in shorthand,
$\max_{j\le p}|\widehat\Psi_j-\Psi_j|
=O_{\Pp}(\sqrt{\log(pn)/(n\alpha h)})+o(\Delta_{\min})$.
\end{theorem}

The stochastic term is the exponential-spacings fluctuation of
$O(n)$ local Hill statistics per coordinate, each built from at least
$c\,n\alpha h$ spacings, controlled simultaneously over $p$ coordinates;
$\sqrt{\log(pn)/(n\alpha h)}$ is the resulting effective-tail-sample-size
rate. The $o(\Delta_{\min})$ term collects the tail bias (E2), the spatial
bias (E3$^{\uparrow}$) evaluated at $C_0h$ --- the first condition in
\eqref{eq:rate-cond} makes the rank displacement smaller than $(C_0-1)h$,
so windows never leave the range covered by (E3$^{\uparrow}$) --- and the
quadrature
error $\omega_\xi(1/n)+O(1/n)$ from (E1).

\begin{theorem}[Sure screening]\label{thm:sis}
Assume the conditions of Theorem~\ref{thm:score} (including
$\Delta_{\min}>0$, i.e.\ every active coordinate detectable in the
sense of Proposition~\ref{prop:geometry}) and
\begin{equation}\label{eq:sis-cond}
 \log(pn)=o\bigl(n\alpha h\,\Delta_{\min}^2\bigr).
\end{equation}
Then, for every $d\ge s$,
\[
 \Pp\bigl(\A\subseteq\widehat\A_d\bigr)\longrightarrow1,
\]
and with probability tending to one every active score is smaller than
every inactive score.
\end{theorem}

\begin{corollary}[Exact recovery]\label{cor:exact}
Under the conditions of Theorem~\ref{thm:sis}, if $d=s$ then
$\Pp(\widehat\A_s=\A)\to1$.
\end{corollary}

The proofs are in Appendix~\ref{app:proofs}, supported by the auxiliary
results of Appendix~\ref{app:aux}, whose centrepiece is a finite-sample
uniform profile bound (Proposition~\ref{thm:app-rank-profile}). The
roadmap: (i) a Dvoretzky--Kiefer--Wolfowitz event per coordinate controls
all empirical-rank displacements simultaneously, contributing the
$\log(pn)=o(nh^2)$ requirement; (ii) rank windows have deterministic
local counts and pass through at most $4n+1$ membership states per
coordinate; (iii) conditionally on a covariate column, the R\'enyi
representation turns each local Hill statistic into a mean of
i.i.d.\ exponentials built from about $n\alpha h$ extremes, plus a
slowly varying remainder bounded through (E2); (iv) a union bound over
the $p(4n+1)$ states yields the $\log(pn)$ factor; (v) uniform score
accuracy below $\Delta_{\min}/2$ gives deterministic separation of
active from inactive scores. \citet[Theorem~3]{gardesPodgorny2025}
establishes the corresponding uniform consistency in fixed dimension
over a compact class of projection matrices; their covering-number
argument and density assumptions are replaced here by the rank
construction, and their result does not imply the growing-$p$ bound.

Condition \eqref{eq:sis-cond} quantifies the cost of dimension:
screening $p$ coordinates costs $\log p$ against the effective local
extreme count $n\alpha h$, deflated by the squared gap. With
$\alpha=n^{-a}$ and $h=n^{-b}/2$ and $b\le a$, dimensions up to
$\log p=o(n^{1-a-b}\Delta_{\min}^2)$ are admissible, a polynomial range
of $p$ for fixed gaps; for $b>a$ the rank-displacement condition
$\log(pn)=o(nh^2)=o(n^{1-2b})$ binds first. These ranges treat (E3$^{\uparrow}$) as given; a spatial modulus of order
$h\log(1/\tau_n)$, the natural benchmark, adds the constraint
$\log(pn)=o(n^{b}\Delta_{\min})$.
\section{Simulation study}\label{sec:simulation}

Each coordinate is scored in $O(n\log n)$ time: one sort by empirical
rank, then a sliding window over the rank grid whose local Hill statistics
are maintained by a Fenwick tree indexed by global response ranks.
Covariates are generated and processed in blocks, so the total cost is
$O(pn\log n)$ time and $O(n)$ working memory beyond the scores. At
$(n,p)=(2000,1000)$ one complete screen takes under a second on a laptop
core. Scores are sorted increasingly.

\subsection{Simulation models}\label{sec:models}

\paragraph{Family A: tail-index structure.}
Models A1, A2 and A3 share the covariate design, the slowly varying
factor, and the tail-index active set $\A_\gamma=\{1,2,3,4\}$; they
differ only in the tail-index function. Covariates are latent Gaussian
AR(1),
\[
 Z_{i1}\sim N(0,1),\qquad
 Z_{ij}=\rho Z_{i,j-1}+\sqrt{1-\rho^2}\,\epsilon_{ij},
 \qquad j=2,\ldots,p,
\]
transformed to $U_{ij}=\Phi(Z_{ij})$; the tail-index functions are
evaluated at the population transforms $U$, while the estimator sees
only $Z$ and computes its own empirical ranks. Responses are generated
by quantile inversion, $Y=V^{-\gamma(U)}\ell(1/V,\cdot)$ with
$V\sim\mathrm{Unif}(0,1)$, the mechanism of
\citet{gardesPodgorny2025}, and all three models carry their slowly
varying factor
\[
 \ell_1(s,u)=\bigl[1+\exp\{v(u)-s\}\bigr]^{-1},
 \qquad
 v(u)=2.5\,(u_5+u_6+u_7+u_8),
\]
which satisfies $\ell_1(s,u)\to1$ as $s\to\infty$: its effect vanishes
in the limit while perturbing estimation at finite tail levels. The
tail-index functions are
\begin{align*}
 \gamma_{A1}(u)&=0.5\exp\bigl(-\textstyle\sum_{j=1}^4u_j\bigr),\\
 \gamma_{A2}(u)&=0.5\exp\bigl(-0.5\textstyle\sum_{j=1}^4u_j\bigr),\\
 \gamma_{A3}(u)&=0.5\exp\bigl(-0.80\textstyle\sum_{j=1}^4u_j
 -0.20\,u_1u_2-0.15\,u_3u_4\bigr).
\end{align*}
A1 is the monotone reference, with trimmed population gap $0.186$
(Section~\ref{sec:population}); its tail-index function and its Gaussian
AR(1) covariate design are those of \citet{yoshidaUmezu2026}, adopted
deliberately so that the comparison of Section~\ref{sec:comparison} runs
on a benchmark chosen by a competing method rather than on one of our
own. The response mechanism differs: they invert a Burr-type conditional
survival function, whereas A1--A3 carry the slowly varying factor
$\ell_1$ above. A2 weakens the signal, with gap
$0.107$. A3 adds interactions among the active coordinates while
keeping a reasonably strong signal, with gap $0.158$. All three share
the maximum $\gmax=0.5$.

\paragraph{Family B: finite-quantile scale contamination.}
Model B1 has the same tail-index function as A1,
$\gamma_{B1}=\gamma_{A1}$, hence the same active set
$\A_\gamma=\{1,2,3,4\}$ and the same gap $0.186$, and \emph{replaces}
the slowly varying factor: responses are generated with $\ell_2$ below
in place of $\ell_1$, which does not appear in B1. Let
$F_i\sim N(0,1)$ be independent of the active coordinates and of the
base AR(1) innovations, and set
\[
 \log\ell_2(s,F)=-\frac1{2s}+0.5\,\bigl\{\Phi(F)-0.5\bigr\},
\]
so the persistent scale factor lies in $[e^{-0.25},e^{0.25}]$, strictly
positive and uniformly bounded. The factor $F$ is not observed. Instead,
a block of $q_{\rm scale}=20$ proxies occupies the coordinates
$\A_{\rm scale}=\{5,\ldots,24\}$:
\[
 Z_{ij}=0.7\,F_i+\sqrt{1-0.7^2}\,\epsilon_{ij},
 \qquad j\in\A_{\rm scale},
\]
with independent standard normal proxy noises, so each proxy is
marginally standard normal, loads on $F$ at $0.7$, and the pairwise
correlation within the block is $0.49$. The B1 covariate vector is
built from a base AR(1) vector of length $p$: coordinates $1$--$4$ and
$25$--$p$ are kept from the base vector, coordinates $5$--$24$ are
replaced by the proxies, and the AR(1) chain beyond coordinate $24$
continues from the base coordinate $24$, not from the overwritten
proxy. The observed B1 covariance is therefore deliberately not AR(1)
over all coordinates; the AR(1) background is intact outside the proxy
block. Because B1 changes both the slowly varying factor and the
covariate design relative to A1, A1--B1 differences measure the joint
effect of the scale cluster and of the exchange of $\ell_1$ for
$\ell_2$, not the cluster alone.

In B1, variables $1$--$4$ determine the limiting tail index; variables
$5$--$24$ are informative about the finite conditional scale through
their common information about $F$, and do not change the limiting
exponent. A screen built on a fixed high quantile can therefore
legitimately rank the proxy cluster highly --- those variables move its
target --- even though the cluster is inactive for tail-index
screening. B1 measures the practical consequence of this target
difference. Section~\ref{sec:estimation} places the mechanism in the
theory: a bounded scale nuisance enters (E3) through a term of order
$L_ch$, which vanishes with the bandwidth, so the screen is
asymptotically immune to it; B1 measures how much sample that immunity
costs at a finite bandwidth.

\paragraph{Monte Carlo design.}
The common design is
\[
 n=2000,\qquad p\in\{500,\,1000,\,2000\},\qquad \rho=0.25,
 \qquad d=20,\qquad \varepsilon=0.05,
\]
with $1{,}000$ replications per cell and per model in the method
comparison (binomial standard errors at most $0.016$) and $200$
common-random-number replications per model in the tuning study
(standard errors at most $0.036$); Sure-$d$ denotes
$\Pp(R_{\max}\le d)$ with $R_{\max}=\max_{j\in\A_\gamma}\mathrm{rank}(j)$,
and all methods within a replication are evaluated on the same dataset.

\paragraph{Relation to the asymptotic theory.}
The simulation models instantiate the theory rather than sit beside it.
Corollary~\ref{cor:app-models} verifies the primitive conditions of
Appendix~\ref{app:primitive} for A1--A3 and B1, and shows that (E1),
(E2) and (E3$^\uparrow$) hold along the tuning sequences
$\alpha=n^{-a}$, $h=n^{-b}/2$ used throughout, with bias orders
$O(n^{-1})$, $O(1/\log n)$ and $O(n^{-b}\log n)$. The effective local
extreme count is $n\alpha h=n^{1-a-b}/2$ while the population gaps are
constants, so conditions \eqref{eq:rate-cond} and \eqref{eq:sis-cond}
reduce to $b<1/2$ and $a+b<1$: the selected $(a,b)=(0.30,0.15)$ and all
nine settings of the aggregation block satisfy them, and
Theorems~\ref{thm:score} and~\ref{thm:sis} apply. The same corollary
covers the aggregated screen, by intersecting the nine separation
events.

What the finite design does not deliver is a numerically tight bound:
at $n=2000$ and $p=1000$, $\sqrt{\log(pn)/(n\alpha h)}\approx0.67$
against $\Delta_{\min}=0.186$, and the lower endpoint entering (E2) is
only $\log\{(3\alpha)^{-1}\}\approx1.18$, so the asymptotic bounds are
not yet informative at this sample size, as is usual for rate results at
moderate $n$ (Remark~\ref{rem:app-finite}). This is why the recovery rates below are reported as
direct finite-sample evidence, and it is the practically relevant
question: the theory says the screen eventually separates active from
inactive coordinates, and the simulations measure how much sample the
separation actually needs.

\subsection{Tuning and the aggregated screen}\label{sec:tuning}

The two smoothing parameters are parameterized as
\[
 \alpha=n^{-a},\qquad h=n^{-b}/2,
\]
following the convention of \citet{gardesPodgorny2025}. Sure-20 is
evaluated for each model over the crossed grid
$a\in\{.25,.30,.35,.40,.45,.50\}$,
$b\in\{0,.05,.10,.15,.20,.30,.40\}$ at $n=2000$, $p=1000$,
$\rho=0.25$, on $200$ common-random-number replications per model: the
$200$ datasets are generated once per model and every cell of the grid
is evaluated on the same datasets, so that differences between cells
are not inflated by simulation noise. Figure~\ref{fig:tuning} shows the
four heatmaps.

The column $a=0.25$ collapses for the three A-family models (Sure-20 at
most $.125$, $.080$ and $.225$): at $\alpha=n^{-0.25}\approx0.15$ the
local Hill statistics reach into quantile regions where $\ell_1$, whose
argument $v(u)-s$ is dominated by the $U_5$--$U_8$ combination at
moderate $s$, deforms the projected distributions. B1, whose bounded
scale factor produces no such deformation, degrades only mildly there
(up to $.580$). Performance also falls toward large $a$ and large $b$
together, where the local extreme count collapses, and the degenerate
bandwidth $b=0$ is never competitive. In between, the block
$a\in\{0.30,0.35\}$, $b\in\{0.10,0.15,0.20\}$ is the best-behaved
region of the grid, but it is not a flat plateau. The four-model mean
of Sure-20 is maximized at
\[
 (a^\star,b^\star)=(0.30,\ 0.15),
\]
at $0.780$; the five other cells of the block lie between $0.718$ and
$0.773$, and the per-model values move by more than their paired
standard errors across the block --- A2, the weakest signal, falls
from $.555$ at the argmax to $.360$ at $(0.35,0.20)$. The residual
sensitivity within an otherwise good region is precisely what the
aggregated screen below is designed to absorb. We keep the argmax as
the single selected setting, with per-model Sure-20 $.935$, $.555$,
$.955$ and $.675$; because it is selected on the same $200$
replications shown in Figure~\ref{fig:tuning}, the comparison of
Section~\ref{sec:comparison} reruns everything on $1{,}000$
independent replications. At $n=2000$ the selected pair gives
$\alpha\approx0.10$, $h\approx0.16$, and about
$2n\alpha h=n^{1-a-b}\approx65$ effective local extremes per interior
window. No model is given its own tuning.

\begin{figure}[t]
\centering
\includegraphics[width=.95\linewidth]{figures/tuning_heatmaps.pdf}
\caption{Sure-20 over the $(a,b)$ tuning grid for the four models at
$n=2000$, $p=1000$, $\rho=0.25$ ($200$ common-random-number
replications per model; $\alpha=n^{-a}$, $h=n^{-b}/2$; standard errors
at most $0.036$). The dashed rectangle marks the aggregation block
$\mathcal N_9$.}
\label{fig:tuning}
\end{figure}

\paragraph{The aggregated screen.}
A single pair $(a,b)$ makes the final ranking depend entirely on one
tuning. The aggregated version of the screen uses the $3\times3$ block
of settings
\[
 \mathcal N_9=\{0.30,0.35,0.40\}\times\{0.10,0.15,0.20\},
\]
the neighbourhood of $(a^\star,b^\star)$ inside the well-behaved region
of the grid. For each setting $(a,b)\in\mathcal N_9$ and each
coordinate $j$, compute the score $\widehat\Psi_j(a,b)$ and its rank
$r_j(a,b)$ among the $p$ coordinates, ranks increasing with the score;
aggregate by
\[
 r_j^{\rm agg}=\min_{(a,b)\in\mathcal N_9}r_j(a,b),
\]
and rank the coordinates by increasing $r_j^{\rm agg}$; coordinates
sharing the same minimum are ordered by the median, then the mean, then
the maximum of their nine ranks, and any remaining ties by a seeded
auxiliary randomization. A coordinate is thus retained early whenever at least
one reasonable tuning finds evidence for it. This direction of
conservatism matches the purpose of a sure-screening first stage: false
exclusion is final, while false retention is corrected by the
second-stage multivariate analysis that the screen feeds. We make no
optimality claim for the construction; it is the procedure used in the
experiments and the application, and its effect is measured against the
single-setting screen in Section~\ref{sec:comparison}. Its cost is nine
score passes instead of one, leaving the $O(pn\log n)$ complexity
unchanged up to the constant.

The two versions are reported side by side throughout:
\emph{Tail-index SIS, selected $(a^\star,b^\star)$} uses the central
pair alone; \emph{Tail-index SIS, aggregated} uses $r_j^{\rm agg}$.

\subsection{Comparison with competing screens}\label{sec:comparison}

Six rules are compared on the same $1{,}000$ simulated datasets per
model and per dimension $p\in\{500,1000,2000\}$, at $n=2000$,
$\rho=0.25$: the two versions of the tail-index screen; the
conditional-Pickands screen of \citet{yoshidaUmezu2026}, implemented
from their equations (2)--(4) with the Epanechnikov kernel and the
tuning of their numerical work transplanted proportionally --- their
bandwidth $h=1$ on standardized covariates is kept, and their
extreme-value count is carried over as the same fraction of the sample
size, $k=\lfloor0.072\,n\rfloor$ --- on 25 equally spaced evaluation
points in $[.02,.98]$; and the quantile-adaptive SIS of \citet{heWangHong2013} at
the levels $\tau\in\{.90,.95,.99\}$, using cubic B-spline marginal
quantile regression with three degrees of freedom, identical at every
level. The competitors rank large utilities first; ours ranks small
scores first. Tables~\ref{tab:cmp500}--\ref{tab:cmp2000} report, for
each model and method, Sure-4 (exact recovery), Sure-20, and the mean
and median of the worst active rank. Both competitors run at their
authors' settings, with no counterpart of the grid search of
Section~\ref{sec:tuning}; in this respect the comparison favours the
tail screen, and Section~\ref{sec:tuning} shows how much its own
performance moves across reasonable tunings.

\input{tables/cmp500.tex}
\input{tables/cmp1000.tex}
\input{tables/cmp2000.tex}
\input{tables/b1comp.tex}

\paragraph{Sensitivity to the dimension.}
Every method degrades as $p$ grows, at very different rates, and the
mean--median contrast locates the degradation. For the selected-pair
tail screen on A1, $\E(R_{\max})$ grows from $5.7$ to $12.9$ between
$p=500$ and $p=2000$ while $\mathrm{Med}(R_{\max})$ stays at $4$: the
typical replication is unaffected, and the loss is carried by a
minority of replications in which an active coordinate falls into the
null range, where its rank scales with $p$. Aggregation compresses that
tail on every model ($\E(R_{\max})$ $4.6\to6.2$ on A1, $9.8\to22.6$
against $29.8\to108.7$ on B1, $5.1\to8.1$ against $7.0\to22.0$ on A3).
The quantile screen at $\tau=.95$ is the least dimension-sensitive rule
on the A family ($\E(R_{\max})$ $4.2\to4.8$ on A1). A2 is the most
dimension-sensitive model for every method (aggregated Sure-20
$.797\to.561$, quantile $.95$ $.930\to.782$): a weak signal loses the
most from added competition. Yoshida--Umezu degrades from a low base
throughout ($\mathrm{Med}(R_{\max})$ $29\to168$ on A1).

\paragraph{Sensitivity to the quantile level.}
On the A family the quantile screen is non-monotone in $\tau$ in every
table, in opposite directions on either side of $\tau=.95$; on B1 its
behaviour tracks its target, and Table~\ref{tab:b1comp} reports the
resulting composition of its leading positions. At $\tau=.90$ its top
four contain on average $3.8$ of the twenty scale proxies and all
twenty fill its top 24; exact recovery is zero to three decimals at
every $p$. This is coherent behaviour for a $.90$-quantile target,
which the proxies genuinely move, and it is fatal for tail-index
recovery. At $\tau=.95$ the contamination remains substantial
($2.4$--$2.5$ proxies in the top four, $19.1$--$19.5$ in the top 24,
median worst active rank $20$--$21$). At $\tau=.99$ the contamination
largely clears ($0.2$--$0.4$ proxies in the top four) because the
quantile target approaches the tail index, but the variance cost is
severe everywhere: on the A family, Sure-4 falls from $.234$--$.891$
at $\tau=.95$ to at most $.125$ at $\tau=.99$, and on B1 the level
$.99$ never reaches Sure-4 $.04$. Yoshida--Umezu, built on
Pickands-type quantile differences at its own fixed levels, is
dominated on all models and dimensions here, with $6.8$--$9.8$ of the
twenty B1 proxies in its top 24 (Table~\ref{tab:b1comp}).

\paragraph{Effect of the aggregation.}
Aggregation improves the tail screen most where failures are worst. On
B1 it lifts exact recovery from $.308/.249/.234$ to $.448/.368/.339$
across the three dimensions, divides $\E(R_{\max})$ by three to five,
and brings $\mathrm{Med}(R_{\max})$ to $5$--$6$; on A1 and A3 it moves
Sure-4 by a few points in either direction ($+3.5$, $+1.4$, $+5.1$ on
A1; $+1.5$, $-4.1$, $+0.7$ on A3) while compressing the
rank tail; Sure-20 improves in all twelve model--dimension cells. On the weak-signal A2 it
trades exact recovery ($.256\to.196$ at $p=500$,
$.149\to.088$ at $p=2000$) for substantially better inclusion
(Sure-20 $.463\to.561$ at $p=2000$; $\E(R_{\max})$ $144.8\to76.9$).
The trade is the intended one: promotion by a single reasonable tuning
occasionally admits a lucky null into the top four when the signal is
weak, and reliably rescues active coordinates from single-setting
failures further down the list. Aggregation stabilizes inclusion ---
the sure-screening objective --- rather than exact selection.

The quantile screen at a moderate level remains the strongest detector
when high quantiles and the limiting exponent share their covariates
and the signal is weak (A2); the aggregated tail screen is competitive
or better on A1 and A3, and on B1 it is the only rule whose leading
positions stay on the active set --- Sure-4 $.448/.368/.339$ and
Sure-20 $.907/.849/.805$ across the three dimensions, against at most
$.038$ and $.665$ for any fixed-quantile level, with $0.12$--$0.19$
proxies in its top four (Table~\ref{tab:b1comp}). It is degraded, not
immune: the same rule reaches Sure-4 $.884/.823/.802$ on A1, so the
scale cluster, jointly with the exchange of slowly varying factors,
costs it about forty-five points of exact recovery. For a question
about the limiting exponent itself, the two rankings are best read
side by side, as in Section~\ref{sec:realdata}: variables high on one
list and low on the other are where the two targets part ways.
\section{Application: violent crime across U.S. communities}
\label{sec:realdata}

The analysis contains $n=1{,}993$ communities and $p=99$ candidate
predictors. Under the tuning selected in Section~\ref{sec:tuning},
$\alpha=n^{-0.3}\approx0.102$ and $h=n^{-0.15}/2\approx0.160$, so the
sample contains $n\alpha\approx204$ upper-tail observations in total and
about $2n\alpha h\approx65$ within a typical interior rank window. The
dimension is therefore comparable to the number of observations
available for local tail estimation. Unrestricted or nonparametric joint
tail modelling with $99$ covariates is poorly supported by roughly $200$
upper-tail observations, and by fewer still within any local covariate
neighbourhood; a marginal screen can rank the $99$ candidates from
exactly this information.

\subsection{Data and protocol}

The Communities and Crime data set
\citep{redmond2002}, distributed by the UCI Machine Learning Repository
in its unnormalized version, links socio-economic census attributes,
law-enforcement survey attributes and FBI crime counts for $2{,}215$ U.S.
communities. The response is the violent-crime rate per $100{,}000$
inhabitants. We keep the communities with an observed positive response
($n=1{,}993$). \citet{yoshidaUmezu2026} analyse the same data set with
their own screen, taking the robbery rate as response; we deliberately
stay on the same data so that the two screens can be read side by side,
and take the violent-crime rate, whose larger positive count leaves a
usable upper tail after the response-blind filtering described below. Preprocessing of the $124$ documented predictive
attributes is response-blind: we retain numeric predictors with fewer
than $5\%$ missing values, which removes the $22$ law-enforcement
variables reported by only $343$ communities, and with at least $100$
distinct observed values, which removes $3$ low-cardinality variables,
keeping the application close to the continuous-marginal framework of
the theory. This leaves $p=99$ continuous predictors --- demographic
composition, income and poverty, education, employment, family
structure, immigration, language, housing, and population density. The
at most one missing value per retained predictor is median-imputed, and
every method below works from the same imputed matrix.

The response is strongly right-skewed: median $374$, $99$th percentile
$2{,}943$, maximum $4{,}877$. The left panel of Figure~\ref{fig:crime}
shows the Hill curve of the response over $k\in[20,600]$ upper order
statistics; it is positive and moves slowly, between $0.283$ and
$0.374$ for $k\in[100,300]$, with value $0.328$ at the working level
$k=n\alpha=204$. The middle panel shows the Pareto QQ plot of the $600$
largest responses, approximately linear over the working range with
least-squares slope $0.27$. The gap between the QQ slope and the Hill
value at the same level is mild evidence that regular variation is only
approximate over the working range. The response exhibits a
sufficiently heavy and stable upper tail to make a tail-index analysis
plausible; neither diagnostic establishes regular variation beyond
doubt, and no single threshold carries a structural interpretation.

The screening protocol is the one fixed in Section~\ref{sec:tuning},
applied unchanged: the aggregated tail-index screen computes the score
at each of the nine settings
$\mathcal N_9=\{0.30,0.35,0.40\}\times\{0.10,0.15,0.20\}$, converts
each pass into ranks, and orders the predictors by the minimum rank
$r_j^{\rm agg}$ across the nine settings; the selected-pair screen at
$(a^\star,b^\star)=(0.30,0.15)$ is reported alongside. The two
competing screens of Section~\ref{sec:comparison} are run on the same
imputed matrix: Yoshida--Umezu at its transplanted tuning and the
quantile screen at $\tau\in\{.90,.95,.99\}$.

\subsection{Results}

Poverty and family structure dominate the aggregated ranking, and they
dominate every competing ranking as well. The two-parent household
variable \texttt{pctKids2Par} is first (per-setting ranks $1$--$4$
across the nine settings), the poverty rate \texttt{pctPoverty} second
(ranks $1$--$6$) and \texttt{pct2Par} third (ranks $1$--$5$); the two
family-structure variables are also leading variables for the quantile
screen at every level (\texttt{pctKids2Par} is first at all three
$\tau$) and for Yoshida--Umezu (ranks $3$--$7$). The related
\texttt{pctKids4w2Par} and \texttt{pct12to17w2Par} complete the
consensus block at aggregated ranks $7$--$8$. The high-quantile and
tail-index rankings share these leading predictors.

The informative disagreements are the tail-specific candidates. The
share of workers using public transportation
(\texttt{pctUsePubTrans}, aggregated rank $5$) and the farm-employment
rate (\texttt{pctWfarm}, rank $6$) rank far lower under every
fixed-quantile method: ranks $49$--$73$ for the quantile screen across
its three levels and $57$--$60$ for Yoshida--Umezu. Their low
tail-index scores are consistent with marginal variation in upper-tail
heaviness across their values, of a kind that fixed-quantile screens are
not designed to detect; Section~\ref{sec:comparison} produces the same
dissociation under controlled conditions, where the target is known.
The two variables differ in stability across the aggregation block:
\texttt{pctUsePubTrans} ranks within the top $15$ at every one of the
nine settings, while \texttt{pctWfarm} moves between $2$ and $47$
(right panel of Figure~\ref{fig:crime}). The never-married male rate
(\texttt{pctMaleNevMar}, aggregated rank $4$, per-setting ranks
$1$--$9$, quantile ranks $28$--$35$) sits between the consensus block
and the tail-specific pair. We report all three as candidate
tail-specific predictors for a second-stage analysis.

Minimum-rank aggregation is deliberately conservative: a predictor is
promoted whenever at least one tuning in the prespecified neighbourhood
ranks it highly, which suits a first-stage sure screen, where false
exclusion is more costly than retaining additional candidates.
\texttt{pctHousWOplumb} (aggregated rank $9$) and \texttt{pctLargHous}
(rank $10$) illustrate the resulting behaviour: their positions rest on
one favourable setting each (per-setting ranks $6$--$55$ and $6$--$74$;
selected-pair ranks $36$ and $46$), so their rankings are markedly more
tuning-sensitive than those of the consensus block or of the
tail-specific candidates, which hold their positions across the whole
block. The conservative screen retains them intentionally; the
per-setting ranks describe this sensitivity, and discrimination among
the retained predictors is left to the subsequent multivariate
analysis.

We do not interpret any ranking causally: the screen identifies
covariates along which the conditional tail index varies marginally,
and confounding among the $99$ predictors is certain.

\begin{figure}[t]
\centering
\includegraphics[width=.98\linewidth]{figures/crime_application.pdf}
\caption{Communities and Crime application ($n=1{,}993$, $p=99$,
$\alpha=n^{-0.3}$, $h=n^{-0.15}/2$). Left: Hill estimates for the
violent-crime rate; the dashed line marks the working level
$k=n\alpha=204$. Middle: Pareto QQ plot of the $600$ largest responses;
the dashed line is fitted over the working range $i\le204$. Right: for
the ten leading predictors of the aggregated screen, the nine
per-setting ranks (open circles) and the aggregated rank (filled
circles); the dotted line is the default retained size $d=20$.}
\label{fig:crime}
\end{figure}
\section{Auxiliary results}
\label{app:aux}
\renewcommand{\theequation}{A.\arabic{equation}}
\setcounter{equation}{0}

This appendix collects the technical tools used by the proofs of
Appendix~\ref{app:proofs}: elementary probabilistic lemmas
(\ref{app:elementary}), the projection lemmas behind
Proposition~\ref{prop:projection} (\ref{app:projlem}), the deterministic
geometry of empirical-rank windows (\ref{app:rankgeom}), the local Hill
representation and the finite-sample uniform profile bound
(\ref{app:hill}), and the passage from profiles to scores, ending with a
deterministic separation lemma (\ref{app:scores}).

Throughout, the dimension is written $p=p_n$, and two auxiliary moduli
quantify the content of (E2)--(E3):
\[
 b_n(a)=\max_{j\le p}\sup_{u\in\I}\sup_{s\ge a^{-1},\,t>1}
 \frac{|\log\{\ell_j(ts,u)/\ell_j(s,u)\}|}{\log t},
 \qquad
 \omega_n(r,\tau)=\max_{j\le p}
 \sup_{u\in\I}\;
 \sup_{\substack{v\in[0,1]\\|v-u|\le r}}\;
 \sup_{w\in[\tau,1-\tau]}
 \left|\log\frac{Q_j(w^{-1},v)}{Q_j(w^{-1},u)}\right|,
\]
so that (E2) reads $b_n(3\alpha)=o(\Delta_{\min})$ and (E3) reads
$\omega_n(C_0h,\tau_n)=o(\Delta_{\min})$ with $\tau_n=(pn)^{-2}$;
$\omega_n(\cdot,\tau)$ is nondecreasing in its first argument, so (E3)
bounds $\omega_n(r,\tau_n)$ for every $r\le C_0h$.

The proofs below never use the whole range $w\in[\tau,1-\tau]$: only the
$k+1$ largest responses of a window enter a Hill statistic, and
Lemma~\ref{lem:app-tailfreeze} shows that on the working event these
correspond to $w\le3\alpha$.  We therefore also write, for
$0<\tau\le3\alpha<1$,
\begin{equation}\label{eq:omega-up}
 \omega_n^{\uparrow}(r,\tau,\alpha)=\max_{j\le p}
 \sup_{u\in\I}\;
 \sup_{\substack{v\in[0,1]\\|v-u|\le r}}\;
 \sup_{w\in[\tau,3\alpha]}
 \left|\log\frac{Q_j(w^{-1},v)}{Q_j(w^{-1},u)}\right|
 \;\le\;\omega_n(r,\tau),
\end{equation}
and record the tail-only condition
\begin{equation}\label{eq:E3up}
 \omega_n^{\uparrow}(C_0h,\tau_n,\alpha)=o(\Delta_{\min}),
 \qquad\tau_n=(pn)^{-2},
 \tag{E3$^{\uparrow}$}
\end{equation}
which is implied by (E3) and suffices everywhere below.  Every result of
this appendix and of Appendix~\ref{app:proofs} is stated with
$\omega_n^{\uparrow}$; Appendix~\ref{app:primitive} gives primitive
sufficient conditions for (E1), (E2) and (E3$^{\uparrow}$).  Put
$\delta_n=(n+1)^{-1}$,
\[
 m_n^-=\bigl\lfloor\{h+\min(h,\varepsilon)\}(n+1)\bigr\rfloor-1,
 \qquad
 K_n=\lfloor\alpha m_n^-\rfloor,
 \qquad
 S_n=4n+1.
\]
All conditional distributions are the pointwise compatible versions of
Section~\ref{sec:population}, and constants are uniform over the
triangular array.

\subsection{Elementary probabilistic lemmas}\label{app:elementary}

The first lemma transfers pointwise perturbations to order statistics, as
in \citet[Appendix, Lemma~3]{gardesPodgorny2025}.

\begin{lemma}[Sorting is sup-norm nonexpansive]\label{lem:app-sorting}
Let $a_1,\ldots,a_m$ and $b_1,\ldots,b_m$ satisfy $a_i\le b_i$ for all $i$.
Then $a_{(i)}\le b_{(i)}$ for all $i$.  Consequently
$|a_i-b_i|\le\epsilon$ for all $i$ implies
$|a_{(i)}-b_{(i)}|\le\epsilon$ for all $i$.
\end{lemma}

\begin{proof}
Fix $i$.  At least $m-i+1$ indices $r$ satisfy $a_r\ge a_{(i)}$, and for
each such $r$, $b_r\ge a_r\ge a_{(i)}$; hence $b_{(i)}\ge a_{(i)}$.  For
the second claim apply the first to $(a_i-\epsilon,b_i)$ and
$(b_i-\epsilon,a_i)$.
\end{proof}

\begin{lemma}[Chernoff bound for exponential means]\label{lem:app-chernoff}
Let $L_k$ be the mean of $k$ i.i.d.\ standard exponentials.  For
$0<t\le1$,
$\Pp(|L_k-1|>t)\le2e^{-kt^2/4}$.
\end{lemma}

\begin{proof}
Markov's inequality applied to $e^{\lambda kL_k}$, optimized at
$\lambda=t/(1+t)$, gives
$\Pp(L_k>1+t)\le e^{-k\{t-\log(1+t)\}}$.  The function
$\varphi(t)=t-\log(1+t)-t^2/4$ has $\varphi(0)=0$ and
\[
 \varphi'(t)=1-\frac1{1+t}-\frac t2=\frac{t(1-t)}{2(1+t)}\ge0
 \quad\text{on }[0,1],
\]
so $t-\log(1+t)\ge t^2/4$ on $(0,1]$.  For the lower tail,
$\Pp(L_k<1-t)\le e^{-k\{-t-\log(1-t)\}}$ and
$-t-\log(1-t)\ge t^2/2\ge t^2/4$.
\end{proof}

\begin{lemma}[Threshold order statistic of uniforms]\label{lem:app-threshold}
Let $V_{(1)}<\cdots<V_{(m)}$ be the order statistics of $m$ i.i.d.\
standard uniforms, $0<\alpha\le1/3$ and
$k=\lfloor\alpha m\rfloor\ge1$.  Then
$\Pp\{V_{(k+1)}>3\alpha\}\le e^{-(k+1)/4}$.
\end{lemma}

\begin{proof}
If $3\alpha\ge1$ the probability is zero.  Otherwise
$\{V_{(k+1)}>3\alpha\}=\{N\le k\}$ with
$N\sim\mathrm{Bin}(m,3\alpha)$ and mean $\mu=3\alpha m\ge3k$.  The
multiplicative Chernoff bound
$\Pp\{N\le(1-\delta)\mu\}\le e^{-\delta^2\mu/2}$
\citep[Section~2.2]{boucheronLugosiMassart2013} with $\delta=2/3$ gives
$\Pp(N\le\mu/3)\le e^{-2\mu/9}\le e^{-2k/3}\le e^{-(k+1)/4}$, the last
step because $8k\ge3(k+1)$ for $k\ge1$.
\end{proof}

\subsection{Projection lemmas}\label{app:projlem}

Let $\mathcal G\subseteq\sigma(U)$, $T=\gamma(U)$, and
$\xi=\esssup(T\mid\mathcal G)$, the smallest $\mathcal G$-measurable random
variable dominating $T$ almost surely.

\begin{lemma}[Uniform remainder replacement]\label{lem:app-remainder}
Define $N(y)=\E\{c(U)y^{-1/T}\mid\mathcal G\}$.  Under (C1)--(C2), almost
surely,
\begin{equation}
 \left|\frac{\Pp(Y>y\mid\mathcal G)}{N(y)}-1\right|
 \le\rho(y):=\sup_u\left|\frac{L(y,u)}{c(u)}-1\right|\longrightarrow0.
 \label{eq:A1}
\end{equation}
\end{lemma}

\begin{proof}
By the tower property,
$\Pp(Y>y\mid\mathcal G)=\E\{y^{-1/T}L(y,U)\mid\mathcal G\}$.  Writing
$L(y,U)=c(U)\{1+r(y,U)\}$ with $|r(y,U)|\le\rho(y)$,
\[
 |\Pp(Y>y\mid\mathcal G)-N(y)|
 =\bigl|\E\{c(U)y^{-1/T}r(y,U)\mid\mathcal G\}\bigr|
 \le\rho(y)N(y),
\]
and $\rho(y)\to0$ is (C2).
\end{proof}

The next lemma is the point at which the argument departs from the proof
of Proposition~1 in \citet{gardesPodgorny2025}: their condition (C.1)
assumes the law of $\gamma(X)$ absolutely continuous, which identifies
$\{\gamma>t\}$ with $\{\gamma\ge t\}$ in their computation.  The tilting
argument below needs no such assumption and allows purely atomic
conditional laws of $T$ --- the relevant case for inactive coordinates
and constant-index null models.

\begin{lemma}[Tilted-measure concentration]\label{lem:app-tilt}
For $y>1$ define the random probability measure
\[
 K_y^\star(A)=
 \frac{\E\{c(U)y^{-1/T}\mathbf 1_{\{T\in A\}}\mid\mathcal G\}}{N(y)},
 \qquad A\in\mathcal B([\underline\gamma,\Gamma]).
\]
On a common probability-one event, $K_y^\star\Rightarrow\delta_\xi$ weakly
as $y\to\infty$.
\end{lemma}

\begin{proof}
Since $K_y^\star\{T>\xi+\delta\}=0$ for every $\delta>0$ by definition of
$\xi$, it suffices to prove, for fixed
$\delta\in(0,\underline\gamma)$,
\begin{equation}
 K_y^\star\{T\le\xi-\delta\}\longrightarrow0\qquad\text{almost surely.}
 \label{eq:A2}
\end{equation}
Split by the $\mathcal G$-measurable event
$G_\delta=\{\xi\ge\underline\gamma+\delta\}$.  On $G_\delta^c$,
$\xi-\delta<\underline\gamma\le T$, so $\{T\le\xi-\delta\}$ is empty and
\eqref{eq:A2} is trivial.  On $G_\delta$, set $r=\xi-\delta$,
$s=\xi-\delta/2$ and $\pi=\Pp(T>s\mid\mathcal G)$.  Then $\pi>0$ almost
surely on $G_\delta$: otherwise
$\xi'=\xi-(\delta/2)\mathbf 1_{G_\delta\cap\{\pi=0\}}$ would be a
$\mathcal G$-measurable almost-sure upper bound for $T$, strictly smaller
than $\xi$ with positive probability, contradicting minimality of the
conditional essential supremum.  Bounding numerator and denominator
separately, using only $c_-\le c\le c_+$ and monotonicity of
$g\mapsto y^{-1/g}$,
\[
 K_y^\star\{T\le r\}
 \le\frac{c_+y^{-1/r}}{c_-\pi y^{-1/s}}
 =\frac{c_+}{c_-\pi}\,y^{-(1/r-1/s)}\longrightarrow0.
\]
A rational sequence $\delta\downarrow0$ and a countable intersection give
a common probability-one event.
\end{proof}

\begin{lemma}[Continuity of the envelope]\label{lem:app-envelope-cont}
Under (C1) and (S), so that $\xi_j(u)$ equals the fibre maximum
$\max_v\gamma(u,v)$, each $\xi_j$ is continuous on $[0,1]$, with
modulus bounded by the modulus of continuity of $\gamma$.
\end{lemma}

\begin{proof}
$\gamma$ is uniformly continuous on the compact cube; if
$|u-u'|\le\delta$ then $|\gamma(u,v)-\gamma(u',v)|\le\varpi(\delta)$ for
every $v$, and taking maxima preserves the bound.
\end{proof}

\subsection{Empirical-rank geometry}\label{app:rankgeom}

The three lemmas of this subsection are purely deterministic or
rank-based; no property of the response enters.  They replace the
covariate-density conditions (H.1)--(H.2) of
\citet{gardesPodgorny2025} in the present empirical-rank setting.

\begin{lemma}[Rank displacement]\label{lem:app-dkw}
Let $\widehat U_{ij}=R_{ij}/(n+1)$ and
\[
 D_n(\eta)=\Bigl\{\max_{j\le p}\max_{i\le n}
 |U_{ij}-\widehat U_{ij}|\le\eta+\delta_n\Bigr\}.
\]
Then $\Pp\{D_n(\eta)^c\}\le2pe^{-2n\eta^2}$.
\end{lemma}

\begin{proof}
Fix $j$ and let $F_{n,j}$ be the empirical distribution function of the
$U_{ij}$, $i\le n$.  By the Dvoretzky--Kiefer--Wolfowitz inequality with
Massart's constant \citep{massart1990},
$\sup_t|F_{n,j}(t)-t|\le\eta$ outside probability
$2e^{-2n\eta^2}$.  On that event, continuity gives
$F_{n,j}(U_{(r)j})=r/n$, hence $|U_{(r)j}-r/n|\le\eta$ and
\[
 |U_{(r)j}-r/(n+1)|\le\eta+\frac{r}{n(n+1)}\le\eta+\delta_n
\]
for every rank $r$ simultaneously.  A union bound over $j$ finishes; no
union over ranks is needed.
\end{proof}

\begin{lemma}[Deterministic window counts]\label{lem:app-counts}
For every $j$, every realization and every $u\in\I$,
$M_j(u)\ge m_n^-$ and $k_j(u)\ge K_n$; moreover $M_j(u)$ is a nonrandom
function of $u$.
\end{lemma}

\begin{proof}
The ranks form a permutation of $1,\ldots,n$, so
$M_j(u)=\#\{r\in\{1,\ldots,n\}:|r/(n+1)-u|\le h\}$, a deterministic
function of $u$: the number of integers in
$[x(u),y(u)]$ with $x(u)=\max\{(u-h)(n+1),1\}$ and
$y(u)=\min\{(u+h)(n+1),n\}$.  Fix $u\in\I$ and write
$m^-=m_n^-=\lfloor\{h+\min(h,\varepsilon)\}(n+1)\rfloor-1$; note
$h+\min(h,\varepsilon)\le\min(2h,\,h+\varepsilon)$.  Four cases.

If neither endpoint is clipped, the interval has length $2h(n+1)$ and
contains at least $\lfloor2h(n+1)\rfloor-1\ge m^-$ integers.

If only the left endpoint is clipped, then $x(u)=1$ and the count is
$\lfloor y(u)\rfloor$.  Either $y(u)=(u+h)(n+1)\ge(\varepsilon+h)(n+1)$,
so the count is at least
$\lfloor(h+\varepsilon)(n+1)\rfloor\ge m^-+1$; or $y(u)=n$ and the count
is $n\ge\lfloor(h+\varepsilon)(n+1)\rfloor-1\ge m^-$, using
$h\le1-\varepsilon$, i.e.\ $(h+\varepsilon)(n+1)\le n+1$.

If only the right endpoint is clipped, then $y(u)=n$ and the count is
$n-\lceil x(u)\rceil+1\ge n-x(u)\ge n-(1-\varepsilon-h)(n+1)
=(h+\varepsilon)(n+1)-1$, hence at least
$\lceil(h+\varepsilon)(n+1)\rceil-1\ge m^-$.

If both endpoints are clipped, the count is $n$, which is at least
$m^-$ as in the second case.  In every case $M_j(u)\ge m^-$, and
$k_j(u)=\lfloor\alpha M_j(u)\rfloor\ge\lfloor\alpha m_n^-\rfloor=K_n$.
\end{proof}

\begin{lemma}[Finite-state structure]\label{lem:app-states}
For fixed $j$ and fixed ranks, $u\mapsto\mathcal W_j(u)$ takes at most
$S_n=4n+1$ distinct values on $[0,1]$.
\end{lemma}

\begin{proof}
The membership indicator of observation $i$ changes only at
$\widehat U_{ij}\pm h$: at most $2n$ breakpoints, hence at most $2n+1$
constancy intervals plus the $2n$ breakpoint states.
\end{proof}

\begin{lemma}[Local extreme count]\label{lem:app-kn}
If $nh\ge4$ and $\alpha m_n^-\ge2$, then $m_n^-\ge nh/2$ and
\[
 \frac{n\alpha h}{4}\;\le\;K_n\;\le\;2n\alpha h+1 .
\]
\end{lemma}

\begin{proof}
$m_n^-\ge\lfloor h(n+1)\rfloor-1\ge h n-2\ge nh/2$ for $nh\ge4$, and
$K_n=\lfloor\alpha m_n^-\rfloor\ge\alpha m_n^-/2\ge n\alpha h/4$ when
$\alpha m_n^-\ge2$ (since $\lfloor x\rfloor\ge x/2$ for $x\ge2$).  For the upper bound, $m_n^-\le2h(n+1)-1\le2hn+2h-1\le2hn+1$, so
$K_n\le\alpha m_n^-\le2n\alpha h+\alpha\le2n\alpha h+1$.
\end{proof}

The rate conditions \eqref{eq:rate-cond} force $nh\to\infty$ and
$K_n\to\infty$, so Lemma~\ref{lem:app-kn} applies for $n$ large; it is
the only conversion between $K_n$ and $n\alpha h$ used in
Appendix~\ref{app:proofs}.

\subsection{Local Hill representation and uniform profile control}
\label{app:hill}

The scheme parallels the proof of Theorem~3 in
\citet{gardesPodgorny2025} --- oscillation control,
exponential-spacings representation, slowly varying remainder --- with
the covering net over their compact matrix class replaced by exact union
bounds over $p$ coordinates and $S_n$ window states.

\begin{lemma}[Conditional PIT and R\'enyi representation]\label{lem:app-pit}
Fix $j$ and condition on the column $\bm U_j=(U_{1j},\ldots,U_{nj})$.
Under (E1):
\begin{enumerate}
\item $V_{ij}=1-F_{Y|j}(Y_i\mid U_{ij})$, $i\le n$, are i.i.d.\ standard
uniform, and $Y_i=Q_j(1/V_{ij},U_{ij})$ almost surely;
\item for any index set $\mathcal W$ of size $m$, deterministic given
$\bm U_j$, with $V_{(1)}<\cdots<V_{(m)}$ the ordered PIT values of its
members and $k=\lfloor\alpha m\rfloor$,
\[
 L_k=\frac1k\sum_{r=1}^k\log\frac{V_{(k+1)}}{V_{(r)}}
\]
is distributed as the mean of $k$ i.i.d.\ standard exponentials,
independent of $V_{(k+1)}$.
\end{enumerate}
\end{lemma}

\begin{proof}
(i) Given $\bm U_j$, the pairs $(U_{ij},Y_i)$ are independent in $i$ and
$Y_i$ has continuous conditional distribution function
$F_{Y|j}(\cdot\mid U_{ij})$, so $V_{ij}$ is standard uniform.  For the
representation, $Q_j(1/V_{ij},U_{ij})=F_{Y|j}^{-1}(F_{Y|j}(Y_i\mid
U_{ij})\mid U_{ij})\le Y_i$ always, with strict inequality only when
$Y_i$ lies in the interior of an interval on which
$F_{Y|j}(\cdot\mid U_{ij})$ is constant; each such interval has
conditional probability zero and there are countably many, so the
representation holds almost surely.  Strict monotonicity is not needed.
Conditioning on the whole column, rather than on window counts as in
\citet[Appendix, Lemma~4]{gardesPodgorny2025}, makes the window sets
deterministic (Lemma~\ref{lem:app-counts}); only coordinate $j$ is
conditioned on, and no independence across coordinates is used.

(ii) With $W_i=-\log V_i$ standard exponentials and order statistics
$W_{(1)}<\cdots<W_{(m)}$, $W_{(m-r+1)}=-\log V_{(r)}$, Abel summation
gives
\[
 kL_k=\sum_{r=1}^k\{W_{(m-r+1)}-W_{(m-k)}\}
 =\sum_{i=1}^k i\,\{W_{(m-i+1)}-W_{(m-i)}\},
\]
and by the R\'enyi representation of exponential order statistics
\citep{renyi1953} the normalized top spacings
$i\{W_{(m-i+1)}-W_{(m-i)}\}$, $i\le k$, are i.i.d.\ standard exponential
and independent of $W_{(m-k)}$, hence of $V_{(k+1)}$.
\end{proof}

The next lemma is the reason no control of central quantiles is needed:
only the $k+1$ largest responses enter a Hill statistic, and on the
working event they are generated by PIT values below $3\alpha$.

\begin{lemma}[Tail-only freezing of the upper order statistics]
\label{lem:app-tailfreeze}
Fix $0<\tau\le3\alpha<1$ and a window state $\mathcal W$ of size $m$
realized at $u$, with $k=\lfloor\alpha m\rfloor$.  Write
$Z_i=\log Q_j(V_{ij}^{-1},U_{ij})$ and
$\check Z_i=\log Q_j(V_{ij}^{-1},u)$, and let $Z_{[1]}\ge\cdots\ge Z_{[m]}$
and $\check Z_{[1]}\ge\cdots\ge\check Z_{[m]}$ be the descending order
statistics.  On the event $\{\min_{i\in\mathcal W}V_{ij}\ge\tau\}\cap
\{V_{(k+1)}\le q:=3\alpha\}$, if every member satisfies
$|U_{ij}-u|\le r$, then
\begin{equation}\label{eq:tailfreeze}
 \max_{1\le s\le k+1}\bigl|Z_{[s]}-\check Z_{[s]}\bigr|
 \le\omega:=\omega_n^{\uparrow}(r,\tau,\alpha),
\end{equation}
and consequently the corresponding Hill statistics satisfy
$|H-\check H|\le2\omega$.
\end{lemma}

\begin{proof}
Let $T=\{i\in\mathcal W:V_{ij}\le q\}$.  For $i\in T$ we have
$V_{ij}\in[\tau,q]$, so \eqref{eq:omega-up} gives
$|Z_i-\check Z_i|\le\omega$.  For $i\notin T$, monotonicity of
$s\mapsto Q_j(s,\cdot)$ gives
\begin{equation}\label{eq:tailfreeze-out}
 Z_i\le\log Q_j(q^{-1},U_{ij})\le\log Q_j(q^{-1},u)+\omega
 \le\check Z_{[k+1]}+\omega,
\end{equation}
the last inequality because $V_{(k+1)}\le q$ forces at least $k+1$ members
to have $\check Z_i\ge\log Q_j(q^{-1},u)$.  Fix $s\le k+1$.  The indices
carrying the $s$ smallest PIT values lie in $T$, and their frozen values
are at least $\check Z_{[s]}$ by monotonicity of $Q_j(\cdot,u)$; their
actual values are therefore at least $\check Z_{[s]}-\omega$, whence
$Z_{[s]}\ge\check Z_{[s]}-\omega$.  Conversely at most $s-1$ indices of
$T$ have $\check Z_i>\check Z_{[s]}$; every remaining index of $T$ has
$Z_i\le\check Z_{[s]}+\omega$, while \eqref{eq:tailfreeze-out} and
$\check Z_{[k+1]}\le\check Z_{[s]}$ give the same bound outside $T$.
Hence $Z_{[s]}\le\check Z_{[s]}+\omega$, proving
\eqref{eq:tailfreeze}.  The Hill statistic is an average of differences
of the first $k$ order statistics and the $(k+1)$st, so its perturbation
is at most $2\omega$.
\end{proof}

\begin{lemma}[Frozen-location Hill expansion]\label{lem:app-hill}
Fix $j$, condition on $\bm U_j$, and fix a window state $\mathcal W$ of
size $m$ realized at $u\in\I$, with $k=\lfloor\alpha m\rfloor$ and Hill
statistic $H$.  On the event $D_n(\eta)$ intersected with
$\{V_{ij}\ge\tau\ \forall i\}$ and $\{V_{(k+1)}\le3\alpha\}$,
\begin{equation}
 |H-\xi_j(u)|
 \le2\omega_n^{\uparrow}(h+\eta+\delta_n,\tau,\alpha)
 +\xi_j(u)|L_k-1|+b_n(3\alpha)L_k.
 \label{eq:A3app}
\end{equation}
\end{lemma}

\begin{proof}
\emph{Spatial freezing.}  Each $i\in\mathcal W$ has
$|\widehat U_{ij}-u|\le h$, hence $|U_{ij}-u|\le r_n:=h+\eta+\delta_n$ on
$D_n(\eta)$.  Lemma~\ref{lem:app-tailfreeze}, applied with $r=r_n$, gives
$|H-\check H|\le2\omega_n^{\uparrow}(r_n,\tau,\alpha)$, where $\check H$
is the Hill statistic of the frozen responses
$\check Y_i=Q_j(1/V_{ij},u)$.

\emph{Exact decomposition.}  Since $s\mapsto Q_j(s,u)$ is nondecreasing,
$\check Y_{(m-r+1)}=Q_j(1/V_{(r)},u)$, and
$\log Q_j(s,u)=\xi_j(u)\log s+\log\ell_j(s,u)$ yields
\[
 \check H=\xi_j(u)L_k+R,
 \qquad
 R=\frac1k\sum_{r=1}^k
 \log\frac{\ell_j(1/V_{(r)},u)}{\ell_j(1/V_{(k+1)},u)}.
\]

\emph{Remainder.}  On $\{V_{(k+1)}\le3\alpha\}$, each ratio is
$\ell_j(ts,u)/\ell_j(s,u)$ with $s=1/V_{(k+1)}\ge(3\alpha)^{-1}$ and
$t=V_{(k+1)}/V_{(r)}>1$, so
$|\log(\cdot)|\le b_n(3\alpha)\log(V_{(k+1)}/V_{(r)})$; averaging gives
$|R|\le b_n(3\alpha)L_k$.
\end{proof}

The decomposition behind \eqref{eq:A3app} is the shape of the whole
argument:
\[
 \widehat\xi_j(u)-\xi_j(u)
 =\text{rank/spatial error}
 +\text{exponential fluctuation}
 +\text{tail remainder}.
\]

\begin{proposition}[Finite-sample uniform profile bound]
\label{thm:app-rank-profile}
Assume (C1)--(C2), (S), the structural part of (E1) (continuity of
every $F_{Y|j}(\cdot\mid u)$, hence the quantile representation), that
$b_n(3\alpha)$ and $\omega_n^{\uparrow}(\cdot,\cdot,\cdot)$ are finite,
$0<h\le1-\varepsilon$, $0<\alpha\le1/3$, and $K_n\ge2$.  For every
$\eta>0$, $0<t\le1$, $0<\tau\le3\alpha$,
\begin{align}
&\Pp\left\{\max_{j\le p}\sup_{u\in\I}
 |\widehat\xi_j(u)-\xi_j(u)|
 >2\omega_n^{\uparrow}(h+\eta+\delta_n,\tau,\alpha)
 +\Gamma t+(1+t)b_n(3\alpha)\right\}
 \nonumber\\
&\qquad\le
 \underbrace{2pe^{-2n\eta^2}}_{\text{ranks}}
 +\underbrace{pn\tau}_{\text{PIT truncation}}
 +\underbrace{2pS_ne^{-K_nt^2/4}}_{\text{spacings}}
 +\underbrace{pS_ne^{-(K_n+1)/4}}_{\text{threshold}}.
 \label{eq:A4app}
\end{align}
\end{proposition}

\begin{proof}
Call a window state \emph{admissible} for coordinate $j$ if it equals
$\mathcal W_j(u)$ for some $u\in\I$; by
Lemmas~\ref{lem:app-counts} and~\ref{lem:app-states} there are at most
$S_n$ admissible states per coordinate and each has
$k\ge K_n\ge2$, so $L_k$ and $V_{(k+1)}$ are well defined on every one of
them.  (States realized only outside $\I$ may have smaller windows; they
are not needed.)  Define the exceptional events
$T_1=D_n(\eta)^c$;
$T_2=\{\exists\,i,j:V_{ij}<\tau\}$;
$T_3=\{\exists\,j,\text{ admissible state }\mathcal W:|L_k-1|>t\}$;
$T_4=\{\exists\,j,\text{ admissible state }\mathcal W:
V_{(k+1)}>3\alpha\}$.

Outside $T_1\cup T_2\cup T_3\cup T_4$, fix $j$ and $u\in\I$: the realized
window is an admissible state, and
Lemma~\ref{lem:app-hill} with $\xi_j\le\Gamma$, $|L_k-1|\le t$,
$L_k\le1+t$ gives the displayed deviation bound; the modulus
$\omega_n^{\uparrow}(h+\eta+\delta_n,\tau,\alpha)$ covers every selected
observation because its definition allows the displaced point
$v=U_{ij}$ to range over all of $[0,1]$, and, by
Lemma~\ref{lem:app-tailfreeze}, only PIT values in $[\tau,3\alpha]$ are
involved.  The upper truncation $V_{ij}\le1-\tau$ is never used.

Probabilities: $\Pp(T_1)\le2pe^{-2n\eta^2}$ is
Lemma~\ref{lem:app-dkw}; $\Pp(T_2)\le pn\tau$ since each $V_{ij}$ is
uniform.  For $T_3$, $T_4$, condition on $\bm U_j$; the admissible states
are then deterministic with $k\ge K_n$,
Lemma~\ref{lem:app-pit} applies to each, and
Lemmas~\ref{lem:app-chernoff} and~\ref{lem:app-threshold} give conditional
bounds $2e^{-kt^2/4}\le2e^{-K_nt^2/4}$ and
$e^{-(k+1)/4}\le e^{-(K_n+1)/4}$, both decreasing in $k$.  A union bound
over the at most $pS_n$ admissible states and
integration over $\bm U_j$ complete the proof; overlapping states and
dependent coordinates cost nothing beyond the union bound.
\end{proof}

\subsection{From profiles to scores}\label{app:scores}

Recall the score grid $G=\{r\delta_n:r=1,\ldots,n\}\cap\I$, of size
$L_n=|G|$, and that $K_n\ge2$ makes every grid point usable, so
$\widehat\Psi_j=L_n^{-1}\sum_{u\in G}\widehat\xi_j(u)$.

\begin{lemma}[Quadrature]\label{lem:app-quadrature}
If $L_n\ge1$, then for every $j$,
\begin{equation}
 \Bigl|L_n^{-1}\sum_{u\in G}\xi_j(u)-\Psi_j\Bigr|
 \le q_n:=\omega_\xi(\delta_n)+\frac{4\Gamma\delta_n}{1-2\varepsilon}.
 \label{eq:A5app}
\end{equation}
\end{lemma}

\begin{proof}
Write the grid points $u_1<\cdots<u_{L_n}$, consecutive multiples of
$\delta_n$.  If $L_n=1$, then $\I$ contains a single multiple of
$\delta_n$, which forces $|\I|<2\delta_n$; the left side is at most
$\Gamma$, while $4\Gamma\delta_n/|\I|>2\Gamma$, so \eqref{eq:A5app}
holds trivially.  Assume $L_n\ge2$ and partition $\I$
exactly: $C_1=[\varepsilon,u_1+\delta_n/2]$,
$C_{L_n}=[u_{L_n}-\delta_n/2,1-\varepsilon]$, and
$C_r=[u_r-\delta_n/2,u_r+\delta_n/2]$ otherwise.  Interior cells have
length $\delta_n$, boundary cells length in $[\delta_n/2,3\delta_n/2)$,
and every point of $C_r$ is within $\delta_n$ of $u_r$.  Then
\[
 \Bigl|\frac1{L_n}\sum_r\xi_j(u_r)-\frac1{|\I|}\int_{\I}\xi_j\Bigr|
 \le\frac1{|\I|}\sum_r\int_{C_r}|\xi_j-\xi_j(u_r)|
 +\Gamma\sum_r\Bigl|\frac{\lambda(C_r)}{|\I|}-\frac1{L_n}\Bigr|.
\]
The first term is at most $\omega_\xi(\delta_n)$.  For the second, both
weight vectors sum to one;
$|\I|=\sum_r\lambda(C_r)\in[(L_n-1)\delta_n,(L_n+1)\delta_n]$ gives
$||\I|/L_n-\delta_n|\le\delta_n/L_n$, so interior cells contribute at most
$\delta_n$ in total and the two boundary cells at most
$\delta_n/2+\delta_n/L_n\le\delta_n$ each; the sum is at most $3\delta_n$,
and the second term is at most $3\Gamma\delta_n/|\I|\le
4\Gamma\delta_n/|\I|$.
\end{proof}

\begin{proposition}[Finite-sample score bound]\label{thm:app-score}
Under the conditions of Proposition~\ref{thm:app-rank-profile} (with in
addition $\omega_\xi(\delta_n)<\infty$) and $L_n\ge1$, for every
$\eta>0$, $0<t\le1$, $0<\tau<1/2$,
\begin{align}
 \Pp\Bigl\{\max_{j\le p}|\widehat\Psi_j-\Psi_j|
 >B_n(\eta,t,\tau)+q_n\Bigr\}
 \le
 2pe^{-2n\eta^2}+pn\tau+2pS_ne^{-K_nt^2/4}+pS_ne^{-(K_n+1)/4},
 \label{eq:A6app}
\end{align}
where
$B_n(\eta,t,\tau)=2\omega_n^{\uparrow}(h+\eta+\delta_n,\tau,\alpha)+\Gamma t
+(1+t)b_n(3\alpha)$.
\end{proposition}

\begin{proof}
Pathwise,
\[
 |\widehat\Psi_j-\Psi_j|
 \le\sup_{u\in\I}|\widehat\xi_j(u)-\xi_j(u)|
 +\Bigl|L_n^{-1}\sum_{u\in G}\xi_j(u)-\Psi_j\Bigr|,
\]
and Lemma~\ref{lem:app-quadrature} bounds the second term by $q_n$ on
every sample path; reusing the same observations for estimation and
averaging therefore costs no extra probability.
Proposition~\ref{thm:app-rank-profile} gives \eqref{eq:A6app}.
\end{proof}

\begin{lemma}[Deterministic score separation]\label{lem:app-separation}
Suppose $\Delta_{\min}>0$ and, on some event,
\[
 \max_{j\le p}\,|\widehat\Psi_j-\Psi_j|<\frac{\Delta_{\min}}{2}.
\]
Then, on that event, every active coordinate has a strictly smaller
score than every inactive coordinate; consequently
$\A\subseteq\widehat\A_d$ for every $d\ge s$, and
$\widehat\A_s=\A$.
\end{lemma}

\begin{proof}
Inactive coordinates have $\Psi_j=\gmax$
(Proposition~\ref{prop:geometry}, proved independently of this lemma
in Appendix~B), hence
$\widehat\Psi_j>\gmax-\Delta_{\min}/2$; active coordinates have
$\Psi_j\le\gmax-\Delta_{\min}$, hence
$\widehat\Psi_j<\gmax-\Delta_{\min}/2$.  Every active score therefore
precedes every inactive score in the increasing order, and the screening
and recovery statements follow from the definition of
$\widehat\A_d$.
\end{proof}
\section{Proofs of the main results}
\label{app:proofs}
\renewcommand{\theequation}{B.\arabic{equation}}
\setcounter{equation}{0}

This appendix proves the propositions of Section~\ref{sec:population}
and the theorems of Section~\ref{sec:estimation} from the auxiliary
results of Appendix~\ref{app:aux}.

\subsection{Proof of Proposition~\ref{prop:projection}}

For fixed $t>0$, with $N$ and $K_y^\star$ as in
Lemmas~\ref{lem:app-remainder} and~\ref{lem:app-tilt} applied to
$\mathcal G$,
\[
 \frac{N(ty)}{N(y)}
 =\int_{[\underline\gamma,\Gamma]}t^{-1/g}\,K_y^\star(dg)
 \longrightarrow t^{-1/\xi}
\]
almost surely, by Lemma~\ref{lem:app-tilt} and boundedness and
continuity of $g\mapsto t^{-1/g}$ on $[\underline\gamma,\Gamma]$.
Lemma~\ref{lem:app-remainder} transfers the limit to
$\Pp(Y>ty\mid\mathcal G)/\Pp(Y>y\mid\mathcal G)$.  Convergence
simultaneously for all $t>0$ follows from rational $t$, monotonicity of
the ratio in $t$, and continuity of $t\mapsto t^{-1/\xi}$.

For the fibre identity \eqref{eq:fibre}, fix $u\in\I$ and apply the above
on the probability space $([0,1]^{p-1},K_j(u,\cdot))$ with the trivial
$\sigma$-field; $\xi_j(u)$ is by definition the essential supremum of
$v\mapsto\gamma(u,v)$ under $K_j(u,\cdot)$.  Let
$m=\max_v\gamma(u,v)$, attained by compactness and continuity.  The
essential supremum is at most $m$.  Conversely, for $\eta>0$ the set
$\{v:\gamma(u,v)>m-\eta\}$ is open and nonempty, so has positive
$K_j(u,\cdot)$-probability under (S); the essential supremum is at least
$m-\eta$ for every $\eta>0$.  \qed

\subsection{Proof of Proposition~\ref{prop:geometry}}

\emph{Claim 1.}  For $j\notin\A$, $\gamma(u,v)=g\{(u,v)_{\A}\}$ does not
involve $u$, and as $v$ ranges over the fibre the active subvector ranges
over all of $[0,1]^{s}$; by \eqref{eq:fibre},
$\xi_j(u)=\max_{w}g(w)=\gmax$ for every $u$.

\emph{Claim 2.}  By \eqref{eq:fibre}, $\xi_j(u)=\gmax$ exactly when the
fibre through $u$ intersects $\M$, i.e.\ when $u\in\pi_j(\M)$.

\emph{Claim 3.}  $u\mapsto\gmax-\xi_j(u)$ is continuous
(Lemma~\ref{lem:app-envelope-cont}), nonnegative, and by Claim~2 vanishes
on $\I$ exactly on $\I\cap\pi_j(\M)$.  A continuous nonnegative function
has positive integral over $\I$ if and only if its positivity set has
positive Lebesgue measure.

\emph{Final statement.}  Write $\M^{\A}$ for the maximizer set of $g$ in
$[0,1]^{s}$; then $\M=\{u:u_{\A}\in\M^{\A}\}$ and, for $j\in\A$,
$\pi_j(\M)$ is the projection of $\M^{\A}$, a finite set when $\M^{\A}$ is
finite; its intersection with $\I$ is Lebesgue-null and Claim~3 applies.
\qed

\subsection{Proof of Proposition~\ref{prop:margin}}

Every $x$ on the fibre $x_j=u$ has
$\|x_{\A}-u^\star_{\A}\|\ge|u-u^\star_j|$, so every point of the fibre
satisfies
$\gamma(x)\le\gmax-c_1\,|u-u^\star_j|^\nu$; taking the fibre maximum
gives $\xi_j(u)\le\gmax-c_1\,|u-u^\star_j|^\nu$.  Averaging
over $u\in\I$ gives the score bound.  For $\nu=2$,
$\E(W-u_j^\star)^2=\mathrm{Var}(W)+(1/2-u_j^\star)^2
=(1-2\varepsilon)^2/12+(1/2-u_j^\star)^2$.  \qed

\subsection{Proof of Theorem~\ref{thm:score}}

Choose $\eta_n=\sqrt{\log(pn)/n}$, $\tau=\tau_n=(pn)^{-2}$, and
$t_n=\sqrt{(4/K_n)\log\{pS_n(pn)\}}$.  The four probability terms in
\eqref{eq:A6app} are then each $o(1)$: the first is
$2pe^{-2\log(pn)}=2p(pn)^{-2}\to0$; the second is $2pn(pn)^{-2}\to0$; the
third is $2pS_n\{pS_n(pn)\}^{-1}=2(pn)^{-1}\to0$; the fourth is
$pS_ne^{-(K_n+1)/4}\to0$ because $\log(pS_n)\le C\log(pn)=o(K_n)$ by the
second condition in \eqref{eq:rate-cond} and
Lemma~\ref{lem:app-kn}.  The same condition makes $t_n\le1$ for large
$n$, as required by Lemma~\ref{lem:app-chernoff}, and
$\Gamma t_n=O(\sqrt{\log(pn)/(n\alpha h)})$, again by
Lemma~\ref{lem:app-kn}.

It remains to check that the bias part of $B_n(\eta_n,t_n,\tau_n)+q_n$ in
Proposition~\ref{thm:app-score} is
$o(\Delta_{\min})$.  First, $\log(pn)=o(nh^2)$ gives
$\eta_n+\delta_n=o(h)$, so $h+\eta_n+\delta_n\le C_0h$ eventually and
$2\omega_n^{\uparrow}(h+\eta_n+\delta_n,\tau_n,\alpha)\le2\omega_n^{\uparrow}(C_0h,\tau_n,\alpha)
=o(\Delta_{\min})$ by (E3) and monotonicity of
$\omega_n^{\uparrow}(\cdot,\tau,\alpha)$.  Second,
$(1+t_n)b_n(3\alpha)\le2b_n(3\alpha)=o(\Delta_{\min})$ by (E2).  Third,
$q_n=\omega_\xi(\delta_n)+O(1/n)$: the first part is
$o(\Delta_{\min})$ by (E1), and
$1/n\le\sqrt{\log(pn)/(n\alpha h)}$ since $\alpha h\le1\le n\log(pn)$, so
the second part is absorbed into the stochastic rate.  Collecting terms,
outside a vanishing probability,
\[
 \max_{j\le p}|\widehat\Psi_j-\Psi_j|
 \le a_n+\beta_n,\qquad
 a_n=O\Bigl(\sqrt{\log(pn)/(n\alpha h)}\Bigr),\quad
 \beta_n=o(\Delta_{\min}),
\]
which is the assertion.  \qed

\subsection{Proof of Theorem~\ref{thm:sis} and
Corollary~\ref{cor:exact}}

Under \eqref{eq:sis-cond},
$\sqrt{\log(pn)/(n\alpha h)}=o(\Delta_{\min})$, so the bound of
Theorem~\ref{thm:score} is below $\Delta_{\min}/2$ for $n$ large,
outside a vanishing probability.  On the complementary event,
Lemma~\ref{lem:app-separation} places every active score strictly before
every inactive score, so $\A\subseteq\widehat\A_d$ for every $d\ge s$,
and $\widehat\A_s=\A$ when $d=s$.  \qed
\section{Primitive sufficient conditions}
\label{app:primitive}
\renewcommand{\theequation}{C.\arabic{equation}}
\setcounter{equation}{0}

Conditions (E1)--(E3) are stated on the projected quantiles. This
appendix gives primitive sufficient conditions --- conditions on
$\gamma$, on the full conditional survival function, and on the fibre
kernels $K_j(u,\cdot)$ --- that imply (E1), (E2) and the tail-only
condition (E3$^\uparrow$) of Appendix~\ref{app:aux}, with explicit
rates, and verifies them for the simulation models of
Section~\ref{sec:simulation}.

All constants and envelopes below are uniform in the triangular array,
in $j\le p_n$, and on the enlarged interval
$I_{\varepsilon/2}=[\varepsilon/2,1-\varepsilon/2]$; scores continue to
be integrated over $\I$. For $x=\iota_j(u,v)$ put
\begin{equation}\label{eq:deficit}
 \xi_j(u)=\esssup_{v\sim K_j(u,\cdot)}\gamma\{\iota_j(u,v)\},
 \qquad
 D_{j,u}(v)=\frac1{\gamma\{\iota_j(u,v)\}}-\frac1{\xi_j(u)}\ge0.
\end{equation}
The reciprocal-index deficit $D_{j,u}$, rather than
$\xi_j(u)-\gamma$, is what makes the Laplace representation below exact.

\subsection{Why (C1)--(C2) do not suffice}

\begin{proposition}[Necessary repairs]\label{prop:app-necessary}
The assumptions (C1)--(C2), even together with (S), imply none of (E1),
(E2), (E3).
\end{proposition}

\begin{proof}
(C1) constrains only the upper tail and does not exclude atoms, so it
does not imply the continuity of the projected conditional distribution
required by (E1).

For (E2), a level bound on a Hall remainder is insufficient. Let
$\bar F(e^t)=\exp\{-t+r(t)\}$, where $r$ is a sum of disjoint smooth
negative bumps centred at $t_m\to\infty$, the $m$th of height
$a_m=e^{-\beta t_m}$, with a decreasing side of width $o(a_m)$ and a
recovery side of slope at most $1/2$. Then $-1+r'(t)<0$, so $\bar F$ is
a valid survival function and $|r(t)|\le e^{-\beta t}$, yet
$r'(t_m)\to-\infty$. With $T(x)=\log Q(e^x)$, quantile inversion gives
\[
 \frac{\mathrm d}{\mathrm dx}\log\{Q(e^x)e^{-x}\}
 =\frac{r'\{T(x)\}}{1-r'\{T(x)\}},
\]
whose absolute value stays away from zero along a subsequence, so the
normalized increment in (E2) need not vanish.

Finally, make all conditional survival functions identical above a fixed
level but let the conditional median below it vary discontinuously with
the conditioning coordinate. This changes none of (C1), (C2), (S), and
none of the tail quantities, yet the printed (E3) fails because its
range $w\in[\tau_n,1-\tau_n]$ contains $w=1/2$. This last obstruction is
an artefact of the statement, not of the method: by
Lemma~\ref{lem:app-tailfreeze} the proofs use only $w\le3\alpha_n$, so
the correct repair is to weaken (E3) to (E3$^\uparrow$) rather than to
impose regularity on the body of the distribution.
\end{proof}

\subsection{The primitive conditions}

\begin{itemize}
\item[(P1)] \emph{Geometry and atomlessness.} For constants
 $0<\underline\gamma\le\Gamma<\infty$ and $L_\gamma<\infty$,
 \[
  \underline\gamma\le\gamma(x)\le\Gamma,
  \qquad
  |\gamma(x)-\gamma(x')|\le L_\gamma\|x-x'\|_\infty .
 \]
 Full fibre support (S) holds, and every conditional law $Y\mid U=x$ is
 atomless.
\item[(P2)] \emph{Differentiable uniform tail approximation.} For
 $t\ge t_0$,
 \[
  \Pp(Y>e^t\mid U=x)=c(x)e^{-t/\gamma(x)}\{1+r(t,x)\},
  \qquad 0<c_-\le c(x)\le c_+<\infty,
 \]
 with $\sup_x|r(t,x)|\le R_0(t)$ and
 $\sup_x|\partial_tr(t,x)|\le R_1(t)$, where $R_0(t)\to0$ and
 $\overline R_1(T)=\sup_{t\ge T}R_1(t)\to0$; eventually
 $R_0(t)\le1/2$. The maps $r(\cdot,x)$ are locally absolutely
 continuous.
\item[(P3)] \emph{Uniform doubling of the fibre-deficit mass.} With
 \[
  H_{j,u}(z)=\int\mathbf 1_{\{D_{j,u}(v)\le z\}}\,
  c\{\iota_j(u,v)\}\,K_j(u,\mathrm dv),
 \]
 there are $z_0,h_0>0$ and $C_D<\infty$ with
 \[
  H_{j,u}(z_0)\ge h_0,
  \qquad
  H_{j,u}(2z)\le C_DH_{j,u}(z)
  \quad(0<z\le z_0/2)
 \]
 uniformly. This allows an atom at $D=0$, includes the degenerate case
 $D\equiv0$ of an inactive coordinate, and rules out only fibres on
 which the probability of coming near the maximal index becomes
 arbitrarily thinner under a fixed rescaling.
\item[(P4)] \emph{Horizontal coupling of the primitive tail.} Whenever
 $u\in\I$, $v\in I_{\varepsilon/2}$ and $|u-v|\le\varepsilon/2$, the
 kernels $K_j(u,\cdot)$ and $K_j(v,\cdot)$ admit a coupling
 $X_u=\iota_j(u,V_u)$, $X_v=\iota_j(v,V_v)$ such that, almost surely and
 for every $t\ge t_0$,
 \[
  \left|\frac1{\gamma(X_v)}-\frac1{\gamma(X_u)}\right|\le L_D|v-u|,
 \]
 \[
  |\log c(X_v)-\log c(X_u)|
  +|\log\{1+r(t,X_v)\}-\log\{1+r(t,X_u)\}|
  \le L_H|v-u|+\zeta_n,
 \]
 with $\zeta_n\ge0$ and $\zeta_n=o(\Delta_{\min})$. This concerns only
 the primitive conditional tail of (P2), never the body.
\end{itemize}

A sufficient construction for (P4) is a coupling with
$\|X_v-X_u\|_\infty\le C_K|v-u|$ together with Lipschitz $1/\gamma$,
$\log c$ and $\log(1+r)$ on the coupled states; (P4) is weaker because it
records only the three quantities actually used.

\subsection{The transfer theorem}

\begin{lemma}[Tilted deficit bound]\label{lem:app-tilted}
Under (P3), for all sufficiently large $t$ the quantity
\[
 B_{j,u}(t)=\int e^{-tD_{j,u}(v)}c\{\iota_j(u,v)\}K_j(u,\mathrm dv)
\]
satisfies $B_{j,u}(t/2)\le C_BB_{j,u}(t)$ for a uniform constant $C_B$,
and consequently
\begin{equation}\label{eq:tiltedmoment}
 \frac{\int D_{j,u}(v)e^{-tD_{j,u}(v)}c\{\iota_j(u,v)\}K_j(u,\mathrm dv)}
      {B_{j,u}(t)}
 \;\le\;\frac{2C_B}{e\,t}
\end{equation}
uniformly in $n,j,u$.
\end{lemma}

\begin{proof}
Suppress $(j,u)$. For $t\ge z_0^{-1}$, restricting the integral to
$\{D\le1/t\}$ gives $B(t)\ge e^{-1}H(1/t)$. Split $B(t/2)$ over
$[0,1/t]$ and the dyadic bands $(2^m/t,2^{m+1}/t]$. Before the bands
reach $z_0$ their total contribution is at most
\[
 H(1/t)\Bigl[1+\sum_{m\ge0}e^{-2^{m-1}}C_D^{m+1}\Bigr],
\]
a convergent series. Iterating the doubling bound backwards from $z_0$
gives $H(1/t)\ge ct^{-\log_2C_D}$ for a uniform $c>0$, so the remaining
part, which starts no lower than $z_0/2$ and is at most
$c_+e^{-tz_0/4}$, is also bounded by a constant times $H(1/t)$.
Combining the two displays proves the doubling property of $B$. For
$z\ge0$, $ze^{-tz}\le(2/e)t^{-1}e^{-tz/2}$; integrating and using the
doubling property gives \eqref{eq:tiltedmoment}.
\end{proof}

\begin{theorem}[Primitive sufficient conditions]\label{thm:app-transfer}
Assume (P1)--(P4), let $\alpha_n\to0$, and suppose eventually
$C_0h_n\le\varepsilon/2$ and $\tau_n=(p_nn)^{-2}\le3\alpha_n<1$. Put
$x_\alpha=\log\{(3\alpha_n)^{-1}\}$ and $c_\gamma=\underline\gamma/2$.
Then, after a common deterministic onset,
\begin{align}
 \omega_\xi(r)&\le L_\gamma r, \label{eq:transfer-E1}\\
 b_n(3\alpha_n)&\le\frac{C_E}{x_\alpha}
 +4\Gamma^2\overline R_1(c_\gamma x_\alpha), \label{eq:transfer-E2}\\
 \omega_n^{\uparrow}(C_0h_n,\tau_n,\alpha_n)
 &\le2\Gamma C_0h_n\bigl[L_H+L_D\Gamma\{2\log(p_nn)+C_Q\}\bigr]
 +2\Gamma\zeta_n, \label{eq:transfer-E3}
\end{align}
where $C_E<\infty$ depends only on $\underline\gamma$, $\Gamma$, $C_B$
and the onset constants, and $C_Q=\max\{0,\log(2c_+)\}$. Consequently
(E1), (E2) and (E3$^\uparrow$) hold as soon as
\begin{equation}\label{eq:transfer-rates}
 n^{-1}=o(\Delta_{\min}),
 \quad
 x_\alpha^{-1}+\overline R_1(c_\gamma x_\alpha)=o(\Delta_{\min}),
 \quad
 h_n\log(p_nn)+\zeta_n=o(\Delta_{\min}).
\end{equation}
For a fixed gap and a Hall derivative envelope of order $O(t^{-1})$ or
smaller, the substantive requirement from (E2) is simply
$\log(1/\alpha_n)\to\infty$.
\end{theorem}

\begin{proof}
\emph{Step 1: (E1).} Under (S), $\xi_j(u)=\max_v\gamma\{\iota_j(u,v)\}$;
comparing the two maxima at the same $v$ in both directions gives
\eqref{eq:transfer-E1}, hence $\omega_\xi(1/n)\le L_\gamma/n$. A mixture
of atomless laws is atomless, because
$\Pp(Y=y\mid U_j=u)=\int\Pp\{Y=y\mid U=\iota_j(u,v)\}K_j(u,\mathrm dv)=0$;
this gives the continuity and the conditional quantile representation
required by (E1).

\emph{Step 2: exact projected factorization.} All derivative identities
below hold almost everywhere, and local absolute continuity justifies
their integrated forms. Let
\[
 \mathcal A_{j,u}(t)=\int e^{-tD_{j,u}(v)}c\{\iota_j(u,v)\}
 \{1+r(t,\iota_j(u,v))\}K_j(u,\mathrm dv),
\]
so that $\bar F_j(e^t\mid u)=e^{-t/\xi_j(u)}\mathcal A_{j,u}(t)$
exactly. With $M_{j,u}(t)=\int D_{j,u}e^{-tD_{j,u}}c\,K_j(u,\mathrm dv)$,
(P2), Lemma~\ref{lem:app-tilted} and $R_0\le1/2$ give
\begin{equation}\label{eq:app-g}
 \left|\frac{\mathcal A_{j,u}'(t)}{\mathcal A_{j,u}(t)}\right|
 \le\frac{\{1+R_0(t)\}M_{j,u}(t)+R_1(t)B_{j,u}(t)}
         {\{1-R_0(t)\}B_{j,u}(t)}
 \le\frac{C_A}{t}+2R_1(t)=:g(t),
\end{equation}
and $g(t)\to0$. Moreover
$(1-R_0)B_{j,u}\le\mathcal A_{j,u}\le(1+R_0)B_{j,u}$, while the same
tilted bound gives $-(\log B_{j,u})'(t)\le C/t$. Integrating from a
fixed $t_1$ at which $R_0(t_1)\le1/2$, and using $D\le\underline\gamma^{-1}$
and $c\ge c_-$ to bound $B_{j,u}(t_1)$ away from zero, yields uniform
positive constants with $ct^{-C}\le\mathcal A_{j,u}(t)\le2c_+$, hence
\begin{equation}\label{eq:app-logA}
 |\log\mathcal A_{j,u}(t)|\le C\log t
\end{equation}
uniformly for large $t$.

\emph{Step 3: inversion and (E2).} Let
$T_{j,u}(x)=\log Q_j(e^x,u)$. Once $g(t)<1/(2\Gamma)$ the map
$t\mapsto\bar F_j(e^t\mid u)$ is strictly decreasing and invertible, and
\begin{equation}\label{eq:app-inv}
 x=\frac{T_{j,u}(x)}{\xi_j(u)}-\log\mathcal A_{j,u}\{T_{j,u}(x)\}.
\end{equation}
By \eqref{eq:app-logA}, $T_{j,u}(x)\ge c_\gamma x$ for all sufficiently
large $x$. Writing $Q_j(e^x,u)=e^{\xi_j(u)x}\ell_j(e^x,u)$ and
differentiating \eqref{eq:app-inv},
\[
 \frac{\mathrm d}{\mathrm dx}\log\ell_j(e^x,u)
 =\frac{\xi_j(u)^2a_{j,u}\{T_{j,u}(x)\}}
        {1-\xi_j(u)a_{j,u}\{T_{j,u}(x)\}},
 \qquad
 a_{j,u}=\frac{\mathcal A_{j,u}'}{\mathcal A_{j,u}},
\]
so \eqref{eq:app-g} and $T_{j,u}(x)\ge c_\gamma x$ give
\[
 \sup_{j,u}\left|\frac{\mathrm d}{\mathrm dx}\log\ell_j(e^x,u)\right|
 \le\frac{C_E}{x}+4\Gamma^2\overline R_1(c_\gamma x).
\]
For $q>1$ the ratio appearing in (E2) is an average of this derivative
over $[\log s,\log(qs)]$; taking $s\ge(3\alpha_n)^{-1}$ and then all
suprema proves \eqref{eq:transfer-E2}.

\emph{Step 4: tail coupling and (E3$^\uparrow$).} For $t\ge t_0$ set
$f_t(x)=c(x)e^{-t/\gamma(x)}\{1+r(t,x)\}$. Under the coupling of (P4),
with $d=|u-v|$,
\[
 e^{-\delta_t}f_t(X_u)\le f_t(X_v)\le e^{\delta_t}f_t(X_u),
 \qquad
 \delta_t=(L_H+tL_D)d+\zeta_n,
\]
so taking expectations gives
$|\log\bar F_j(e^t\mid v)-\log\bar F_j(e^t\mid u)|\le\delta_t$. By
\eqref{eq:app-g}, after increasing the onset if necessary,
$-\frac{\mathrm d}{\mathrm dt}\log\bar F_j(e^t\mid u)\ge1/(2\Gamma)$,
while $\mathcal A_{j,u}\le2c_+$ and \eqref{eq:app-inv} give
$T_{j,u}(x)\le\Gamma(x+C_Q)$. Evaluating the survival comparison at
$t=T_{j,u}(x)$ and using the inverse-slope bound for the $v$-law,
\[
 |T_{j,v}(x)-T_{j,u}(x)|
 \le2\Gamma\bigl[\{L_H+L_D\Gamma(x+C_Q)\}d+\zeta_n\bigr].
\]
Taking $x=\log(1/w)$ with $w\in[\tau_n,3\alpha_n]$: the lower endpoint
$x\ge x_\alpha\to\infty$ places both quantiles beyond the common onset,
while $x\le2\log(p_nn)$; for $d\le C_0h_n\le\varepsilon/2$ the displaced
point lies in $I_{\varepsilon/2}$. Taking suprema proves
\eqref{eq:transfer-E3}.
\end{proof}

\begin{remark}[The body is never used]
Step~4 controls only $w\le3\alpha_n$. This is exactly the range that
Lemma~\ref{lem:app-tailfreeze} shows to be relevant, and it is why (P4)
constrains the primitive tail alone. No coupling of central quantiles is
required anywhere.
\end{remark}

\subsection{Verifying deficit-mass doubling}

\begin{proposition}[A one-step ratio bound implies (P3)]
\label{prop:app-doubling}
Suppose that, uniformly in $n,j,u$, either $H_{j,u}(0)\ge h_0>0$, or
\begin{equation}\label{eq:app-margin}
 H_{j,u}(z)=a_{j,u}z^{\kappa_{j,n}}L_{j,u}(1/z)\{1+o(1)\},
 \qquad z\downarrow0,
\end{equation}
with $a_{j,u}$ bounded above and away from zero, $0\le\kappa_{j,n}\le
\kappa_+<\infty$, the $o(1)$ uniform, $L_{j,u}$ bounded above and away
from zero at one common point beyond the uniform onset, and
\begin{equation}\label{eq:app-onestep}
 \sup_{n,j,u}\frac{L_{j,u}(t/2)}{L_{j,u}(t)}\le C_L
\end{equation}
for all sufficiently large $t$. Then (P3) holds.
\end{proposition}

\begin{proof}
In the atomic case $H(z)\ge h_0$ for every $z>0$, so
$H(2z)\le c_+\le(c_+/h_0)H(z)$. Otherwise
\[
 \frac{H_{j,u}(2z)}{H_{j,u}(z)}
 =2^{\kappa_{j,n}}\frac{L_{j,u}\{1/(2z)\}}{L_{j,u}(1/z)}\{1+o(1)\},
\]
uniformly bounded by \eqref{eq:app-onestep} and $\kappa_{j,n}\le\kappa_+$.
Taking a fixed small $z_0$, the expansion and the uniform lower bound on
$a_{j,u}L_{j,u}(1/z_0)$ give $H_{j,u}(z_0)\ge h_0'>0$.
\end{proof}

A pure two-sided power margin is not necessary, and is in fact false for
the models below: their corner margins carry powers of $\log(1/z)$ and
factors $\exp\{b\sqrt{\log(1/z)}\}$, which are slowly varying and
therefore satisfy \eqref{eq:app-onestep}.

\subsection{Verification for A1--A3 and B1}

The active set in all four simulation models is $\{1,2,3,4\}$ and the
fibre maximum is attained when every free active PIT coordinate equals
zero. The following boundary calculation is the only model-specific
ingredient; its dimension is at most four, independently of $p_n$.

\begin{lemma}[Gaussian corner]\label{lem:app-corner}
Let $X\sim N_d(\mu,\Sigma)$ with $1\le d\le d_0$, let $\Omega=\Sigma^{-1}$,
and let the parameters range over a compact family of positive definite
covariance matrices and bounded means with
$\min_i(\Omega\mathbf 1)_i\ge\nu_0>0$. Suppose that on a fixed
neighbourhood $N_0$ of the lower corner $D$ is uniformly $C^2$ with
\[
 D\{\Phi(X)\}=\sum_{i=1}^da_i\Phi(X_i)+O(\|\Phi(X)\|^2),
 \qquad 0<a_-\le a_i\le a_+<\infty,
\]
and that $c_*\sum_iu_i\le D(u)\le C_*\sum_iu_i$ on $N_0$ with
$\inf_{u\notin N_0}D(u)>0$. Write $\kappa=\mathbf 1^\top\Omega\mathbf 1$,
$\nu=\Omega\mathbf 1$, $b=\mathbf 1^\top\Omega\mu$ and
$r_z=-\Phi^{-1}(z)$. Then, uniformly over the family,
\begin{equation}\label{eq:app-corner}
 \Pp[D\{\Phi(X)\}\le z]
 =C(\mu,\Sigma,a)\,z^{\kappa}\{\sqrt{2\pi}r_z\}^{\kappa-d}e^{-br_z}
 \{1+O(r_z^{-1})\},
\end{equation}
with
$C(\mu,\Sigma,a)=|\Sigma|^{-1/2}e^{-\mu^\top\Omega\mu/2}
\prod_i\Gamma(\nu_i)a_i^{-\nu_i}/\Gamma(\kappa+1)$.
The conclusion persists after multiplication by a positive weight that
is uniformly Lipschitz at the corner and bounded there above and away
from zero. In particular (P3) holds.
\end{lemma}

\begin{proof}
The density of $U=\Phi(X)$ is
\[
 f_U(u)=|\Sigma|^{-1/2}
 \exp\bigl\{-\tfrac12z(u)^\top(\Omega-I)z(u)+z(u)^\top\Omega\mu
 -\tfrac12\mu^\top\Omega\mu\bigr\},
 \quad z_i(u)=\Phi^{-1}(u_i).
\]
Set $u_i=\delta x_i$ and $r=r_\delta$. Uniformly on compact subsets of
$(0,\infty)^d$,
$\Phi^{-1}(\delta x_i)=-r+r^{-1}\log x_i+o(r^{-1})$ and
$e^{-r^2/2}=\sqrt{2\pi}r\delta\{1+O(r^{-2})\}$, so with $\nu=\Omega\mathbf 1$
\[
 f_U(\delta x)\delta^d=|\Sigma|^{-1/2}e^{-\mu^\top\Omega\mu/2}
 \delta^{\kappa}\{\sqrt{2\pi}r\}^{\kappa-d}e^{-br}
 \prod_ix_i^{\nu_i-1}\{1+o(1)\}.
\]
The two-sided bound on $D$ confines $\{D(\delta x)\le\delta\}$ to a
fixed rescaled simplex and excludes points outside $N_0$. Choose
$0<\theta<\tfrac12\min\{1,\nu_0,\underline\lambda\}$ with
$\underline\lambda=\inf\lambda_{\min}(\Omega)$, and split the simplex into
$\mathcal B_r=\{x:\min_ix_i<e^{-\sqrt r}\}$ and its complement. On
$\mathcal B_r$, Mills bounds and convexity of the Gaussian quadratic form
give the envelope $C\prod_ix_i^{\theta-1}$, and a one-dimensional beta
integral bounds the normalized mass of $\mathcal B_r$ by
$Ce^{-\theta\sqrt r}=o(r^{-1})$. On $\mathcal B_r^c$ we have
$|\log x_i|\le\sqrt r$ from below and a uniform bound from above, and the
uniform normal-quantile expansion
\[
 \Phi^{-1}(\delta x_i)=-r+\frac{\log x_i}{r}
 +O\!\left(\frac{1+|\log x_i|+(\log x_i)^2+\log r}{r^3}\right)
\]
shows that the quadratic Gaussian term contributes a relative error
$r^{-1}$ times a fixed polynomial in $1+\sum_i|\log x_i|$, whose
Dirichlet log moments are finite uniformly because $\nu_i\ge\nu_0$. The
quadratic remainder in $D$ moves the limiting simplex boundary by
$O(\delta)=o(r^{-1})$. Integrating the density expansion and using
\[
 \int_{\{\sum_ia_ix_i\le1\}}\prod_ix_i^{\nu_i-1}\,\mathrm dx
 =\frac{\prod_i\Gamma(\nu_i)a_i^{-\nu_i}}{\Gamma(\kappa+1)}
\]
gives \eqref{eq:app-corner}. The factor multiplying $z^\kappa$ is
slowly varying and obeys \eqref{eq:app-onestep}, so
Proposition~\ref{prop:app-doubling} applies. A positive weight
$w(u)=w(0)+O(\|u\|_1)$ changes the rescaled integrand by a relative
$O(\delta)$ and therefore only multiplies the leading constant.
\end{proof}

For the AR(1) covariance of the four active Gaussian coordinates at
$\rho=1/4$, the corner exponents are
\begin{equation}\label{eq:app-kappas}
 \kappa_j=
 \begin{cases}
  34/15, & j=1,4,\\
  41/15, & j=2,3,\\
  14/5+\lambda_j^2/(1-\lambda_j^2), & j\notin\{1,2,3,4\},
 \end{cases}
 \qquad \lambda_j=\rho^{\,j-4}\ (j\ge5),
\end{equation}
so $34/15\le\kappa_j\le43/15$ uniformly; in B1 a proxy coordinate is
independent of the active block and has $\kappa_j=14/5$. All components
of $\Omega\mathbf 1$ are bounded away from zero. The reciprocal deficit
has the required expansion: in A1 and A2 this follows by expanding the
exponential index surface at the active corner, and in A3 the relevant
first derivatives of its exponent are $0.8+0.20u$, $0.8+0.15u$ or $0.8$,
hence uniformly at least $0.8$. In the A family the weight $c$ is
identically one; in B1 the bounded latent scale factor, integrated over
the conditional Gaussian law of $F$, has a uniform expansion
$w_{j,u}(z)=w_{j,u}(0)+O(\|z\|_1)$ with $0<w_-\le w_{j,u}(0)\le w_+$, so
the weighted clause of Lemma~\ref{lem:app-corner} applies.

The remaining primitive assumptions hold as follows. The conditional
generators are strictly increasing and atomless, which gives (P1). For
A1--A3, $q(s,u)=s^{\gamma(u)}/[1+\exp\{v(u)-s\}]$, and both the survival
remainder and its log-level derivative decay faster than every
exponential in $t=\log y$; for B1,
$q(s,u,F)=s^{\gamma(u)}\exp\{A(F)-1/(2s)\}$ with
$A(F)=\{\Phi(F)-1/2\}/2$, and inversion gives
$\Pp(Y>y\mid U=x)=c_B(x)y^{-1/\gamma(x)}\{1+O(y^{-2})\}$ with
$c_B(x)=\E[e^{A(F)/\gamma}\mid U=x]$, the same order holding for the
log-level derivative. This is (P2). A common-residual coupling of the
conditional Gaussian laws moves the other PIT coordinates by
$O(|u-v|)$ on $I_{\varepsilon/2}$, which gives the $1/\gamma$ clause of
(P4); the $c$ and $\log(1+r)$ clauses follow by differentiating the
bounded posterior exponential moments, with $\zeta_n=0$. The projected
profiles are exactly
\[
 \xi_j(u)=
 \begin{cases}
  \tfrac12e^{-u}, & \text{A1, B1},\ j\le4,\\
  \tfrac12e^{-u/2}, & \text{A2},\ j\le4,\\
  \tfrac12e^{-0.8u}, & \text{A3},\ j\le4,\\
  \tfrac12, & j>4.
 \end{cases}
\]

\begin{corollary}[The simulation models satisfy the estimation
conditions]\label{cor:app-models}
For A1--A3 and B1, conditions (E1), (E2) and (E3$^\uparrow$) hold along
every sequence
\[
 \alpha_n=n^{-a},\qquad h_n=\tfrac12n^{-b},\qquad p_n\le n^{C},
\]
with fixed $0<a<1$ and $b,C>0$, with deterministic bias orders
\[
 \omega_\xi(1/n)=O(n^{-1}),
 \qquad
 b_n(3\alpha_n)=O(1/\log n),
 \qquad
 \omega_n^{\uparrow}(C_0h_n,\tau_n,\alpha_n)=O(n^{-b}\log n).
\]
If in addition $b<1/2$ and $a+b<1$, then
$\log(p_nn)=o(nh_n^2)$ and $\log(p_nn)=o(n\alpha_nh_n)$, so the
population gaps of A1--A3 and B1 being fixed and positive,
Theorems~\ref{thm:score} and~\ref{thm:sis} apply. This covers all nine
exponent pairs used by the aggregated screen: intersecting the nine
separation events, every active coordinate has rank at most $s$ and
every inactive coordinate rank at least $s+1$ at every setting, so
taking the minimum of the nine ranks preserves the same strict
separation, and a finite union bound makes the intersection probability
tend to one.
\end{corollary}

\begin{remark}[Finite-design honesty]\label{rem:app-finite}
Corollary~\ref{cor:app-models} is asymptotic. It does not certify that
the $o(\Delta_{\min})$ biases are numerically small at the design used
in Section~\ref{sec:simulation}: there $\log(1/\alpha)\simeq2.28$, the
lower endpoint entering (E2) is only
$x_\alpha=\log\{(3\alpha)^{-1}\}\simeq1.18$, and $h\log(pn)$ is much
larger than the population gaps under the available uniform upper
bounds. The finite-sample tables are empirical evidence about that
regime; they are not a numerical verification of (E2), (E3$^\uparrow$)
or of the concentration bound.
\end{remark}
\section{Discussion}\label{sec:discussion}

The screen is a first stage. Its intended
pipeline reduces $p$ coordinates to a manageable $d$, after which a joint
method --- the central tail-index subspace of \citet{gardesPodgorny2025}, a
structured tail regression \citep{wangTsai2009,sasakiTaoWang2026}, or
multivariate local estimation --- can resolve interactions and correlated
neighbours among the retained variables. The two stages solve different
problems: rapid elimination of evidently irrelevant coordinates scales
to large $p$, while accurate joint tail inference does not.

The scope of a marginal screen for this target is delimited exactly
rather than assumed. An active coordinate every one of whose fibres
meets the global argmax set leaves no trace in the score, at any sample
size; Proposition~\ref{prop:geometry} characterizes this case exactly,
in terms of the argmax projections, where a signal-strength condition
would only have hidden it. Knowing the boundary is what allows a user
to tell a screen that has failed from a target that a marginal screen
cannot see, and it marks where iterative or grouped extensions, with
their own identification arguments, would be required. Two further
features of the target are worth stating in the same spirit. Projected
mixture bias is logarithmic, so (E2) covers every gap with
$\Delta_{\min}\log(1/\alpha)\to\infty$ and in particular every fixed
gap; whether gaps of the excluded order $1/\log n$ are detectable by
any estimator is a question about the target, not about this screen.
And the envelope near its maximum is estimated from observations whose
remaining active coordinates are simultaneously near their maximizing
values, so the sample cost grows with the number of active coordinates:
the weak-signal model A2 in Section~\ref{sec:simulation} measures it,
and aggregation over neighbouring tunings recovers most of it.

Appendix~\ref{app:primitive} derives the estimation conditions from
primitive assumptions that are verifiable on concrete models; what
remains is to establish whether the logarithmic projected bias, which
restricts the covered gaps to $\Delta_{\min}\log(1/\alpha)\to\infty$, is a
property of the target or of this estimator. The score functional also
invites variants: a minimum or a low quantile of the profile would
emphasize depressions confined to short intervals, at the price of new
uniformity arguments for empirical profile extrema.

In practice, tuning matters more than the theory prescribes. The tail
fraction trades Hill bias against variance; the bandwidth additionally
determines how often near-maximizing configurations of the other active
coordinates enter a window. We recommend the workflow of
Sections~\ref{sec:simulation} and~\ref{sec:realdata}: run the screen on
the block of nine neighbouring tunings and rank by the aggregated
minimum rank. The aggregation is deliberately conservative --- a
predictor is retained when at least one reasonable tuning finds
evidence for it --- which is the right direction of error for a
first-stage sure screen: false exclusion is final, while false
retention is corrected by the second-stage multivariate analysis. The
per-setting ranks describe tuning sensitivity and should accompany the
ranking as a diagnostic, not act as a second screening threshold.
Dependence deserves caution of its own: it can help, by concentrating
active configurations, or hurt, by promoting inactive neighbours, and
the direction is design-specific.

Projection turns a tail-index surface into an upper envelope. Variation
of that envelope is a scalable marginal signal for tail-index activity,
and its absence, characterized by the argmax projections, tells the
analyst exactly what a marginal screen can and cannot see.


%% ===== FIN DU MANUSCRIT — DEBUT DU JOURNAL DE CAMPAGNE =====

# Journal — renforcement mathématique (v2)

Objet : obtenir des preuves plus fortes / plus générales que celles de manuscript/main-v2.tex,
avec moins d'hypothèses. Chaque itération est menée par deux sessions GPT en parallèle ;
la vague suivante reçoit ce journal (obstacles + acquis) pour ne pas refaire le même chemin.

Fichier transmis à chaque session : maths/paper-v2-math.tex
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
   détectabilité Δ_j > 0 comme équivalence, sous des hypothèses minimales (sans (S), avec le
   support conditionnel réel).
2. Le cas non détectable (π_j(M) = I_ε, l'exemple diagonal) est réglé : soit un théorème
   d'impossibilité pour toute statistique marginale, soit un score d'ordre supérieur qui le voit,
   avec sa propre théorie.
3. La cible est stable sous dépendance : ce que devient l'enveloppe quand les fibres se déforment,
   énoncé quantitativement (continuité en la copule, pas seulement « ça change »).

**II. L'estimation est optimale, pas seulement suffisante.**
4. Borne inférieure minimax sur une classe contenant les modèles du papier, **appariée** à la borne
   supérieure : le régime log(pn) ≍ nαh Δ_min² est le vrai seuil, constantes comprises si possible,
   sinon à facteur log près explicité. Si le seuil actuel n'est pas optimal, on le remplace.
5. Un seuil de transition net (sharp threshold) : en deçà, aucune procédure ne réussit ; au-delà,
   celle du papier réussit.
6. Les hypothèses de nuisance tombent : (E2), (E3↑) remplacées par une condition de second ordre
   standard (ρ, A(·)) avec compromis biais-variance explicite, et (S) éliminée.
7. **Adaptativité** : choix de h, α, d sans connaître Δ_min, ρ, ni s, avec garantie
   d'oracle. L'agrégation par rangs sur le bloc de neuf tunings reçoit enfin une théorie, au lieu
   d'être une recette.

**III. L'inférence existe.**
8. Loi limite pour Ψ̂_j, et pour le maximum sur p coordonnées (extrême de champ), donnant un test
   de Δ_j = 0 valide uniformément sur la classe.
9. Contrôle d'erreur multiple sur le jeu retenu : FDR ou FWER, non asymptotique de préférence.
10. Correction de biais explicite (Hill à biais réduit local), avec le gain de taux correspondant.

**IV. Le domaine de validité est large.**
11. Dépendance temporelle : mélange (β-mixing) ou champ dépendant, mêmes conclusions, taux dégradé
    explicitement quantifié.
12. γ non continu / non borné inférieurement, ou support de Y avec indice nul : ce qui survit.
13. Robustesse aux contaminations : une fraction ε_n d'observations arbitraires, seuil de rupture.

**V. Tout est vérifié.**
14. Chaque théorème retenu a été confié à une session indépendante dont le seul but est de le casser
    (contre-exemple, hypothèse cachée, circularité) et a survécu.
15. Les énoncés sont compilables en LaTeX, numérotés, cohérents entre eux, et cohérents avec les
    simulations et l'application déjà dans le papier ; les constantes annoncées sont vérifiées
    numériquement au moins une fois.
16. **Saturation** : deux vagues consécutives n'apportent plus ni théorème nouveau, ni hypothèse
    affaiblie, ni contre-exemple.

Tant qu'une ligne de I–IV n'est ni démontrée ni démontrée impossible, la campagne continue.

**Règle de conduite.** Ne pas s'arrêter avant que le cahier soit rempli. Pas d'arrêt pour cause de
temps écoulé, de vague décevante, de session qui rend une réponse partielle, ni de « c'est déjà
bien ». Une vague qui échoue produit au minimum une entrée « obstacle » dans ce journal, et la
vague suivante repart de là. Le seul arrêt légitime est le critère lui-même — ou une instruction
explicite de Luc.

## Itérations

### Vague 0 — 2026-08-16
Connexion navigateur testée (Brave, compte Pro). Aucune session lancée pour l'instant.

### Vague 1 — 2026-08-16 (GPT-5.6 Sol, effort Pro)
Deux sessions parallèles, mêmes pièces jointes (paper-v2-math.tex), consignes volontairement larges.

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
   C_j(u) = supp K_j(u,·) compact et γ continue,
   ξ_j(u) = ess sup_{v∼K_j(u,·)} γ(u,v) = max_{v∈C_j(u)} γ(u,v).
   Preuve (identique des deux côtés) : le max m est atteint sur le compact ; tout voisinage d'un
   maximiseur a masse conditionnelle > 0 ; continuité ⟹ ess sup ≥ m−η pour tout η.
   Puis ξ_j(u) = γ* ⟺ C_j(u) ∩ M_{j,u} ≠ ∅, et
   Δ_j > 0 ⟺ λ{u ∈ I_ε : C_j(u) ∩ M_{j,u} = ∅} > 0, **sans continuité de u ↦ ξ_j(u)**.
   Hypothèse nulle minimale (A : « (NA) », B : « (NS) ») : pour j ∉ A, C_j(u) ∩ M_{j,u} ≠ ∅
   p.p. (S) l'implique et est strictement plus forte.
   → **critère I.1 atteint** (sous réserve de vérification).

2. **L'atomlessness / la continuité de (E1) tombent complètement.** PIT conditionnel *randomisé*
   (transformée distributionnelle) : avec Z_ij ∼ U(0,1) auxiliaires,
   W_ij = F(Y_i−|U_ij) + Z_ij{F(Y_i|U_ij) − F(Y_i−|U_ij)}, V_ij = 1 − W_ij sont i.i.d.
   uniformes conditionnellement à la colonne, et Y_i = Q_j(1/V_ij, U_ij) p.s. Les auxiliaires
   n'entrent jamais dans l'estimateur : pur artefact de preuve. (P1) perd sa clause d'atomlessness.

3. **La cible doit être D = {j : Δ_j > 0}, pas A.** Les théorèmes se réécrivent sans supposer
   que toute coordonnée active est détectable ; D = A devient une condition d'identification
   séparée. Plus honnête et strictement plus fort que la version actuelle.

4. **log(pn) = o(nh²) n'est pas nécessaire.** A : l'union DKW ne porte que sur les p colonnes,
   d'où log(2p) = o(nh²). B : la condition disparaît — on ancre la fenêtre sur la statistique
   d'ordre U_(r)j et on contrôle la largeur par concentration des espacements uniformes
   (U_(b)−U_(a) ∼ Beta(d, N−d), Chernoff), d'où ρ_n(h,x) = h + 2√(hx/n) + 2x/n = h{1+o(1)}
   dès que log(pn) = o(nαh). **B est plus fort ici.**

5. **L'équicontinuité de ξ_j tombe.** A : remplacée par variation bornée (q°_n ≤ C_ε(V_n+Γ)/n).
   B : supprimée purement — la moyenne des ξ_j(U_(r)j) est une intégrale empirique tronquée,
   contrôlée par Hoeffding (Lemme 6), sans aucun module. **B est plus fort ici.**

6. **Le doublement uniforme (P3) tombe**, remplacé par la même quantité des deux côtés : la moyenne
   de déficit *inclinée* m̄(T) = sup_{j,u} sup_{t≥T} M_{j,u}(t)/B_{j,u}(t) → 0
   (A : « (P3-L) », B : « (TD) »). (P3) implique m̄(T) = O(1/T), donc rien n'est perdu côté
   simulations, mais le doublement n'est plus présenté comme nécessaire.

#### Apports propres à A

7. **Nouveau score par blocs de rangs disjoints.** B ≍ 1/h blocs disjoints de m ≍ nh rangs, un Hill
   par bloc sur k ≍ αm extrêmes, score = moyenne des B statistiques. Conditionnellement à la colonne,
   les termes de Rényi sont **indépendants** ⟹ MGF exacte ⟹
   max_j |Ψ̂_j − Ψ_j| = O_P(√(log 2p/(nα)) + log 2p/(nα)), soit **le facteur h gagné**, avec des
   modules seulement *en moyenne* sur les blocs (plus faible que le sup). Condition de screening :
   log(2p) = o(nα Δ_D²) au lieu de log(pn) = o(nαh Δ_min²). Subsiste une condition locale
   log(p/h) = o(nαh) (chaque bloc a besoin d'assez d'extrêmes). A signale honnêtement que les
   simulations actuelles utilisent les fenêtres glissantes : le taux nα ne peut pas leur être
   attribué sans refaire les campagnes.

8. **Recouvrement exact sans connaître s.** γ̂max = max_j Ψ̂_j, Δ̂_j = γ̂max − Ψ̂_j ; tout seuil
   2e < λ < Δ_D − 2e donne {j : Δ̂_j > λ} = D. Supprime l'exigence d = s.

9. **Impossibilité I — la barrière 1/log n est informationnelle sous (C1)–(C2) nues.** Construction
   à deux points : T_n = n^{2γ0}, δ_n = a/(2γ0 log n), lois identiques sous T_n, Pareto
   d'indices différents au-delà. Alors Δ_1,n = a/(4γ0 log n) > 0 mais
   ‖P_0^{⊗n} − P_1^{⊗n}‖_TV ≤ 1/n, donc **aucune procédure** (même multivariée, même non
   marginale) ne détecte uniformément les écarts d'ordre 1/log n. Un seuil d'apparition
   quantifié (type (P2)) est donc *nécessaire*, pas commode.

10. **Impossibilité II — l'invisibilité marginale vaut pour la loi marginale entière.** Avec
    g(t) = γ0 + a cos 2πt, comparer γ0(u) = g(u_2) et γ1(u) = g((u_1+u_2) mod 1) :
    la coordonnée 1 est inactive dans l'un, active dans l'autre, et pourtant la loi de (U_1, Y)
    est *identique*. Toute règle marginale a erreur ≥ 1/2 à tout n. Ce n'est donc pas un défaut du
    score moyenné, ni même de l'indice de queue.

11. **Ordre d'interaction exact.** r*(j) = min{|J| : j ∈ J, Δ_J > 0} ≤ s, borne **atteinte** :
    pour γ(u) = g((u_1+⋯+u_s) mod 1), r*(j) = s pour toute coordonnée active. Donc aucun screen
    de groupe d'ordre fixe ne répare l'invisibilité marginale.
    → **critère I.2 : le cas non détectable est tranché par l'impossibilité** (9–11).

#### Apports propres à B

12. **Impossibilité III — sous dépendance non restreinte, A n'est pas identifiable du tout.**
    U_1 = U_2 = Z ; γ^(1) = 1+u_1 et γ^(2) = 1+(u_1+u_2)/2 coïncident sur la diagonale, donc
    les lois jointes observées sont identiques alors que A^(1) = {1} ≠ {1,2} = A^(2).
    La phrase « unrestricted cross-coordinate dependence » du papier est **incompatible** avec la
    récupération de A. Il faut soit (S)/(NS), soit redéfinir la cible (le classement des Ψ_j),
    soit définir l'activité modulo le support de U.

13. **(C2) est essentielle, mais s'affaiblit en condition de reste *moyennée* (AR).** Preuve de
    projection refaite sous (AR) via la mesure aléatoire inclinée K_y* ⟹ δ_ξ.
    Contre-exemple montrant que la variation régulière *ponctuelle* ne suffit pas : Y = B(u) = u^{-1/β}
    avec prob. 1/2, sinon Pareto(1) ; chaque loi conditionnelle a indice 1 et L ≡ 1 à terme, mais
    le mélange a indice 1/β > 1. La cause est l'absence d'un seuil d'apparition commun.

14. **Le biais spatial passe de h log(pn) à h log(1/α).** Condition horizontale affine en log
    (HQ) : |log{Q_j(s,v)/Q_j(s,u)}| ≤ A(ρ) + B(ρ) log s pour |u−v| ≤ ρ, avec B(ρ) < ξ̲ − b.
    Plus faible que (E3↑) : aucune borne uniforme en s → ∞ n'est exigée, la croissance naturelle
    |ξ_j(v)−ξ_j(u)| log s est permise. Le gel spatial (Lemme 8) donne
    |H − Ȟ| ≤ 2A(ρ) + B(ρ){L_k + 2 log(3/α)}.

15. **La troncature globale min_{i,j} V_ij ≥ (pn)^{-2} disparaît**, remplacée par un contrôle
    bilatéral P{V_(k+1) ∉ [α/3, 3α]} ≤ 2e^{-k/4} (Chernoff multiplicatif dans les deux sens).

16. **Ce que (P3) cachait : la codimension effective.** Si H_{j,u}(z) ≍ a z^κ, alors
    H(2z)/H(z) → 2^κ, donc un doublement uniforme **borne κ uniformément**. Or après projection
    sur une coordonnée active, κ ≈ le nombre d'autres coordonnées actives devant approcher
    simultanément leur maximiseur. Transformée de Laplace : −B'(t)/B(t) ∼ κ/t, donc le biais de
    variation lente est d'ordre **κ/log(1/α)** et la condition honnête est
    κ_n/log(1/α_n) = o(g_n), pas 1/log(1/α_n). **C'est ici que la taille de l'ensemble actif
    entre dans la théorie** — le papier la faisait disparaître par une constante de doublement uniforme.

17. **Théorème fini non asymptotique complet** (Théorème 10 de B) avec les quatre termes de
    probabilité 2pe^{−nε²/32} + (6p+2pL)e^{−x} + 2pL e^{−Kt²/4} + 2pL e^{−K/4} et
    E_n(x,t) = 2A(ρ) + B(ρ){1+t+2log(3/α)} + Γt + (1+t)b + C_ε Γ{√(x/n) + 1/n}.

#### Critiques d'audit (B)

18. Les définitions des modèles de simulation (A1–A3, B1) **n'étaient pas dans le fichier transmis**
    — erreur de ma part, computation_simulation.tex et real_data.tex avaient été omis.
    Corrigé : maths/paper-v2-full.tex les inclut désormais. La vérification de
    Corollary~\ref{cor:app-models} n'a donc pas pu être auditée.
19. Le lemme du coin gaussien est une **esquisse, pas une preuve** : enveloppe intégrable uniforme,
    contrôle uniforme de l'inverse gaussien sur la région mobile, erreur intégrée O(1/r_z),
    stabilité sous la frontière non linéaire, extension pondérée, uniformité sur les familles
    conditionnelles — tout cela est affirmé, pas établi.
20. **L'abstract surinterprète** : il annonce que la dimension polynomiale est abordable dès que le
    compte local d'extrêmes domine log p, alors que le théorème imprimé exige en plus le
    déplacement des rangs et les biais de quantiles projetés, et que le transfert primitif ajoute
    h log(pn).

#### Tension à trancher (le point le plus intéressant de la vague)

A obtient √(log 2p / (nα)), B conserve √(log(pn) / (nαh)). Le facteur h sépare les deux.
A change l'estimateur (blocs disjoints) et gagne l'indépendance entre blocs ; B garde les fenêtres
glissantes. Question ouverte : le facteur h est-il un artefact du chevauchement, ou une vraie
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

- Aucune borne inférieure ne concerne encore le régime principal log p vs nαh : celle de A porte
  sur la barrière 1/log n sous (C1)–(C2) non quantifiées.
- Le score par blocs de A n'a aucune validation numérique.
- La codimension κ_n (point 16) doit être calculée pour les modèles du papier ; si κ_n croît avec
  s, plusieurs énoncés du papier sont à corriger, pas seulement à affaiblir.

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