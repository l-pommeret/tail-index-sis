"""Amas de nuisance d'echelle a verite EXACTE, greffe sur des donnees reelles.

La longueur est une nuisance d'echelle dont le statut vient de l'approximation
asymptotique P(max > y) ~ theta T c y^(-1/gamma). On construit ici une seconde
nuisance dont le statut est exact.

Lemme de Breiman. Si Y est a variation reguliere d'indice -1/gamma et si V > 0
est bornee et INDEPENDANTE de Y, alors Y*V est a variation reguliere du MEME
indice, et P(YV > y) ~ E[V^(1/gamma)] P(Y > y) : seule l'echelle change. Aucune
approximation, aucune condition asymptotique supplementaire.

Construction, calquee sur le M2 de l'article :

    F ~ N(0,1) latent, independant de tout
    V = exp{ kappa (Phi(F) - 1/2) }              bornee dans [e^-k/2, e^k/2]
    Y' = Y * V                                    meme gamma, echelle deplacee
    z_j = lambda F + sqrt(1-lambda^2) eps_j       s proxys observes de F

Les s coordonnees z_j deplacent donc tous les quantiles finis de Y' et n'entrent
pas dans son indice de queue -- exactement, et sur des donnees reelles ou Y
conserve sa vraie queue et sa vraie structure de dependance.

Un balayage en kappa donne une courbe dose-reponse : a partir de quelle force de
nuisance chaque ecran se laisse-t-il prendre ?

usage: python code/py/inject_scale_cluster.py IN.parquet OUTDIR [S] [LAMBDA]
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

IN = sys.argv[1] if len(sys.argv) > 1 else "results/llm/main_joint.parquet"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "results/llm/cluster"
S = int(sys.argv[3]) if len(sys.argv) > 3 else 20
LAMBDA = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
KAPPAS = [0.0, 0.5, 1.0, 2.0, 4.0]
SEED = 424242

os.makedirs(OUTDIR, exist_ok=True)
df = pd.read_parquet(IN)
n = len(df)
M = df["max_S_small"].values.astype(float)      # M = log Y
rng = np.random.default_rng(SEED)

# facteur latent et ses proxys : tires UNE fois, partages par tous les kappa,
# pour que seule la force de la nuisance change d'un bras a l'autre
f = rng.standard_normal(n)
Z = LAMBDA * f[:, None] + np.sqrt(1 - LAMBDA ** 2) * rng.standard_normal((n, S))
from scipy.stats import norm
u_f = norm.cdf(f)                                # Phi(F), vectorise

meta = {"doc_id", "url", "dump", "date", "argmax_tok", "top_S", "n_tok_doc",
        "T_scored", "max_S_small", "max_S_large", "max_delta",
        "mean_S_small", "mean_S_large"}
keep = [c for c in df.columns if c not in meta and
        pd.api.types.is_numeric_dtype(df[c])]
X = df[keep].copy()
for j in range(S):
    X[f"scale_{j:02d}"] = Z[:, j]
print(f"n = {n}, {len(keep)} covariables reelles + {S} proxys d'echelle "
      f"= {X.shape[1]}")
print(f"lambda = {LAMBDA} (correlation intra-amas {LAMBDA**2:.2f})")

for kap in KAPPAS:
    log_v = kap * (u_f - 0.5)
    Mp = M + log_v                    # log(Y * V) = log Y + log V
    tag = f"k{int(kap*10):03d}"
    out = X.copy()
    out.insert(0, "max_S_small", Mp)
    out.to_parquet(os.path.join(OUTDIR, f"{tag}.parquet"), index=False)
    rho = pd.Series(Z[:, 0]).corr(pd.Series(Mp), method="spearman")
    print(f"  kappa = {kap:4.1f} : sd(log V) = {np.std(log_v):.4f}, "
          f"deplacement median de M = {np.median(log_v):+.4f}, "
          f"rho(proxy_0, M') = {rho:+.3f}  -> {tag}.parquet")

print(f"\ngamma est identique dans les cinq bras par le lemme de Breiman :")
print(f"V est bornee dans [exp(-kappa/2), exp(kappa/2)] et independante de Y.")
