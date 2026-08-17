"""Covariables de surface d'un document, denses et continues par construction.

Contrainte de conception. Le score de l'article ne depend des covariables que
par leurs RANGS, et les fenetres de rangs degenerent des qu'une colonne porte
un atome de masse. Le pretraitement de l'article ecarte deja les variables a
moins de 100 valeurs distinctes ou a atome > 2 %. Toute covariable produite ici
doit donc etre continue sur la plage des documents reels : ce sont des RATIOS
et des MOYENNES, jamais des comptages bruts ni des indicatrices.

C'est ce qui exclut TF-IDF du fond de nuisance : une colonne a 95 % de zeros
produit des ex aequo massifs et casse les comptages deterministes des fenetres.
Les trigrammes de caracteres ne sont retenus que s'ils apparaissent dans une
fraction elevee des documents (df >= DF_MIN), ce qui les rend denses.

Trois blocs :
  surface     ratios de classes de caracteres, ponctuation, mots, phrases,
              repetition, compression, entropie -- interpretables
  trigrammes  frequences de trigrammes de caracteres a forte df -- denses,
              semi-interpretables
  (l'encodeur dense est produit separement par encode_docs.py)

usage: python code/py/doc_covariates.py IN.jsonl OUT.parquet [NPROC]
"""
from __future__ import annotations

import collections
import gzip
import math
import re
import sys
import unicodedata

import numpy as np

DF_MIN = 0.30          # un trigramme doit apparaitre dans >= 30 % des documents
NTRI_MAX = 3000        # plafond sur le nombre de trigrammes retenus
_WORD = re.compile(r"[^\W\d_]+", re.UNICODE)
_SENT = re.compile(r"[.!?]+[\s\"')\]]*")

PUNCT = ".,;:!?'\"()[]{}-–—/\\&%$#@*+=<>|~^`"


def _ratios(t: str) -> dict:
    """Ratios de classes de caracteres. Denominateur = longueur, donc continus."""
    n = len(t)
    if n == 0:
        return {}
    cnt = collections.Counter(t)
    alpha = digit = upper = space = punct = nonascii = ctrl = 0
    for ch, c in cnt.items():
        o = ord(ch)
        if ch.isalpha():
            alpha += c
            if ch.isupper():
                upper += c
        elif ch.isdigit():
            digit += c
        if ch.isspace():
            space += c
        if ch in PUNCT:
            punct += c
        if o > 127:
            nonascii += c
        if o < 32 and ch not in "\n\t\r":
            ctrl += c
    out = {
        "r_alpha": alpha / n, "r_digit": digit / n, "r_upper": upper / max(alpha, 1),
        "r_space": space / n, "r_punct": punct / n, "r_nonascii": nonascii / n,
        "r_ctrl": ctrl / n,
    }
    # ponctuation detaillee : chacune rapportee au total de caracteres
    for name, ch in [("period", "."), ("comma", ","), ("quote", '"'),
                     ("apos", "'"), ("lparen", "("), ("hyphen", "-"),
                     ("colon", ":"), ("semi", ";"), ("excl", "!"),
                     ("quest", "?"), ("slash", "/"), ("amp", "&"),
                     ("pct", "%"), ("star", "*"), ("eq", "="), ("gt", ">")]:
        out["p_" + name] = cnt.get(ch, 0) / n
    # categories unicode agregees : capte mojibake et OCR casse
    cats = collections.Counter(unicodedata.category(ch) for ch in t)
    for c in ("Ll", "Lu", "Nd", "Po", "Zs", "So", "Cf", "Mn", "Lo"):
        out["u_" + c] = cats.get(c, 0) / n
    return out


def _entropy(counts, total) -> float:
    if total <= 0:
        return 0.0
    return -sum((c / total) * math.log(c / total) for c in counts if c > 0)


def _words(t: str) -> dict:
    w = _WORD.findall(t.lower())
    nw = len(w)
    if nw < 2:
        return {"w_meanlen": 0.0, "w_sdlen": 0.0, "w_ttr": 0.0, "w_hapax": 0.0,
                "w_long": 0.0, "w_short": 0.0, "w_entropy": 0.0}
    L = np.fromiter((len(x) for x in w), dtype=np.float64, count=nw)
    c = collections.Counter(w)
    # TTR brut depend de la longueur ; on le corrige par sqrt (indice de Guiraud)
    return {
        "w_meanlen": float(L.mean()), "w_sdlen": float(L.std()),
        "w_ttr": len(c) / math.sqrt(nw),
        "w_hapax": sum(1 for v in c.values() if v == 1) / nw,
        "w_long": float((L >= 10).mean()), "w_short": float((L <= 3).mean()),
        "w_entropy": _entropy(c.values(), nw),
    }


def _struct(t: str) -> dict:
    lines = t.split("\n")
    nl = len(lines)
    LL = np.fromiter((len(x) for x in lines), dtype=np.float64, count=nl)
    sents = [s for s in _SENT.split(t) if s.strip()]
    ns = len(sents)
    SL = (np.fromiter((len(s.split()) for s in sents), dtype=np.float64, count=ns)
          if ns else np.zeros(1))
    return {
        "s_meanline": float(LL.mean()), "s_sdline": float(LL.std()),
        "s_shortline": float((LL < 40).mean()),
        "s_emptyline": float((LL == 0).mean()),
        "s_meansent": float(SL.mean()), "s_sdsent": float(SL.std()),
        "s_sentper1k": 1000.0 * ns / max(len(t), 1),
    }


def _repetition(t: str) -> dict:
    """Fraction de n-grammes de mots dupliques : capte boilerplate et boucles."""
    w = t.lower().split()
    out = {}
    for k in (2, 3, 5, 10):
        if len(w) <= k:
            out[f"rep_{k}gram"] = 0.0
            continue
        g = [" ".join(w[i:i + k]) for i in range(len(w) - k + 1)]
        c = collections.Counter(g)
        out[f"rep_{k}gram"] = 1.0 - len(c) / len(g)
        # part de masse dans le n-gramme le plus frequent
        out[f"rep_{k}top"] = c.most_common(1)[0][1] / len(g)
    return out


def _compress(t: str) -> dict:
    b = t.encode("utf-8", "ignore")
    if not b:
        return {"z_ratio1": 0.0, "z_ratio9": 0.0, "z_gain": 0.0}
    r1 = len(gzip.compress(b, 1)) / len(b)
    r9 = len(gzip.compress(b, 9)) / len(b)
    return {"z_ratio1": r1, "z_ratio9": r9, "z_gain": r1 - r9}


def _char_entropy(t: str) -> dict:
    n = len(t)
    if n < 3:
        return {"h_uni": 0.0, "h_bi": 0.0, "h_cond": 0.0}
    c1 = collections.Counter(t)
    c2 = collections.Counter(t[i:i + 2] for i in range(n - 1))
    h1, h2 = _entropy(c1.values(), n), _entropy(c2.values(), n - 1)
    return {"h_uni": h1, "h_bi": h2, "h_cond": h2 - h1}


def surface_features(t: str) -> dict:
    """Toutes les covariables de surface d'un document."""
    f = {}
    f.update(_ratios(t)); f.update(_words(t)); f.update(_struct(t))
    f.update(_repetition(t)); f.update(_compress(t)); f.update(_char_entropy(t))
    # la longueur est conservee mais c'est une NUISANCE D'ECHELLE connue :
    # elle deplace les quantiles de Y sans entrer dans gamma (cf. AMAS_ECHELLE)
    f["len_chars"] = float(len(t))
    f["len_words"] = float(len(t.split()))
    f["len_lines"] = float(t.count("\n") + 1)
    return f


def trigram_counter(t: str) -> collections.Counter:
    s = re.sub(r"\s+", " ", t.lower())
    return collections.Counter(s[i:i + 3] for i in range(len(s) - 2))


def trigram_features(t: str, keys: list) -> np.ndarray:
    """Frequences relatives des trigrammes retenus. Denses car df elevee."""
    c = trigram_counter(t)
    tot = max(sum(c.values()), 1)
    return np.fromiter((c.get(k, 0) / tot for k in keys),
                       dtype=np.float32, count=len(keys))
