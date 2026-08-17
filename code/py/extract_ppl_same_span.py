"""X = profil d'un PETIT modele, Y = queue d'un GRAND modele, MEME span.

C'est la question de curation posee directement : le profil de surprisal d'un
petit modele permet-il d'identifier les textes ou un grand modele va echouer
catastrophiquement ? Le decoupage prefixe/suffixe repondait a une question plus
dure (predire l'avenir du document) et donnait un resultat nul.

Non-degenerescence : X ne fait intervenir QUE pythia-70m et pythia-160m, Y que
pythia-410m. Aucun quantile de X n'est une fonction du vecteur de NLL qui
definit Y. L'ecart 70m-160m est admis pour la meme raison ; l'ecart avec 410m
serait exclu, puisqu'il contient la NLL dont Y est le maximum.

usage: python code/py/extract_ppl_same_span.py DOCS.parquet OUTDIR [N] [GPU]
"""
import os, sys, time, warnings; warnings.filterwarnings("ignore")
DOCS = sys.argv[1] if len(sys.argv) > 1 else "results/wild/main_docs.parquet"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "results/wild/ppl_same"
NDOC = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
os.environ["CUDA_VISIBLE_DEVICES"] = sys.argv[4] if len(sys.argv) > 4 else "1"
os.environ.setdefault("HF_HOME", "/people/pommeret/.cache/huggingface")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
import numpy as np, pandas as pd, torch
from transformers import AutoModelForCausalLM, AutoTokenizer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_ppl_covariates import profile

BURN, NSPAN = 64, 512
NEED = BURN + NSPAN
BATCH = 8
MODEL_Y = "EleutherAI/pythia-410m"
MODELS_X = ["EleutherAI/pythia-70m", "EleutherAI/pythia-160m"]

os.makedirs(OUTDIR, exist_ok=True)
df = pd.read_parquet(DOCS, columns=["doc_id", "text"])
tok = AutoTokenizer.from_pretrained(MODEL_Y)
ids, toks = [], []
for s in range(0, len(df), 2000):
    part = df.iloc[s:s+2000]
    for e, did in zip(tok(part["text"].tolist(), truncation=True,
                          max_length=NEED)["input_ids"], part["doc_id"].values):
        if len(e) >= NEED:
            toks.append(np.asarray(e[:NEED], dtype=np.int64)); ids.append(int(did))
    if len(toks) >= NDOC: break
toks, ids = toks[:NDOC], ids[:NDOC]
n = len(toks); Tk = np.stack(toks)
print(f"{n} fenetres de {NEED} jetons, span score = {NSPAN}")

def nll_all(mn):
    m = AutoModelForCausalLM.from_pretrained(mn, torch_dtype=torch.float32).cuda().eval()
    out = np.empty((n, NEED-1), dtype=np.float32); t = time.time()
    for s in range(0, n, BATCH):
        x = torch.tensor(Tk[s:s+BATCH], dtype=torch.long, device="cuda")
        with torch.no_grad():
            lp = torch.log_softmax(m(x).logits.float()[:, :-1], dim=-1)
            out[s:s+BATCH] = (-lp.gather(2, x[:, 1:, None]).squeeze(2)).cpu().numpy()
    del m; torch.cuda.empty_cache(); print(f"  {mn} : {time.time()-t:.0f} s")
    return out[:, BURN-1:]

SY = nll_all(MODEL_Y)
Ylog = SY.max(axis=1)
del SY
tid = Tk[:, BURN:BURN+NSPAN]   # out[:,BURN-1:] couvre x[BURN..NEED-1]

# accumulation dans un tableau prealloue : n dictionnaires de 270 entrees
# coutent plusieurs Go a n = 1e5, un tableau float32 coute 100 Mo
store, cols, blocks_out = {}, None, []
for mn in MODELS_X:
    store[mn] = nll_all(mn)
D = store[MODELS_X[0]] - store[MODELS_X[1]]
sources = [("s" + MODELS_X[0].split("-")[-1], store[MODELS_X[0]]),
           ("s" + MODELS_X[1].split("-")[-1], store[MODELS_X[1]]),
           ("d70_160", D)]
mats = []
for tag, S in sources:
    f0 = profile(S[0].astype(np.float64), tid[0], tag)
    ks = list(f0.keys())
    A = np.empty((n, len(ks)), dtype=np.float32)
    A[0] = np.fromiter(f0.values(), dtype=np.float64, count=len(ks))
    for i in range(1, n):
        d_ = profile(S[i].astype(np.float64), tid[i], tag)
        A[i] = np.fromiter((d_[k] for k in ks), dtype=np.float64, count=len(ks))
    mats.append(A); blocks_out += ks
    print(f"  profils {tag} : {A.shape[1]} covariables", flush=True)
F = pd.DataFrame(np.concatenate(mats, axis=1), columns=blocks_out)
F.insert(0, "doc_id", ids); F["Ylog"] = Ylog
bad = F.isna().any(axis=0)
if bad.any(): print(f"  colonnes a NaN ecartees : {int(bad.sum())}"); F = F.loc[:, ~bad]
F.to_parquet(os.path.join(OUTDIR, "ppl_cov.parquet"), index=False)
print(f"ECRIT {OUTDIR}/ppl_cov.parquet : {n} x {F.shape[1]-2} covariables")
print(f"  correlation de rang Y vs s70m_mean : "
      f"{pd.Series(Ylog).corr(F['s70m_mean'], method='spearman'):+.3f}")
print(f"  correlation de rang Y vs s70m_max  : "
      f"{pd.Series(Ylog).corr(F['s70m_max'], method='spearman'):+.3f}")
