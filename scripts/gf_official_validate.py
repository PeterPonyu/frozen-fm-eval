#!/usr/bin/env python3
# R2 / loader validation for Geneformer-V2: run the OFFICIAL geneformer tokenizer code
# (rank_genes) as a hacked module, re-embed at the official 4096 budget vs my loader's
# 2048, and compare to the cached self-written embedding. Tests whether the only
# difference (token budget) changes the verdict; my 2048 truncation is pessimistic-for-FM.
import os, sys, types, pickle, importlib.util, numpy as np, anndata as ad, scipy.sparse as sp, torch
from pathlib import Path
warn = __import__("warnings"); warn.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
# Official geneformer source (code only; weights/dicts are already local). Fetch once with:
#   GIT_LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/ctheodoris/Geneformer $GF_SRC
GFREPO = os.environ.get("GF_SRC", os.path.expanduser("~/geneformer_src"))
GFDICTS = ".../data/models/Geneformer/geneformer"   # real LFS dicts
MDL = ".../data/models/Geneformer/Geneformer-V2-104M"

# ---- module-hack: stand up minimal geneformer package so official tokenizer.py imports ----
sys.modules["loompy"] = types.ModuleType("loompy")                  # not used (h5ad path)
pkg = types.ModuleType("geneformer"); pkg.__path__ = [GFREPO + "/geneformer"]
for nm, fn in [("GENE_MEDIAN_FILE", "gene_median_dictionary_gc104M.pkl"),
               ("TOKEN_DICTIONARY_FILE", "token_dictionary_gc104M.pkl"),
               ("ENSEMBL_DICTIONARY_FILE", "gene_name_id_dict_gc104M.pkl"),
               ("ENSEMBL_MAPPING_FILE", "ensembl_mapping_dict_gc104M.pkl")]:
    setattr(pkg, nm, Path(GFDICTS) / fn)
for nm in ["GENE_MEDIAN_FILE_30M", "TOKEN_DICTIONARY_FILE_30M", "ENSEMBL_DICTIONARY_FILE_30M", "ENSEMBL_MAPPING_FILE_30M"]:
    setattr(pkg, nm, Path(GFDICTS) / "gene_median_dictionary_gc104M.pkl")  # V1-only; unused for V2
sys.modules["geneformer"] = pkg
spec = importlib.util.spec_from_file_location("geneformer.tokenizer", GFREPO + "/geneformer/tokenizer.py")
tk = importlib.util.module_from_spec(spec); sys.modules["geneformer.tokenizer"] = tk; spec.loader.exec_module(tk)
official_rank_genes = tk.rank_genes
print("imported OFFICIAL geneformer.tokenizer.rank_genes:", official_rank_genes)

# ---- official dicts + model ----
tok = pickle.load(open(f"{GFDICTS}/token_dictionary_gc104M.pkl", "rb"))
med = pickle.load(open(f"{GFDICTS}/gene_median_dictionary_gc104M.pkl", "rb"))
n2i = pickle.load(open(f"{GFDICTS}/gene_name_id_dict_gc104M.pkl", "rb"))
CLS, EOS, PAD = tok["<cls>"], tok["<eos>"], tok["<pad>"]
sys.modules["torchvision"] = None
from transformers.models.bert.modeling_bert import BertModel
dev = "cuda" if torch.cuda.is_available() else "cpu"
model = BertModel.from_pretrained(MDL, output_hidden_states=True, add_pooling_layer=False)
model = (model.to(dev).half() if dev == "cuda" else model).eval()
HID = int(model.config.hidden_size)
print("model loaded", dev, "hidden", HID, flush=True)

def embed_atlas(adata, budget):
    a = adata
    syms = [str(s) for s in a.var_names]; ens = [n2i.get(s) for s in syms]
    gi = [i for i, e in enumerate(ens) if e is not None and e in tok and e in med]
    toks_k = np.array([tok[ens[i]] for i in gi]); med_k = np.array([med[ens[i]] for i in gi], np.float32)
    X = a.X[:, gi]; X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    n = a.n_obs; embs = np.zeros((n, HID), np.float32); B = 16; GB = budget - 2
    for s in range(0, n, B):
        rows, lens = [], []
        for r in range(s, min(s + B, n)):
            x = X[r].toarray().ravel(); totc = x.sum()
            if totc <= 0: rows.append([CLS, EOS]); lens.append(2); continue
            val = (x / totc) * 1e4 / med_k                       # official norm: X/n_counts*1e4/median
            nz = np.nonzero(val)[0]
            ranked = official_rank_genes(val[nz], toks_k[nz])     # OFFICIAL ranking
            ids = [CLS] + ranked[:GB].tolist() + [EOS]; rows.append(ids); lens.append(len(ids))
        L = max(lens); ic = np.full((len(rows), L), PAD, np.int64); am = np.zeros((len(rows), L), np.int64)
        for i, ids in enumerate(rows): ic[i, :len(ids)] = ids; am[i, :len(ids)] = 1
        with torch.no_grad():
            o = model(input_ids=torch.tensor(ic, device=dev), attention_mask=torch.tensor(am, device=dev))
            h = o.hidden_states[-1].float()
        m = torch.tensor(am, device=dev).float().unsqueeze(-1).clone(); m[:, 0, :] = 0   # drop <cls>
        for i, ln in enumerate(lens): m[i, ln - 1, :] = 0                                 # drop <eos>
        embs[s:s + len(rows)] = ((h * m).sum(1) / m.sum(1).clamp(min=1)).cpu().numpy()
    return embs

from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import roc_auc_score
def knn_auroc(Z, y, b):
    bv, bc = np.unique(b, return_counts=True); tb = bv[np.argmax(bc)]; te = b == tb; tr = ~te
    cls = np.unique(y); yi = np.array([np.where(cls == v)[0][0] for v in y]); NC = len(cls)
    Zs = StandardScaler().fit(Z[tr]).transform(Z)
    c = KNeighborsClassifier(15, weights="distance").fit(Zs[tr], yi[tr])
    P = np.zeros((te.sum(), NC)); P[:, c.classes_] = c.predict_proba(Zs[te])
    yb = label_binarize(yi[te], classes=range(NC)); pr = [j for j in range(NC) if 0 < yb[:, j].sum() < len(yb)]
    return float(roc_auc_score(yb[:, pr], P[:, pr], average="macro"))

ATLAS = "bm_all"   # the atlas where FMs trail -> most stringent test
A = ad.read_h5ad(f"expand_results/labeled_raw/{ATLAS}.h5ad")
import collections
y0 = A.obs["cell_type"].astype(str).values; b0 = A.obs["batch"].astype(str).values
cnt = collections.Counter(y0); keep = np.array([cnt[v] >= 10 for v in y0])
A = A[keep].copy(); y0, b0 = y0[keep], b0[keep]
print(f"atlas {ATLAS}: {A.n_obs} cells, {len(np.unique(y0))} types", flush=True)

e2048 = embed_atlas(A, 2048); print("embedded budget=2048", flush=True)
e4096 = embed_atlas(A, 4096); print("embedded budget=4096", flush=True)
cache = np.load(f"expand_results/fm_emb/gf_lr_{ATLAS}.npz", allow_pickle=True)["X"][keep]

def cos(U, V):
    Un = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-9); Vn = V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-9)
    return float((Un * Vn).sum(1).mean())

# PCA reference
Xl = A.X.toarray() if sp.issparse(A.X) else np.asarray(A.X); Xl = np.log1p(Xl / Xl.sum(1, keepdims=True).clip(1) * 1e4)
hv = np.argsort(-Xl.var(0))[:2000]; pca = PCA(50, random_state=0).fit_transform(StandardScaler().fit_transform(Xl[:, hv]))

res = {
    "atlas": ATLAS, "n_cells": int(A.n_obs),
    "cosine_cache_vs_recompute2048": round(cos(cache, e2048), 4),
    "cosine_2048_vs_4096": round(cos(e2048, e4096), 4),
    "knn_auroc_cache2048": round(knn_auroc(cache, y0, b0), 4),
    "knn_auroc_recompute2048": round(knn_auroc(e2048, y0, b0), 4),
    "knn_auroc_official4096": round(knn_auroc(e4096, y0, b0), 4),
    "knn_auroc_PCA": round(knn_auroc(pca, y0, b0), 4),
}
import json; json.dump(res, open("expand_results/gf_official_validate.json", "w"), indent=1)
print("\n=== GENEFORMER-V2 OFFICIAL-CODE VALIDATION ===")
for k, v in res.items(): print(f"  {k:32s} {v}")
print("\n  cache≈recompute (cosine ~1) => my loader reproduces; 4096≥2048 kNN => official budget no worse (my 2048 pessimistic)")
