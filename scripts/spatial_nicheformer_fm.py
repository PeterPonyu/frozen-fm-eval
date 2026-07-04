#!/usr/bin/env python3
"""Second FM for cluster-H: Nicheformer ZERO-SHOT on the SAME CosMx lymph node
data (no new dataset). HF path (aletlvl/Nicheformer, trust_remote_code) + the
benchmark's committed CosMx technology-mean. Genes mapped symbol->ENSEMBL via HGNC
(6070/20310 vocab overlap, validated). Embeddings -> KMeans(4) -> ARI vs niche.

Env: nfspatial (GPU).  Run: conda run -n nfspatial python scripts/spatial_nicheformer_fm.py
"""
import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
for _v in ("ALL_PROXY","all_proxy","HTTP_PROXY","http_proxy","HTTPS_PROXY","https_proxy"):
    _u = os.environ.get(_v)
    if _u and _u.startswith("socks://"):
        os.environ[_v] = _u.replace("socks://", "socks5://", 1)
import sys, numpy as np, pandas as pd, anndata as ad, torch, scipy.sparse as sp
from transformers import AutoModelForMaskedLM
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SP = os.path.join(ROOT, "raw_pulls/spatial/nicheid")
BM = os.path.join(SP, "spatial-niche-benchmark/benchmark/nicheformer_data/model_means")
MEANS = os.path.join(BM, "cosmx_mean_script.npy")
REF = os.path.join(BM, "model.h5ad")
sys.path.insert(0, os.path.join(SP, "nf_hf"))
from tokenization_nicheformer import NicheformerTokenizer  # custom, bypass AutoTokenizer

# --- load + map genes symbol->ENSEMBL (HGNC) ---
a = ad.read_h5ad(os.path.join(SP, "lymph.h5ad"))
y = a.obs["niche"].astype("category").cat.codes.values
K = int(a.obs["niche"].nunique())
h = pd.read_csv(os.path.join(SP, "hgnc.tsv"), sep="\t", dtype=str, low_memory=False)
sym2ens = {}
for _, r in h.iterrows():
    eid = r.get("ensembl_gene_id")
    if isinstance(eid, str) and eid.startswith("ENSG"):
        for col in ("symbol", "alias_symbol", "prev_symbol"):
            v = r.get(col)
            if isinstance(v, str):
                for s in v.split("|"):
                    sym2ens.setdefault(s, eid)
# Nicheformer tokenizer expects RAW COUNTS (it does sf_normalize+log1p internally).
# lymph.h5ad .X is log-normalized -> use a.raw (raw integer counts, 6519 genes).
raw = a.raw
Xraw = raw.X.toarray() if sp.issparse(raw.X) else np.asarray(raw.X)
raw_genes = list(raw.var_names)
ens = [sym2ens.get(g) for g in raw_genes]
keep = [i for i, e in enumerate(ens) if e is not None]
in_ens = [ens[i] for i in keep]
Xin = Xraw[:, keep].astype(np.float32)

# --- pre-align to the EXACT 20,310 reference genes (zero-fill missing), matching the
#     benchmark's run_nicheformer subsetting so technology_mean (20310) aligns ---
ref_genes = list(ad.read_h5ad(REF).var_names)
col = {g: j for j, g in enumerate(ref_genes)}
Xa = np.zeros((Xin.shape[0], len(ref_genes)), dtype=np.float32)
n_hit = 0
for j, g in enumerate(in_ens):
    c = col.get(g)
    if c is not None:
        Xa[:, c] += Xin[:, j]; n_hit += 1
a = ad.AnnData(X=Xa, obs=a.obs.copy())
a.var_names = ref_genes
a.obs["modality"] = 4; a.obs["specie"] = 5; a.obs["assay"] = 8   # spatial, human, CosMx
print(f"cells={a.n_obs} ref_genes={a.n_vars} genes_filled={n_hit} K={K}")

# --- model (HF) + custom tokenizer (direct instantiation; AutoTokenizer routing is broken) ---
model = AutoModelForMaskedLM.from_pretrained("aletlvl/Nicheformer", trust_remote_code=True)
tok = NicheformerTokenizer(vocab_file=os.path.join(SP, "nf_hf/vocab.json"),
                           technology_mean=np.load(MEANS), max_length=1500, aux_tokens=30)
tok.name_or_path = "aletlvl/Nicheformer"
# our X is already aligned to ref genes in ref order (== technology_mean order);
# disable the tokenizer's internal reference-concat so it does NOT reorder columns.
tok._load_reference_model = lambda: None
dev = "cuda" if torch.cuda.is_available() else "cpu"
model.eval().to(dev)
print(f"[Nicheformer] device={dev}; tokenizing...")
inputs = tok(a)
ids = torch.as_tensor(np.asarray(inputs["input_ids"]))
att = torch.as_tensor(np.asarray(inputs["attention_mask"]))
print(f"[Nicheformer] input_ids {tuple(ids.shape)}; batched embedding...")

embs = []
B = 8
with torch.no_grad():
    for i in range(0, ids.shape[0], B):
        e = model.get_embeddings(input_ids=ids[i:i+B].to(dev),
                                 attention_mask=att[i:i+B].to(dev),
                                 layer=-1, with_context=False)
        embs.append(np.asarray(e.detach().cpu().numpy()))
        del e; torch.cuda.empty_cache()
emb = np.concatenate(embs, axis=0)
if emb.ndim == 3:                       # (cells, tokens, dim) -> mean pool
    emb = emb.mean(axis=1)
import os as _os; _os.makedirs("expand_results/spatial_emb",exist_ok=True)
np.savez("expand_results/spatial_emb/nicheformer.npz", X=np.asarray(emb), niche=a.obs["niche"].astype(str).values, celltype=a.obs["cell_type"].astype(str).values)
print(f"[Nicheformer] embeddings {emb.shape}; saved nicheformer.npz", flush=True)

lab = KMeans(n_clusters=K, n_init=10, random_state=20260623).fit_predict(emb)
ari, nmi = round(float(ARI(y, lab)), 4), round(float(NMI(y, lab)), 4)
print(f"[Nicheformer zero-shot] ARI={ari} NMI={nmi}")

out = os.path.join(ROOT, "cluster_H_spatial_selfcomputed.csv")
df = pd.read_csv(out)
df = df[df.method != "Nicheformer (zero-shot FM, self-computed)"]
df = pd.concat([df, pd.DataFrame([dict(method="Nicheformer (zero-shot FM, self-computed)",
        method_type="spatial-FM", dataset="CosMx_LymphNode", ARI=ari, NMI=nmi)])],
        ignore_index=True).sort_values("ARI", ascending=False)
df.to_csv(out, index=False)
print("\n=== cluster-H table updated ===")
print(df.to_string(index=False))
