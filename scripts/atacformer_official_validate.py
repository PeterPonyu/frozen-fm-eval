#!/usr/bin/env python3
# R2 / loader validation for Atacformer: run the OFFICIAL geniml AtacformerModel (a plain
# nn.TransformerEncoder, activation="relu", no positional embeddings) as a hacked module and
# compare to my self-written reconstruction (which used activation="gelu" -- a real bug found
# while reading the official code). Same tokenization; swap only the model. Quantify whether
# the gelu/relu error changed the embedding and the cross-sample cell-type accuracy.
#   Official code: GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/databio/geniml $GENIML_SRC
import os, sys, json, gzip, types, importlib.util, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, torch, torch.nn as nn
import warnings; warnings.filterwarnings("ignore")
os.chdir(os.path.join(os.path.dirname(__file__), ".."))
GENIML = os.environ.get("GENIML_SRC", os.path.expanduser("~/geniml_src"))
CK = ".../data/models/atacformer-base-hg38"
SEED = 20260623; N_CELLS = 20000; MAXLEN = 4096; B = 16; torch.manual_seed(SEED)
dev = "cuda" if torch.cuda.is_available() else "cpu"
c = json.load(open(f"{CK}/config.json")); D, H, L, FF = c["hidden_size"], c["num_attention_heads"], c["num_hidden_layers"], c["intermediate_size"]
V, PAD, CLS = c["vocab_size"], c["pad_token_id"], c["cls_token_id"]

# ---- module-hack: import official AtacformerModel + AtacformerConfig in isolation ----
def load_iso(modname, path, pkg):
    spec = importlib.util.spec_from_file_location(modname, path); m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m; spec.loader.exec_module(m); return m
AF = GENIML + "/geniml/atacformer"
geniml_pkg = types.ModuleType("geniml"); geniml_pkg.__path__ = [GENIML + "/geniml"]; sys.modules["geniml"] = geniml_pkg
af_pkg = types.ModuleType("geniml.atacformer"); af_pkg.__path__ = [AF]; sys.modules["geniml.atacformer"] = af_pkg
load_iso("geniml.atacformer.configuration_atacformer", AF + "/configuration_atacformer.py", "geniml.atacformer")
load_iso("geniml.atacformer.modeling_utils", AF + "/modeling_utils.py", "geniml.atacformer")
load_iso("geniml.atacformer.functional", AF + "/functional.py", "geniml.atacformer")
mod = load_iso("geniml.atacformer.modeling_atacformer", AF + "/modeling_atacformer.py", "geniml.atacformer")
AtacformerModel = mod.AtacformerModel
AtacformerConfig = sys.modules["geniml.atacformer.configuration_atacformer"].AtacformerConfig
from safetensors.torch import load_file
sd = load_file(f"{CK}/model.safetensors")
sd_base = {k[len("atacformer."):]: v for k, v in sd.items() if k.startswith("atacformer.")}  # strip base-model prefix
def build_official():
    # The official class targets transformers 4.51; from_pretrained's finalize step breaks on
    # transformers 5.3 (all_tied_weights_keys). Hack-as-module: instantiate the official class
    # directly (running its __init__ -> relu encoder, no PE) and load weights manually; if even
    # __init__ breaks on the modern API, fall back to a faithful relu replica of the same arch.
    try:
        cfg = AtacformerConfig(**{k: v for k, v in c.items() if k not in ("architectures", "torch_dtype", "transformers_version")})
        mdl = AtacformerModel(cfg); info = mdl.load_state_dict(sd_base, strict=False)
        print("OFFICIAL AtacformerModel class instantiated; missing", len(info.missing_keys), "unexpected", len(info.unexpected_keys), flush=True)
        return mdl
    except Exception as e:
        print("official class init failed (", str(e)[:70], ") -> faithful relu replica", flush=True)
        class Rep(nn.Module):
            def __init__(s):
                super().__init__(); s.embeddings = nn.Module(); s.embeddings.token_embeddings = nn.Embedding(V, D, padding_idx=PAD)
                s.encoder = nn.TransformerEncoder(nn.TransformerEncoderLayer(D, H, FF, dropout=0.0, activation="relu", batch_first=True, norm_first=False, layer_norm_eps=c["norm_eps"]), L)
            def forward(s, ids, mask): return s.encoder(s.embeddings.token_embeddings(ids), src_key_padding_mask=mask)
        r = Rep(); r.load_state_dict(sd_base, strict=False); return r
official = build_official().to(dev).eval()
if dev == "cuda": official = official.half()

# ---- my reconstruction (gelu) ----
class MyModel(nn.Module):
    def __init__(s):
        super().__init__(); s.tok = nn.Embedding(V, D, padding_idx=PAD)
        layer = nn.TransformerEncoderLayer(D, H, FF, dropout=0.0, activation="gelu", batch_first=True, norm_first=False, layer_norm_eps=c["norm_eps"])
        s.encoder = nn.TransformerEncoder(layer, L)
    def forward(s, ids, mask): return s.encoder(s.tok(ids), src_key_padding_mask=mask)
from safetensors.torch import load_file
sd = load_file(f"{CK}/model.safetensors"); remap = {}
for k, v in sd.items():
    if k == "atacformer.embeddings.token_embeddings.weight": remap["tok.weight"] = v
    elif k.startswith("atacformer.encoder."): remap[k.replace("atacformer.encoder.", "encoder.")] = v
mine = MyModel(); mine.load_state_dict(remap, strict=False); mine = mine.to(dev).eval()
if dev == "cuda": mine = mine.half()
print("MY reconstruction loaded (gelu)", flush=True)

# ---- data + peak->universe tokenization (identical to atacformer_embed.py) ----
chrom = []; us = []; ue = []
with gzip.open(f"{CK}/universe.bed.gz", "rt") as f:
    for ln in f:
        p = ln.split()
        if len(p) < 3: continue
        chrom.append(p[0]); us.append(int(p[1])); ue.append(int(p[2]))
chrom = np.array(chrom); us = np.array(us); ue = np.array(ue); utok = np.arange(len(chrom))
chr_idx = {}
for ch in np.unique(chrom):
    mk = np.where(chrom == ch)[0]; o = mk[np.argsort(us[mk])]; chr_idx[ch] = (us[o], ue[o], utok[o])
a = ad.read_h5ad(os.path.expanduser("~/Desktop/data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"))
m = pd.read_csv("raw_pulls/scatac/atac_cell_meta.csv.gz"); m["Barcode"] = m["Barcode"].astype(str); m = m.drop_duplicates("Barcode").set_index("Barcode")
obs = a.obs_names.astype(str); keep = obs.isin(m.index); a = a[np.where(keep)[0]]; md = m.loc[obs[keep]]
X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
df = md.reset_index(); df["_i"] = np.arange(len(df)); frac = min(1.0, N_CELLS / len(df))
sel = df.groupby(["Sample.ID", "Cell.Type"], group_keys=False).apply(lambda g: g.sample(max(1, int(round(len(g) * frac))), random_state=SEED))["_i"].values
sel = np.sort(sel); X = X[sel]; md = md.iloc[sel]
pchr = a.var["chr"].astype(str).values; pst = a.var["start"].astype(np.int64).values; pen = a.var["end"].astype(np.int64).values
peak2tok = np.full(len(pchr), -1, np.int64)
for ch in np.unique(pchr):
    if ch not in chr_idx: continue
    cs, ce, ct = chr_idx[ch]; pm = np.where(pchr == ch)[0]
    for pi in pm:
        ps, pe = pst[pi], pen[pi]; lo = np.searchsorted(cs, ps - 2000); hi = np.searchsorted(cs, pe)
        if hi <= lo: continue
        ov = np.minimum(pe, ce[lo:hi]) - np.maximum(ps, cs[lo:hi]); k = int(np.argmax(ov))
        if ov[k] > 0: peak2tok[pi] = ct[lo + k]
Xc = X.tocsr(); n = Xc.shape[0]
print("matched cells", n, flush=True)

def embed(model):
    emb = np.zeros((n, D), np.float32)
    for s0 in range(0, n, B):
        rows = []
        for r in range(s0, min(s0 + B, n)):
            tk = peak2tok[Xc[r].indices]; tk = tk[tk >= 0]
            if len(tk) == 0: rows.append([CLS]); continue
            tk = np.unique(tk)
            if len(tk) > MAXLEN - 1: tk = tk[:MAXLEN - 1]
            rows.append([CLS] + tk.tolist())
        Lm = max(len(x) for x in rows); ids = np.full((len(rows), Lm), PAD, np.int64); am = np.ones((len(rows), Lm), bool)
        for i, rr in enumerate(rows): ids[i, :len(rr)] = rr; am[i, :len(rr)] = False
        with torch.no_grad():
            h = model(torch.tensor(ids, device=dev), torch.tensor(am, device=dev)).float()
        mm = torch.tensor(~am, device=dev).float().unsqueeze(-1).clone(); mm[:, 0, :] = 0  # drop CLS (both, isolate activation)
        emb[s0:s0 + len(rows)] = ((h * mm).sum(1) / mm.sum(1).clamp(min=1)).cpu().numpy()
    return emb

e_off = embed(official); print("embedded OFFICIAL", flush=True)
e_mine = embed(mine); print("embedded MINE", flush=True)

# ---- cross-sample cell-type probe accuracy (held-out sample) ----
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
y = md["Cell.Type"].astype(str).values; samp = md["Sample.ID"].astype(str).values
def probe_acc(Z):
    sv, sc = np.unique(samp, return_counts=True); ts = sv[np.argmax(sc)]; te = samp == ts; tr = ~te
    Zs = StandardScaler().fit(Z[tr]).transform(Z)
    clf = LogisticRegression(max_iter=300).fit(Zs[tr], y[tr])
    return float((clf.predict(Zs[te]) == y[te]).mean())
def cos(U, Vv):
    Un = U / (np.linalg.norm(U, axis=1, keepdims=True) + 1e-9); Vn = Vv / (np.linalg.norm(Vv, axis=1, keepdims=True) + 1e-9)
    return float((Un * Vn).sum(1).mean())
res = {"n_cells": int(n), "cosine_official_vs_mine": round(cos(e_off, e_mine), 4),
       "probe_acc_official_relu": round(probe_acc(e_off), 4), "probe_acc_mine_gelu": round(probe_acc(e_mine), 4),
       "my_loader_bug": "activation gelu (mine) vs relu (official)"}
json.dump(res, open("expand_results/atacformer_official_validate.json", "w"), indent=1)
print("\n=== ATACFORMER OFFICIAL-CODE VALIDATION ===")
for k, v in res.items(): print(f"  {k:28s} {v}")
