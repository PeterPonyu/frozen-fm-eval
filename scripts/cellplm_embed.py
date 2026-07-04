# Self-driven CellPLM (OmicsFormer) zero-shot cell embedder via the vendored CellPLM package.
# Uses CellEmbeddingPipeline.predict (GMVAE latent, 512-d). Atlas gene SYMBOLS are pre-mapped to
# Ensembl IDs (CellPLM's vocab) with Geneformer's gene_name_id dict, so we set ensembl_auto_conversion=
# False (avoids the mygene network dep). numba is stubbed (env numpy 2.5 > numba cap) — only pure-numpy
# scanpy paths are used. Output mirrors the other FMs: fm_emb/cellplm_{atlas}.npz with X/y/batch.
import os, sys, types, numpy as np, anndata as ad, scipy.sparse as sp, pickle, glob as _glob
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
# --- stub numba (scanpy import; numpy 2.5 incompatible with installed numba) ---
nb=types.ModuleType("numba")
def _id(*a,**k): return a[0] if (a and callable(a[0])) else (lambda f: f)
nb.njit=nb.jit=_id; nb.prange=range; nb.vectorize=_id; nb.guvectorize=_id; nb.generated_jit=_id
nb.config=types.SimpleNamespace(NUMBA_NUM_THREADS=1); nb.set_num_threads=lambda *a,**k:None
for sub in ["core","types","typed","extending","core.types"]: sys.modules["numba."+sub]=types.ModuleType("numba."+sub)
nb.types=sys.modules["numba.types"]; sys.modules["numba"]=nb
sys.path.insert(0, os.path.join(os.path.dirname(__file__),"cellplm_repo"))
import torch
from CellPLM.pipeline.cell_embedding import CellEmbeddingPipeline

DIR=os.environ.get("CELLPLM_DIR",".../data/models/CellPLM"); PREFIX="20230926_85M"
GF=".../data/models/Geneformer/geneformer"
N2I=pickle.load(open(f"{GF}/gene_name_id_dict_gc104M.pkl","rb"))  # gene SYMBOL -> Ensembl ID
dev="cuda" if torch.cuda.is_available() else "cpu"
PFX=os.environ.get("CELLPLM_PREFIX","cellplm")

pipe=CellEmbeddingPipeline(pretrain_prefix=PREFIX, pretrain_directory=DIR)
GSET=set(pipe.model.gene_set)
print(f"CellPLM loaded | gene_set={len(GSET)} | dev={dev}",flush=True)

@torch.no_grad()
def embed(adata, ctcol, bcol, out):
    a=adata.copy(); syms=[str(s) for s in a.var_names]
    ens=[N2I.get(s) for s in syms]
    keep=[i for i,e in enumerate(ens) if e is not None and e in GSET]
    if len(keep) < 2000:
        print("SKIP (CellPLM vocab overlap",len(keep),")",out,flush=True); return
    a=a[:, keep].copy(); a.var_names=[ens[i] for i in keep]; a.var_names_make_unique()
    print(out,"mapped genes",len(keep),"of",len(syms),"cells",a.n_obs,flush=True)
    emb=pipe.predict(a, device=dev, ensembl_auto_conversion=False)  # [n,512]
    X=emb.detach().cpu().numpy().astype(np.float32)
    np.savez(out, X=X, y=adata.obs[ctcol].astype(str).values, batch=adata.obs[bcol].astype(str).values)
    print("SAVED",out,X.shape,flush=True)

ATL=[(".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident",f"{PFX}_GSE130148_lung"),
     (".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch",f"{PFX}_GSE165784_retina"),
     (".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch",f"{PFX}_lung24k")]
os.makedirs("expand_results/fm_emb",exist_ok=True)
if len(sys.argv)>1 and sys.argv[1]=="labeled_raw":
    ATL=[(f,"cell_type","batch",f"{PFX}_lr_"+os.path.basename(f)[:-5]) for f in sorted(_glob.glob("expand_results/labeled_raw/*.h5ad"))]
    print("labeled_raw atlases:",len(ATL),flush=True)
for f,ct,bt,name in ATL:
    out=f"expand_results/fm_emb/{name}.npz"
    if os.path.exists(out): print("skip exists",out,flush=True); continue
    try:
        A=ad.read_h5ad(f); embed(A,ct,bt,out)
    except Exception as e:
        import traceback; print("FAIL",name,repr(e)[:120],flush=True); traceback.print_exc()
