# Self-driven UCE (Universal Cell Embedding) zero-shot embedder via the vendored snap-stanford/UCE pipeline.
# Drives eval_single_anndata.main per atlas (4-layer model), reads obsm['X_uce'] (1280-d) -> fm_emb/uce_{atlas}.npz.
# Pre-downloaded figshare model files are symlinked into uce_repo/model_files/ (figshare_download skips existing).
# numba stubbed (env numpy 2.5 > numba cap). Only the 11 human-symbol atlases are embedded (species='human');
# UCE protein embeddings are keyed by human symbol, so mouse-symbol / ENSG atlases are skipped (vocab guard).
import os, sys, types, argparse, numpy as np, anndata as ad, pickle, glob as _glob, shutil
os.environ.pop("ALL_PROXY",None); os.environ.pop("all_proxy",None)
nb=types.ModuleType("numba")
def _id(*a,**k): return a[0] if (a and callable(a[0])) else (lambda f: f)
nb.njit=nb.jit=_id; nb.prange=range; nb.vectorize=_id; nb.guvectorize=_id; nb.generated_jit=_id
nb.config=types.SimpleNamespace(NUMBA_NUM_THREADS=1); nb.set_num_threads=lambda *a,**k:None
for sub in ["core","types","typed","extending","core.types"]: sys.modules["numba."+sub]=types.ModuleType("numba."+sub)
nb.types=sys.modules["numba.types"]; sys.modules["numba"]=nb

ROOT="..."
UCEREPO=os.path.join(ROOT,"scripts","uce_repo"); DL=".../data/models/UCE"
MF=os.path.join(UCEREPO,"model_files"); os.makedirs(MF,exist_ok=True)
# symlink pre-downloaded files into model_files/ (UCE's figshare_download skips existing paths)
for fn in ["4layer_model.torch","all_tokens.torch","species_chrom.csv","species_offsets.pkl"]:
    dst=os.path.join(MF,fn); src=os.path.join(DL,fn)
    if os.path.exists(src) and not os.path.exists(dst): os.symlink(src,dst)
# extract protein embeddings once
if not os.path.isdir(os.path.join(MF,"protein_embeddings")):
    import tarfile; tp=os.path.join(DL,"protein_embeddings.tar.gz")
    print("extracting protein embeddings...",flush=True)
    with tarfile.open(tp) as t: t.extractall(MF)
    # tar may nest under model_files/protein_embeddings or a subdir; normalize
    if not os.path.isdir(os.path.join(MF,"protein_embeddings")):
        for d in _glob.glob(os.path.join(MF,"*protein*")):
            if os.path.isdir(d): os.symlink(d, os.path.join(MF,"protein_embeddings")); break

os.chdir(UCEREPO); sys.path.insert(0, UCEREPO)
GF=".../data/models/Geneformer/geneformer"
HSYM=set(pickle.load(open(f"{GF}/gene_name_id_dict_gc104M.pkl","rb")).keys())  # human gene symbols
OUT=os.path.join(ROOT,"expand_results","fm_emb"); WORK=os.path.join(ROOT,"expand_results","uce_work"); os.makedirs(WORK,exist_ok=True)
from accelerate import Accelerator
from eval_single_anndata import main

def args_for(adata_path):
    return argparse.Namespace(adata_path=adata_path, dir=WORK+"/", species="human", filter=True, skip=True,
        model_loc=os.path.join(MF,"4layer_model.torch"), batch_size=25, pad_length=1536, pad_token_idx=0,
        chrom_token_left_idx=1, chrom_token_right_idx=2, cls_token_idx=3, CHROM_TOKEN_OFFSET=143574,
        sample_size=1024, CXG=True, nlayers=4, output_dim=1280, d_hid=5120, token_dim=5120, multi_gpu=False,
        spec_chrom_csv_path=os.path.join(MF,"species_chrom.csv"), token_file=os.path.join(MF,"all_tokens.torch"),
        protein_embeddings_dir=os.path.join(MF,"protein_embeddings")+"/", offset_pkl_path=os.path.join(MF,"species_offsets.pkl"))

def run(src_h5ad, ctcol, bcol, name):
    out=os.path.join(OUT,f"uce_{name}.npz")
    if os.path.exists(out): print("skip exists",out,flush=True); return
    A=ad.read_h5ad(src_h5ad)
    ov=len(set(str(s) for s in A.var_names)&HSYM)
    if ov<2000: print(f"SKIP (human-symbol overlap {ov})",name,flush=True); return
    work_h5=os.path.join(WORK,f"{name}.h5ad"); A.write(work_h5)
    print(name,"cells",A.n_obs,"overlap",ov,flush=True)
    acc=Accelerator(project_dir=WORK)
    main(args_for(work_h5), acc)
    emb_path=os.path.join(WORK,f"{name}_uce_adata.h5ad")
    E=ad.read_h5ad(emb_path); X=np.asarray(E.obsm["X_uce"],dtype=np.float32)
    np.savez(out, X=X, y=A.obs[ctcol].astype(str).values, batch=A.obs[bcol].astype(str).values)
    print("SAVED",out,X.shape,flush=True)

ATL=[(".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident","GSE130148_lung"),
     (".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch","GSE165784_retina"),
     (".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch","lung24k")]
if len(sys.argv)>1 and sys.argv[1]=="labeled_raw":
    ATL=[(f,"cell_type","batch","lr_"+os.path.basename(f)[:-5]) for f in sorted(_glob.glob(os.path.join(ROOT,"expand_results/labeled_raw/*.h5ad")))]
    print("labeled_raw atlases:",len(ATL),flush=True)
for f,ct,bt,name in ATL:
    try: run(f,ct,bt,name)
    except Exception as e:
        import traceback; print("FAIL",name,repr(e)[:140],flush=True); traceback.print_exc()
