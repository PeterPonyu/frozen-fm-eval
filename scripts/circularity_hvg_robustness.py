# Gap #1: is the metric-circularity result ("human labels explain LESS expression variance than the PCA
# clustering, 20/20 atlases") an artifact of measuring expr_R2 in a 2000-HVG space (which favors PCA, a
# linear method living in that space)? Test: hold the partitions fixed (PCA-50 KMeans / FM KMeans / labels,
# all built as in the main analysis) and recompute expr_R2 in DIFFERENT evaluation spaces:
#   top-1000 HVG, top-2000 HVG (the paper's setting), top-4000 HVG, and ALL genes (no HVG selection = no bias).
# Count atlases where labels' expr_R2 < PCA-clustering's (and < BOTH FM and PCA clusterings) at each space.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
rng=np.random.RandomState(20260625)
LR={os.path.basename(f)[:-5]:f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
NATIVE={"GSE130148_lung":(".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident"),
        "GSE165784_retina":(".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch"),
        "lung24k":(".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch")}
PICK=["GSE130148_lung","GSE165784_retina","lung24k"]+["lr_"+os.path.basename(f)[:-5] for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad"))]
def load(name):
    if name in NATIVE: f,ct,bt=NATIVE[name]
    else: f=LR[name[3:]]; ct,bt="cell_type","batch"
    A=ad.read_h5ad(f); X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; X=np.log1p(X/tot*1e4)
    return X, A.obs[ct].astype(str).values, A.obs[bt].astype(str).values
def fmemb(name, fm):
    key=f"{fm}_{name}" if name in NATIVE else f"{fm}_lr_{name[3:]}"
    p=f"expand_results/fm_emb/{key}.npz"
    return np.load(p,allow_pickle=True)["X"] if os.path.exists(p) else None
def expr_R2(M, part):
    tot=((M-M.mean(0))**2).sum()
    if tot<=0: return float("nan")
    wit=sum(((M[part==c]-M[part==c].mean(0))**2).sum() for c in np.unique(part))
    return float(1-wit/tot)
SPACES=["hvg1000","hvg2000","hvg4000","allgenes"]
rows=[]
for name in PICK:
    try: X,y0,b=load(name)
    except Exception as e: print("skip",name,str(e)[:40]); continue
    cnt=collections.Counter(y0); keep=np.array([cnt[v]>=10 for v in y0]); X,y0,b=X[keep],y0[keep],b[keep]
    cls=np.unique(y0); y=np.array([np.where(cls==v)[0][0] for v in y0]); NC=len(cls)
    if NC<3: continue
    bv,bc=np.unique(b,return_counts=True); tb=bv[np.argmax(bc)]; te=b==tb; tr=~te
    if te.sum()<150 or tr.sum()<400: continue
    var=X.var(0); order=np.argsort(-var)
    # partitions built exactly as in the main analysis: PCA-50 on top-2000 HVG, KMeans=NC; FM KMeans; labels
    Xh2=X[:,order[:2000]]; sc=StandardScaler().fit(Xh2[tr]); Xhs=sc.transform(Xh2)
    Xp=PCA(min(50,Xhs.shape[1]-1),random_state=0).fit(Xhs[tr]).transform(Xhs)
    parts={"label":y[te], "PCA":KMeans(NC,n_init=5,random_state=0).fit_predict(Xp[te])}
    for fm,nm in (("gf","Geneformer-V2"),("scgpt","scGPT")):
        e=fmemb(name,fm)
        if e is not None and len(e)==len(keep):
            Z=StandardScaler().fit(e[keep][tr]).transform(e[keep])
            parts[nm]=KMeans(NC,n_init=5,random_state=0).fit_predict(Z[te])
    row={"atlas":name,"NC":NC}
    for spc in SPACES:
        if spc=="allgenes": M=X[te]
        else: k=int(spc[3:]); M=X[:,order[:k]][te]
        r2={p:expr_R2(M,parts[p]) for p in parts}
        row[spc]=dict(label=round(r2["label"],4), PCA=round(r2["PCA"],4),
                      lab_lt_PCA=bool(r2["label"]<r2["PCA"]),
                      lab_lt_both=bool(r2["label"]<r2["PCA"] and all(r2["label"]<r2[k] for k in("Geneformer-V2","scGPT") if k in r2)))
    rows.append(row); print(name,"done",flush=True)
# summary counts at each space
N=len(rows)
summ={}
for spc in SPACES:
    lt_pca=sum(r[spc]["lab_lt_PCA"] for r in rows)
    lt_both=sum(r[spc]["lab_lt_both"] for r in rows)
    summ[spc]=dict(labels_lt_PCA=f"{lt_pca}/{N}", labels_lt_both=f"{lt_both}/{N}")
json.dump({"rows":rows,"summary":summ,"N":N},open("expand_results/circularity_hvg_robustness.json","w"),indent=1)
print(f"\n=== CIRCULARITY ROBUSTNESS to evaluation space (N={N} atlases) ===")
print(f"{'eval space':10s} {'labels<PCA-clust':>16} {'labels<both':>12}")
for spc in SPACES: print(f"{spc:10s} {summ[spc]['labels_lt_PCA']:>16} {summ[spc]['labels_lt_both']:>12}")
print("\nallgenes = no HVG selection (removes the disclosed HVG-space bias toward PCA).")
print("DONE")
