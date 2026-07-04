
import json as _json, math as _math
try:
    import numpy as _np; _FL=(float,_np.floating)
except Exception: _FL=(float,)
_origdump=_json.dump
def _san(o):
    if isinstance(o,_FL): o=float(o); return None if not _math.isfinite(o) else o
    if isinstance(o,dict): return {k:_san(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)): return [_san(v) for v in o]
    return o
def _safedump(o,f,**k): k.pop("allow_nan",None); return _origdump(_san(o),f,allow_nan=False,**k)
_json.dump=_safedump  # NaN/Infinity -> null (valid JSON for R jsonlite)

# Reference-free + non-linear-probe re-check of "baseline beats FM", testing the circularity (labels
# derived from expression space) and the linear-probe bias (geometry!=structure).
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score
rng=np.random.RandomState(20260623)
LR={os.path.basename(f)[:-5]:f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
NATIVE={"GSE130148_lung":(".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident"),
        "GSE165784_retina":(".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch"),
        "lung24k":(".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch")}
MOUSE={"lr_lps_mm","lr_lsk_batch","lr_progastin","lr_urine","lr_astrocytes_sci"}
# representative HUMAN atlases (avoid mouse cross-species + circular louvain)
PICK=["GSE130148_lung","GSE165784_retina","lung24k"]+["lr_"+os.path.basename(f)[:-5] for f in sorted(glob.glob("expand_results/labeled_raw/*.h5ad"))]
def load(name):
    if name in NATIVE: f,ct,bt=NATIVE[name]
    else: f=LR[name[3:]]; ct,bt="cell_type","batch"
    A=ad.read_h5ad(f); X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; X=np.log1p(X/tot*1e4)
    return X, A.obs[ct].astype(str).values, A.obs[bt].astype(str).values
def fmemb(name, fm):  # fm in gf/scgpt
    key = f"{fm}_{name}" if name in NATIVE else f"{fm}_lr_{name[3:]}"
    p=f"expand_results/fm_emb/{key}.npz"
    return np.load(p,allow_pickle=True)["X"] if os.path.exists(p) else None
def expr_R2(Xhvg, part):  # fraction of HVG expression variance explained by a partition
    tot=((Xhvg-Xhvg.mean(0))**2).sum()
    wit=0.0
    for c in np.unique(part):
        m=part==c
        if m.sum()<1: continue
        wit+=((Xhvg[m]-Xhvg[m].mean(0))**2).sum()
    return float(1-wit/tot)
def auroc(Xtr,ytr,Xte,yte,NC,probe):
    if probe=="logreg": c=LogisticRegression(max_iter=300,multi_class="multinomial")
    else: c=KNeighborsClassifier(n_neighbors=15,weights="distance")
    c.fit(Xtr,ytr); P=np.zeros((len(Xte),NC)); P[:,c.classes_]=c.predict_proba(Xte)
    yb=label_binarize(yte,classes=range(NC)); pr=[k for k in range(NC) if 0<yb[:,k].sum()<len(yb)]
    return float(roc_auc_score(yb[:,pr],P[:,pr],average="macro")) if pr else float("nan")
out=[]
for name in PICK:
    try: X,y0,b=load(name)
    except Exception as e: print("skip",name,str(e)[:50]); continue
    cnt=collections.Counter(y0); keep=np.array([cnt[v]>=10 for v in y0]); X,y0,b=X[keep],y0[keep],b[keep]
    cls=np.unique(y0); y=np.array([np.where(cls==v)[0][0] for v in y0]); NC=len(cls)
    if NC<3: continue
    bv,bc=np.unique(b,return_counts=True); tb=bv[np.argmax(bc)]; te=b==tb; tr=~te
    if te.sum()<150 or tr.sum()<400: continue
    var=X.var(0); hv=np.argsort(-var)[:2000]; Xh=X[:,hv]
    sc=StandardScaler().fit(Xh[tr]); Xhs=sc.transform(Xh)
    pca=PCA(min(50,Xhs.shape[1]-1),random_state=0).fit(Xhs[tr]); Xp=pca.transform(Xhs)
    reps={"PCA50":Xp}
    gf=fmemb(name,"gf"); sg=fmemb(name,"scgpt")
    if gf is not None and len(gf)==len(keep): reps["Geneformer-V2"]=StandardScaler().fit(gf[keep][tr]).transform(gf[keep])
    if sg is not None and len(sg)==len(keep): reps["scGPT"]=StandardScaler().fit(sg[keep][tr]).transform(sg[keep])
    spec="mouse" if name in MOUSE else "human"; circ = name=="lung24k"
    row={"atlas":name,"NC":NC,"n":int(len(y)),"species":spec,"circular_louvain":circ}
    # held-out-batch HVG for reference-free expr_R2 (use same Xh on test cells)
    Xh_te=Xh[te]
    row["exprR2_truelabel"]=round(expr_R2(Xh_te, y[te]),4)
    for rn,Z in reps.items():
        row[f"lin_{rn}"]=round(auroc(Z[tr],y[tr],Z[te],y[te],NC,"logreg"),3)
        row[f"knn_{rn}"]=round(auroc(Z[tr],y[tr],Z[te],y[te],NC,"knn"),3)
        km=KMeans(NC,n_init=5,random_state=0).fit_predict(Z[te])
        row[f"exprR2_{rn}"]=round(expr_R2(Xh_te, km),4)
    out.append(row); print(name,"done",flush=True)
json.dump(out,open("expand_results/fair_recheck.json","w"),indent=1)
print("\n=== FAIR RE-CHECK (human atlases) ===")
print(f"{'atlas':18s} {'NC':>3} | {'lin:PCA/GF/scGPT':>22} | {'knn:PCA/GF/scGPT':>22} | {'exprR2 true/PCA/GF/scGPT':>30}")
for r in out:
    lin=f"{r.get('lin_PCA50',0):.2f}/{r.get('lin_Geneformer-V2',float('nan')):.2f}/{r.get('lin_scGPT',float('nan')):.2f}"
    knn=f"{r.get('knn_PCA50',0):.2f}/{r.get('knn_Geneformer-V2',float('nan')):.2f}/{r.get('knn_scGPT',float('nan')):.2f}"
    er=f"{r.get('exprR2_truelabel',0):.3f}/{r.get('exprR2_PCA50',0):.3f}/{r.get('exprR2_Geneformer-V2',float('nan')):.3f}/{r.get('exprR2_scGPT',float('nan')):.3f}"
    print(f"{r['atlas'][:18]:18s} {r['NC']:3d} | {lin:>22} | {knn:>22} | {er:>30}")
