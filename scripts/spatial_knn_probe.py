
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

# Cluster-H fair probe: supervised kNN niche-prediction (spatial-blocked CV) + reference-free
# compo_R2/expr_R2/contiguity, for FM embeddings vs nbhd-composition & PCA baselines.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score
rng=np.random.RandomState(20260623)
a=ad.read_h5ad("raw_pulls/spatial/nicheid/lymph.h5ad")
xy=np.asarray(a.obsm["spatial"]); niche=a.obs["niche"].astype(str).values; ct=a.obs["cell_type"].astype(str).values
X=a.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4); var=Xln.var(0); hv=np.argsort(-var)[:2000]; Xh=Xln[:,hv]
cts=np.unique(ct); ct2i={c:i for i,c in enumerate(cts)}
nn=NearestNeighbors(n_neighbors=20).fit(xy); _,idx=nn.kneighbors(xy)
onehot=np.zeros((len(ct),len(cts)),np.float32); onehot[np.arange(len(ct)),[ct2i[c] for c in ct]]=1
compo=onehot[idx].mean(1)
nlev=np.unique(niche); y=np.array([np.where(nlev==v)[0][0] for v in niche]); K=len(nlev)
def R2(F,part):
    t=((F-F.mean(0))**2).sum(); w=sum(((F[part==c]-F[part==c].mean(0))**2).sum() for c in np.unique(part)); return float(1-w/t)
def contig(part): return float((part[idx[:,1:]]==part[:,None]).mean())
# spatial-blocked 4-fold: grid the tissue
gx=np.digitize(xy[:,0],np.quantile(xy[:,0],[.25,.5,.75])); gy=np.digitize(xy[:,1],np.quantile(xy[:,1],[.25,.5,.75]))
folds=(gx//2)*2+(gy//2)  # 4 contiguous spatial quadrants (prev gx*4+gy %% 4 collapsed to y-stripes)
def knn_probe(Z):
    Zs=StandardScaler().fit_transform(Z); aucs=[]
    for f in range(4):
        te=folds==f; tr=~te
        if len(np.unique(y[tr]))<K or te.sum()<50: continue
        c=KNeighborsClassifier(n_neighbors=30,weights="distance").fit(Zs[tr],y[tr])
        P=np.zeros((te.sum(),K)); P[:,c.classes_]=c.predict_proba(Zs[te])
        yb=label_binarize(y[te],classes=range(K)); pr=[k for k in range(K) if 0<yb[:,k].sum()<len(yb)]
        if pr: aucs.append(roc_auc_score(yb[:,pr],P[:,pr],average="macro"))
    return float(np.mean(aucs)) if aucs else float("nan")
pca=PCA(50,random_state=0).fit_transform(StandardScaler().fit_transform(Xh)); spca=pca[idx].mean(1)
reps={"nbhd-composition":compo,"PCA-expression":pca,"spatial-smoothed-PCA":spca}
for f in sorted(glob.glob("expand_results/spatial_emb/*.npz")):
    nm=os.path.basename(f)[:-4]; d=np.load(f,allow_pickle=True); Z=d["X"]
    if Z.shape[0]==len(y): reps[nm+" (FM)"]=np.asarray(Z,np.float32)
    else: print("WARN shape mismatch",nm,Z.shape)
rows=[dict(method="TRUE niche label",knn_niche_AUROC=None,compo_R2=round(R2(compo,y),4),expr_R2=round(R2(Xh,y),4),contig=round(contig(y),4))]
for nm,Z in reps.items():
    part=KMeans(K,n_init=10,random_state=0).fit_predict(StandardScaler().fit_transform(Z))
    rows.append(dict(method=nm,knn_niche_AUROC=round(knn_probe(Z),3),compo_R2=round(R2(compo,part),4),expr_R2=round(R2(Xh,part),4),contig=round(contig(part),4)))
json.dump(rows,open("expand_results/spatial_knn_probe.json","w"),indent=2)
print(f"{'method':30s} {'kNN-niche-AUROC':>15} {'compo_R2':>9} {'expr_R2':>8} {'contig':>7}")
for r in rows: print(f"{r['method'][:30]:30s} {str(r['knn_niche_AUROC']):>15} {r['compo_R2']:9.4f} {r['expr_R2']:8.4f} {r['contig']:7.4f}")
