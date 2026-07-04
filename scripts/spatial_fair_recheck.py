
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

# Cluster-H fair re-check (reference-free). Tests the circularity: the niche labels are
# neighborhood-composition-defined anatomical zones, so compo-R2 (variance of local cell-type
# composition explained) vs expr-R2 (single-cell expression variance) reveals the niche label is a
# COMPOSITION cut, not an expression cut -> the nbhd-composition baseline winning ARI is near-tautological.
import anndata as ad, numpy as np, json, warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
rng=np.random.RandomState(20260623)
a=ad.read_h5ad("raw_pulls/spatial/nicheid/lymph.h5ad")
xy=np.asarray(a.obsm["spatial"]); niche=a.obs["niche"].astype(str).values; ct=a.obs["cell_type"].astype(str).values
import scipy.sparse as sp
X=a.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4)
var=Xln.var(0); hv=np.argsort(-var)[:2000]; Xh=Xln[:,hv]
# neighborhood cell-type composition (k=20 spatial neighbors)
cts=np.unique(ct); ct2i={c:i for i,c in enumerate(cts)}
nn=NearestNeighbors(n_neighbors=20).fit(xy); _,idx=nn.kneighbors(xy)
onehot=np.zeros((len(ct),len(cts)),np.float32); onehot[np.arange(len(ct)),[ct2i[c] for c in ct]]=1
compo=onehot[idx].mean(1)  # n x n_celltypes neighborhood composition
def R2(F, part):
    tot=((F-F.mean(0))**2).sum(); wit=0.0
    for c in np.unique(part):
        m=part==c
        if m.sum(): wit+=((F[m]-F[m].mean(0))**2).sum()
    return float(1-wit/tot)
def contig(part):  # fraction of spatial neighbors sharing the cell's label
    return float((part[idx[:,1:]]==part[:,None]).mean())
K=len(np.unique(niche)); nlab=np.array([np.where(np.unique(niche)==v)[0][0] for v in niche])
# representations -> KMeans(K) partitions
pca=PCA(50,random_state=0).fit_transform(StandardScaler().fit_transform(Xh))
# spatial-smoothed PCA (mean over neighbors)
spca=pca[idx].mean(1)
reps={"nbhd-composition":compo, "PCA-expression":pca, "spatial-smoothed-PCA":spca}
rows=[]
rows.append(dict(method="TRUE niche label", compo_R2=round(R2(compo,nlab),4), expr_R2=round(R2(Xh,nlab),4), contiguity=round(contig(nlab),4)))
for name,Z in reps.items():
    part=KMeans(K,n_init=10,random_state=0).fit_predict(Z)
    rows.append(dict(method=name+" KMeans", compo_R2=round(R2(compo,part),4), expr_R2=round(R2(Xh,part),4), contiguity=round(contig(part),4)))
json.dump(rows, open("expand_results/spatial_fair_reffree.json","w"), indent=2)
print(f"{'method':28s} {'compo_R2':>9} {'expr_R2':>8} {'contig':>7}")
for r in rows: print(f"{r['method'][:28]:28s} {r['compo_R2']:9.4f} {r['expr_R2']:8.4f} {r['contiguity']:7.4f}")
print("\nINTERPRETATION: if TRUE niche label has HIGH compo_R2 but LOW expr_R2, the niches are composition-defined")
print("-> nbhd-composition baseline reconstructs the definition (circular); expression FMs aren't expected to.")
