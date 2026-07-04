# Push: put cluster H (spatial niche-ID) on a dose-response ruler. The niche task is composition-defined
# (R2_compo 0.234 >> R2_expr 0.007), so a per-cell EXPRESSION model (the FMs) is handicapped by exactly
# one axis: spatial context. We dose that axis -- progressively spatial-smooth a per-cell representation
# (PCA, and one FM) over k spatial neighbours -- and measure niche prediction vs k, with the 3 raw spatial
# FMs and the nbhd-composition baseline as references. Two supervised metrics (spatial-blocked kNN AUROC
# and macro-F1). Prediction: niche quality rises monotonically with the spatial dose and overtakes the FMs.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score, f1_score
from scipy.stats import spearmanr
a=ad.read_h5ad("raw_pulls/spatial/nicheid/lymph.h5ad")
xy=np.asarray(a.obsm["spatial"]); niche=a.obs["niche"].astype(str).values; ct=a.obs["cell_type"].astype(str).values
X=a.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4); var=Xln.var(0); hv=np.argsort(-var)[:2000]; Xh=Xln[:,hv]
nlev=np.unique(niche); y=np.array([np.where(nlev==v)[0][0] for v in niche]); K=len(nlev)
# spatial-blocked 4-fold (identical to spatial_knn_probe.py)
gx=np.digitize(xy[:,0],np.quantile(xy[:,0],[.25,.5,.75])); gy=np.digitize(xy[:,1],np.quantile(xy[:,1],[.25,.5,.75]))
folds=(gx//2)*2+(gy//2)
def probe(Z):
    Zs=StandardScaler().fit_transform(Z); aucs=[]; f1s=[]
    for f in range(4):
        te=folds==f; tr=~te
        if len(np.unique(y[tr]))<K or te.sum()<50: continue
        c=KNeighborsClassifier(n_neighbors=30,weights="distance").fit(Zs[tr],y[tr])
        P=np.zeros((te.sum(),K)); P[:,c.classes_]=c.predict_proba(Zs[te]); pred=P.argmax(1)
        yb=label_binarize(y[te],classes=range(K)); pr=[k for k in range(K) if 0<yb[:,k].sum()<len(yb)]
        if pr: aucs.append(roc_auc_score(yb[:,pr],P[:,pr],average="macro"))
        f1s.append(f1_score(y[te],pred,average="macro"))
    return (float(np.mean(aucs)) if aucs else float("nan"),
            float(np.mean(f1s)) if f1s else float("nan"))
# max-k spatial neighbour index (include self)
KMAX=200; nn=NearestNeighbors(n_neighbors=KMAX+1).fit(xy); _,idx=nn.kneighbors(xy)
def smooth(Z,k):
    if k==0: return Z
    return Z[idx[:,:k+1]].mean(1)   # self + k nearest spatial neighbours
pca=PCA(50,random_state=0).fit_transform(StandardScaler().fit_transform(Xh))
# pick scGPT-spatial as the FM to also dose (a per-cell spatial FM); keep all 3 as raw references
FMZ={}
for f in glob.glob("expand_results/spatial_emb/*.npz"):
    nm=os.path.basename(f)[:-4]; Z=np.load(f,allow_pickle=True)["X"]
    if Z.shape[0]==len(y): FMZ[nm]=np.asarray(Z,np.float32)
# cell-type composition baseline (the niche-definition recipe), as a reference horizontal
cts=np.unique(ct); ct2i={c:i for i,c in enumerate(cts)}
onehot=np.zeros((len(ct),len(cts)),np.float32); onehot[np.arange(len(ct)),[ct2i[c] for c in ct]]=1
compo=onehot[idx[:,:21]].mean(1)  # 20-NN composition (matches original baseline)

KS=[0,3,6,12,25,50,100,200]
curves={"spatial-smoothed PCA":pca}
if "scgpt_spatial" in FMZ: curves["spatial-smoothed scGPT-spatial (FM)"]=FMZ["scgpt_spatial"]
dose=[]
for cname,Z in curves.items():
    for k in KS:
        au,f1=probe(smooth(Z,k))
        dose.append(dict(curve=cname,k=int(k),auroc=round(au,4),f1=round(f1,4)))
        print(f"{cname[:34]:34s} k={k:3d}  AUROC={au:.3f}  F1={f1:.3f}",flush=True)
# references (raw, no dosing)
refs=[]
for nm,Z in FMZ.items():
    au,f1=probe(Z); refs.append(dict(name=nm+" (raw FM)",auroc=round(au,4),f1=round(f1,4)))
cau,cf1=probe(compo); refs.append(dict(name="nbhd-composition baseline",auroc=round(cau,4),f1=round(cf1,4)))
# monotonicity / dose stats per curve (Spearman k vs metric)
stats={}
for cname in curves:
    pts=[d for d in dose if d["curve"]==cname]
    ks=[d["k"] for d in pts]
    stats[cname]={m:dict(spearman=round(float(spearmanr(ks,[d[m] for d in pts])[0]),3),
                         k0=pts[0][m], kmax=pts[-1][m], gain=round(pts[-1][m]-pts[0][m],4)) for m in("auroc","f1")}
json.dump({"dose":dose,"refs":refs,"stats":stats,"K_niches":int(K)},
          open("expand_results/spatial_dose_response.json","w"),indent=1)
print("\n=== references (raw) ===")
for r in refs: print(f"  {r['name'][:34]:34s} AUROC={r['auroc']:.3f}  F1={r['f1']:.3f}")
print("\n=== dose monotonicity (Spearman k vs metric) ===")
for c,s in stats.items():
    print(f"  {c}")
    for m in("auroc","f1"): print(f"     {m}: rho={s[m]['spearman']:+.2f}  k0={s[m]['k0']:.3f} -> kmax={s[m]['kmax']:.3f}  (+{s[m]['gain']:.3f})")
print("DONE")
