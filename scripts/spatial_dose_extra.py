# Probe #3 completeness: (A) DIFFERENTIAL -- does spatial smoothing help a purely per-cell representation
# (PCA) MORE than it helps the already-spatial FMs (Novae/Nicheformer/scGPT-spatial)? If so, the FMs already
# carry partial spatial context, and the residual gap is the remaining aggregation they lack. (B) MECHANISM --
# as the spatial dose k rises, does the smoothed-PCA partition's composition-alignment (R2_compo) rise toward
# the niche-definition structure? That shows the dose literally adds composition signal, not just smooths noise.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score, f1_score
from scipy.stats import spearmanr
a=ad.read_h5ad("raw_pulls/spatial/nicheid/lymph.h5ad")
xy=np.asarray(a.obsm["spatial"]); niche=a.obs["niche"].astype(str).values; ct=a.obs["cell_type"].astype(str).values
X=a.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4); var=Xln.var(0); hv=np.argsort(-var)[:2000]
nlev=np.unique(niche); y=np.array([np.where(nlev==v)[0][0] for v in niche]); K=len(nlev)
gx=np.digitize(xy[:,0],np.quantile(xy[:,0],[.25,.5,.75])); gy=np.digitize(xy[:,1],np.quantile(xy[:,1],[.25,.5,.75]))
folds=(gx//2)*2+(gy//2)
KMAX=200; nn=NearestNeighbors(n_neighbors=KMAX+1).fit(xy); _,idx=nn.kneighbors(xy)
def smooth(Z,k): return Z if k==0 else Z[idx[:,:k+1]].mean(1)
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
    return (float(np.mean(aucs)) if aucs else float("nan"), float(np.mean(f1s)) if f1s else float("nan"))
def R2(F,part):
    t=((F-F.mean(0))**2).sum(); w=sum(((F[part==c]-F[part==c].mean(0))**2).sum() for c in np.unique(part)); return float(1-w/t)
pca=PCA(50,random_state=0).fit_transform(StandardScaler().fit_transform(Xln[:,hv]))
cts=np.unique(ct); ct2i={c:i for i,c in enumerate(cts)}
onehot=np.zeros((len(ct),len(cts)),np.float32); onehot[np.arange(len(ct)),[ct2i[c] for c in ct]]=1
compo=onehot[idx[:,:21]].mean(1)  # neighborhood cell-type composition (the niche-definition substrate)
FMZ={}
for f in glob.glob("expand_results/spatial_emb/*.npz"):
    nm=os.path.basename(f)[:-4]; Z=np.load(f,allow_pickle=True)["X"]
    if Z.shape[0]==len(y): FMZ[nm]=np.asarray(Z,np.float32)

# (A) DIFFERENTIAL: F1 gain from smoothing, per representation
reps={"PCA (per-cell)":pca}
for k in ("scgpt_spatial","novae","nicheformer"):
    if k in FMZ: reps[k+" (FM)"]=FMZ[k]
KS=[0,3,6,12,25,50,100,200]
diff=[]
for rn,Z in reps.items():
    f1s={k:probe(smooth(Z,k))[1] for k in KS}
    best_k=max(f1s,key=f1s.get);
    diff.append(dict(rep=rn, f1_percell=round(f1s[0],4), f1_best=round(f1s[best_k],4), best_k=best_k,
                     gain=round(f1s[best_k]-f1s[0],4)))
    print(f"{rn:30s} f1(k=0)={f1s[0]:.3f}  best f1={f1s[best_k]:.3f} @k={best_k}  gain=+{f1s[best_k]-f1s[0]:.3f}",flush=True)
# (B) MECHANISM: R2_compo of the smoothed-PCA partition vs k (does the dose add composition structure?)
mech=[]
true_compo=R2(compo,y)  # composition variance explained by the TRUE niche labels (the target)
for k in KS:
    part=KMeans(K,n_init=5,random_state=0).fit_predict(StandardScaler().fit_transform(smooth(pca,k)))
    mech.append(dict(k=int(k), compo_R2=round(R2(compo,part),4)))
    print(f"  smoothed-PCA k={k:3d}  partition R2_compo={mech[-1]['compo_R2']:.3f}  (true-niche R2_compo={true_compo:.3f})",flush=True)
rho_mech=spearmanr([m["k"] for m in mech],[m["compo_R2"] for m in mech])[0]
json.dump({"differential":diff,"mechanism":mech,"true_niche_compo_R2":round(true_compo,4),
           "mechanism_spearman_k_vs_compoR2":round(float(rho_mech),3)},
          open("expand_results/spatial_dose_extra.json","w"),indent=1)
print("\n=== DIFFERENTIAL: per-cell PCA should gain MORE than already-spatial FMs ===")
pg=[d["gain"] for d in diff if d["rep"].startswith("PCA")][0]
fmg=[d["gain"] for d in diff if "FM" in d["rep"]]
print(f"PCA gain=+{pg:.3f} vs FM gains {['+%.3f'%g for g in fmg]} (mean +{np.mean(fmg):.3f})")
print(f"=== MECHANISM: smoothed-PCA partition R2_compo rises with dose, Spearman={rho_mech:+.3f} (toward true-niche {true_compo:.3f}) ===")
print("DONE")
