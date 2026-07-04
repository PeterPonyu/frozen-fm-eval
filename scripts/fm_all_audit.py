# Generic multi-FM audit: auto-discovers EVERY FM embedding in fm_emb/ by prefix and scores it
# (kNN fair probe + reference-free expr_R2) against PCA/HVG baselines on the 20 raw atlases.
# Add an FM = drop its {prefix}_{suffix}.npz files + re-run; this picks it up automatically.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score
rng=np.random.RandomState(20260623); EMB="expand_results/fm_emb"
FMNAME={"gf":"Geneformer-V2-104M","gf316":"Geneformer-V2-316M","scgpt":"scGPT","scf":"scFoundation","cellplm":"CellPLM","uce":"UCE"}
LR={os.path.basename(f)[:-5]:f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
ATL=[("GSE130148_lung",".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident"),
     ("GSE165784_retina",".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch"),
     ("lung24k",".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch")]
for nm,f in sorted(LR.items()): ATL.append(("lr_"+nm,f,"cell_type","batch"))
# Data-driven gene-naming classification (replaces the earlier hand-curated MOUSE set, which
# mislabelled 3 mouse-symbol atlases -- lung24k, lr_breast_hm, lr_tcell_cancer -- as "human").
# Human FMs (scGPT/Geneformer/scFoundation) tokenize human symbols; an atlas counts as a valid
# human test only if its var_names ARE human symbols. Mouse-symbol => cross-species token-map
# failure (grouped with mouse); ENSG-only => gene-ID-format incompatible with the symbol vocab.
HSYM=set(pd.read_csv("expand_results/scf_gene_index.tsv",sep="\t")["gene_name"])
def naming_of(A):
    vn=[str(x) for x in A.var_names]
    if len(set(vn)&HSYM)>=2000: return "human-symbol"
    if sum(1 for x in vn[:300] if x.startswith("ENSG"))>50: return "human-ensembl"
    return "mouse"
SPECIES={"human-symbol":"human","mouse":"mouse","human-ensembl":"human-ensembl"}
def R2(F,part):
    t=((F-F.mean(0))**2).sum(); w=sum(((F[part==c]-F[part==c].mean(0))**2).sum() for c in np.unique(part)); return float(1-w/t)
def knn_auroc(Z,y,tr,te,NC):
    Zs=StandardScaler().fit(Z[tr]).transform(Z)
    c=KNeighborsClassifier(n_neighbors=15,weights="distance").fit(Zs[tr],y[tr])
    P=np.zeros((te.sum(),NC)); P[:,c.classes_]=c.predict_proba(Zs[te])
    yb=label_binarize(y[te],classes=range(NC)); pr=[k for k in range(NC) if 0<yb[:,k].sum()<len(yb)]
    return float(roc_auc_score(yb[:,pr],P[:,pr],average="macro")) if pr else float("nan")
out=[]
for suffix,f,ct,bt in ATL:
    try: A=ad.read_h5ad(f)
    except Exception as e: print("skip",suffix,str(e)[:40]); continue
    y0=A.obs[ct].astype(str).values; b=A.obs[bt].astype(str).values
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4); var=Xln.var(0); hv=np.argsort(-var)[:2000]; Xh=Xln[:,hv]
    cnt=collections.Counter(y0); keep=np.array([cnt[v]>=10 for v in y0])
    cls=np.unique(y0[keep]); 
    if len(cls)<3: continue
    pca=PCA(min(50,Xh.shape[1]-1),random_state=0).fit_transform(StandardScaler().fit_transform(Xh))
    reps={"PCA":pca,"HVG":Xh}
    for ef in glob.glob(f"{EMB}/*_{suffix}.npz"):
        pref=os.path.basename(ef)[:-4]; pref=pref[:len(pref)-len(suffix)-1]  # strip _suffix
        d=np.load(ef,allow_pickle=True); Z=d["X"]
        if Z.shape[0]==A.n_obs: reps[FMNAME.get(pref,pref)]=np.asarray(Z,np.float32)
    yi=np.array([np.where(cls==v)[0][0] if v in cls else -1 for v in y0]); m=yi>=0
    Xh,pca,b2=Xh[m],pca[m],b[m]; yi=yi[m]; NC=len(cls)
    for k in list(reps): reps[k]=reps[k][m]
    bv,bc=np.unique(b2,return_counts=True); tb=bv[np.argmax(bc)]; te=b2==tb; tr=~te
    if te.sum()<100 or tr.sum()<300 or len(np.unique(yi[tr]))<3: continue
    nmcls=naming_of(A); row={"atlas":suffix,"naming":nmcls,"species":SPECIES[nmcls],"NC":int(NC),"reps":{}}
    for rn,Z in reps.items():
        try:
            km=KMeans(NC,n_init=5,random_state=0).fit_predict(StandardScaler().fit_transform(Z))
            row["reps"][rn]=dict(knn_auroc=round(knn_auroc(Z,yi,tr,te,NC),4), expr_R2=round(R2(Xh,km),4))
        except Exception as e: row["reps"][rn]=dict(error=str(e)[:60])
    out.append(row); print(suffix, {k:row['reps'][k].get('knn_auroc') for k in row['reps']}, flush=True)
json.dump(out, open("expand_results/fm_all_audit.json","w"), indent=1)
print("DONE", len(out), "atlases")
