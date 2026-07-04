# Depth probe #2, COMPREHENSIVE extension: does the batch-shift -> coverage-collapse dose-response
# hold (a) when foundation-model representations are added, (b) for a second reliability metric, and
# (c) for a second, classifier-free shift measure? Addresses the "single metric / few weak methods"
# critique: 5 representations (2 classical + 3 FM families), 2 degradation metrics, 2 shift measures.
# Key test: FM points must lie ON the same curve as PCA -> FMs do not repair the exchangeability failure.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp, collections, pandas as pd
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr, pearsonr
EMB="expand_results/fm_emb"
FMNAME={"gf":"Geneformer-V2-104M","scgpt":"scGPT","scf":"scFoundation"}  # 3 FM families (one per arch)
FMSET=set(FMNAME.values())
HSYM=set(pd.read_csv("expand_results/scf_gene_index.tsv",sep="\t")["gene_name"])
def naming_of(A):
    vn=[str(x) for x in A.var_names]
    if len(set(vn)&HSYM)>=2000: return "human-symbol"
    if sum(1 for x in vn[:300] if x.startswith("ENSG"))>50: return "human-ensembl"
    return "mouse"
LR={os.path.basename(f)[:-5]:f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
ATL=[("GSE130148_lung",".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident"),
     ("GSE165784_retina",".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch"),
     ("lung24k",".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch")]
for nm,f in sorted(LR.items()): ATL.append(("lr_"+nm,f,"cell_type","batch"))

def ece(P,pred,y,B=15):
    conf=P.max(1); acc=(pred==y).astype(float); bins=np.linspace(0,1,B+1); e=0.
    for i in range(B):
        m=(conf>bins[i])&(conf<=bins[i+1])
        if m.sum(): e+=m.mean()*abs(acc[m].mean()-conf[m].mean())
    return float(e)
def cov(calP,caly,teP,tey,alpha=0.1):
    s=1-calP[np.arange(len(caly)),caly]; n=len(s)
    q=np.quantile(s,min(1.,np.ceil((n+1)*(1-alpha))/n),method="higher")
    sets=(1-teP)<=q; return float(sets[np.arange(len(tey)),tey].mean())
def fullP(clf,X,NC):
    Pp=clf.predict_proba(X); P=np.zeros((X.shape[0],NC)); P[:,clf.classes_]=Pp; return P
def disc(P,y,NC):  # NaN-safe present-class macro-AUROC (discrimination), unlike multi_class=ovr
    yb=label_binarize(y,classes=range(NC)); pr=[k for k in range(NC) if 0<yb[:,k].sum()<len(yb)]
    return float(roc_auc_score(yb[:,pr],P[:,pr],average="macro")) if pr else float("nan")
def shift_auroc(Z,is_test):
    if is_test.sum()<20 or (~is_test).sum()<20: return float("nan")
    cv=StratifiedKFold(5,shuffle=True,random_state=0)
    P=cross_val_predict(LogisticRegression(max_iter=200),StandardScaler().fit_transform(Z),
                        is_test.astype(int),cv=cv,method="predict_proba")[:,1]
    return float(roc_auc_score(is_test.astype(int),P))
def shift_cdist(Z,is_test):
    # parameter-free: standardized centroid distance between test-batch and training cells
    Zs=StandardScaler().fit_transform(Z); a=Zs[is_test].mean(0); b0=Zs[~is_test].mean(0)
    sd=Zs.std(0).mean()+1e-9; return float(np.linalg.norm(a-b0)/(np.sqrt(Zs.shape[1])*sd))

rng=np.random.RandomState(20260625); rows=[]
for suffix,f,ct,bt in ATL:
    try: A=ad.read_h5ad(f)
    except Exception as e: print("skip",suffix,str(e)[:40]); continue
    nm=naming_of(A); y0=A.obs[ct].astype(str).values; b=A.obs[bt].astype(str).values
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; Xln=np.log1p(X/tot*1e4)
    var=Xln.var(0); hv=np.argsort(-var)[:2000]; Xh=Xln[:,hv]
    cnt=collections.Counter(y0); keep=np.array([cnt[v]>=10 for v in y0]); cls=np.unique(y0[keep])
    if len(cls)<3: continue
    pca=PCA(min(50,Xh.shape[1]-1),random_state=0).fit_transform(StandardScaler().fit_transform(Xh))
    reps={"PCA":pca,"HVG":Xh}
    for ef in glob.glob(f"{EMB}/*_{suffix}.npz"):
        pref=os.path.basename(ef)[:-4]; pref=pref[:len(pref)-len(suffix)-1]
        if pref in FMNAME and nm=="human-symbol":
            Z=np.load(ef,allow_pickle=True)["X"]
            if Z.shape[0]==A.n_obs: reps[FMNAME[pref]]=np.asarray(Z,np.float32)
    yi=np.array([np.where(cls==v)[0][0] if v in cls else -1 for v in y0]); mk=yi>=0
    yi=yi[mk]; b2=b[mk]; NC=len(cls)
    for k in list(reps): reps[k]=reps[k][mk]
    # pick the largest batch with enough test/train cells and >=3 train classes
    bv,bc=np.unique(b2,return_counts=True); tb=None
    for bb in bv[np.argsort(-bc)]:
        te=b2==bb
        if te.sum()>=150 and (~te).sum()>=450 and len(np.unique(yi[~te]))>=3: tb=bb;break
    if tb is None: continue
    te=b2==tb; nt=np.where(~te)[0]; rng.shuffle(nt)
    a1=int(.5*len(nt)); a2=int(.75*len(nt)); fit_i=nt[:a1]; cal_i=nt[a1:a2]; rnd_i=nt[a2:]
    if len(cal_i)<30 or len(rnd_i)<30: continue
    for rn,Z in reps.items():
        try:
            ss=StandardScaler().fit(Z[fit_i]); Zs=ss.transform(Z)
            clf=LogisticRegression(max_iter=300).fit(Zs[fit_i],yi[fit_i])
            Pcal=fullP(clf,Zs[cal_i],NC); Pxb=fullP(clf,Zs[te],NC); Prnd=fullP(clf,Zs[rnd_i],NC)
            cov_rnd=cov(Pcal,yi[cal_i],Prnd,yi[rnd_i]); cov_xb=cov(Pcal,yi[cal_i],Pxb,yi[te])
            ece_rnd=ece(Prnd,Prnd.argmax(1),yi[rnd_i]); ece_xb=ece(Pxb,Pxb.argmax(1),yi[te])
            disc_rnd=disc(Prnd,yi[rnd_i],NC); disc_xb=disc(Pxb,yi[te],NC)  # discrimination (ranking) in/cross-batch
            rows.append(dict(atlas=suffix,naming=nm,rep=rn,is_fm=rn in FMSET,
                shift_auroc=round(shift_auroc(Z,te),4), shift_cdist=round(shift_cdist(Z,te),4),
                cov_gap=round(cov_rnd-cov_xb,4), ece_gap=round(ece_xb-ece_rnd,4),
                cov_xb=round(cov_xb,4), ece_xb=round(ece_xb,4),
                disc_xb=round(disc_xb,4), disc_drop=round(disc_rnd-disc_xb,4)))
        except Exception as e: print("  err",suffix,rn,str(e)[:50])
    print(f"{suffix[:22]:22s} nm={nm[:6]} reps={list(reps)}",flush=True)

def corr(rs,xk,yk):
    x=np.array([r[xk] for r in rs]); y=np.array([r[yk] for r in rs]); ok=np.isfinite(x)&np.isfinite(y)
    x,y=x[ok],y[ok]
    if len(x)<4: return None
    sr,sp_=spearmanr(x,y); pr,pp=pearsonr(x,y)
    return dict(n=int(len(x)),spearman=round(float(sr),3),spearman_p=float(sp_),pearson=round(float(pr),3),pearson_p=float(pp))
fm=[r for r in rows if r["is_fm"]]; cl=[r for r in rows if not r["is_fm"]]
stats={}
for yk in("cov_gap","ece_gap"):
    for xk in("shift_auroc","shift_cdist"):
        stats[f"ALL:{yk}~{xk}"]=corr(rows,xk,yk)
    stats[f"FM-only:{yk}~shift_auroc"]=corr(fm,"shift_auroc",yk)
    stats[f"classical-only:{yk}~shift_auroc"]=corr(cl,"shift_auroc",yk)
json.dump({"rows":rows,"stats":stats},open("expand_results/batch_shift_fm_probe.json","w"),indent=1)
print("\n=== COMPREHENSIVE BATCH-SHIFT (FM + classical) ===")
print(f"points: {len(rows)}  (FM={len(fm)}, classical={len(cl)})  reps:",sorted(set(r['rep'] for r in rows)))
for k,v in stats.items(): print(f"  {k:34s} {v}")
# do FM points lie on the same line as classical? compare mean cov_gap in low/high shift bins per group
for grp,rs in (("FM",fm),("classical",cl)):
    lo=[r["cov_gap"] for r in rs if np.isfinite(r["shift_auroc"]) and r["shift_auroc"]<0.9]
    hi=[r["cov_gap"] for r in rs if np.isfinite(r["shift_auroc"]) and r["shift_auroc"]>=0.9]
    print(f"  {grp:9s} cov_gap: weak-shift mean={np.mean(lo) if lo else float('nan'):+.3f} (n={len(lo)}) | strong mean={np.mean(hi) if hi else float('nan'):+.3f} (n={len(hi)})")
print("DONE")
