
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

# FM-vs-baseline reliability on the 3 raw+labeled atlases.
# Representations: Geneformer-V2 (self-embedded), scGPT (cached, where available), PCA50, HVG-logreg, kNN, centroid.
# Same cross-batch (largest held-out batch) protocol + ECE + conformal coverage gap + abstention.
import anndata as ad, numpy as np, json, os, warnings, scipy.sparse as sp, collections
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score, accuracy_score
from scipy.special import softmax
rng=np.random.RandomState(20260623)
def ece(P,pred,y,B=15):
    conf=P.max(1); acc=(pred==y).astype(float); bins=np.linspace(0,1,B+1); e=0.
    for i in range(B):
        m=(conf>bins[i])&(conf<=bins[i+1])
        if m.sum(): e+=m.mean()*abs(acc[m].mean()-conf[m].mean())
    return float(e)
def conf_cov(calP,caly,teP,tey,alpha=0.1):
    s=1-calP[np.arange(len(caly)),caly]; n=len(s)
    q=np.quantile(s,min(1.,np.ceil((n+1)*(1-alpha))/n),method="higher")
    sets=(1-teP)<=q; return float(sets[np.arange(len(tey)),tey].mean())
def probe(Xtr,ytr,Xte,nclasses):
    c=LogisticRegression(max_iter=300,multi_class="multinomial"); c.fit(Xtr,ytr)
    Pp=c.predict_proba(Xte); P=np.zeros((Xte.shape[0],nclasses),dtype=float)
    P[:,c.classes_]=Pp  # scatter into full class space (train split may miss classes)
    return P, c
def audit_rep(Z,y_raw,b,label):
    classes=np.unique(y_raw); cls={c:i for i,c in enumerate(classes)}; y=np.array([cls[v] for v in y_raw])
    cnt=collections.Counter(y); keep=np.array([cnt[v]>=10 for v in y]); Z,y,b=Z[keep],y[keep],b[keep]
    classes=np.unique(y); y=np.array([np.where(classes==v)[0][0] for v in y])
    bvals,bc=np.unique(b,return_counts=True); order=bvals[np.argsort(-bc)]; tb=None
    for bb in order:
        te=b==bb
        if te.sum()>=200 and (~te).sum()>=500 and len(np.unique(y[~te]))>=3: tb=bb;break
    if tb is None: return None
    te=b==tb; tr=~te
    Zs=StandardScaler().fit(Z[tr]).transform(Z)
    yb=label_binarize(y[te],classes=range(len(classes)))
    NC=len(classes)
    P,_=probe(Zs[tr],y[tr],Zs[te],NC); pred=P.argmax(1)
    # cal split for conformal
    tri=np.where(tr)[0]; rng.shuffle(tri); ci=tri[:int(.4*len(tri))]; fi=tri[int(.4*len(tri)):]
    Pcal,_=probe(Zs[fi],y[fi],Zs[ci],NC)
    nte=te.sum(); perm=rng.permutation(len(y)); rte=np.zeros(len(y),bool); rte[perm[:nte]]=True
    Prnd,_=probe(Zs[~rte],y[~rte],Zs[rte],NC)
    cov_xb=conf_cov(Pcal,y[ci],P,y[te]); cov_rnd=conf_cov(Pcal,y[ci],Prnd,y[rte])
    present=[c for c in range(NC) if 0<yb[:,c].sum()<len(yb)]
    try: auroc=float(roc_auc_score(yb[:,present],P[:,present],average="macro")) if present else float("nan")
    except Exception: auroc=float("nan")
    conf=P.max(1); thr=np.quantile(conf,0.2); kp=conf>=thr
    return dict(rep=label,n=int(len(y)),n_ct=int(len(classes)),test_batch=str(tb),
        xb_auroc=auroc,ece=ece(P,pred,y[te]),cov_xb=cov_xb,cov_rnd=cov_rnd,cov_gap=float(cov_rnd-cov_xb),
        acc_full=float(accuracy_score(y[te],pred)),acc_at80=float(accuracy_score(y[te][kp],pred[kp])))
ATL=[("GSE130148_lung",".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad","celltype","orig.ident",
      "expand_results/fm_emb/gf_GSE130148_lung.npz","expand_results/fm_emb/scgpt_GSE130148_lung.npz"),
     ("GSE165784_retina",".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad","cell_type","batch",
      "expand_results/fm_emb/gf_GSE165784_retina.npz","expand_results/fm_emb/scgpt_GSE165784_retina.npz"),
     ("lung24k",".../data/datasets/DevelopmentDatasets/lung.h5ad","louvain","batch",
      "expand_results/fm_emb/gf_lung24k.npz","expand_results/fm_emb/scgpt_lung24k.npz")]
import glob as _glob
for _f in sorted(_glob.glob("expand_results/labeled_raw/*.h5ad")):
    _nm=os.path.basename(_f)[:-5]
    ATL.append(("lr_"+_nm,_f,"cell_type","batch",f"expand_results/fm_emb/gf_lr_{_nm}.npz",f"expand_results/fm_emb/scgpt_lr_{_nm}.npz"))
print("total atlases in FM audit:",len(ATL),flush=True)
out=[]
for name,f,ct,bt,gfp,scp in ATL:
    A=ad.read_h5ad(f); y=A.obs[ct].astype(str).values; b=A.obs[bt].astype(str).values
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    tot=X.sum(1,keepdims=True); tot[tot==0]=1; X=np.log1p(X/tot*1e4)  # fair baseline normalization
    var=X.var(0); hv=np.argsort(-var)[:2000]; Xh=X[:,hv]
    pca=PCA(n_components=min(50,Xh.shape[1]-1),random_state=0).fit(StandardScaler().fit_transform(Xh))
    Xp=pca.transform(StandardScaler().fit_transform(Xh))
    reps=[("PCA50",Xp,y,b),("HVG2000",Xh,y,b)]
    if os.path.exists(gfp):
        g=np.load(gfp,allow_pickle=True); reps.append(("Geneformer-V2(FM)",g["X"],g["y"].astype(str),g["batch"].astype(str)))
    if scp and os.path.exists(scp):
        s=np.load(scp,allow_pickle=True); reps.append(("scGPT(FM)",s["X"],s["y"].astype(str),s["batch"].astype(str)))
    res=[audit_rep(Z,yy,bb,lab) for lab,Z,yy,bb in reps]; res=[r for r in res if r]
    out.append(dict(atlas=name,reps=res))
    print(f"=== {name} ===")
    for r in res: print(f"  {r['rep']:20s} auroc={r['xb_auroc']:.3f} ece={r['ece']:.3f} cov_gap={r['cov_gap']:.3f} acc@80-full={r['acc_at80']-r['acc_full']:+.3f}")
json.dump(out,open("expand_results/fm_vs_baseline_raw.json","w"),indent=1)
print("[fm-vs-baseline-done]")
