
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

# LEAN multi-atlas baseline reliability sweep (clusters J/K). 2 fits/method, RF capped.
import anndata as ad, numpy as np, json, os, warnings, scipy.sparse as sp, collections, time
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score, accuracy_score
from scipy.special import softmax
rng=np.random.RandomState(20260623)
man=json.load(open("expand_results/atlas_manifest.json"))
DIRS={os.path.basename(d):d for d in [
 ".../data/datasets/extra_preprocessed",".../data/datasets/CancerDatasets",
 ".../data/datasets/CancerDatasets2",".../data/datasets/DevelopmentDatasets",
 ".../data/datasets/DevelopmentDatasets2"]}
usable=[r for r in man if r.get("usable")]; NCELL=6000
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
def full(clf,Xte,NC):
    Pp=clf.predict_proba(Xte); P=np.zeros((Xte.shape[0],NC)); P[:,clf.classes_]=Pp; return P
def fit_predict(m,Xf,yf,Xs,NC):
    if m=="logreg": c=LogisticRegression(max_iter=300,multi_class="multinomial")
    elif m=="knn": c=KNeighborsClassifier(n_neighbors=15,weights="distance")
    elif m=="rf": c=RandomForestClassifier(n_estimators=60,n_jobs=4,random_state=0)
    elif m=="centroid":
        c=NearestCentroid(); c.fit(Xf,yf)
        def cp(X):
            d=np.stack([np.linalg.norm(X-cc,axis=1) for cc in c.centroids_],1); P=softmax(-d,1)
            FP=np.zeros((X.shape[0],NC)); FP[:,c.classes_]=P; return FP
        return [cp(x) for x in Xs], c
    c.fit(Xf,yf); return [full(c,x,NC) for x in Xs], c
METHODS=[("pca-logreg","pca","logreg"),("hvg-logreg","hvg","logreg"),
         ("knn","pca","knn"),("centroid","pca","centroid"),("rf","pca","rf")]
results=[]; t0=time.time()
for r in usable:
    f=os.path.join(DIRS[r["dir"]],r["file"])
    try: A=ad.read_h5ad(f)
    except Exception as e: print("SKIP",r["file"],str(e)[:40],flush=True); continue
    y0=A.obs[r["ct"]].astype(str).values; b=A.obs[r["batch"]].astype(str).values
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    cnt=collections.Counter(y0); keep=np.array([cnt[v]>=10 for v in y0]); X,y0,b=X[keep],y0[keep],b[keep]
    if len(X)>NCELL: idx=rng.choice(len(X),NCELL,replace=False); X,y0,b=X[idx],y0[idx],b[idx]
    cls=np.unique(y0); 
    if len(cls)<3: continue
    y=np.array([np.where(cls==v)[0][0] for v in y0]); NC=len(cls)
    bvals,bc=np.unique(b,return_counts=True); tb=None
    for bb in bvals[np.argsort(-bc)]:
        te=b==bb
        if te.sum()>=200 and (~te).sum()>=600 and len(np.unique(y[~te]))>=3: tb=bb;break
    if tb is None: continue
    te=b==tb; nt=np.where(~te)[0]; rng.shuffle(nt)
    a1=int(.5*len(nt)); a2=int(.75*len(nt)); fit_i=nt[:a1]; cal_i=nt[a1:a2]; rnd_i=nt[a2:]
    var=X.var(0); hv=np.argsort(-var)[:2000]; Xh=X[:,hv]
    ssh=StandardScaler().fit(Xh[fit_i]); Xhs=ssh.transform(Xh)
    ncomp=min(50,Xhs.shape[1]-1,len(fit_i)-1); pca=PCA(ncomp,random_state=0).fit(Xhs[fit_i]); Xp=pca.transform(Xhs)
    REP={"pca":Xp,"hvg":Xhs}; yb_te=label_binarize(y[te],classes=range(NC)); yb_rd=label_binarize(y[rnd_i],classes=range(NC))
    arow={"atlas":r["file"].replace("_prepped.h5ad","").replace(".h5ad",""),"n":int(len(y)),"n_ct":NC,
          "n_batch":int(len(bvals)),"test_batch":str(tb),"methods":{}}
    for mname,rep,algo in METHODS:
        Z=REP[rep]
        try:
            (Pcal,Pxb,Prnd),clf=fit_predict(algo,Z[fit_i],y[fit_i],[Z[cal_i],Z[te],Z[rnd_i]],NC)
            pred=Pxb.argmax(1)
            xb=float(roc_auc_score(yb_te,Pxb,average="macro",multi_class="ovr"))
            rd=float(roc_auc_score(yb_rd,Prnd,average="macro",multi_class="ovr"))
            d=dict(xb_auroc=xb,rnd_auroc=rd,leak_gap=rd-xb,ece_raw=ece(Pxb,pred,y[te]),
                   cov_rnd=cov(Pcal,y[cal_i],Prnd,y[rnd_i]),cov_xb=cov(Pcal,y[cal_i],Pxb,y[te]),
                   acc_full=float(accuracy_score(y[te],pred)))
            d["cov_gap"]=d["cov_rnd"]-d["cov_xb"]
            conf=Pxb.max(1); thr=np.quantile(conf,0.2); kp=conf>=thr
            d["acc_at80"]=float(accuracy_score(y[te][kp],pred[kp])) if kp.sum() else float("nan")
            if algo=="logreg":
                Ts=np.linspace(0.5,5,46); dfc=clf.decision_function(Z[cal_i]); best=(1,1e9)
                for T in Ts:
                    p=softmax(dfc/T,1); vi=[(i,list(clf.classes_).index(v)) for i,v in enumerate(y[cal_i]) if v in clf.classes_]
                    if not vi: continue
                    ri=np.array([x[0] for x in vi]); ci=np.array([x[1] for x in vi])
                    nll=-np.log(np.clip(p[ri,ci],1e-9,1)).mean()
                    if nll<best[1]: best=(T,nll)
                T=best[0]; dft=clf.decision_function(Z[te]); Pt=np.zeros((te.sum(),NC)); Pt[:,clf.classes_]=softmax(dft/T,1)
                d["ece_temp"]=ece(Pt,pred,y[te]); d["T"]=float(T)
            ys=y[fit_i].copy(); rng.shuffle(ys); (_,Ps,_) ,_=fit_predict(algo,Z[fit_i],ys,[Z[cal_i],Z[te],Z[cal_i][:1]],NC)
            d["shuffle_auroc"]=float(roc_auc_score(yb_te,Ps,average="macro",multi_class="ovr"))
            arow["methods"][mname]=d
        except Exception as e: arow["methods"][mname]=dict(error=str(e)[:90])
    results.append(arow)
    print(f"[{time.time()-t0:5.0f}s] {arow['atlas'][:26]:26s} ct={NC:2d} | "+" ".join(f"{m}:gap={arow['methods'][m].get('cov_gap',float('nan')):.3f}" for m in ['pca-logreg','rf']),flush=True)
json.dump(results,open("expand_results/multiatlas_baseline.json","w"),indent=1)
print("DONE",len(results),"atlases in",int(time.time()-t0),"s",flush=True)
