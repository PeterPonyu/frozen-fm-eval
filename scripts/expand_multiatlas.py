# Multi-atlas baseline reliability expansion (clusters J/K, baseline side).
# No FM weights, no downloads. 24 local atlases x expanded simple-baseline panel.
# Per atlas, one cross-batch holdout (largest batch) + matched random split.
import anndata as ad, numpy as np, json, os, warnings, scipy.sparse as sp
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.calibration import CalibratedClassifierCV
rng = np.random.RandomState(20260623)
man = json.load(open("expand_results/atlas_manifest.json"))
DIRS = {os.path.basename(d): d for d in [
 ".../data/datasets/extra_preprocessed",
 ".../data/datasets/CancerDatasets",".../data/datasets/CancerDatasets2",
 ".../data/datasets/DevelopmentDatasets",".../data/datasets/DevelopmentDatasets2"]}
usable = [r for r in man if r.get("usable")]
NCELL=8000  # subsample cap for speed
def ece(p, yhat, y, B=15):
    conf=p.max(1); pred=yhat; acc=(pred==y).astype(float)
    bins=np.linspace(0,1,B+1); e=0.0
    for i in range(B):
        m=(conf>bins[i])&(conf<=bins[i+1])
        if m.sum(): e+=m.mean()*abs(acc[m].mean()-conf[m].mean())
    return float(e)
def temp_scale(logits, y):
    # 1D grid search T on NLL
    from scipy.special import softmax
    Ts=np.linspace(0.5,5,46); best=(1.0,1e9)
    for T in Ts:
        p=softmax(logits/T,axis=1); nll=-np.log(np.clip(p[np.arange(len(y)),y],1e-9,1)).mean()
        if nll<best[1]: best=(T,nll)
    return best[0]
def macro_auroc(yb, score):
    try: return float(roc_auc_score(yb, score, average="macro", multi_class="ovr"))
    except Exception: return np.nan
def conformal_cov(cal_p, cal_y, test_p, test_y, alpha=0.1):
    s=1-cal_p[np.arange(len(cal_y)),cal_y]; n=len(s)
    q=np.quantile(s, min(1.0,np.ceil((n+1)*(1-alpha))/n), method="higher")
    sets=(1-test_p)<=q
    cov=sets[np.arange(len(test_y)),test_y].mean()
    return float(cov), float(sets.sum(1).mean())
def reps_and_probe(Xtr,ytr,Xte, method):
    if method=="pca-logreg":
        clf=LogisticRegression(max_iter=300,C=1.0,multi_class="multinomial")
        clf.fit(Xtr,ytr); P=clf.predict_proba(Xte); L=clf.decision_function(Xte); return P,L,clf
    if method=="knn":
        clf=KNeighborsClassifier(n_neighbors=15,weights="distance"); clf.fit(Xtr,ytr)
        P=clf.predict_proba(Xte); return P,np.log(np.clip(P,1e-6,1)),clf
    if method=="centroid":
        clf=NearestCentroid(metric="euclidean"); clf.fit(Xtr,ytr)
        # pseudo-proba via softmax of neg distance to centroids
        d=np.stack([np.linalg.norm(Xte-c,axis=1) for c in clf.centroids_],1)
        from scipy.special import softmax; P=softmax(-d,axis=1); return P,-d,clf
    if method=="svm-linear":
        base=LinearSVC(C=0.5,max_iter=3000)
        clf=CalibratedClassifierCV(base,cv=3); clf.fit(Xtr,ytr)
        P=clf.predict_proba(Xte); return P,np.log(np.clip(P,1e-6,1)),clf
    if method=="rf":
        clf=RandomForestClassifier(n_estimators=120,n_jobs=-1,random_state=0); clf.fit(Xtr,ytr)
        P=clf.predict_proba(Xte); return P,np.log(np.clip(P,1e-6,1)),clf
METHODS=["pca-logreg","hvg-logreg","knn","centroid","svm-linear","rf"]
results=[]
for r in usable:
    f=os.path.join(DIRS[r["dir"]], r["file"])
    try:
        A=ad.read_h5ad(f)
    except Exception as e:
        print("SKIP",r["file"],str(e)[:50]); continue
    y_raw=A.obs[r["ct"]].astype(str).values; b=A.obs[r["batch"]].astype(str).values
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,dtype=np.float32)
    # drop tiny classes (<10)
    import collections; cnt=collections.Counter(y_raw); keep=np.array([cnt[v]>=10 for v in y_raw])
    X,y_raw,b=X[keep],y_raw[keep],b[keep]
    if len(X)>NCELL:
        idx=rng.choice(len(X),NCELL,replace=False); X,y_raw,b=X[idx],y_raw[idx],b[idx]
    classes=np.unique(y_raw); cls2i={c:i for i,c in enumerate(classes)}; y=np.array([cls2i[v] for v in y_raw])
    if len(classes)<3: continue
    # held-out batch = largest batch with >=200 cells and not all classes missing
    bvals,bc=np.unique(b,return_counts=True); order=bvals[np.argsort(-bc)]
    test_batch=None
    for bb in order:
        te=b==bb
        if te.sum()>=200 and (~te).sum()>=500 and len(np.unique(y[~te]))>=3: test_batch=bb;break
    if test_batch is None: continue
    te=b==test_batch; trv=~te
    # representation: HVG2000 (by variance) + PCA50 on standardized HVG
    var=X.var(0); hv=np.argsort(-var)[:2000]; Xh=X[:,hv]
    sc_=StandardScaler().fit(Xh[trv]); Xhs=sc_.transform(Xh)
    ncomp=min(50,Xhs.shape[1]-1,trv.sum()-1)
    pca=PCA(n_components=ncomp,random_state=0).fit(Xhs[trv]); Xp=pca.transform(Xhs)
    # random split (matched sizes) for leakage gap + conformal random coverage
    nte=te.sum(); perm=rng.permutation(len(y)); rnd_te=np.zeros(len(y),bool); rnd_te[perm[:nte]]=True; rnd_tr=~rnd_te
    yb_te=label_binarize(y[te],classes=range(len(classes)))
    yb_rte=label_binarize(y[rnd_te],classes=range(len(classes)))
    arow={"atlas":r["file"].replace("_prepped.h5ad","").replace(".h5ad",""),"n":int(len(y)),
          "n_ct":int(len(classes)),"n_batch":int(len(bvals)),"test_batch":str(test_batch),"methods":{}}
    for m in METHODS:
        feat = Xh if m=="hvg-logreg" else Xp
        mm = "pca-logreg" if m=="hvg-logreg" else m
        try:
            # cross-batch
            P,L,clf = reps_and_probe(feat[trv],y[trv],feat[te], mm)
            xb_auroc=macro_auroc(yb_te,P)
            pred=P.argmax(1); acc_full=accuracy_score(y[te],pred)
            ece_raw=ece(P,pred,y[te])
            # temp scale (only meaningful for logit-bearing); fit T on train-as-cal via cv proxy: use random split logits
            # calibration split: take 40% of train as cal
            ntr=trv.sum(); calmask=np.zeros(len(y),bool)
            tri=np.where(trv)[0]; rng.shuffle(tri); cal_idx=tri[:int(.4*ntr)]; fit_idx=tri[int(.4*ntr):]
            Pc,Lc,clf2=reps_and_probe(feat[fit_idx],y[fit_idx],feat[cal_idx],mm)
            if mm=="pca-logreg":
                T=temp_scale(clf2.decision_function(feat[cal_idx]),y[cal_idx])
                from scipy.special import softmax
                ece_temp=ece(softmax(clf.decision_function(feat[te])/T,axis=1),pred,y[te])
            else: T=np.nan; ece_temp=np.nan
            # conformal: cal on held-in cal_idx, test on cross-batch te and random rnd_te
            Pte=P
            Pcal,_,_=reps_and_probe(feat[fit_idx],y[fit_idx],feat[cal_idx],mm)
            cov_xb,sz_xb=conformal_cov(Pcal,y[cal_idx],Pte,y[te])
            Prnd,_,_=reps_and_probe(feat[rnd_tr],y[rnd_tr],feat[rnd_te],mm)
            cov_rnd,_=conformal_cov(Pcal,y[cal_idx],Prnd,y[rnd_te])
            # random-split auroc for leakage gap
            rnd_auroc=macro_auroc(yb_rte,Prnd)
            # abstention acc@80 coverage
            conf=Pte.max(1); thr=np.quantile(conf,0.2); kp=conf>=thr
            acc80=accuracy_score(y[te][kp],pred[kp]) if kp.sum() else np.nan
            # shuffle control
            ys=y[trv].copy(); rng.shuffle(ys); Ps,_,_=reps_and_probe(feat[trv],ys,feat[te],mm)
            sh_auroc=macro_auroc(yb_te,Ps)
            arow["methods"][m]=dict(xb_auroc=xb_auroc,rnd_auroc=rnd_auroc,leak_gap=float(rnd_auroc-xb_auroc) if rnd_auroc==rnd_auroc else np.nan,
                ece_raw=ece_raw,ece_temp=ece_temp,cov_rnd=cov_rnd,cov_xb=cov_xb,cov_gap=float(cov_rnd-cov_xb),
                acc_full=float(acc_full),acc_at80=float(acc80),shuffle_auroc=sh_auroc)
        except Exception as e:
            arow["methods"][m]=dict(error=str(e)[:80])
    results.append(arow)
    print(f"{arow['atlas'][:28]:28s} n={arow['n']:5d} ct={arow['n_ct']:2d} | "
          + " ".join(f"{m}:cov_gap={arow['methods'][m].get('cov_gap',float('nan')):.3f}" for m in ["pca-logreg","knn","rf"] if 'cov_gap' in arow['methods'][m]), flush=True)
json.dump(results, open("expand_results/multiatlas_baseline.json","w"), indent=1)
print("DONE atlases:",len(results))
