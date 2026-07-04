#!/usr/bin/env python3
"""G005 fairness fix (per verifier): MATCHED-condition comparison on the SAME 20k cells
+ same top-2048 peaks. Compares, with identical audit logic:
  - ChromFound FM embedding (from npz)
  - raw top-2048 log-TF-IDF  (= the FM's literal input, no SVD)  [matched-dim control]
  - LSI = SVD(raw-TFIDF, 50)  [matched-peak, low-dim control]
  - (reference) full-peak LSI-50 on the same 20k cells
Shows whether the FM's apparent acc/ECE deficit is a real FM failure or a subset/dim artifact.
Env: dl. Run: conda run -n dl python3 scripts/scatac_fm_matched.py
"""
import os, json, numpy as np, anndata as ad, pandas as pd, scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=os.path.join(ROOT,"scatac_results")
SEED=20260623; rng=np.random.default_rng(SEED); L_PEAKS,N_CELLS=2048,20000

# --- reproduce EXACT 20k subsample + top-2048 peaks + raw log-TF-IDF (mirror embed script) ---
a=ad.read_h5ad(os.path.expanduser("~/Desktop/data/datasets/ATAC_data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5ad"))
m=pd.read_csv(os.path.join(ROOT,"raw_pulls/scatac/atac_cell_meta.csv.gz")); m["Barcode"]=m["Barcode"].astype(str); m=m.drop_duplicates("Barcode").set_index("Barcode")
obs=a.obs_names.astype(str); keep=obs.isin(m.index); a=a[np.where(keep)[0]]; md=m.loc[obs[keep]]
X=a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
df=md.reset_index(); df["_i"]=np.arange(len(df)); frac=min(1.0,N_CELLS/len(df))
sel=df.groupby(["Sample.ID","Cell.Type"],group_keys=False).apply(lambda g:g.sample(max(1,int(round(len(g)*frac))),random_state=SEED))["_i"].values
sel=np.sort(sel); Xs=X[sel]; md=md.iloc[sel]
acc_pk=np.asarray((Xs>0).sum(0)).ravel(); top=np.sort(np.argsort(-acc_pk)[:L_PEAKS])
Xb=Xs[:,top].tocsr().astype(np.float32); Xb.data[:]=1.0
rs=np.asarray(Xb.sum(1)).ravel(); rs[rs==0]=1; cs=np.asarray(Xb.sum(0)).ravel(); cs[cs==0]=1
val=np.log1p(np.asarray(Xb.multiply(1/rs[:,None]).multiply(np.log(1+Xb.shape[0]/cs)[None,:]).todense())*1e4).astype(np.float32)  # raw top-2048 TF-IDF (FM input)
# full-peak LSI-50 (same cells)
Xfull=Xs.tocsr().astype(np.float32); Xfb=Xfull.copy(); Xfb.data[:]=1.0
kp=np.asarray(Xfb.sum(0)).ravel()>=0.01*Xfb.shape[0]; Xfb=Xfb[:,kp]
rs2=np.asarray(Xfb.sum(1)).ravel(); rs2[rs2==0]=1; cs2=np.asarray(Xfb.sum(0)).ravel(); cs2[cs2==0]=1
tf=np.log1p(np.asarray(Xfb.multiply(1/rs2[:,None]).multiply(np.log(1+Xfb.shape[0]/cs2)[None,:]).todense())*1e4)
lsi_full=StandardScaler().fit_transform(TruncatedSVD(51,random_state=SEED).fit_transform(tf)[:,1:])
lsi_top=StandardScaler().fit_transform(TruncatedSVD(51,random_state=SEED).fit_transform(val)[:,1:])

d=np.load(os.path.join(OUT,"scatac_chromfound_emb.npz"),allow_pickle=True)
fm=(d["emb"]-d["emb"].mean(0))/(d["emb"].std(0)+1e-8)
y=d["y"].astype(str); sample=d["sample"].astype(str)
assert (md["Cell.Type"].values==y).all(), "cell order mismatch"  # confirm same cells/order
cls={c:i for i,c in enumerate(sorted(set(y)))}; yi=np.array([cls[c] for c in y])
samples=np.array(sorted(set(sample))); rng.shuffle(samples); tr=np.array([s in set(samples[:len(samples)//2]) for s in sample])
ALPHA=0.10
def ece(P,yt,b=15):
    cf=P.max(1); co=(P.argmax(1)==yt).astype(float); ed=np.linspace(0,1,b+1); e=0.
    for i in range(b):
        mm=(cf>ed[i])&(cf<=ed[i+1])
        if mm.sum(): e+=mm.mean()*abs(co[mm].mean()-cf[mm].mean())
    return e
def audit(F):
    clf=LogisticRegression(max_iter=400,n_jobs=-1).fit(F[tr],yi[tr])
    hi=np.where(~tr)[0]; P=clf.predict_proba(F[hi]); yh=yi[hi]; sh=sample[hi]
    a_=accuracy_score(yh,P.argmax(1)); e_=ece(P,yh)
    hs=np.array(sorted(set(sh))); rng2=np.random.default_rng(1); rng2.shuffle(hs)
    cg=np.array([s in set(hs[:len(hs)//2]) for s in sh]); tg=~cg
    s=1-P[cg][np.arange(cg.sum()),yh[cg]]; n=len(s); q=np.quantile(s,min(np.ceil((n+1)*(1-ALPHA))/n,1),method="higher")
    cov=(P[tg]>=1-q)[np.arange(tg.sum()),yh[tg]].mean()
    return round(a_,4),round(float(e_),4),round(float(cov),4)
rows={}
for name,F in [("ChromFound FM (2048-d, faithful recipe)",fm),
               ("raw top-2048 log-TFIDF (FM input, 2048-d)",val),
               ("LSI top-2048 -> SVD-50",lsi_top),
               ("LSI full-peak -> SVD-50 (ref)",lsi_full)]:
    rows[name]=audit(F)
print(f"{'substrate':42s} {'acc':>7} {'ECE':>8} {'x-samp cov':>10}")
for k,(a_,e_,c_) in rows.items(): print(f"{k:42s} {a_:>7} {e_:>8} {c_:>10}")
# FM vs its own input correlation (passthrough check)
corr=np.corrcoef(d["emb"].mean(0), val.mean(0))[0,1]
print(f"\nFM-output vs FM-input per-peak mean corr: {corr:.3f}")
json.dump({k:dict(acc=v[0],ECE=v[1],xsample_cov=v[2]) for k,v in rows.items()}|{"fm_vs_input_corr":round(float(corr),3)},
          open(os.path.join(OUT,"scatac_fm_matched.json"),"w"),indent=2)
print("saved scatac_fm_matched.json")
