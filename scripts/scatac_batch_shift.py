# Push: put scATAC on the SAME batch-shift dose-response ruler as scRNA, to quantify the modality
# contrast (cluster I). For each held-out sample we measure batch-shift strength and the conformal
# coverage gap (+ ECE gap), in 3 representations: peak-LSI (full 130k) and two scATAC FMs (20k subset).
# Prediction: unlike scRNA (steep gap-vs-shift slope), scATAC stays near zero gap even at high shift.
import numpy as np, json, warnings, collections
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
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
def shift_auroc(Z,is_test):
    if is_test.sum()<20 or (~is_test).sum()<20: return float("nan")
    cv=StratifiedKFold(5,shuffle=True,random_state=0)
    P=cross_val_predict(LogisticRegression(max_iter=200),StandardScaler().fit_transform(Z),
                        is_test.astype(int),cv=cv,method="predict_proba")[:,1]
    return float(roc_auc_score(is_test.astype(int),P))

REPS={}
d=np.load("scatac_results/scatac_lsi.npz",allow_pickle=True)
REPS["peak-LSI (full)"]=(np.asarray(d["lsi"]),np.asarray(d["y"]).astype(str),np.asarray(d["sample"]).astype(str))
for nm,fn in (("ChromFound (FM)","scatac_chromfound_emb"),("Atacformer (FM)","scatac_atacformer_emb")):
    e=np.load(f"scatac_results/{fn}.npz",allow_pickle=True)
    REPS[nm]=(np.asarray(e["emb"]),np.asarray(e["y"]).astype(str),np.asarray(e["sample"]).astype(str))

rng=np.random.RandomState(20260625); rows=[]
for rn,(Z,y0,samp) in REPS.items():
    cls=np.unique(y0); yi=np.array([np.where(cls==v)[0][0] for v in y0]); NC=len(cls)
    for s in np.unique(samp):
        te=samp==s
        if te.sum()<150 or (~te).sum()<600 or len(np.unique(yi[~te]))<3: continue
        nt=np.where(~te)[0]; rng.shuffle(nt); a1=int(.5*len(nt)); a2=int(.75*len(nt))
        fit_i,cal_i,rnd_i=nt[:a1],nt[a1:a2],nt[a2:]
        if len(cal_i)<40 or len(rnd_i)<40: continue
        try:
            ss=StandardScaler().fit(Z[fit_i]); Zs=ss.transform(Z)
            clf=LogisticRegression(max_iter=300).fit(Zs[fit_i],yi[fit_i])
            Pcal=fullP(clf,Zs[cal_i],NC); Pxb=fullP(clf,Zs[te],NC); Prnd=fullP(clf,Zs[rnd_i],NC)
            cov_rnd=cov(Pcal,yi[cal_i],Prnd,yi[rnd_i]); cov_xb=cov(Pcal,yi[cal_i],Pxb,yi[te])
            ece_rnd=ece(Prnd,Prnd.argmax(1),yi[rnd_i]); ece_xb=ece(Pxb,Pxb.argmax(1),yi[te])
            rows.append(dict(modality="scATAC",rep=rn,is_fm="FM" in rn,sample=str(s),
                shift_auroc=round(shift_auroc(Z,te),4),cov_gap=round(cov_rnd-cov_xb,4),
                ece_gap=round(ece_xb-ece_rnd,4),cov_xb=round(cov_xb,4)))
        except Exception as ex: print("err",rn,s,str(ex)[:50])
    pts=[r for r in rows if r["rep"]==rn and np.isfinite(r["shift_auroc"])]
    print(f"{rn:18s} n={len(pts):2d} shift[min..max]={min(r['shift_auroc'] for r in pts):.2f}..{max(r['shift_auroc'] for r in pts):.2f} "
          f"cov_gap mean={np.mean([r['cov_gap'] for r in pts]):+.3f}",flush=True)

def corr(rs,xk,yk):
    x=np.array([r[xk] for r in rs]); y=np.array([r[yk] for r in rs]); ok=np.isfinite(x)&np.isfinite(y); x,y=x[ok],y[ok]
    if len(x)<4: return None
    sr,sp_=spearmanr(x,y); return dict(n=int(len(x)),spearman=round(float(sr),3),spearman_p=float(sp_))
allp=[r for r in rows if np.isfinite(r["shift_auroc"])]
hi=[r for r in allp if r["shift_auroc"]>=0.9]
stats=dict(scatac_all_covgap_vs_shift=corr(allp,"shift_auroc","cov_gap"),
           scatac_all_ecegap_vs_shift=corr(allp,"shift_auroc","ece_gap"),
           scatac_highshift_mean_covgap=round(float(np.mean([r["cov_gap"] for r in hi])),4) if hi else None,
           scatac_highshift_n=len(hi),
           scatac_mean_covgap=round(float(np.mean([r["cov_gap"] for r in allp])),4),
           scatac_max_covgap=round(float(np.max([r["cov_gap"] for r in allp])),4))
json.dump({"rows":rows,"stats":stats},open("expand_results/scatac_batch_shift.json","w"),indent=1)
print("\n=== scATAC batch-shift (modality contrast) ===")
for k,v in stats.items(): print(f"  {k}: {v}")
print("CONTRAST: scRNA high-shift mean cov_gap +0.48..+0.51 (steep) vs scATAC above.")
print("DONE")
