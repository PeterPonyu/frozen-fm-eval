# Depth probe #2: batch-shift dose-response.
# The paper claims cross-batch conformal coverage collapse is a GENERAL exchangeability failure whose
# magnitude "scales with batch-effect strength" -- but states that scaling only qualitatively. Here we
# make it quantitative: for each of the 24 atlases we measure how strongly the held-out test batch
# violates exchangeability (= how discriminable test-batch cells are from training cells, a CV AUROC in
# the same PCA space the coverage gap was computed in), and correlate it with the stored coverage gap.
# Prediction: bigger batch shift -> bigger coverage collapse, across atlases AND across methods.
import anndata as ad, numpy as np, json, os, warnings, scipy.sparse as sp, collections
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr, pearsonr
man=json.load(open("expand_results/atlas_manifest.json"))
DIRS={os.path.basename(d):d for d in [
 ".../data/datasets/extra_preprocessed",".../data/datasets/CancerDatasets",
 ".../data/datasets/CancerDatasets2",".../data/datasets/DevelopmentDatasets",
 ".../data/datasets/DevelopmentDatasets2"]}
usable=[r for r in man if r.get("usable")]; NCELL=6000
GAP={r["atlas"]:r for r in json.load(open("expand_results/multiatlas_baseline.json"))}

def shift_auroc(Z,is_test):
    # 5-fold CV AUROC of predicting test-batch membership from representation Z (exchangeability violation)
    if is_test.sum()<20 or (~is_test).sum()<20: return float("nan")
    cv=StratifiedKFold(5,shuffle=True,random_state=0)
    P=cross_val_predict(LogisticRegression(max_iter=200),StandardScaler().fit_transform(Z),
                        is_test.astype(int),cv=cv,method="predict_proba")[:,1]
    return float(roc_auc_score(is_test.astype(int),P))

rows=[]
for r in usable:
    f=os.path.join(DIRS[r["dir"]],r["file"]); atlas=r["file"].replace("_prepped.h5ad","").replace(".h5ad","")
    try: A=ad.read_h5ad(f)
    except Exception as e: print("SKIP",r["file"],str(e)[:40],flush=True); continue
    y0=A.obs[r["ct"]].astype(str).values; b=A.obs[r["batch"]].astype(str).values
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float32)
    cnt=collections.Counter(y0); keep=np.array([cnt[v]>=10 for v in y0]); X,y0,b=X[keep],y0[keep],b[keep]
    rng=np.random.RandomState(20260623)  # per-atlas deterministic subsample (independent of loop order)
    if len(X)>NCELL: idx=rng.choice(len(X),NCELL,replace=False); X,y0,b=X[idx],y0[idx],b[idx]
    cls=np.unique(y0)
    if len(cls)<3: continue
    y=np.array([np.where(cls==v)[0][0] for v in y0]); NC=len(cls)
    if atlas not in GAP: print("no gap row for",atlas,flush=True); continue
    tb=str(GAP[atlas]["test_batch"])             # force the SAME held-out batch the coverage gap used
    te=b==tb
    if te.sum()<20 or (~te).sum()<100: print("skip tiny test batch",atlas,flush=True); continue
    nt=np.where(~te)[0]; rng.shuffle(nt); a1=int(.5*len(nt)); fit_i=nt[:a1]
    var=X.var(0); hv=np.argsort(-var)[:2000]; Xh=X[:,hv]
    ssh=StandardScaler().fit(Xh[fit_i]); Xhs=ssh.transform(Xh)
    ncomp=min(50,Xhs.shape[1]-1,len(fit_i)-1); pca=PCA(ncomp,random_state=0).fit(Xhs[fit_i]); Xp=pca.transform(Xhs)
    bs_pca=shift_auroc(Xp,te); bs_hvg=shift_auroc(Xhs,te)
    gm=GAP[atlas]["methods"]
    row=dict(atlas=atlas,n_batch=int(GAP[atlas]["n_batch"]),test_batch=str(tb),
             batch_shift_pca=round(bs_pca,4),batch_shift_hvg=round(bs_hvg,4),
             cov_gap={m:gm[m].get("cov_gap") for m in gm})
    rows.append(row); print(f"{atlas[:24]:24s} shift_pca={bs_pca:.3f} cov_gap[pca-logreg]={gm['pca-logreg'].get('cov_gap'):+.3f}",flush=True)

# ---- correlations ----
def corr(xs,ys):
    xs=np.array(xs); ys=np.array(ys); ok=np.isfinite(xs)&np.isfinite(ys); xs,ys=xs[ok],ys[ok]
    if len(xs)<4: return None
    sr,sp_=spearmanr(xs,ys); pr,pp=pearsonr(xs,ys)
    return dict(n=int(len(xs)),spearman=round(float(sr),3),spearman_p=float(sp_),pearson=round(float(pr),3),pearson_p=float(pp))
# (a) primary: pca-logreg gap vs PCA batch-shift, across atlases
prim=corr([r["batch_shift_pca"] for r in rows],[r["cov_gap"]["pca-logreg"] for r in rows])
# (b) pooled across the 4 PCA-based methods (shared PCA shift x)
PCAM=["pca-logreg","knn","centroid","rf"]
px=[r["batch_shift_pca"] for r in rows for m in PCAM]; py=[r["cov_gap"][m] for r in rows for m in PCAM]
pool=corr(px,py)
# (c) hvg-logreg gap vs HVG batch-shift (independent representation)
hv=corr([r["batch_shift_hvg"] for r in rows],[r["cov_gap"]["hvg-logreg"] for r in rows])
stats=dict(primary_pca_logreg=prim,pooled_4pca_methods=pool,hvg_logreg=hv)
json.dump({"rows":rows,"stats":stats},open("expand_results/batch_shift_dose_response.json","w"),indent=1)
print("\n=== BATCH-SHIFT DOSE-RESPONSE ===")
print("atlases:",len(rows))
print("primary  (pca-logreg gap vs PCA batch-shift):",prim)
print("pooled   (4 PCA methods):",pool)
print("hvg-logreg (independent rep):",hv)
lo=[r["cov_gap"]["pca-logreg"] for r in rows if r["batch_shift_pca"]<0.9]
hi=[r["cov_gap"]["pca-logreg"] for r in rows if r["batch_shift_pca"]>=0.9]
print(f"\nweak shift (<0.9): n={len(lo)} mean cov_gap={np.mean(lo):+.3f}")
print(f"strong shift(>=0.9): n={len(hi)} mean cov_gap={np.mean(hi):+.3f}")
print("DONE")
