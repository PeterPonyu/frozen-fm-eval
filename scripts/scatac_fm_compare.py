import os, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT=os.path.join(ROOT,"scatac_results"); FIG=os.path.join(OUT,"figures")
M=json.load(open(os.path.join(OUT,"scatac_fm_matched.json")))
labels=["ChromFound\nFM (2048d)","raw TF-IDF\n(FM input,2048d)","LSI top-2048\n(SVD-50)","LSI full-peak\n(SVD-50)"]
keys=["ChromFound FM (2048-d, faithful recipe)","raw top-2048 log-TFIDF (FM input, 2048-d)","LSI top-2048 -> SVD-50","LSI full-peak -> SVD-50 (ref)"]
acc=[M[k]["acc"] for k in keys]; ece=[M[k]["ECE"] for k in keys]; cov=[M[k]["xsample_cov"] for k in keys]
col=["#d95f0e","#fdae6b","#9ecae1","#2c7fb8"]
fig,ax=plt.subplots(1,3,figsize=(13,4.3))
for i,(vals,title,fmt) in enumerate([(acc,"cell-type accuracy",".3f"),(ece,"ECE (lower=better)",".3f"),(cov,"x-sample conformal cov (target .90)",".3f")]):
    ax[i].bar(labels,vals,color=col); ax[i].set_title(title)
    for j,v in enumerate(vals): ax[i].text(j,v+(max(vals)*0.01),format(v,fmt),ha="center",fontsize=8)
    ax[i].tick_params(axis="x",labelsize=7)
ax[2].axhline(0.90,color="k",ls="--",lw=1)
plt.suptitle("scATAC MATCHED comparison (same 20k cells): FM ≈ its own raw input; calibration is a DIMENSIONALITY effect (2048d≈0.09 vs SVD-50≈0.01); coverage robust for all")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig4_FM_vs_FMfree.png"),dpi=150); plt.close(); print("fig4 regenerated (matched)")
