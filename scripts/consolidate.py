"""Cross-cluster consolidated summary. Clusters kept SEPARATE (different tasks/metrics)."""
import json, pandas as pd, os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ws=pd.read_csv(os.path.join(ROOT,"winrate_strata.csv")).set_index("stratum")
H=json.load(open(os.path.join(ROOT,"cluster_H_summary.json")))[0]
def row(cl,desc,key=None,wr=None,lo=None,hi=None,k=None,note=""):
    if key: r=ws.loc[key]; wr,lo,hi,k=r.win_rate_studyEW,r.ci_lo,r.ci_hi,int(r.n_studies)
    return dict(cluster=cl,claim=desc,baseline_winrate=wr,ci_lo=lo,ci_hi=hi,n_studies=k,note=note)
rows=[
 row("A","zero-shot scFM <= simple baselines (representation/annotation/integration)","cluster A vs FM",note="decisive; classical-DR 0.967"),
 row("B","linear/additive >= DL/FM (perturbation)","cluster B vs FM",note="metric-contingent; k>=3 -> 0.518"),
 row("H","zero-shot spatial FM <= simple baselines (niche-ID)",wr=H["winrate_variant_level"],k=1,
     note="self-computed, 1 dataset, 3 FMs; best baseline "+str(H["headline_ratio"])+"x best FM"),
 row("E/F","raw scFM signal <= gene-level baselines (GRN)",wr=1.0,k=1,
     note="Geneformer attention 0.70 < Variance0.88 (3/3); scFM-emb-in-trained-GNN win (usage-dependent)"),
 row("G","trajectory/velocity: NO scFM evaluated (open gap)",k=0,note="veloBench 15 methods x17 datasets, zero FMs"),
 row("D","integration metrics scIB/kBET/LISI unreliable at scale (caveat)",k=0,note="RBET kBET/LISI CV->0; scIB gameable; caveats cluster-A scIB"),
 row("C","ensemble/calibration give reliable UQ; FM UQ untested",k=0,note="popV consensus score8->95%/<=3-><50%; HCE OOD 24-32pt"),
 row("I-ATAC","scATAC cell-type calibration ROBUST to shift; FM no gain over input (first audit)",k=1,
     note="FM-free peak-LSI ECE 0.006, cov-gap ~0 (no collapse) vs scRNA SC18 +0.123. FM arm (ChromFound zero-shot, FIRST scATAC-FM calib audit): MATCHED 20k/top-2048 -> FM 0.834/ECE0.094 does NOT beat its own raw TF-IDF input 0.850/0.088; ECE gap is dimensionality (2048d~0.09 vs SVD50~0.01) not FM-vs-baseline; coverage robust all -> FM adds nothing over input (simple>=FM, fairly framed)"),
]
df=pd.DataFrame(rows); df.to_csv(os.path.join(ROOT,"ALL_CLUSTERS_SUMMARY.csv"),index=False)
print(df[["cluster","claim"]].to_string(index=False))
