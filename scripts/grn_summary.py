import pandas as pd, json, os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G=os.path.join(ROOT,"raw_pulls/grn")
# scRegNet BEELINE (TFs+500): classical vs DL vs scFM-GNN, AUROC across 8 datasets
s=pd.read_csv(os.path.join(G,"scregnet/scRegNet_results_TFs500.csv"))
auroc_cols=[c for c in s.columns if c.endswith("_AUROC")]
s["meanAUROC"]=s[auroc_cols].mean(axis=1)
byc=s.groupby("method_category").meanAUROC.agg(["mean","min","max"]).round(3)
# Kendiukhov: gene-level baselines vs raw scFM attention (K562 CRISPRi, n=151)
k=pd.read_csv(os.path.join(G,"kendiukhov/kendiukhov_summary_aurocs_trivial_baselines.csv"))
attn=0.704  # Geneformer L13 raw attention (paper-reported, in kendiukhov_paper_reported)
genelevel=k[k.predictor.isin(["Mean Expression","1 - Dropout Rate","Variance"])]
base_beats_rawFM=int((genelevel.mean_auroc>attn).sum()); n=len(genelevel)
summary=dict(
 contrast_A_raw_scFM_signal_vs_gene_level_baseline=dict(
   gene_level_AUROC={r.predictor:round(r.mean_auroc,3) for r in genelevel.itertuples()},
   raw_Geneformer_attention_AUROC=attn,
   baseline_beats_rawFM=f"{base_beats_rawFM}/{n}", verdict="gene-level baselines BEAT raw scFM attention"),
 contrast_B_trained_GNN_on_scFM_embeddings_vs_classical=dict(
   classical_meanAUROC=float(byc.loc["classical","mean"]) if "classical" in byc.index else None,
   scFM_GNN_meanAUROC=float(byc.loc["scFM_GNN","mean"]) if "scFM_GNN" in byc.index else None,
   verdict="scFM-embeddings fed to a SUPERVISED GNN BEAT classical (but so do non-FM GNNs ~0.85)"))
json.dump(summary, open(os.path.join(ROOT,"grn_summary.json"),"w"), indent=2)
print("=== scRegNet meanAUROC by category ==="); print(byc.to_string())
print("\n=== Kendiukhov raw-signal contrast ==="); print(genelevel[["predictor","mean_auroc"]].to_string(index=False), f"\n  raw Geneformer attention={attn}  -> baseline beats raw-FM {base_beats_rawFM}/{n}")
