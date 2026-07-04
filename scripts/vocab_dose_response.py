# Depth probe #1: vocabulary dose-response curve.
# X = fraction of an atlas's genes (and expression mass) that map into each FM's gene vocabulary;
# Y = that FM's zero-shot kNN-AUROC on the atlas (read from fm_all_audit.json).
# Tests artifact #2 (gene-vocabulary mismatch): if "FM failures" on mouse/ENSG atlases are a
# token-mapping artifact rather than a representational deficit, AUROC should rise monotonically
# with vocabulary coverage and plateau at PCA-parity once coverage is high.
import anndata as ad, numpy as np, json, os, glob, warnings, pickle, scipy.sparse as sp, pandas as pd
warnings.filterwarnings("ignore")
from scipy.stats import spearmanr, pearsonr
M=".../data/models"
# ---- vocabularies (per FM) ----
GF=set(pickle.load(open(f"{M}/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl","rb")).keys())  # human symbols
SYM2ENS=pickle.load(open(f"{M}/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl","rb"))          # symbol->ENSG
SCGPT=set(json.load(open(f"{M}/scGPT-human/vocab.json")).keys())                                   # human symbols
SCF=set(pd.read_csv("expand_results/scf_gene_index.tsv",sep="\t")["gene_name"])                    # human symbols
CELLPLM=set(json.load(open(f"{M}/CellPLM/20230926_85M.config.json"))["gene_list"])                 # ENSG ids
import torch
UCE=set(torch.load(f"scripts/uce_repo/model_files/protein_embeddings/Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt",weights_only=False).keys())  # human symbols
# membership test per FM: returns boolean mask over var_names
def maps(fm, vn):
    if fm in("Geneformer-V2-104M","Geneformer-V2-316M"): V=GF
    elif fm=="scGPT": V=SCGPT
    elif fm=="scFoundation": V=SCF
    elif fm=="UCE": V=UCE
    elif fm=="CellPLM":  # symbols -> ENSG -> membership in CellPLM ENSG list
        return np.array([ (vn[i] in SYM2ENS and SYM2ENS[vn[i]] in CELLPLM) for i in range(len(vn)) ])
    else: return None
    return np.array([x in V for x in vn])

# ---- atlas list (identical to fm_all_audit.py) ----
LR={os.path.basename(f)[:-5]:f for f in glob.glob("expand_results/labeled_raw/*.h5ad")}
ATL=[("GSE130148_lung",".../data/datasets/DevelopmentDatasets2/GSE130148_LungHmDev.h5ad"),
     ("GSE165784_retina",".../data/datasets/DevelopmentDatasets2/GSE165784_RetinaHmDev.h5ad"),
     ("lung24k",".../data/datasets/DevelopmentDatasets/lung.h5ad")]
for nm,f in sorted(LR.items()): ATL.append(("lr_"+nm,f))

AUD={r["atlas"]:r for r in json.load(open("expand_results/fm_all_audit.json"))}
FMS=["scGPT","Geneformer-V2-104M","Geneformer-V2-316M","scFoundation","CellPLM","UCE"]
rows=[]
for suffix,f in ATL:
    if suffix not in AUD: continue
    try: A=ad.read_h5ad(f)
    except Exception as e: print("skip",suffix,str(e)[:40]); continue
    vn=[str(x) for x in A.var_names]
    X=A.X; X=X.toarray() if sp.issparse(X) else np.asarray(X); X=np.asarray(X,np.float64)
    gmass=X.sum(0); total=gmass.sum()  # per-gene total expression
    naming=AUD[suffix].get("naming")
    for fm in FMS:
        rep=AUD[suffix]["reps"].get(fm)
        if not rep or "knn_auroc" not in rep or rep["knn_auroc"]!=rep["knn_auroc"]: continue  # need finite AUROC
        m=maps(fm,vn)
        frac_gene=float(m.mean())
        frac_expr=float(gmass[m].sum()/total) if total>0 else 0.0
        rows.append(dict(atlas=suffix,fm=fm,naming=naming,ngene=int(len(vn)),
                         n_mapped=int(m.sum()),frac_gene=round(frac_gene,4),
                         frac_expr=round(frac_expr,4),knn_auroc=rep["knn_auroc"],
                         expr_R2=rep.get("expr_R2")))
    print(suffix,"done",flush=True)

# ---- correlations ----
def corr(rs,xkey):
    x=np.array([r[xkey] for r in rs]); y=np.array([r["knn_auroc"] for r in rs])
    if len(rs)<4: return None
    sr,sp_=spearmanr(x,y); pr,pp=pearsonr(x,y)
    return dict(n=len(rs),spearman=round(float(sr),3),spearman_p=float(sp_),
                pearson=round(float(pr),3),pearson_p=float(pp))
stats={"pooled_all":{},"universal3":{}}  # universal3 = the 3 FMs run on every atlas (full x-range)
UNIV={"scGPT","Geneformer-V2-104M","Geneformer-V2-316M"}
for xk in("frac_gene","frac_expr"):
    stats["pooled_all"][xk]=corr(rows,xk)
    stats["universal3"][xk]=corr([r for r in rows if r["fm"] in UNIV],xk)
# per-FM (only the universal 3 span enough range to be meaningful)
stats["per_fm"]={}
for fm in FMS:
    rs=[r for r in rows if r["fm"]==fm]
    stats["per_fm"][fm]=corr(rs,"frac_expr")
out={"rows":rows,"stats":stats}
json.dump(out,open("expand_results/vocab_dose_response.json","w"),indent=1)
print("\n=== DOSE-RESPONSE ===")
print(f"points: {len(rows)}  (universal-3 span full range; scf/cellplm/uce high-frac only)")
print("pooled-all   frac_expr:",stats["pooled_all"]["frac_expr"])
print("universal-3  frac_expr:",stats["universal3"]["frac_expr"])
print("universal-3  frac_gene:",stats["universal3"]["frac_gene"])
print("\nlow-coverage (frac_expr<0.5) AUROC vs high (>=0.5), universal-3:")
u=[r for r in rows if r["fm"] in UNIV]
lo=[r["knn_auroc"] for r in u if r["frac_expr"]<0.5]; hi=[r["knn_auroc"] for r in u if r["frac_expr"]>=0.5]
print(f"  low  n={len(lo)} mean AUROC={np.mean(lo):.3f}")
print(f"  high n={len(hi)} mean AUROC={np.mean(hi):.3f}")
print("DONE")
