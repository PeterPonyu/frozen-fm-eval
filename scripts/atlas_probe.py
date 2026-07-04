import anndata as ad, numpy as np, glob, os, json, warnings
warnings.filterwarnings("ignore")
dirs = [".../data/datasets/extra_preprocessed",
        ".../data/datasets/CancerDatasets",
        ".../data/datasets/CancerDatasets2",
        ".../data/datasets/DevelopmentDatasets",
        ".../data/datasets/DevelopmentDatasets2"]
CT = ["cell_type","celltype","CellType","cell.type","cluster","leiden","louvain","annotation","labels","cell_ontology_class","majority_voting","predicted_labels"]
BATCH = ["batch","donor","sample","Sample","patient","Patient","orig.ident","dataset","study","donor_id","sample_id","batchlb","condition"]
def pick(cols, cands):
    low = {c.lower(): c for c in cols}
    for k in cands:
        if k.lower() in low: return low[k.lower()]
    return None
rows=[]
for d in dirs:
    for f in sorted(glob.glob(os.path.join(d,"*.h5ad"))):
        try:
            a = ad.read_h5ad(f, backed="r")
        except Exception as e:
            rows.append(dict(file=os.path.basename(f), dir=os.path.basename(d), error=str(e)[:60])); continue
        obs=a.obs; cols=list(obs.columns)
        ctcol=pick(cols,CT); bcol=pick(cols,BATCH)
        # find ANY categorical col with 2..n levels as fallback batch
        nct = int(obs[ctcol].nunique()) if ctcol else 0
        nb = int(obs[bcol].nunique()) if bcol else 0
        # raw counts? check a small slice max/integer
        israw=None
        try:
            import scipy.sparse as sp
            X=a.X[:200]; X=X.toarray() if sp.issparse(X) else np.asarray(X)
            israw = bool(np.allclose(X, np.round(X)) and X.max()>30)
        except Exception: israw=None
        rows.append(dict(file=os.path.basename(f), dir=os.path.basename(d), n=int(a.n_obs), g=int(a.n_vars),
                         ct=ctcol, nct=nct, batch=bcol, nb=nb, raw=israw,
                         usable=bool(ctcol and bcol and nb>=2 and nct>=3 and a.n_obs>=800)))
        a.file.close() if hasattr(a,"file") else None
json.dump(rows, open("expand_results/atlas_manifest.json","w"), indent=1)
us=[r for r in rows if r.get("usable")]
print(f"TOTAL {len(rows)} | USABLE {len(us)}")
print(f"{'file':40s} {'n':>7} {'ct':>14} {'nct':>4} {'batch':>12} {'nb':>4} raw")
for r in sorted(us, key=lambda x:-x['n']):
    print(f"{r['file'][:40]:40s} {r['n']:7d} {str(r['ct'])[:14]:>14} {r['nct']:4d} {str(r['batch'])[:12]:>12} {r['nb']:4d} {r['raw']}")
print("--- NOT usable (reason) ---")
for r in rows:
    if not r.get("usable"):
        why = r.get("error") or f"ct={r.get('ct')} nct={r.get('nct')} batch={r.get('batch')} nb={r.get('nb')} n={r.get('n')}"
        print(f"{r['file'][:40]:40s} {why}")
