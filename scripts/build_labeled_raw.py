# Transfer cell_type+batch labels from *_prepped (labeled, normalized) to raw-count GSE twins by barcode.
# Produces FM-embeddable labeled raw-count atlases.
import anndata as ad, numpy as np, json, os, glob, warnings, scipy.sparse as sp
warnings.filterwarnings("ignore")
man=json.load(open("expand_results/atlas_manifest.json"))
DIRS={os.path.basename(d):d for d in [
 ".../data/datasets/extra_preprocessed",".../data/datasets/CancerDatasets",
 ".../data/datasets/CancerDatasets2",".../data/datasets/DevelopmentDatasets",
 ".../data/datasets/DevelopmentDatasets2"]}
# labeled prepped atlases (have cell_type + batch)
prepped=[r for r in man if r.get("ct") and r.get("nct",0)>=3 and r.get("nb",0)>=2 and "_prepped" in r["file"]]
# raw candidate pool: integer counts, NO cell_type (the unlabeled raw twins) + the few already-raw labeled handled separately
raw_pool=[r for r in man if r.get("raw") and not r.get("ct")]
def names(r): 
    a=ad.read_h5ad(os.path.join(DIRS[r["dir"]],r["file"]),backed="r"); return [str(x) for x in a.obs_names]
print("prepped labeled:",len(prepped),"| raw unlabeled candidates:",len(raw_pool),flush=True)
rawnames={r["file"]:set(names(r)) for r in raw_pool}
os.makedirs("expand_results/labeled_raw",exist_ok=True)
out=[]
for p in prepped:
    pa=ad.read_h5ad(os.path.join(DIRS[p["dir"]],p["file"]))
    pn=[str(x) for x in pa.obs_names]; pset=set(pn)
    best=max(raw_pool,key=lambda r: len(pset & rawnames[r["file"]]))
    ov=len(pset & rawnames[best["file"]]); frac=ov/max(1,len(pset))
    nm=p["file"].replace("_prepped.h5ad","")
    if frac<0.5:
        print(f"  SKIP {nm}: best overlap {frac:.0%} ({best['file']})",flush=True); continue
    ra=ad.read_h5ad(os.path.join(DIRS[best['dir']],best["file"]))
    ra.obs_names=[str(x) for x in ra.obs_names]; pa.obs_names=pn
    # dedupe
    ra=ra[~ra.obs_names.duplicated()].copy() if hasattr(ra.obs_names,'duplicated') else ra
    import pandas as pd
    pmap=pd.Series(pa.obs[p["ct"]].astype(str).values, index=pn)
    pbat=pd.Series(pa.obs[p["batch"]].astype(str).values, index=pn)
    common=[x for x in ra.obs_names if x in pset]
    ra=ra[common].copy(); ra.obs["cell_type"]=pmap.reindex(ra.obs_names).values; ra.obs["batch"]=pbat.reindex(ra.obs_names).values
    ra=ra[~ra.obs["cell_type"].isna()].copy()
    op=f"expand_results/labeled_raw/{nm}.h5ad"; ra.write(op)
    out.append(dict(name=nm,raw=best["file"],n=int(ra.n_obs),nct=int(ra.obs['cell_type'].nunique()),nb=int(ra.obs['batch'].nunique()),overlap=round(frac,3)))
    print(f"  OK {nm:20s} <- {best['file'][:34]:34s} n={ra.n_obs} ct={ra.obs['cell_type'].nunique()} batch={ra.obs['batch'].nunique()}",flush=True)
json.dump(out,open("expand_results/labeled_raw_manifest.json","w"),indent=1)
print("BUILT",len(out),"labeled raw atlases",flush=True)
