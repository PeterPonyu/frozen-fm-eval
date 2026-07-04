#!/usr/bin/env python3
# Provenance manifest: per-result + per-script SHA-256 hashes plus the runtime environment.
# Closes the "no per-result provenance manifest" gap noted in the paper's Reproducibility section.
# Output: expand_results/provenance_manifest.json (machine) + provenance_manifest.md (human).
import hashlib, json, os, sys, glob, platform, datetime

os.chdir(os.path.join(os.path.dirname(__file__), ".."))  # research/sc-fm-benchmark

def sha256(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b: break
            h.update(b)
    return h.hexdigest()

def entry(path):
    st = os.stat(path)
    return dict(path=path, sha256=sha256(path), bytes=st.st_size,
                mtime=datetime.datetime.utcfromtimestamp(st.st_mtime).isoformat() + "Z")

# result tables consumed by the paper (tracked JSONs) + scATAC result JSONs
results = sorted(glob.glob("expand_results/*.json")) + sorted(glob.glob("scatac_results/*.json"))
results = [p for p in results if os.path.basename(p) != "provenance_manifest.json"]
scripts = sorted(glob.glob("scripts/*.py"))

def pkg(name):
    try:
        import importlib.metadata as m
        return m.version(name)
    except Exception:
        return None

env = dict(python=sys.version.split()[0], platform=platform.platform(),
           packages={k: pkg(k) for k in ["numpy", "scipy", "scikit-learn", "anndata", "pandas", "torch"]})

manifest = dict(
    generated_utc=datetime.datetime.utcnow().isoformat() + "Z",
    note="SHA-256 of every paper-backing result table and every first-party script, plus the runtime env. "
         "Re-run a script and re-hash to confirm a result table is reproduced bit-for-bit.",
    environment=env,
    n_results=len(results), n_scripts=len(scripts),
    results=[entry(p) for p in results],
    scripts=[entry(p) for p in scripts],
)
json.dump(manifest, open("expand_results/provenance_manifest.json", "w"), indent=1)

with open("expand_results/provenance_manifest.md", "w") as f:
    f.write(f"# Provenance manifest\n\nGenerated {manifest['generated_utc']}. "
            f"{len(results)} result tables, {len(scripts)} scripts.\n\n")
    f.write(f"**Environment:** python {env['python']}; "
            + ", ".join(f"{k} {v}" for k, v in env['packages'].items() if v) + "\n\n")
    f.write("## Result tables (paper-backing)\n\n| file | sha256 (first 16) | bytes |\n|---|---|---|\n")
    for e in manifest["results"]:
        f.write(f"| `{os.path.basename(e['path'])}` | `{e['sha256'][:16]}` | {e['bytes']} |\n")
    f.write("\n## First-party scripts\n\n| file | sha256 (first 16) | bytes |\n|---|---|---|\n")
    for e in manifest["scripts"]:
        f.write(f"| `{os.path.basename(e['path'])}` | `{e['sha256'][:16]}` | {e['bytes']} |\n")

print(f"manifest: {len(results)} result tables, {len(scripts)} scripts hashed")
print(f"env: python {env['python']}; " + ", ".join(f"{k} {v}" for k, v in env['packages'].items() if v))
print("wrote expand_results/provenance_manifest.json + .md")
