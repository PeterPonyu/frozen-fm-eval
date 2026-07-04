#!/usr/bin/env python3
# Formal test that the coverage-gap-vs-batch-shift slope differs between scRNA and scATAC,
# upgrading the qualitative "rho=0.82 vs -0.24" contrast to a tested interaction.
# scRNA points = batch_shift_dose_response.json (headline = PCA+logreg). scATAC = scatac_batch_shift.json.
# Headline test: permutation test on Delta(slope) = slope_scRNA - slope_scATAC (modality labels permuted).
import json, os, numpy as np
from scipy import stats
os.chdir(os.path.join(os.path.dirname(__file__), ".."))

bsd = json.load(open("expand_results/batch_shift_dose_response.json"))["rows"]
rna = [(float(r["batch_shift_pca"]), float(r["cov_gap"]["pca-logreg"])) for r in bsd
       if r.get("batch_shift_pca") is not None and r["cov_gap"].get("pca-logreg") is not None]
rna = np.array([p for p in rna if np.isfinite(p[0]) and np.isfinite(p[1])])

sba = json.load(open("expand_results/scatac_batch_shift.json"))["rows"]
atac = np.array([(float(r["shift_auroc"]), float(r["cov_gap"])) for r in sba
                 if r.get("shift_auroc") is not None and np.isfinite(float(r["shift_auroc"]))])

def slope(xy):  # OLS slope of gap on shift
    x, y = xy[:, 0], xy[:, 1]
    return np.polyfit(x, y, 1)[0]

s_rna, s_atac = slope(rna), slope(atac)
rho_rna = stats.spearmanr(rna[:, 0], rna[:, 1]).correlation
rho_atac = stats.spearmanr(atac[:, 0], atac[:, 1]).correlation
d_slope_obs = s_rna - s_atac

# pooled permutation test on the slope difference (H0: same slope; modality is exchangeable)
X = np.vstack([rna, atac]); n_rna = len(rna)
rng = np.random.default_rng(0); B = 20000; cnt = 0
for _ in range(B):
    idx = rng.permutation(len(X)); a, b = X[idx[:n_rna]], X[idx[n_rna:]]
    if abs(slope(a) - slope(b)) >= abs(d_slope_obs): cnt += 1
p_perm = (cnt + 1) / (B + 1)

# interaction-coefficient view: gap ~ b0 + b1*shift + b2*atac + b3*(shift*atac); test b3 (= slope diff)
mod = np.r_[np.zeros(len(rna)), np.ones(len(atac))]  # 1 = scATAC
sh = X[:, 0]; D = np.c_[np.ones(len(X)), sh, mod, sh * mod]; yv = X[:, 1]
beta, *_ = np.linalg.lstsq(D, yv, rcond=None)
resid = yv - D @ beta; dof = len(X) - D.shape[1]
sigma2 = (resid @ resid) / dof; cov = sigma2 * np.linalg.inv(D.T @ D)
se_b3 = np.sqrt(cov[3, 3]); t_b3 = beta[3] / se_b3; p_b3 = 2 * stats.t.sf(abs(t_b3), dof)

# restricted to the overlapping shift range [0.85, 0.96]
ov = lambda xy: xy[(xy[:, 0] >= 0.85) & (xy[:, 0] <= 0.96)]
rna_o, atac_o = ov(rna), ov(atac)

res = dict(
    n_scRNA=len(rna), n_scATAC=len(atac),
    slope_scRNA=round(s_rna, 4), slope_scATAC=round(s_atac, 4), slope_diff=round(d_slope_obs, 4),
    spearman_scRNA=round(float(rho_rna), 3), spearman_scATAC=round(float(rho_atac), 3),
    perm_p_slopediff=round(p_perm, 5), interaction_b3=round(float(beta[3]), 4),
    interaction_t=round(float(t_b3), 2), interaction_p=float(f"{p_b3:.2e}"),
    overlap_0p85_0p96=dict(n_scRNA=len(rna_o), n_scATAC=len(atac_o),
                           mean_gap_scRNA=round(float(rna_o[:, 1].mean()), 4) if len(rna_o) else None,
                           mean_gap_scATAC=round(float(atac_o[:, 1].mean()), 4) if len(atac_o) else None),
)
json.dump(res, open("expand_results/modality_contrast_test.json", "w"), indent=1)
print(json.dumps(res, indent=1))
print(f"\n=> scRNA slope {s_rna:+.3f} vs scATAC slope {s_atac:+.3f}; difference {d_slope_obs:+.3f}")
print(f"   permutation p (slope diff) = {p_perm:.1e};  interaction t={t_b3:.1f}, p={p_b3:.1e}")
if len(rna_o) and len(atac_o):
    print(f"   overlapping shift 0.85-0.96: scRNA mean gap {rna_o[:,1].mean():+.3f} (n={len(rna_o)}) vs scATAC {atac_o[:,1].mean():+.3f} (n={len(atac_o)})")
