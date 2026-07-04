# The Frozen-FM-Eval protocol

A reusable design for **auditing frozen foundation-model embeddings** rather than ranking them. It is written for single-cell genomics but the moves are domain-general; the last section gives a recipe for a new setting.

## Why a protocol, not a leaderboard

A frozen FM produces an embedding; a downstream metric scores it against baselines. Three things routinely make that score mean something other than "representation quality":

1. **Metric–method entanglement (circularity).** When the "ground truth" is itself derived in the same space a baseline occupies (e.g. cell-type labels from clustering + marker scoring in expression space), agreement metrics (ARI, NMI, label-AUROC) reward proximity to the *label-generating process*, not to biology. A linear baseline that lives in that space is flattered; an FM that departs from it is penalised — regardless of which is more faithful.
2. **Coverage–quality confound.** A pretrained tokenizer reads only the genes in its vocabulary. Cross-species orthologs, symbol-vs-Ensembl mismatches, or novel panels drop most tokens, so an atlas can score the model's *tokenizer coverage* while appearing to score its *representation*.
3. **Broken exchangeability.** Split-conformal coverage and calibration are only valid under exchangeability, which cross-batch/donor annotation violates. Reported coverage then silently fails to hold out of distribution.

The protocol neutralises each, and — crucially — **does not stop at "the metric is confounded."** It converts every apparent win or loss into a *manipulable, causal mechanism* by placing the criticized method class on a dose–response curve.

## The four deconfounding moves

**(i) Pair every agreement metric with a circularity-free re-check.**
For each headline agreement score, also compute:
- a **non-linear (kNN) probe** of the same labels — a linear-vs-non-linear gap localizes an apparent baseline advantage to *probe geometry*, not representation;
- a **reference-free, label-independent structure metric** (`R²_expr`: how much held-out HVG-expression variance a partition explains). If the ground-truth labels explain *less* structure than the clustering they grade, the agreement metric is scoring a coarser cut than the data supports — the circularity is quantified, not merely asserted.

Scripts: `fair_recheck.py`, `parity_probe_types.py`, `direct_circularity.py`, `r2_null_control.py`, `circularity_hvg_robustness.py`.

**(ii) Treat tokenizer/vocabulary coverage as a measured variable.**
Record, per atlas × model, the fraction of features that map into the model's vocabulary (and the expression mass they carry). Report every quality number *against* coverage rather than marginalising over it, so a "representational deficit" that is really a coverage artifact is visible.

Scripts: `vocab_dose_response.py`.

**(iii) Argue parity by equivalence, not by an absent p-value.**
"No significant difference" from a non-significant test is not evidence of equivalence. Use a two-one-sided-test (TOST) equivalence procedure against a pre-declared margin (here ±0.02 AUROC), per family, and report the outcome honestly — a **boundary** equivalence result is reported as "no significant difference," *not* as a clean parity pass.

Scripts: `parity_pooled_cluster_ci.py`, `parity_robustness.py`.

**(iv) Score the reliability axis explicitly.**
Discrimination is not reliability. For each embedding + readout, report:
- **split-conformal coverage** under the deployment shift (cross-batch / cross-sample / cross-disease), not just in-distribution;
- **calibration** (ECE), including how it moves with readout capacity and temperature scaling;
- **selective abstention** (risk–coverage, AURC), and whether abstention is fair to rare classes.

Scripts: `scatac_audit.py`, `scatac_fm_matched.py`, `fm_all_audit.py`, `scatac_batch_shift.py`.

## The signature move: dose–response with the criticized class on the curve

A binary verdict ("A beats B") is replaced by a **graded curve** whose x-axis is a manipulable *dose* and on which the criticized method class is explicitly plotted. This turns a leaderboard rank into a mechanism that is either **attributable to a named cause** or **dissolved**:

- **Vocabulary dose** — sweep the fraction of readable genes (across atlases, and causally *within a fixed atlas* by renaming genes so the tokenizer cannot read them). An apparent cross-species "FM failure" becomes a tokenization *cliff* to chance once the vocabulary reads almost no genes, at parity otherwise (`vocab_ablation.py`, `vocab_ablation_scgpt.py`).
- **Batch-shift dose** — score coverage/ECE degradation against measured shift strength (held-out-batch discriminability). A "calibration collapse" attributed to FMs is shown to afflict classical baselines by the same amount, i.e. it is a property of the *regime* (`batch_shift_dose_response.py`, `batch_shift_fm_probe.py`, `b2_fm_shift_intrinsic.py`).
- **Spatial-aggregation dose** — smooth per-cell representations over k spatial neighbours. A "spatial-FM gap" closes once a per-cell baseline is given the one operation it lacked (`spatial_dose_response.py`, `spatial_dose_foldblocked.py`).

The rule: **whenever you would report a rank, instead find the dose that produces it, put the loser on the same curve, and report the mechanism.**

## Inputs and outputs

- **Input:** a frozen embedding matrix per (dataset, model) plus one or more baselines (PCA / HVG / scVI, or the analogue in your domain), on identical splits, with labels and a shift axis (batch/donor/disease).
- **Output per stage:** a JSON summary carrying the probe scores, the reference-free metric, coverage, the equivalence outcome, the reliability triple, and — where applicable — the fitted dose–response and the criticized class's position on it.

## Reproducibility contract

Each script is standalone and deterministic with a fixed seed; it reads inputs and writes one summary at the paths named in its header. Baselines are computed locally on the same splits as the FM. FM embedders are isolated (`scgpt_embed.py`, `geneformer_v2_embed.py`, `scfoundation_embed.py`, `cellplm_embed.py`, `uce_embed.py`, `scatac_chromfound_embed.py`, `atacformer_embed.py`); where a loader is faithful-but-not-official it is disclosed, and where it reproduces the official embeddings exactly (scGPT, Geneformer-V2) that is stated.

## Applying this to a new frozen-embedding setting

The protocol is not specific to genomics. For any frozen encoder (vision, text, audio, graph) scored against baselines:

1. **Name the label-generating process.** If your "ground truth" is derived in a space a baseline occupies, add the non-linear probe + a reference-free structure metric (variance explained on held-out raw features) and check whether the labels are coarser than unsupervised structure.
2. **Measure the tokenizer/preprocessing coverage** that varies across datasets (OOV rate, resolution, channel/modality mismatch) and report quality against it.
3. **Pre-declare an equivalence margin** and use TOST instead of reading equivalence off a non-significant test.
4. **Pick your deployment shift axis** and report conformal coverage + ECE + selective abstention *under that shift*, not in-distribution.
5. **For every rank you are tempted to report, build a dose–response** on the manipulable driver behind it and place the criticized method on the curve — attribute the effect to a cause or dissolve it.

The deliverable is the audit, not the verdict: a frozen embedding is characterised by *where on each curve it sits and why*, not by a single leaderboard cell.
