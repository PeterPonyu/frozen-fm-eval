# Cover Letter

**To:** The Editor-in-Chief, *Neurocomputing*
**Date:** 2026-07-01
**Manuscript:** Auditing Frozen Foundation-Model Embeddings: A Leakage-Controlled, Calibration-Aware Evaluation Protocol Stress-Tested on Single-Cell Genomics
**Author:** Zeyu Fu (Army Medical University, Chongqing, China)

Dear Editor,

I submit the enclosed manuscript for consideration as an original research article in *Neurocomputing*.

The contribution is a method for *evaluating* pretrained neural representations, not a new model or system. Foundation models are increasingly deployed as frozen feature extractors and ranked on leaderboards against cheaper baselines, yet the yardsticks are entangled with the methods they rank: agreement metrics presuppose a ground truth derived in the space they grade (label circularity), tokenizer coverage is confounded with representation quality, and the exchangeability that underwrites conformal guarantees is broken by distribution shift. A leaderboard rank cannot say *why* a method placed where it did, cannot be intervened on, and inherits whatever confounds its metric carries. The paper offers a reusable, deterministic, leakage-controlled, calibration-aware protocol that addresses these three general hazards.

I anticipate the reviewer question "where is the new model or system?" and answer it directly: the deliverable is the evaluation instrument itself. It pairs every agreement metric with a non-linear probe and a reference-free structure metric to expose circularity; promotes tokenizer/vocabulary coverage to a measured variable; argues parity with an equivalence test (TOST) rather than an absent p-value; and adds a reliability axis of split-conformal coverage, ECE calibration, and selective abstention. Its signature move replaces a binary leaderboard verdict with a dose-response analysis that places the criticized method class *on* the curve, converting a rank into a manipulable, causal mechanism (vocabulary coverage, batch-shift strength, spatial-aggregation dose). The framing is deliberately domain-general; single-cell genomics is the stress test that makes circularity, coverage mismatch, and exchangeability failure severe at once, across four modalities, five foundation-model families, and a supporting meta-analysis. I believe this methodological focus fits the journal's scope in the evaluation and reliability of neural learning systems.

The framework code is released at https://github.com/PeterPonyu/frozen-fm-eval, and a full reproducibility archive, including a provenance manifest with SHA-256 hashes of every result table and script and the captured runtime environment, is deposited under Zenodo concept DOI 10.5281/zenodo.21071826.

I confirm that this work is original, has not been published previously, and is not under consideration elsewhere. There is a single author and no competing interests to declare.

Thank you for your consideration.

Sincerely,
Zeyu Fu
Institute of Combined Injury, College of Preventive Medicine, Army Medical University, Chongqing, China
fuzeyu99@126.com | ORCID 0009-0001-8329-0108
