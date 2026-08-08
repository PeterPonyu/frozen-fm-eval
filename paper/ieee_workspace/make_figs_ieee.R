#!/usr/bin/env Rscript
# Figures via tikzDevice: text is typeset by LaTeX (pdfLaTeX), so fonts/sizes match the body 1:1.
# Engine is pdftex to match the IEEEtran build: IEEEtran forces Times (\rmdefault=ptm), which
# resolves natively under pdflatex (T1/OT1) -> real Times in body AND figures. lualatex has no
# TUptm.fd and silently falls back to Latin Modern (wrong IEEE font), so pdftex is the correct engine.
# Labels are written as LaTeX (sanitize=FALSE): math as $...$, percent as \%, arrows as $\to$.
# Dense point/jitter layers are rasterised (ggrastr) so the .tex stays small; axes/text stay vector.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(jsonlite); library(tikzDevice); library(ggrastr)})
options(tikzDefaultEngine = "pdftex")
options(tikzMetricPackages = c(getOption("tikzMetricPackages"), "\\usepackage{amsmath}", "\\usepackage{amssymb}"),
        tikzLatexPackages = c(getOption("tikzLatexPackages"), "\\usepackage{amsmath}", "\\usepackage{amssymb}"))
# Without this, tikzDevice re-measures every unique string's glyph metrics by shelling out to a
# fresh lualatex process every run (no cache persists across Rscript invocations) -- that was the
# actual cost of a full regen here (~3.5min cold, same root cause as ../make_figs.R). Persisting
# the dictionary means only genuinely NEW strings pay the lualatex cost on subsequent runs.
options(tikzMetricsDictionary = ".tikz_metrics_ieee_pdftex")  # cwd (ieee_workspace/); engine-specific (pdftex metrics != luatex), isolated from ../paper/.tikz_metrics
base <- Sys.getenv("SCFM_BASE", "..")  # repo subdir root; run from paper/ or set SCFM_BASE
figs <- file.path(base, "paper", "ieee_workspace", "figs"); dir.create(figs, showWarnings = FALSE, recursive = TRUE)
TW <- 6.5  # \textwidth in inches (letterpaper, 1in margins) -> 1:1, no rescale in LaTeX
theme_set(theme_bw(base_size = 8) + theme(  # base 9->8: the height-trimmed compact panels were crowding text; smaller font relieves it
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 8, face = "plain", hjust = 0.5, margin = margin(b = 1.5)),  # centered; not bold; hug the panel
  strip.text = element_text(size = 7.5, margin = margin(2, 2, 2, 2)),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7, margin = margin(l = 3, r = 8)),  # gap between handle and label
  legend.key.size = unit(8, "pt"), legend.key.spacing.x = unit(4, "pt"),
  legend.margin = margin(0, 1, 0, 1), legend.box.spacing = unit(2, "pt"),  # pull top/bottom legends toward the panel so titles/legends aren't stranded
  axis.title = element_text(size = 7.5),
  axis.title.x = element_text(margin = margin(t = 2)),
  axis.title.y = element_text(margin = margin(r = 3)),  # comfortable constant gap to tick labels; per-panel distance adapts
  axis.text = element_text(size = 6.8), plot.margin = margin(3, 6, 3, 4)))  # 14pp condensation: trimmed per-panel dead margin (was 5,9,5,5)
CB <- c(A="#2c7fb8", B="#d95f0e")
SPP <- c("vocabulary-matched (human)"="#2c7fb8","vocab-mismatch (cross-species/ID)"="#d95f0e")
# Atlases whose gene var_names do NOT match the human-symbol FM vocabulary: 8 mouse-symbol + 1 ENSG-only.
# Human FMs are handicapped by tokenization here (not representation); separated from the 11 matched atlases.
MISMATCH <- c("lr_lps_mm","lr_lsk_batch","lr_progastin","lr_urine","lr_astrocytes_sci",
              "lung24k","lr_breast_hm","lr_tcell_cancer","lr_hesc_hspc_cd8")
ann <- function(t) plot_annotation(tag_levels = "A")  # long descriptions live in the LaTeX caption
# bold (A)/(B)/(C) panel tags on every multi-panel (ann-using) figure
theme_update(plot.tag = element_text(size = 10, face = "bold"), plot.tag.location = "plot")  # panel tags OUTSIDE the panel, at the plot's top-left (journal convention)
pct <- function(x) paste0(round(x*100), "\\%")
RP <- function(g, ...) g  # vector points (few hundred max); rasterise unneeded + fragile w/ tikzDevice
emit <- function(name, p, w, h){ tikz(file.path(figs, paste0(name,".tex")), width=w, height=h, standAlone=FALSE, sanitize=FALSE); print(p); dev.off(); cat(name,"ok\n") }

## ---- Fig 1: per-study baseline win-rate ----
ps <- read.csv(file.path(base, "per_study_winrate.csv")); ps <- ps[ps$cluster %in% c("A","B"), ]
ps$study_id <- gsub("_", " ", ps$study_id)  # underscores break LaTeX (sanitize=FALSE)
ps$study_id <- factor(ps$study_id, levels = ps$study_id[order(ps$cluster, ps$win_rate)])
topLev <- levels(ps$study_id)[nlevels(ps$study_id)]  # anchor for labels: top row of the discrete y-axis (annotate needs a real factor level, not a raw numeric, on a discrete scale)
f1 <- ggplot(ps, aes(win_rate, study_id, color = cluster, size = n)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.8869, linetype = "dotted", color = CB["A"]) +
  geom_vline(xintercept = 0.6533, linetype = "dotted", color = CB["B"]) +
  annotate("text", x = 0.5, y = topLev, label = "chance 0.5", vjust = -2.0, hjust = 0.5, size = 2.0, color = "grey45") +
  annotate("text", x = 0.8869, y = topLev, label = "cluster A mean", vjust = -2.0, hjust = 0.5, size = 2.0, color = CB[["A"]]) +
  annotate("text", x = 0.6533, y = topLev, label = "cluster B mean", vjust = -3.6, hjust = 0.5, size = 2.0, color = CB[["B"]]) +
  RP(geom_point(alpha = .85)) + scale_color_manual(values = CB) +
  scale_size_continuous(range = c(1.5, 5), trans = "sqrt", guide = "none") +
  scale_x_continuous(limits = c(-.02, 1.02), labels = pct) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.9))) +  # top headroom so the mean/chance labels clear the panel frame
  coord_cartesian(clip = "off") +
  labs(x = "baseline win-rate vs FM/DL (per study)", y = NULL, color = "cluster",
       title = "Per-study baseline win-rate") + theme(legend.position = "top", plot.margin = margin(14, 9, 5, 5))

## ---- Fig 2: meta win-rate by stratum ----
ws <- read.csv(file.path(base, "winrate_strata.csv"))
keep <- c("cluster A vs FM","baseline=classical-DR","cluster B vs FM",
          "cluster B vs FM (>=3 comparisons/study)","baseline=mean-baseline","baseline=linear-baseline")
lab <- c("cluster A vs FM"="A: scFM vs simple (all)","baseline=classical-DR"="A: HVG/PCA/scVI/Harmony",
         "cluster B vs FM"="B: perturbation (all)","cluster B vs FM (>=3 comparisons/study)"="B: $k\\geq3$ studies",
         "baseline=mean-baseline"="mean baseline","baseline=linear-baseline"="linear baseline")
d2 <- ws[ws$stratum %in% keep, ]; d2$lab <- lab[d2$stratum]
d2$lab <- factor(d2$lab, levels = d2$lab[order(d2$win_rate_studyEW)]); d2$grp <- ifelse(grepl("^A", d2$lab), "A", "B")
f2 <- ggplot(d2, aes(win_rate_studyEW, lab, color = grp)) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .25, na.rm = TRUE) +
  geom_point(size = 2.4) + scale_color_manual(values = CB, guide = "none") +
  geom_text(aes(label = sprintf("%.2f ($k{=}%d$)", win_rate_studyEW, n_studies)), vjust = -1, size = 2.6) +
  scale_x_continuous(limits = c(0,1.18), labels = pct, breaks = c(0,.3,.6,.9)) +  # modest right headroom so the "0.97 ($k{=}4$)" label clears the panel edge without dead space
  scale_y_discrete(expand = expansion(add = c(0.5, 1.3))) +  # top headroom so the topmost value label clears the frame
  labs(x = "win-rate (EW per study, 95\\% CI)", y = NULL, title = "Meta-analytic baseline win-rates")
emit("fig_meta", (f1 | f2) + plot_annotation(tag_levels="A"), TW, 3.8)

## ---- Fig 3: scATAC coverage vs scRNA ----
a <- fromJSON(file.path(base, "scatac_results", "scatac_audit_result.json"))
d3 <- data.frame(regime = factor(c("random","x-sample","Ctrl$\\to$AD","scRNA\nx-batch"),
                  levels=c("random","x-sample","Ctrl$\\to$AD","scRNA\nx-batch")),
                 cov = c(a$coverage_random, a$coverage_crosssample, a$coverage_control_to_AD, 0.90 - 0.123),
                 kind = c("scATAC","scATAC","scATAC","scRNA"))
f3 <- ggplot(d3, aes(regime, cov, fill = kind)) + geom_col(width = .65) +
  geom_hline(yintercept = .90, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.3f", cov)), vjust = -.4, size = 2.6) +
  scale_fill_manual(values = c(scATAC = "#1b9e77", scRNA = "#b2182b"), name = NULL) +  # scATAC: teal, distinct from the blue the B/C legend assigns to "baseline"; scRNA collapse: crimson, distinct from the orange B/C assigns to "scATAC FM"
  coord_cartesian(ylim = c(.6, .97)) +
  labs(x = NULL, y = "conformal coverage", title = "conformal coverage") +
  theme(axis.text.x = element_text(angle = 22, hjust = 1))

## ---- Fig 4: scATAC matched comparison (both FMs) ----
m <- fromJSON(file.path(base, "scatac_results", "scatac_fm_matched.json"))
af <- fromJSON(file.path(base, "scatac_results", "scatac_atacformer_audit_result.json"))
keys <- c("ChromFound FM (2048-d, faithful recipe)","raw top-2048 log-TFIDF (FM input, 2048-d)",
          "LSI top-2048 -> SVD-50","LSI full-peak -> SVD-50 (ref)")
short <- c("ChromFound","raw TF-IDF","LSI top-2048","LSI full-peak")
dm <- data.frame(sub = short, acc = sapply(keys, function(k) m[[k]]$acc), ECE = sapply(keys, function(k) m[[k]]$ECE),
                 kind = c("scATAC FM","baseline","baseline","baseline"))
dm <- rbind(dm, data.frame(sub="Atacformer", acc=af$cross_sample_probe_acc, ECE=af$ECE, kind="scATAC FM"))
dm$sub <- factor(dm$sub, levels = dm$sub)
fillk <- scale_fill_manual(values=c("scATAC FM"="#d95f0e","baseline"="#2c7fb8"), name=NULL)
p4a <- ggplot(dm, aes(sub, acc, fill=kind)) + geom_col(width=.7) +
  geom_text(aes(label=sprintf("%.3f",acc)), vjust=-.4, size=2.6) + fillk +
  coord_cartesian(ylim=c(0,1.05)) + labs(x=NULL, y="cell-type accuracy", title="accuracy") +
  theme(axis.text.x=element_text(angle=35, hjust=1, size=7))
p4e <- ggplot(dm, aes(sub, ECE, fill=kind)) + geom_col(width=.7) +
  geom_text(aes(label=sprintf("%.3f",ECE)), vjust=-.4, size=2.6) + fillk +
  coord_cartesian(ylim=c(0, 0.108)) +  # headroom so the tallest bar's label clears the panel-corner C tag
  labs(x=NULL, y="ECE", title="calibration error") +
  theme(axis.text.x=element_text(angle=35, hjust=1, size=7))
# Three equal-width panels (A|B|C), not A|(B|C) which gave A half the width and squeezed B/C.
# Guides collected to one bottom row: A carries its own scATAC/scRNA key, B and C share the
# baseline/scATAC-FM key (distinct colour families, so the two keys read unambiguously).
sc <- (f3 | p4a | p4e) + plot_layout(guides="collect") & theme(legend.position="bottom")
emit("fig_scatac", sc + ann("scATAC (cluster I): conformal coverage holds, and neither scATAC FM beats a simple LSI baseline"), TW, 3.9)

## ---- Fig 5: scRNA cross-batch LOBO ----
v1 <- data.frame(atlas = rep(c("lung (GEO)","bcells","cardiac-lymph","intestine"), each = 5),
  method = rep(c("PCA","scVI","HVG","scGPT (FM)","Geneformer (FM)"), times = 4),
  auroc = c(0.987,0.967,0.976,0.929,NA, 0.989,0.992,0.980,0.975,0.960,
            0.962,0.950,0.930,0.932,0.880, 0.972,0.982,0.963,0.968,0.957))
v1$atlas <- factor(v1$atlas, levels = c("lung (GEO)","bcells","cardiac-lymph","intestine"))
v1$kind <- ifelse(grepl("FM", v1$method), "foundation model", "classical baseline")
v1$method <- factor(v1$method, levels = c("PCA","scVI","HVG","scGPT (FM)","Geneformer (FM)"))
v1miss <- data.frame(atlas = factor("lung (GEO)", levels = levels(v1$atlas)),
                      method = factor("Geneformer (FM)", levels = levels(v1$method)))  # the NA cell (Geneformer failed on lung/GEO): label it explicitly, don't leave a silent gap
f5 <- ggplot(v1[!is.na(v1$auroc),], aes(method, auroc, fill = kind)) + geom_col(width = .7) +
  facet_wrap(~atlas, nrow = 1) + geom_text(aes(label = sprintf("%.2f", auroc)), vjust = -.4, size = 2.1) +
  geom_text(data = v1miss, aes(x = method), y = 0.855, label = "n/a", inherit.aes = FALSE, size = 2.1, color = "grey40") +
  scale_fill_manual(values = c("classical baseline"="#2c7fb8","foundation model"="#d95f0e"), name=NULL) +
  coord_cartesian(ylim = c(0.84, 1.0)) +
  labs(x = NULL, y = "LOBO macro-AUROC", title = "cross-batch LOBO") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
emit("fig5_scrna_lobo", f5, TW, 3.8)

## ---- Fig 6: SC18 coverage gap + V4 ----
sc18 <- data.frame(atlas = factor(c("lung","bcells","cardiac-lymph","intestine"),
                   levels = c("bcells","cardiac-lymph","intestine","lung")), gap = c(0.1233,0.0124,0.0196,0.0755))
p6a <- ggplot(sc18, aes(atlas, gap, fill = atlas)) + geom_col(width = .65) +
  geom_text(aes(label = sprintf("%.3f", gap)), vjust = -.4, size = 2.6) +
  scale_fill_brewer(palette = "Blues", guide = "none") + coord_cartesian(ylim = c(0, 0.14)) +
  labs(x = NULL, y = "cross-batch coverage gap", title = "coverage gap (scGPT)") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
# ECE read from deposited tables (not hardcoded): readout-capacity ECE from V4, temperature scaling from V4ts.
v4rc <- fromJSON(file.path(base, "..", "results", "V4", "v4_result.json"), simplifyVector = FALSE)$pooled_by_capacity
v4ts <- fromJSON(file.path(base, "..", "results", "V4ts", "sc_tempscale_result.json"), simplifyVector = FALSE)$pooled
capece <- setNames(sapply(v4rc, function(r) as.numeric(r$ece_lobo)), sapply(v4rc, function(r) r$cap))
v4 <- data.frame(cond = factor(c("linear","mlp-64","mlp-256","raw (T=1)","temp-scaled"),
                 levels = c("linear","mlp-64","mlp-256","raw (T=1)","temp-scaled")),
                 ece = c(capece[["C0_linear"]], capece[["C1_mlp64"]], capece[["C2_mlp256"]],
                         as.numeric(v4ts$ece_raw), as.numeric(v4ts$ece_tempscaled)),
                 grp = c("readout capacity","readout capacity","readout capacity","temperature","temperature"))
p6b <- ggplot(v4, aes(cond, ece, fill = grp)) + geom_col(width = .7) +
  geom_text(aes(label = sprintf("%.3f", ece)), vjust = -.4, size = 2.6) +
  scale_fill_manual(values = c("readout capacity"="#d95f0e","temperature"="#2c7fb8"), guide = "none") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") + coord_cartesian(ylim = c(0, 0.25)) +
  labs(x = NULL, y = "cross-batch (LOBO) ECE", title = "cross-batch calibration") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
emit("fig6_scrna_calib", (p6a | p6b) + ann("Self-computed scRNA calibration/reliability"), TW, 4.0)

## ---- Fig 7: 24-atlas coverage gap ----
mb <- fromJSON(file.path(base, "expand_results", "multiatlas_baseline.json"), simplifyVector = FALSE)
meth <- c("pca-logreg","hvg-logreg","knn","centroid","rf")
labm <- c("pca-logreg"="PCA+logreg","hvg-logreg"="HVG+logreg","knn"="$k$NN","centroid"="centroid","rf"="random-forest")
rows <- do.call(rbind, lapply(mb, function(a){ do.call(rbind, lapply(meth, function(mm){
  g <- a$methods[[mm]]$cov_gap; if (is.null(g)||is.na(g)) return(NULL)
  data.frame(method=labm[[mm]], cov_gap=as.numeric(g)) })) }))
rows$method <- factor(rows$method, levels=labm[meth]); mn <- aggregate(cov_gap~method, rows, mean)
f7 <- ggplot(rows, aes(method, cov_gap)) + geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  RP(geom_jitter(aes(color="per-atlas"), width=.13, height=0, alpha=.5, size=1.3)) +
  geom_point(data=mn, aes(color="mean (24 atlases)"), size=2.6) +
  geom_text(data=mn, aes(label=sprintf("%.2f", cov_gap)), vjust=-1.95, color="#d95f0e", size=2.8) +
  scale_color_manual(name=NULL, values=c("mean (24 atlases)"="#d95f0e","per-atlas"="#2c7fb8"),
                      breaks=c("mean (24 atlases)","per-atlas"),
                      guide=guide_legend(override.aes=list(size=c(2.6,1.3), alpha=c(1,.5)))) +
  labs(x=NULL, y="conformal coverage gap", title="Simple baselines, 24 leakage-controlled atlases") +
  theme(axis.text.x=element_text(angle=20, hjust=1), legend.position="top")
emit("fig7_multiatlas_covgap", f7, TW*0.86, 4.0)

## ---- Fig 8: two FMs vs best baseline, 20 atlases ----
fb <- fromJSON(file.path(base, "expand_results", "fm_vs_baseline_raw.json"), simplifyVector = FALSE)
gv <- function(at,rep,k){ for(r in at$reps) if(r$rep==rep){ v<-r[[k]]; return(if(is.null(v)) NA else as.numeric(v)) }; NA }
FMS <- c("Geneformer-V2(FM)"="Geneformer-V2", "scGPT(FM)"="scGPT")
df <- do.call(rbind, lapply(fb, function(at){
  ba <- suppressWarnings(max(gv(at,"PCA50","xb_auroc"), gv(at,"HVG2000","xb_auroc"), na.rm=TRUE))
  be <- suppressWarnings(min(gv(at,"PCA50","ece"), gv(at,"HVG2000","ece"), na.rm=TRUE))
  do.call(rbind, lapply(names(FMS), function(fk) data.frame(fm=FMS[[fk]],
    species=ifelse(at$atlas %in% MISMATCH, "vocab-mismatch (cross-species/ID)", "vocabulary-matched (human)"),
    base_auroc=ba, fm_auroc=gv(at,fk,"xb_auroc"), base_ece=be, fm_ece=gv(at,fk,"ece")))) }))
df <- df[is.finite(df$fm_auroc) & is.finite(df$base_auroc),]; df$fm <- factor(df$fm, levels=c("scGPT","Geneformer-V2"))
p8a <- ggplot(df, aes(base_auroc, fm_auroc, color=species)) + geom_abline(linetype="dashed", color="grey50") +
  RP(geom_point(size=1.5, alpha=.85)) + facet_wrap(~fm) + scale_color_manual(values=SPP, name=NULL) +
  coord_cartesian(xlim=c(0.45,1.02), ylim=c(0.45,1.02)) +
  labs(x="best simple baseline macro-AUROC", y="FM macro-AUROC", title="cross-batch discrimination", tag="A") +
  theme(legend.position="top", aspect.ratio=1, plot.tag.position=c(-0.06,1.05), plot.tag.location="panel", plot.tag=element_text(face="bold", size=11))
p8e <- ggplot(df, aes(base_ece, fm_ece, color=species)) + geom_abline(linetype="dashed", color="grey50") +
  RP(geom_point(size=1.5, alpha=.85)) + facet_wrap(~fm) + scale_color_manual(values=SPP, guide="none") +
  coord_cartesian(xlim=c(0,0.82), ylim=c(0,0.82)) +
  labs(x="best baseline ECE (lower = better)", y="FM ECE", title="calibration", tag="B") +
  theme(aspect.ratio=1, plot.tag.position=c(-0.06,1.05), plot.tag.location="panel", plot.tag=element_text(face="bold", size=11))
# Manual per-panel tags at a fixed npc position (NOT patchwork's auto tag_levels, which drifts
# with facet/legend layout); A/B now sit at each panel's top-left corner deterministically.
f8 <- (p8a / p8e) + plot_layout(guides="collect") & theme(legend.position="bottom")
emit("fig8_fm_vs_baseline", f8, TW, 5.9)

## ---- Fig 9: fair re-check ----
gx <- function(r,k){ v<-r[[k]]; if(is.null(v)||is.na(v)) NA else as.numeric(v) }
# Panel A: corrected multi-FM audit (scGPT + Geneformer-104M/316M + scFoundation; best-FM per atlas),
# coloured by the data-driven gene-naming split (11 vocabulary-matched vs 9 vocab-mismatch atlases).
ff <- fromJSON(file.path(base, "expand_results", "fm_fair_fig.json"), simplifyVector = FALSE)
fa <- do.call(rbind, lapply(ff, function(r) data.frame(
  species=ifelse(r$naming=="human-symbol","vocabulary-matched (human)","vocab-mismatch (cross-species/ID)"),
  knn_pca=gx(r,"knn_pca"), knn_fm=gx(r,"knn_fm") )))
p9a <- ggplot(fa, aes(knn_pca, knn_fm, color=species)) + geom_abline(linetype="dashed", color="grey50") +
  RP(geom_point(size=1.7, alpha=.85)) + scale_color_manual(values=SPP, name=NULL) + coord_equal(xlim=c(0.45,1.02), ylim=c(0.45,1.02)) +
  labs(x="best baseline $k$NN macro-AUROC", y="best FM $k$NN macro-AUROC", title="fair non-linear probe ($k$NN)") +
  guides(color=guide_legend(nrow=2)) +  # stack the two long keys so they fit the half-width panel (no left-edge clip)
  theme(legend.position="top")
# Panel B: label-vs-clustering circularity (species-agnostic; unchanged from the original fair re-check).
fr2 <- fromJSON(file.path(base, "expand_results", "fair_recheck.json"), simplifyVector = FALSE)
fd <- do.call(rbind, lapply(fr2, function(r) data.frame(
  er_true=gx(r,"exprR2_truelabel"), er_pca=gx(r,"exprR2_PCA50"),
  er_fm=suppressWarnings(max(gx(r,"exprR2_Geneformer-V2"),gx(r,"exprR2_scGPT"),na.rm=TRUE)) )))
fb2 <- rbind(data.frame(rep="PCA clustering", true=fd$er_true, val=fd$er_pca),
             data.frame(rep="FM clustering", true=fd$er_true, val=fd$er_fm))
p9b <- ggplot(fb2, aes(true, val, color=rep, shape=rep)) + geom_abline(linetype="dashed", color="grey50") +
  RP(geom_point(size=1.6, alpha=.8)) + scale_color_manual(values=c("PCA clustering"="#2c7fb8","FM clustering"="#7b3294"), name=NULL) +
  scale_shape_manual(values=c("PCA clustering"=16,"FM clustering"=17), name=NULL) +  # purple/blue hue alone is colorblind-confusable; distinct shapes fix it
  coord_equal(xlim=c(0,0.30), ylim=c(0,0.30)) +
  labs(x="$R^2_{\\mathrm{expr}}$ of the human labels", y="$R^2_{\\mathrm{expr}}$ of the clustering", title="reference-free structure") +
  theme(legend.position="top", axis.title.y=element_text(margin=margin(r=7)))  # wider gap than the r=4 default so the rotated y-title clears its tick numbers
emit("fig9_fair_recheck", (p9a | p9b) + ann("Fair re-check of the FM-vs-baseline gap (cluster J)"), TW, 4.0)

## ---- Fig 10: cluster-H fair re-check ----
sp <- fromJSON(file.path(base, "expand_results", "spatial_knn_probe.json"), simplifyVector = FALSE)
# ARI values read from the deposited result table (not hardcoded).
ariv <- unlist(fromJSON(file.path(base, "expand_results", "spatial_ari.json"), simplifyVector = TRUE)$ari)
sd2 <- do.call(rbind, lapply(sp, function(r){ if (is.null(r$knn_niche_AUROC)) return(NULL); mm <- r$method
  data.frame(method=mm, kind=ifelse(grepl("FM",mm),"spatial FM","simple baseline"),
             knn=as.numeric(r$knn_niche_AUROC), compo=as.numeric(r$compo_R2), expr=as.numeric(r$expr_R2),
             ari=ifelse(mm %in% names(ariv), ariv[[mm]], NA)) }))
lab2 <- c("nbhd-composition"="nbhd-composition","PCA-expression"="PCA (expression)","spatial-smoothed-PCA"="spatial-smoothed PCA",
          "novae (FM)"="Novae (FM)","scgpt_spatial (FM)"="scGPT-spatial (FM)","nicheformer (FM)"="Nicheformer (FM)")
sd2$ml <- ifelse(sd2$method %in% names(lab2), lab2[sd2$method], sd2$method); sd2$ml <- factor(sd2$ml, levels=sd2$ml[order(sd2$knn)])
mA <- rbind(data.frame(ml=sd2$ml, metric="unsupervised ARI (vs label)", val=sd2$ari),
            data.frame(ml=sd2$ml, metric="supervised $k$NN niche-AUROC", val=sd2$knn))
p10a <- ggplot(mA, aes(val, ml, color=metric, shape=metric)) + geom_line(aes(group=ml), color="grey70") +
  geom_point(size=2.6) + scale_color_manual(values=c("unsupervised ARI (vs label)"="#d95f0e","supervised $k$NN niche-AUROC"="#2c7fb8"), name=NULL) +
  scale_shape_manual(values=c(16,17), name=NULL) + labs(x="score", y=NULL, title="niche info present ($k$NN) vs ARI", tag="A") +
  theme(legend.position="top", plot.tag.position=c(-0.02,1.06), plot.tag.location="panel", plot.tag=element_text(face="bold", size=11))
sb <- rbind(data.frame(ml="TRUE niche label", kind="label", compo=0.234, expr=0.007), sd2[,c("ml","kind","compo","expr")])
p10b <- ggplot(sb, aes(compo, expr, color=kind)) + geom_point(size=2.2, alpha=.85) +
  ggrepel::geom_text_repel(aes(label=ml), size=2.4, max.overlaps=Inf, show.legend=FALSE,
                           seed=1, box.padding=0.55, point.padding=0.35, min.segment.length=0,
                           segment.color="grey60", segment.size=0.3, force=2) +  # deterministic de-collision in the cramped half-width panel
  scale_color_manual(values=c("label"="black","simple baseline"="#2c7fb8","spatial FM"="#7b3294"), name=NULL) +
  labs(x="$R^2_{\\mathrm{compo}}$ (neighborhood composition)", y="$R^2_{\\mathrm{expr}}$ (expression)", title="niches are composition-defined", tag="B") +
  theme(legend.position="top", plot.tag.position=c(-0.02,1.06), plot.tag.location="panel", plot.tag=element_text(face="bold", size=11))
# Manual per-panel tags (NOT patchwork's auto tag_levels, whose plot-location strands A/B
# far left of the long y-axis category labels); A/B now sit at each panel's top-left corner.
# free B's left axis so patchwork does not indent its (short, numeric) tick labels to
# match A's long categorical labels -- otherwise B's ticks sit far from the axis.
emit("fig10_spatial_fair", p10a / free(p10b, type="space", side="l"), TW, 6.3)


## ---- Fig 11: cluster L (perturbation) ----
# Values read from the deposited result tables (not hardcoded): full-gene MSE by combo-novelty
# stratum from V5b_FAIR, and the uncertainty-abstention MSEs from V5c_PARTIAL.
v5bj <- fromJSON(file.path(base, "..", "results", "V5b", "v5b_fair_result.json"), simplifyVector = FALSE)$by_novelty
slab <- c(combo_seen2="combo seen2", combo_seen1="combo seen1", combo_seen0="combo seen0")
v5b <- do.call(rbind, lapply(v5bj, function(r) data.frame(
  stratum=slab[[r$stratum]],
  method=c("GEARS (GNN)","additive prior","no-perturb"),
  mse=c(as.numeric(r$gears_fullMSE), as.numeric(r$additive_fullMSE), as.numeric(r$noperturb_fullMSE)))))
v5b$method <- factor(v5b$method, levels=c("GEARS (GNN)","additive prior","no-perturb"))
v5b$stratum <- factor(v5b$stratum, levels=c("combo seen2","combo seen1","combo seen0"))
pL1 <- ggplot(v5b, aes(stratum, mse, fill=method)) + geom_col(position="dodge", width=.72) +
  scale_fill_manual(values=c("GEARS (GNN)"="#d95f0e","additive prior"="#2c7fb8","no-perturb"="grey75"), name=NULL) +
  labs(x=NULL, y="full-gene MSE", title="full-gene MSE") + theme(legend.position="top")
v5cj <- fromJSON(file.path(base, "..", "results", "V5c", "v5c_partial_result.json"), simplifyVector = FALSE)
v5c <- data.frame(cond=factor(c("retain all","abstain top-30\\%","shuffle ctrl"), levels=c("retain all","abstain top-30\\%","shuffle ctrl")),
  mse=c(as.numeric(v5cj$mse_full), as.numeric(v5cj$mse_retain70), as.numeric(v5cj$mse_retain70_shuffled)))
pL2 <- ggplot(v5c, aes(cond, mse)) + geom_col(width=.6, fill="#2c7fb8") +
  geom_text(aes(label=sprintf("%.4f",mse)), vjust=-.5, size=2.5) + coord_cartesian(ylim=c(0.070,0.0810)) +
  labs(x=NULL, y="retained MSE", title="uncertainty abstention") +
  theme(axis.text.x=element_text(size=7.5))
emit("fig11_perturbation", (pL1 | pL2) + ann("Cluster L: perturbation (additive $\\geq$ GEARS)"), TW, 4.0)

## ---- Fig 12: 24-atlas reliability (abstention + ECE temp) ----
labm2 <- c("pca-logreg"="PCA+logreg","hvg-logreg"="HVG+logreg","knn"="$k$NN","centroid"="centroid","rf"="RF")
ml <- do.call(rbind, lapply(mb, function(a) do.call(rbind, lapply(names(a$methods), function(m){
  x<-a$methods[[m]]; if(is.null(x$acc_at80)||is.na(x$acc_at80)) return(NULL)
  data.frame(method=labm2[[m]], lift=as.numeric(x$acc_at80)-as.numeric(x$acc_full)) }))))
ml$method <- factor(ml$method, levels=labm2)
mn12 <- aggregate(lift~method, ml, mean)
pR1 <- ggplot(ml, aes(method, lift)) + geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_jitter(width=.13, height=0, alpha=.5, color="#2c7fb8", size=1.3) +
  geom_point(data=mn12, color="#d95f0e", size=2.6) +
  geom_label(data=mn12, aes(label=sprintf("%+.3f",lift)), vjust=-0.7, color="#d95f0e", size=2.5, fill="white", label.size=0, label.padding=unit(0.6,"pt")) +
  labs(x=NULL, y="accuracy gain @80\\% cov.", title="abstention recovers accuracy") +
  theme(axis.text.x=element_text(angle=20, hjust=1))
et <- do.call(rbind, lapply(mb, function(a){ x<-a$methods[["pca-logreg"]]
  if(is.null(x$ece_temp)||is.na(x$ece_temp)) return(NULL)
  data.frame(atlas=gsub("_"," ",a$atlas), raw=as.numeric(x$ece_raw), temp=as.numeric(x$ece_temp)) }))
pR2 <- ggplot(et) + geom_segment(aes(x=raw, xend=temp, y=reorder(atlas,raw), yend=reorder(atlas,raw)), color="grey60", arrow=arrow(length=unit(3,"pt"))) +
  geom_point(aes(raw, atlas), color="#d95f0e", size=1.4) + geom_point(aes(temp, atlas), color="#2c7fb8", size=1.4) +
  labs(x="cross-batch ECE: raw $\\to$ temp-scaled", y="24 atlases (worst ECE top)", title="temperature scaling reduces ECE") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())  # per-atlas codes intentionally hidden; ylab states the ordering so the blank axis is self-evident, not truncated
emit("fig12_scrna_reliability", (pR1 | pR2) + ann("Self-computed scRNA reliability across 24 atlases"), TW, 4.4)

## ---- Fig (cluster-K add-ons): FM-general capacity effect + abstention fairness ----
# Panel A: the readout-capacity -> out-of-batch ECE effect (scGPT, V4) reproduces on a second FM
# family (Geneformer-V1-10M, pooled LOBO over 3 atlases) -- the MLP-256 readout is worst for both.
gfr <- fromJSON(file.path(base, "..", "results", "GFREADOUT", "gf_readout_result.json"), simplifyVector = FALSE)$pooled
cm <- c(C0_linear="linear", C1_mlp64="MLP-64", C2_mlp256="MLP-256")
dcap <- rbind(
  data.frame(cap=c("linear","MLP-64","MLP-256"), ece=c(0.183,0.166,0.229), fm="scGPT (lung x-batch)"),
  data.frame(cap=cm[sapply(gfr, function(r) r$cap)], ece=sapply(gfr, function(r) as.numeric(r$ece_lobo)), fm="Geneformer (pooled LOBO)"))
dcap$cap <- factor(dcap$cap, levels=c("linear","MLP-64","MLP-256"))
pKa <- ggplot(dcap, aes(cap, ece, fill=fm)) + geom_col(position="dodge", width=.7) +
  geom_text(aes(label=sprintf("%.3f",ece)), position=position_dodge(.7), vjust=-.4, size=2.4) +
  scale_fill_manual(values=c("scGPT (lung x-batch)"="#2c7fb8","Geneformer (pooled LOBO)"="#d95f0e"), name=NULL) +
  coord_cartesian(ylim=c(0,0.27)) +
  labs(x="readout capacity", y="cross-batch ECE", title="MLP-256 readout worsens ECE on both FMs") +
  theme(legend.position="top", legend.text=element_text(size=7))
# Panel B: selective abstention on the lung scGPT cross-batch probe -- confidence beats random, and is
# fair to rare cell types (retained more, not over-rejected; the rejected rare cells are genuine errors).
ab <- fromJSON(file.path(base, "..", "results", "SC1_SC2", "abstention_subgroup_result.json"), simplifyVector = FALSE)
rc <- do.call(rbind, lapply(ab$risk_coverage_curve, function(p) data.frame(cov=as.numeric(p[[1]]), risk=as.numeric(p[[2]]))))
arnd <- as.numeric(ab$aurc_random_abstention); aconf <- as.numeric(ab$aurc_confidence)
pKb <- ggplot(rc, aes(cov, risk)) +
  geom_hline(yintercept=arnd, linetype="dashed", color="grey55") +
  annotate("text", x=0.04, y=arnd+0.009, label="random abstention", hjust=0, size=2.6, color="grey45") +
  geom_line(color="#1b7837", linewidth=.7) + geom_point(color="#1b7837", size=1.1) +
  annotate("text", x=0.06, y=0.045, label="confidence", hjust=0, size=2.6, color="#1b7837") +
  coord_cartesian(ylim=c(0,0.27)) +
  labs(x="coverage (fraction retained)", y="risk (error among retained)", title="selective abstention: effective, rare-type-fair")
emit("fig_clusterk_addon", (pKa | pKb) + ann("Cluster-K add-ons: FM-general capacity effect and abstention fairness"), TW, 3.8)

## ---- Fig 13: integration mixing vs purity ----
si <- fromJSON(file.path(base,"..","results","SCINT","sc_integration_result.json"))$rows
sv <- fromJSON(file.path(base,"..","results","SCVI","scvi_integration_result.json"))$rows
gi <- function(d) data.frame(atlas=d$atlas, rep=if(!is.null(d$rep)) d$rep else "scvi",
  iLISI=if(!is.null(d$iLISI_batch)) d$iLISI_batch else d$iLISI,
  cLISI=if(!is.null(d$cLISI_celltype)) d$cLISI_celltype else d$cLISI)
li <- rbind(gi(si), gi(sv))
li$method <- c("scgpt"="scGPT (FM)","geneformer"="Geneformer (FM)","pca50"="PCA","scvi"="scVI")[li$rep]
li$atlas <- gsub("_"," ", li$atlas); li <- li[!is.na(li$method),]
f13 <- ggplot(li, aes(iLISI, cLISI, color=method, shape=atlas)) +
  geom_point(size=2.8, alpha=.9, position=position_jitter(width=0.05, height=0.05, seed=42)) +  # intestine scGPT/scVI nearly coincide (1.689,1.235 vs 1.774,1.205); tiny fixed-seed jitter separates them
  scale_color_manual(values=c("scGPT (FM)"="#d95f0e","Geneformer (FM)"="#7b3294","PCA"="#2c7fb8","scVI"="#1b9e77"), name=NULL) +
  scale_shape_discrete(name=NULL) +
  labs(x="iLISI (batch mixing) $\\to$ higher better", y="cLISI (purity) $\\to$ lower better",
       title="mixing vs.\\ purity") +
  guides(color=guide_legend(nrow=1, order=1), shape=guide_legend(nrow=1, order=2)) +
  theme(legend.position="bottom", legend.box="horizontal", legend.margin=margin(t=-4, b=0))  # bottom strip frees the right-hand dead space and widens the scatter
emit("fig13_integration", f13, TW*0.82, 4.0)

## ---- Fig 14: vocabulary dose-response (depth probe #1) ----
# X = fraction of an atlas's genes mappable into each FM's gene vocabulary; Y = that FM's zero-shot
# kNN macro-AUROC. The "FM failures" on mouse/ENSG atlases are a token-mapping artifact: AUROC rises
# monotonically with coverage and plateaus at PCA-parity once the vocabulary is matched.
vd <- fromJSON(file.path(base, "expand_results", "vocab_dose_response.json"), simplifyVector = FALSE)$rows
dv <- do.call(rbind, lapply(vd, function(r) data.frame(
  fm=r$fm, frac=as.numeric(r$frac_gene), auroc=as.numeric(r$knn_auroc),
  match=ifelse(r$naming=="human-symbol","vocabulary-matched (human)","vocab-mismatch (cross-species/ID)"))))
spU <- 0.832  # universal-3 frac_gene Spearman (vocab_dose_response.json stats); annotated in caption
f14 <- ggplot(dv, aes(frac, auroc)) +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey55") +
  geom_hline(yintercept=0.894, linetype="dotted", color="grey40") +
  geom_smooth(method="lm", se=TRUE, color="grey30", fill="grey80", linewidth=.5) +
  geom_point(aes(color=match, shape=fm), size=2, alpha=.9) +
  annotate("text", x=0.06, y=0.455, label="chance", hjust=0, size=2.6, color="grey55") +
  annotate("text", x=0.02, y=0.965, label="PCA parity (0.894)", hjust=0, size=2.6, fontface="bold", color="grey20") +
  scale_color_manual(values=SPP, labels=c("vocab-mismatch (cross-species/ID)"="vocab-mismatch","vocabulary-matched (human)"="vocab-matched"), name=NULL) +
  scale_shape_manual(values=c(16,17,15,18,3,8), name="FM family", guide="none") +  # shapes still separate families visually; a 6-entry key collides with the two color legends at bottom, so no key (caption drops the shape clause)
  scale_x_continuous(limits=c(-0.02,1.02), labels=pct) + coord_cartesian(ylim=c(0.40,1.0)) +
  labs(x="fraction of atlas genes mapped (FM vocabulary)", y="zero-shot $k$NN macro-AUROC") +
  guides(color=guide_legend(override.aes=list(shape=16))) +
  theme(legend.position="bottom")

## ---- Fig 15: batch-shift dose-response (depth probe #2) ----
# X = how discriminable the held-out batch is from training (exchangeability-violation AUROC in PCA space);
# Y = cross-batch conformal coverage gap. Across 24 atlases AND 4 baseline methods the collapse magnitude
# is a graded function of batch shift (Spearman ~0.81) -- making the "general exchangeability failure"
# claim quantitative, not qualitative, and showing it is method-general (not FM- or classifier-specific).
# Modality contrast on ONE ruler: scRNA (steep) vs scATAC (flat), same axes. The 24-atlas x 4-classical
# breadth (rho=0.75) lives in the text + fig7; here we juxtapose the FM+classical scRNA curve against the
# scATAC curve (peak-LSI + 2 scATAC FMs) to show coverage collapse is a scRNA-cross-batch property.
GRPCOL <- c("classical (24 atlases, 4 methods)"="#2c7fb8","foundation model (high shift)"="#d95f0e")
XL <- c(0.1,1.0); YL <- c(-0.15,0.95)
# Panel A: scRNA -- the genuine full-RANGE dose-response is the 24-atlas x 4-classical-method sweep
# (spans shift 0.13-1.0, rho=0.82); the fit is on THAT. FM families (vocab-matched atlases) only reach
# the high-shift end, where they collapse by the same amount -- overlaid as confirmation, not their own fit.
bsd <- fromJSON(file.path(base, "expand_results", "batch_shift_dose_response.json"), simplifyVector = FALSE)$rows
mlab <- c("pca-logreg","knn","centroid","rf")
dcl <- do.call(rbind, lapply(bsd, function(r) do.call(rbind, lapply(mlab, function(m)
  data.frame(shift=as.numeric(r$batch_shift_pca), gap=as.numeric(r$cov_gap[[m]]),
             grp="classical (24 atlases, 4 methods)")))))
dcl <- dcl[is.finite(dcl$shift) & is.finite(dcl$gap), ]
dpca <- do.call(rbind, lapply(bsd, function(r) data.frame(shift=as.numeric(r$batch_shift_pca), gap=as.numeric(r$cov_gap[["pca-logreg"]]))))
dpca <- dpca[is.finite(dpca$shift) & is.finite(dpca$gap), ]  # headline curve = PCA+logreg (rho=0.82); other 3 methods shown as points only
bf <- fromJSON(file.path(base, "expand_results", "batch_shift_fm_probe.json"), simplifyVector = FALSE)$rows
dfm <- do.call(rbind, lapply(bf, function(r) if (isTRUE(r$is_fm) && is.finite(as.numeric(r$shift_auroc))) data.frame(
  shift=as.numeric(r$shift_auroc), gap=as.numeric(r$cov_gap), grp="foundation model (high shift)") else NULL))
dall <- rbind(dcl, dfm)
p15a <- ggplot(mapping=aes(shift, gap)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey55") +
  geom_smooth(data=dpca, method="lm", se=TRUE, color="grey30", fill="grey80", linewidth=.5) +
  geom_point(data=dall, aes(color=grp, shape=grp), size=1.6, alpha=.8) +
  scale_color_manual(values=GRPCOL, name=NULL) + scale_shape_manual(values=c(16,17), name=NULL) +
  coord_cartesian(xlim=XL, ylim=YL) +
  labs(x="batch-shift strength (held-out-batch AUROC)", y="cross-batch coverage gap",
       title="scRNA: collapse scales with shift ($\\rho{=}0.82$)") +
  theme(legend.position="top", legend.text=element_text(size=6.5))
# Panel B: scATAC -- same axes; coverage gap stays near zero across the whole shift range (flat, n.s.).
sba <- fromJSON(file.path(base, "expand_results", "scatac_batch_shift.json"), simplifyVector = FALSE)$rows
dsa <- do.call(rbind, lapply(sba, function(r) if (is.finite(as.numeric(r$shift_auroc))) data.frame(
  shift=as.numeric(r$shift_auroc), gap=as.numeric(r$cov_gap),
  grp=ifelse(isTRUE(r$is_fm),"foundation model","classical")) else NULL))
p15b <- ggplot(dsa, aes(shift, gap)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey55") +
  geom_smooth(aes(group=1), method="lm", se=TRUE, color="grey30", fill="grey80", linewidth=.5) +
  geom_point(aes(color=grp, shape=grp), size=1.8, alpha=.85) +
  annotate("text", x=0.55, y=0.60, label="flat: no cross-batch collapse", hjust=0.5, size=2.7, fontface="italic", color="grey45") +  # the shared y-range is empty here on purpose (scATAC does not collapse); label the void so it reads as signal, not a gap
  annotate("segment", x=0.55, xend=0.72, y=0.50, yend=0.06, linewidth=0.3, color="grey65", arrow=arrow(length=unit(3,"pt"))) +
  scale_color_manual(values=c("classical"="#2c7fb8","foundation model"="#d95f0e"), name=NULL) +
  scale_shape_manual(values=c("classical"=16,"foundation model"=17), name=NULL) +
  coord_cartesian(xlim=XL, ylim=YL) +
  labs(x="batch-shift strength (held-out-sample AUROC)", y="cross-batch coverage gap",
       title="scATAC: coverage holds ($\\rho{=}{-}0.24$, n.s.)") + theme(legend.position="top")
emit("fig15_batch_dose", (p15a | p15b) + ann("Batch-shift dose-response: scRNA vs scATAC"), TW, 4.2)

## ---- Fig 16: spatial dose-response (depth probe #3) ----
# Niche-ID is composition-defined, so per-cell expression models (the FMs) are handicapped on one axis:
# spatial context. Dosing it -- spatial-smoothing per-cell PCA (and one FM) over k neighbours -- lifts
# niche prediction monotonically past all 3 spatial FMs. Two supervised metrics; raw FMs + composition as refs.
sd16 <- fromJSON(file.path(base, "expand_results", "spatial_dose_response.json"), simplifyVector = FALSE)
dd <- do.call(rbind, lapply(sd16$dose, function(r) data.frame(
  curve=r$curve, k=as.numeric(r$k), AUROC=as.numeric(r$auroc), F1=as.numeric(r$f1))))
dd$cv <- ifelse(grepl("PCA", dd$curve), "smoothed PCA", "smoothed scGPT-spatial (FM)")
rf <- do.call(rbind, lapply(sd16$refs, function(r) data.frame(name=r$name, AUROC=as.numeric(r$auroc), F1=as.numeric(r$f1))))
rlab <- c("novae (raw FM)"="Novae","scgpt_spatial (raw FM)"="scGPT-sp.",
          "nicheformer (raw FM)"="Nicheformer","nbhd-composition baseline"="composition")
rf$lab <- rlab[rf$name]
# On the F1 panel, "nbhd-composition baseline" (0.5439) and "novae (raw FM)" (0.5257) sit only 0.018
# apart and their labels collide; nudge just those two labels apart (dotted reference lines stay put).
rf$F1_nudge <- 0
rf$F1_nudge[rf$name == "nbhd-composition baseline"] <- 0.028
rf$F1_nudge[rf$name == "novae (raw FM)"] <- -0.028
DCOL <- c("smoothed PCA"="#2c7fb8","smoothed scGPT-spatial (FM)"="#d95f0e")
mk16 <- function(ycol, ttl){
  rfl <- rf; rfl$ylab <- rf[[ycol]] + (if (ycol == "F1") rf$F1_nudge else 0)
  # Labels used to sit at x=185 (just left of the last point, k=200), right-aligned back over the
  # curve -- at this composite's shrunk row height the curve (which rises well above every
  # reference value by k~12-50) then reads as crossing straight through the label text. Anchoring
  # past the last real point (k=200 -> x=201) with hjust=0 and extending the x-range to give it
  # room puts every label in space no curve ever reaches, independent of panel height.
  ggplot(dd, aes(k+1, .data[[ycol]], color=cv)) +
    geom_hline(data=rf, aes(yintercept=.data[[ycol]]), linetype="dotted", color="grey55", linewidth=.35) +
    geom_label(data=rfl, aes(x=230, y=ylab, label=lab), inherit.aes=FALSE, hjust=0, vjust=0.5, size=1.7, color="grey40", fill="white", label.size=0, label.padding=unit(0.5,"pt")) +
    geom_line(linewidth=.6) + geom_point(size=1.7) +
    scale_color_manual(values=DCOL, name=NULL) +
    scale_x_continuous(trans="log", limits=c(0.9, 2000), breaks=c(1,4,7,13,26,51,101,201), labels=c(0,3,6,12,25,50,100,200)) +
    labs(x="spatial dose: $k$ neighbours smoothed", y=ycol, title=ttl) + theme(legend.position="top")
}
p16a <- mk16("AUROC","niche $k$NN-AUROC")
p16b <- mk16("F1","niche macro-F1 ($\\rho{=}0.98$, peaks at $k{\\approx}100$)")  # rho is the F1 Spearman -- label it on the F1 panel, not AUROC
emit("fig16_spatial_dose", (p16a | p16b) + ann("Spatial dose-response (niche-ID)"), TW, 4.4)

cat("ALL tikz figures written to", figs, "\n")

## ---- Fig: vocabulary ablation (probe #1 CAUSAL, within-atlas) ----
# Within a fixed matched atlas, keep only a fraction of the genes the FM tokenizer can read; AUROC is flat
# (parity) down to ~1-2k genes then cliffs to chance at 0 -> causal threshold, no cross-atlas biology confound.
va <- fromJSON(file.path(base, "expand_results", "vocab_ablation.json"), simplifyVector = FALSE)
dva <- do.call(rbind, lapply(va, function(r) data.frame(
  atlas=gsub("_"," ",r$atlas), genes=as.numeric(r$n_genes), auroc=as.numeric(r$knn_auroc), pca=as.numeric(r$pca_ref))))
dz <- dva[dva$genes>0, ]                       # log-x; the genes=0 -> 0.50 point annotated separately
pcaband <- range(dva$pca)
f17 <- ggplot(dz, aes(genes, auroc, color=atlas)) +
  annotate("rect", xmin=300, xmax=22000, ymin=pcaband[1], ymax=pcaband[2], alpha=.16, fill="grey70") +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey55") +
  geom_vline(xintercept=2048, linetype="dotted", color="grey55") +
  geom_line(linewidth=.6) + geom_point(size=1.8) +
  annotate("text", x=330, y=0.515, label="0 genes $\\to$ chance (0.50)", hjust=0, size=2.6, color="grey25") +
  annotate("text", x=2048, y=0.70, label="FM budget (2048)", hjust=0.5, size=2.4, color="grey35", angle=90) +  # shorter/darker; sits in the clear band between the "0 genes" annotation (below) and the plateau curves (above)
  annotate("text", x=330, y=pcaband[2]+0.014, label="PCA parity (all genes)", hjust=0, size=2.6, fontface="bold", color="grey20") +
  scale_x_log10(breaks=c(500,1000,2000,5000,10000,20000)) +
  scale_color_manual(values=c("GSE130148 lung"="#2c7fb8","lr gastric"="#d95f0e","lr stomach cancer"="#7b3294"),
                     labels=c("GSE130148 lung"="lung","lr gastric"="gastric","lr stomach cancer"="stomach cancer"), name=NULL) +
  coord_cartesian(ylim=c(0.45,1.0)) +
  labs(x="readable genes fed to the FM tokenizer (log scale)", y="zero-shot $k$NN macro-AUROC") +
  theme(legend.position="bottom")
fvoc <- (f14 | f17) + ann("The vocabulary artifact: dose-response (left) and causal within-atlas gene ablation (right)")
emit("fig_vocab", fvoc, TW, 3.9)

## ================= IEEE/JBHI compact variants (two-column, page-limited) =================
## Emitted as separate *_ieee.tex so the shared originals stay 1:1 for main.tex / main_els.tex.
## Full-width panels sit in figure* at the IEEE text width (FW); single-panel dose-response
## figures sit single-column (CW). Heights are trimmed hard to conserve the 14-page budget.
FW <- 7.16   # IEEE journal full text width (two columns + gutter), inches
CW <- 3.44   # IEEE journal single-column width, inches
# Panel-A legend keys were clipping ("...cross-specie") and bleeding under panel B's title at
# half-width; shorten the two keys (full meaning is in the caption) and keep them stacked on one
# top row apiece so neither overruns the panel.
p9a_ieee <- p9a + scale_color_manual(values=SPP, name=NULL,
  labels=c("vocabulary-matched (human)"="vocab-matched","vocab-mismatch (cross-species/ID)"="vocab-mismatch"))
# fig9_fair_recheck_ieee / fig10_spatial_fair_ieee / fig_vocab_ieee are superseded by the
# figfair_ieee / figspatial_ieee composites below (not \input by manuscript.tex any more);
# skipped here to save regen time -- their code is kept in git history if ever needed standalone.
# scATAC composite: push the A/B/C plot-tags out into the left/top margin so panel-A's
# "A" clears the rotated y-axis title "conformal coverage (target 0.90)" (was overprinting the ")").
# Height trimmed 2.5->2.3->2.0in (2026-07 condensation pass) for the 14pp target; single dense
# row, ample tag/title clearance at this margin so the cut is safe without a margin/tag rework.
emit("fig_scatac_ieee",         sc + ann("") &
       theme(plot.margin=margin(8,3,2,12), plot.tag.position=c(-0.05,1.06)), FW, 1.6)
# Height trimmed 2.5->2.15->1.85in for the 14pp target; simple 1x2 row, ample tag/title clearance.
emit("fig15_batch_dose_ieee",   (p15a | p15b) + ann(""), FW, 1.45)
# Height trimmed 2.3->2.0->1.75in for the 14pp target; simple 1x2 row.
emit("fig11_perturbation_ieee", (pL1 | pL2)  + ann(""), FW, 1.48)
cat("IEEE compact variants written\n")

## ================= IEEE/JBHI merged composite floats (page-limit condensation) =================
## Each composite fuses several original panels into ONE full-width figure* so results are retained
## while the main/supplement page budget shrinks. Base panel objects are reused (built above); tags
## are (re)assigned A.. via plot_annotation so every panel keeps a letter. Heights kept conservative.
CT <- function() plot_annotation(tag_levels = "A")
# MAIN fig:fair -- cluster J parity + label circularity + vocabulary artifact (2x2).
# Panels p9a/p9b are coord_equal (square); at full width the square would starve the row's
# height and collide the two centred titles, so free the aspect (coord_cartesian, same limits)
# for the composite only -- the shared originals keep their square coord.
# Shorten the (long) axis titles for the 2x2 composite only -- at half-width they were cramping the
# plot area; drop "macro"/"best"/"of the" and stack the A-panel legend on one row.
p9a_f <- p9a + coord_cartesian(xlim=c(0.45,1.02), ylim=c(0.45,1.02)) +
  labs(x="baseline $k$NN AUROC", y="FM $k$NN AUROC") +
  scale_color_manual(values=SPP, name=NULL,
    labels=c("vocabulary-matched (human)"="vocab-matched","vocab-mismatch (cross-species/ID)"="vocab-mismatch")) +
  guides(color=guide_legend(nrow=1))
p9b_f <- p9b + coord_cartesian(xlim=c(0,0.30), ylim=c(0,0.30)) +
  labs(x="label $R^2_{\\mathrm{expr}}$", y="clustering $R^2_{\\mathrm{expr}}$")
f14_f <- f14 + labs(y="zero-shot $k$NN AUROC")
f17_f <- f17 + labs(y="zero-shot $k$NN AUROC")
emit("figfair_ieee",     (p9a_f | p9b_f) / (f14_f | f17_f) + CT()
       & theme(plot.margin=margin(6,4,2,8), plot.tag.position=c(-0.03,1.05)), FW, 3.7)  # 2.85->3.7: 14pp squeeze made text too small after resizebox; restore legible panel/text ratio
# MAIN fig:spatial -- cluster H fair re-check (p10) + spatial dose-response (p16) (2x2)
# free() the bottom-left panel (C=p16a, numeric AUROC axis) from the left-space reserved for
# row-1's long categorical labels (p10a); otherwise C's plot sits far right of its y-axis title.
# Strip p10a/p10b's manual A/B tags + their "panel" tag-location so the composite CT() tags them like
# every other panel (global "plot" location, above the top legend). Panel-corner tags with a top legend
# overprint the legend text ("A" over "supervised...", "B" over "label").
# Panel-A's title + 2-item legend and panel-B's title were sized for the shared build's full-width
# stacked layout (fig10_spatial_fair); side-by-side at half-width here, the untouched versions overran
# their own cell and bled into the neighbouring panel's title/legend. Shorten both (full wording is in
# the caption) so each fits its ~3.5in half.
p10a_c <- p10a + labs(tag=NULL, title="niche info ($k$NN) vs ARI") +
  theme(plot.tag.location="plot", plot.tag.position=c(-0.03,1.05)) +
  scale_color_manual(values=c("unsupervised ARI (vs label)"="#d95f0e","supervised $k$NN niche-AUROC"="#2c7fb8"),
                      name=NULL, labels=c("unsupervised ARI (vs label)"="unsup. ARI","supervised $k$NN niche-AUROC"="sup. $k$NN AUROC")) +
  scale_shape_manual(values=c(16,17), name=NULL, labels=c("unsupervised ARI (vs label)"="unsup. ARI","supervised $k$NN niche-AUROC"="sup. $k$NN AUROC"))
p10b_c <- p10b + labs(tag=NULL, title="composition-defined niches") +
  theme(plot.tag.location="plot", plot.tag.position=c(-0.03,1.05))
# Top row (A=p10a_c, B=p10b_c): each panel's independent top legend is positioned without a
# reserved column boundary, so at half-width B's "label" key/text lands inside A's legend area and
# overpaints the tail of "unsup. ARI". guides="collect" reserves the two legends' combined width as
# one row instead of letting them float independently (they keep separate keys -- different scales
# -- just laid out together, mirroring the row16 fix below for the identical-scale case).
# guides="collect" alone drops the merged legend to patchwork's default position (right), silently
# discarding each panel's own legend.position="top" -- the resulting right-hand legend column steals
# enough panel width that the A/B titles (and, in row16, the C/D x-axis titles) run into each other.
# The explicit "& theme(legend.position=\"top\")" restores the intended top strip.
row10 <- (p10a_c | p10b_c) + plot_layout(guides="collect") & theme(legend.position="top")
# Bottom row (C=p16a, D=p16b) share the same 2-item cv color scale; local guides="collect" merges
# their identical legends into one instead of duplicating it under both panels.
# C is free()'d (its numeric y-axis needs far less left margin than row1's long ARI/AUROC axis
# title), which shrinks its own reserved left space -- reusing row1's -0.03 tag offset here would
# land "C" in what is now dead space to the left of the panel (see figScontext_ieee's identical
# bug/fix). Give C its own near-zero offset and D the shared -0.03, set directly on each panel
# (not the composite-wide "&" below) so the two rows can never silently clobber each other's tag.
p16a_c <- p16a + theme(plot.tag.position=c(0.01,1.06))
p16b_c <- p16b + theme(plot.tag.position=c(-0.03,1.05))
row16 <- (free(p16a_c, type="space", side="l") | p16b_c) + plot_layout(guides="collect") & theme(legend.position="top")
emit("figspatial_ieee",  row10 / row16 + CT()
       & theme(plot.margin=margin(6,4,2,8)), FW, 3.6)  # 2.7->3.6: compressed height starved ggrepel in panel B (label collisions) and shrank text; restore room
# SUPP figS:context -- meta-analysis (f1,f2) + LOBO (f5) + integration trade-off (f13).
# f1 needs vertical room (~9 study rows) and f5 needs full width (4 facet strips), so this is a
# 3-row stack, not a 2x2: (f1|f2) tall on top, f5 full-width, f13 full-width.
# free() the two full-width rows (C=f5, D=f13) from the left-space patchwork reserves for row-1's
# long categorical study labels (f1); otherwise their short numeric y-axes sit far right of the title.
# f1/f2's top scale expansion was sized for the shared standalone fig_meta render (headroom for the
# annotate() call-out labels above the topmost row); at this composite's shrunk row height that same
# absolute headroom reads as a large dead band between the title/legend and the first data row.
# Composite-only rebuilds (same data/aes, not layered on f1/f2 -- annotate() layers can't be
# swapped in-place) halve the top expansion and pull the labels' vjust in to match; f1 has
# clip="off" so the labels are never cut, they just need less empty panel above them to float in.
f1_c <- ggplot(ps, aes(win_rate, study_id, color = cluster, size = n)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.8869, linetype = "dotted", color = CB["A"]) +
  geom_vline(xintercept = 0.6533, linetype = "dotted", color = CB["B"]) +
  annotate("text", x = 0.5, y = topLev, label = "chance 0.5", vjust = -1.0, hjust = 0.5, size = 1.8, color = "grey45") +
  annotate("text", x = 0.8869, y = topLev, label = "cluster A mean", vjust = -1.0, hjust = 0.5, size = 1.8, color = CB[["A"]]) +
  annotate("text", x = 0.6533, y = topLev, label = "cluster B mean", vjust = -2.6, hjust = 0.5, size = 1.8, color = CB[["B"]]) +
  RP(geom_point(alpha = .85)) + scale_color_manual(values = CB) +
  scale_size_continuous(range = c(1.5, 5), trans = "sqrt", guide = "none") +
  scale_x_continuous(limits = c(-.02, 1.02), labels = pct) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.7))) +
  coord_cartesian(clip = "off") +
  labs(x = "baseline win-rate vs FM/DL (per study)", y = NULL, color = "cluster",
       title = "Per-study baseline win-rate") +
  theme(legend.position = "top", plot.margin = margin(14, 9, 5, 5),
        axis.text.y = element_text(size = 5.6), plot.tag.position = c(-0.03, 1.05))
f2_c <- ggplot(d2, aes(win_rate_studyEW, lab, color = grp)) +
  geom_vline(xintercept = .5, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .25, na.rm = TRUE) +
  geom_point(size = 2.4) + scale_color_manual(values = CB, guide = "none") +
  geom_text(aes(label = sprintf("%.2f ($k{=}%d$)", win_rate_studyEW, n_studies)), vjust = -0.6, size = 2.6) +
  scale_x_continuous(limits = c(0,1.18), labels = pct, breaks = c(0,.3,.6,.9)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.6))) +
  coord_cartesian(clip = "off") +  # w/o this the topmost row's value label gets cut by the panel frame at this composite's shrunk row height (matches f1_c's fix)
  labs(x = "win-rate (EW per study, 95\\% CI)", y = NULL, title = "Meta-analytic baseline win-rates") +
  theme(plot.tag.position = c(-0.03, 1.05))
# Panel-tag position is NOT set via a shared composite-wide "& theme(...)" here (that was the bug):
# free()'d rows (C=f5, D=f13) shrink their own reserved left margin down to what their numeric axis
# actually needs, so the SAME fractional x-offset tuned for row1's wide categorical-label margin
# (-0.03) lands the C/D tags in what is now dead space well to the left of the panel. Each panel
# instead carries its own plot.tag.position, set directly on that panel object before combining, so
# a tweak to one row's tag can never be silently overridden by a later blanket "&" call.
f5_c  <- free(f5  + theme(plot.tag.position = c(0.005, 1.06)), type = "space", side = "l")
f13_c <- free(f13 + theme(plot.tag.position = c(0.005, 1.06)), type = "space", side = "l")
emit("figScontext_ieee", (f1_c | f2_c) / f5_c / f13_c + plot_layout(heights=c(1.9,0.85,0.55)) + CT()
       & theme(plot.margin=margin(6,4,2,8)), FW, 4.1)  # row1's height share increased a lot (1.2->1.9 of the total) plus +0.25in total: the 12-row study-label overlap needed real room, not just reallocation
# SUPP figS:naive -- metric-artifact naive linear-probe view (fig8), laid side-by-side for full width.
# Height 3.5->2.7->2.3in: the panels are square (aspect.ratio=1, width-bound at ~1.5in), so the extra
# height was pure dead space that stranded the titles at top and the shared legend at bottom.
# p8a/p8e are aspect.ratio=1 (2 square facets each) and, at this composite's trimmed height, are
# HEIGHT-bound: the square panels only need ~40% of FW's width, so emitting at full FW left a huge
# empty gap between panel A and panel B (and outer margins) -- a canvas-vs-content mismatch, not a
# real design choice. Narrow the canvas to what the height-bound squares actually need instead of
# stretching a wide canvas around two small blocks (resizebox still scales the result to \textwidth
# in the document, so this only removes dead space, it doesn't change the final print size).
emit("figSnaive_ieee",   (p8a | p8e) + plot_layout(guides="collect") & theme(legend.position="bottom"), 4.5, 1.6)
# MAIN figS:reliabK -- cluster-K reliability, collapse + recovery, unified into ONE 7-panel figure.
# Row 1 (collapse): calibration | capacity/temp | 24-atlas reframe.  Rows 2-3 (recovery/robustness):
# abstention-accuracy | temp-ECE / capacity-ECE | abstention-fairness. Row-height RATIOS unchanged
# (2.5:2.2:2.2); total 6.9->5.9in for the 14pp target scales the whole figure down proportionally.
f7_s <- f7 + labs(title="24 leakage-controlled atlases")
emit("figSreliabK_ieee",
     (p6a | p6b | f7_s) / (pR1 | pR2) / (pKa | pKb) + plot_layout(heights=c(2.5,2.2,2.2)) + CT()
     & theme(plot.margin=margin(8,5,2,6), plot.tag.position=c(-0.03,1.06)), FW, 5.2)  # 3.9->5.2: 7-panel figure was over-squeezed; restore legible text-to-panel ratio
cat("IEEE merged composites written\n")
