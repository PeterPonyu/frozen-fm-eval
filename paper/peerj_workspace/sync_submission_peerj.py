#!/usr/bin/env python3
"""Refresh paper/submission_peerj/ from the live peerj_workspace.

- Export each TikZ panel to source/FigureN.pdf (appearance order)
- Convert manuscript.tex \\input{figs/...} -> \\includegraphics{FigureN.pdf}
- Copy manuscript.pdf + refresh CONTENTS.md
- Strip .omc / aux junk from the upload tree
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_PAPER = HERE.parent
SUB = REPO_PAPER / "submission_peerj"
SRC = SUB / "source"
FIGS = HERE / "figs"

# Appearance order in peerj_workspace/manuscript.tex (Figure1 = combined schematic)
FIG_MAP: list[tuple[str, str]] = [
    ("__schematic__", "Figure1.pdf"),  # fig0_overview + fig0b_roadmap
    ("fig9_fair_recheck.tex", "Figure2.pdf"),
    ("fig_meta.tex", "Figure3.pdf"),
    ("fig10_spatial_fair.tex", "Figure4.pdf"),
    ("fig16_spatial_dose.tex", "Figure5.pdf"),
    ("fig_scatac.tex", "Figure6.pdf"),
    ("fig15_batch_dose.tex", "Figure7.pdf"),
    ("fig13_integration.tex", "Figure8.pdf"),
    ("fig5_scrna_lobo.tex", "Figure9.pdf"),
    ("fig8_fm_vs_baseline.tex", "Figure10.pdf"),
    ("fig_vocab.tex", "Figure11.pdf"),
    ("fig6_scrna_calib.tex", "Figure12.pdf"),
    ("fig12_scrna_reliability.tex", "Figure13.pdf"),
    ("fig7_multiatlas_covgap.tex", "Figure14.pdf"),
    ("fig_clusterk_addon.tex", "Figure15.pdf"),
    ("fig11_perturbation.tex", "Figure16.pdf"),
]


PREAMBLE = r"""
\documentclass[border=2pt]{standalone}
\usepackage[T1]{fontenc}
\usepackage{times}
\usepackage{amsmath,amssymb}
\usepackage{xcolor}
\definecolor{cData}{HTML}{1b4965}\definecolor{cFM}{HTML}{d95f0e}\definecolor{cBase}{HTML}{2c7fb8}
\definecolor{cEval}{HTML}{6a51a3}\definecolor{cFind}{HTML}{1b7837}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepgfplotslibrary{groupplots}
\usetikzlibrary{positioning,fit,arrows.meta,backgrounds,calc,shapes.geometric,shadows}
\begin{document}
"""


def run(cmd, cwd=None, env=None) -> None:
    print("+", " ".join(cmd))
    r = subprocess.run(cmd, cwd=cwd or HERE, env=env, text=True)
    if r.returncode != 0:
        raise SystemExit(r.returncode)


def compile_workspace() -> str:
    run(["latexmk", "-pdf", "-interaction=nonstopmode", "manuscript.tex"])
    info = subprocess.run(
        ["pdfinfo", str(HERE / "manuscript.pdf")], capture_output=True, text=True
    ).stdout
    for line in info.splitlines():
        if line.startswith("Pages"):
            return line.split()[-1]
    return "?"


def export_figures() -> None:
    build = HERE / "_figbuild"
    build.mkdir(exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)

    for src, out_name in FIG_MAP:
        out = SRC / out_name
        if src == "__schematic__":
            tex = build / "Figure1.tex"
            tex.write_text(
                PREAMBLE
                + r"\begin{minipage}{16cm}"
                + "\n"
                + r"\input{../figs/fig0_overview.tex}\\[6pt]"
                + "\n"
                + r"\input{../figs/fig0b_roadmap.tex}"
                + "\n"
                + r"\end{minipage}"
                + "\n\\end{document}\n"
            )
        else:
            n = out_name.replace("Figure", "").replace(".pdf", "")
            tex = build / f"Figure{n}.tex"
            tex.write_text(
                PREAMBLE + f"\\input{{../figs/{src}}}\n\\end{{document}}\n"
            )
        log = build / (tex.stem + ".log")
        r = subprocess.run(
            [
                "pdflatex",
                "-interaction=nonstopmode",
                "-halt-on-error",
                tex.name,
            ],
            cwd=build,
            capture_output=True,
            text=True,
        )
        if r.returncode != 0:
            print(r.stdout[-2000:])
            raise SystemExit(f"FAILED {out_name}")
        pdf = build / (tex.stem + ".pdf")
        subprocess.run(
            ["pdfcrop", "--margins", "2", str(pdf), str(out)],
            check=True,
            capture_output=True,
        )
        print(f"  {out_name} <- {src}")


def convert_manuscript() -> None:
    s = (HERE / "manuscript.tex").read_text()
    # Replace the two-panel schematic body (titles + both TikZ panels) with Figure1.pdf
    s = re.sub(
        r"\\hbox to \\linewidth\{.*?Study architecture.*?\\hfil\}\\vspace\{2pt\}\s*"
        r"\\resizebox\{\\linewidth\}\{!\}\{\\input\{figs/fig0_overview\.tex\}\}\\\\\[6pt\]\s*"
        r"\\hbox to \\linewidth\{.*?Argument roadmap.*?\\hfil\}\\vspace\{2pt\}\s*"
        r"\\resizebox\{\\linewidth\}\{!\}\{\\input\{figs/fig0b_roadmap\.tex\}\}",
        r"\\includegraphics[width=\\linewidth]{Figure1.pdf}",
        s,
        count=1,
        flags=re.S,
    )

    for src, out_name in FIG_MAP:
        if src == "__schematic__":
            continue
        s = s.replace(
            f"\\fitfig{{\\input{{figs/{src}}}}}",
            f"\\fitfig{{\\includegraphics{{{out_name}}}}}",
        )
        s = s.replace(f"\\input{{figs/{src}}}", f"\\includegraphics{{{out_name}}}")

    if "figs/fig0_" in s or "figs/fig0b_" in s:
        raise SystemExit("schematic TikZ inputs survived conversion")
    if re.search(r"\\input\{figs/", s):
        left = re.findall(r"\\input\{figs/[^}]+\}", s)
        raise SystemExit(f"TikZ inputs survived conversion: {left[:5]}")

    s = re.sub(r"\\graphicspath\{\{figs/\}\}\s*", "", s)
    (SRC / "manuscript.tex").write_text(s)


def refresh_bundle(pages: str) -> None:
    shutil.copy2(HERE / "manuscript.pdf", SUB / "manuscript.pdf")
    # keep references.bib / wlpeerj.cls for local test; strip .omc
    omc = SUB / ".omc"
    if omc.exists():
        shutil.rmtree(omc)
    omc_src = SRC / ".omc"
    if omc_src.exists():
        shutil.rmtree(omc_src)
    for p in SRC.glob("manuscript.*"):
        if p.suffix in {".aux", ".log", ".out", ".fls", ".fdb_latexmk"}:
            p.unlink()

    contents = SUB / "CONTENTS.md"
    if contents.exists():
        text = contents.read_text()
        text = re.sub(
            r"(?m)^(\*\*?PeerJ.*)$",
            r"\1",
            text,
        )
        # stamp rebuild note at end
        stamp = (
            f"\n\n---\nRefreshed from `peerj_workspace/` "
            f"({pages} pp workspace PDF; 16 Figure PDFs re-exported).\n"
        )
        if "Refreshed from" in text:
            text = re.sub(r"\n---\nRefreshed from.*", stamp.rstrip() + "\n", text, flags=re.S)
        else:
            text = text.rstrip() + stamp
        contents.write_text(text)


def main() -> int:
    pages = compile_workspace()
    print("workspace pages", pages)
    export_figures()
    convert_manuscript()
    refresh_bundle(pages)
    # local smoke: compile source with flat figures
    env = dict(**{k: v for k, v in __import__("os").environ.items()})
    env["TEXINPUTS"] = str(HERE) + "//:"
    r = subprocess.run(
        [
            "latexmk",
            "-pdf",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "manuscript.tex",
        ],
        cwd=SRC,
        capture_output=True,
        text=True,
        env=env,
    )
    if r.returncode != 0:
        print(r.stdout[-2500:])
        print("WARN: source smoke compile failed; workspace PDF still copied")
    else:
        log = (SRC / "manuscript.log").read_text(errors="ignore")
        m = re.search(r"Output written on manuscript\.pdf \((\d+) pages", log)
        print("source smoke pages", m.group(1) if m else "?")
        shutil.copy2(SRC / "manuscript.pdf", SUB / "manuscript.pdf")
        for ext in (".aux", ".log", ".out", ".fls", ".fdb_latexmk"):
            p = SRC / f"manuscript{ext}"
            if p.exists():
                p.unlink()
    print("submission_peerj refreshed:", SUB)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
