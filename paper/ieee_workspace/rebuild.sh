#!/usr/bin/env bash
# Regenerate IEEE figures + compile manuscript. Usage: ./rebuild.sh [--figs-only|--tex-only]
set -e
cd "$(dirname "$0")"
REPO=/home/zeyufu/Desktop/singlecell-genomics-research
export SCFM_BASE="$REPO/research/sc-fm-benchmark"
MODE="${1:-all}"
if [ "$MODE" != "--tex-only" ]; then
  echo "=== regen figures ==="
  Rscript make_figs_ieee.R
fi
if [ "$MODE" != "--figs-only" ]; then
  echo "=== compile ==="
  latexmk -pdf -interaction=nonstopmode manuscript.tex >/dev/null 2>&1 || true
  latexmk -pdf -interaction=nonstopmode manuscript.tex 2>&1 | tail -3
  echo -n "PAGES: "; pdfinfo manuscript.pdf | awk '/Pages/{print $2}'
fi
