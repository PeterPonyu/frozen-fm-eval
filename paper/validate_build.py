#!/usr/bin/env python3
"""
Build validation for sc-fm-benchmark paper.

Checks:
- LaTeX compilation status (errors, warnings, overfull boxes)
- Font embedding (via pdffonts)
- Page counts vs. expected
- Undefined references
- Figure file existence
- Stale/generated output markers

Usage:
    python validate_build.py [--workspace main|peerj|ieee|all]
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Expected page counts (update as paper evolves)
EXPECTED_PAGES = {
    "main": 39,
    "peerj": 33,
    "ieee_main": 15,
    "ieee_supp": 1,
}

# Build configurations
BUILD_CONFIGS = {
    "main": {
        "dir": Path("."),
        "pdf": "main.pdf",
        "tex": "main.tex",
        "engine": "lualatex",
        "expected_fonts": ["LMRoman", "CMMI", "CMR", "CMSY"],  # Latin Modern text + CM math
    },
    "peerj": {
        "dir": Path("peerj_workspace"),
        "pdf": "manuscript.pdf",
        "tex": "manuscript.tex",
        "engine": "pdflatex",
        "expected_fonts": ["NimbusRomNo9L", "NimbusSanL"],  # Times/Helvetica clones
    },
    "ieee": {
        "dir": Path("ieee_workspace"),
        "pdf": "manuscript.pdf",
        "tex": "manuscript.tex",
        "engine": "pdflatex",
        "expected_fonts": ["CMR", "CMMI", "CMSY"],  # Computer Modern
    },
}


class ValidationResult:
    def __init__(self):
        self.errors = []
        self.warnings = []
        self.info = []

    def add_error(self, msg: str):
        self.errors.append(f"ERROR: {msg}")

    def add_warning(self, msg: str):
        self.warnings.append(f"WARNING: {msg}")

    def add_info(self, msg: str):
        self.info.append(f"INFO: {msg}")

    def is_ok(self) -> bool:
        return len(self.errors) == 0

    def report(self) -> str:
        lines = []
        if self.errors:
            lines.append("ERRORS:")
            lines.extend(f"  {e}" for e in self.errors)
        if self.warnings:
            lines.append("WARNINGS:")
            lines.extend(f"  {w}" for w in self.warnings)
        if self.info:
            lines.append("INFO:")
            lines.extend(f"  {i}" for i in self.info)
        if not lines:
            lines.append("All checks passed ✓")
        return "\n".join(lines)


def check_pdf_exists(pdf_path: Path, result: ValidationResult):
    """Check if PDF exists."""
    if not pdf_path.exists():
        result.add_error(f"PDF not found: {pdf_path}")
        return False
    result.add_info(f"PDF exists: {pdf_path}")
    return True


def check_page_count(pdf_path: Path, expected: int, result: ValidationResult):
    """Check page count via pdfinfo."""
    try:
        output = subprocess.check_output(
            ["pdfinfo", str(pdf_path)],
            stderr=subprocess.STDOUT,
            text=True
        )
        match = re.search(r"Pages:\s+(\d+)", output)
        if match:
            actual = int(match.group(1))
            if actual == expected:
                result.add_info(f"Page count OK: {actual} pages")
            else:
                result.add_warning(
                    f"Page count mismatch: expected {expected}, got {actual}"
                )
        else:
            result.add_error("Could not parse page count from pdfinfo")
    except subprocess.CalledProcessError as e:
        result.add_error(f"pdfinfo failed: {e.output}")
    except FileNotFoundError:
        result.add_warning("pdfinfo not available (install poppler-utils)")


def check_fonts_embedded(pdf_path: Path, expected_fonts: List[str], result: ValidationResult):
    """Check font embedding via pdffonts."""
    try:
        output = subprocess.check_output(
            ["pdffonts", str(pdf_path)],
            stderr=subprocess.STDOUT,
            text=True
        )
        lines = output.strip().split("\n")

        # Skip the header separator. The type field may contain spaces (for
        # example ``CID Type 0C``), so parse the stable trailing columns by
        # regex instead of indexing split tokens.
        if len(lines) < 3:
            result.add_warning("pdffonts returned no font rows")
            return
        font_lines = [l for l in lines[2:] if l.strip()]
        row_pattern = re.compile(
            r"^(\S+).*\s+(yes|no)\s+(yes|no)\s+(yes|no)\s+\d+\s+\d+\s*$"
        )

        unembedded = []
        malformed = 0
        for line in font_lines:
            match = row_pattern.match(line)
            if not match:
                malformed += 1
                continue
            name, emb, _sub, _uni = match.groups()
            if emb.lower() != "yes":
                unembedded.append(name)

        if malformed:
            result.add_error(f"Malformed pdffonts rows: {malformed}")
        if unembedded:
            result.add_error(f"Unembedded fonts: {', '.join(unembedded)}")
        else:
            result.add_info("All fonts embedded ✓")

        # Check for expected font families
        found_fonts = any(
            any(exp in line for exp in expected_fonts)
            for line in font_lines
        )
        if found_fonts:
            result.add_info(f"Expected font families present")
        else:
            result.add_warning(
                f"Expected fonts not clearly present: {expected_fonts}"
            )

    except subprocess.CalledProcessError as e:
        result.add_error(f"pdffonts failed: {e.output}")
    except FileNotFoundError:
        result.add_warning("pdffonts not available (install poppler-utils)")


def check_log_for_issues(log_path: Path, result: ValidationResult):
    """Check LaTeX log for errors, undefined refs, overfull boxes."""
    if not log_path.exists():
        result.add_warning(f"Log file not found: {log_path}")
        return

    content = log_path.read_text(errors="ignore")

    # Check for errors
    error_pattern = r"^! .*$"
    errors = re.findall(error_pattern, content, re.MULTILINE)
    if errors:
        result.add_error(f"LaTeX errors found: {len(errors)} errors")
        for err in errors[:3]:  # Show first 3
            result.add_error(f"  {err}")

    # Check for undefined references
    undef_pattern = r"LaTeX Warning: Reference .* undefined"
    undef_refs = re.findall(undef_pattern, content)
    if undef_refs:
        result.add_error(f"Undefined references: {len(undef_refs)}")

    # Check for overfull hboxes > 20pt
    overfull_pattern = r"Overfull \\hbox \((\d+\.\d+)pt too wide\)"
    overfull = [
        float(m) for m in re.findall(overfull_pattern, content)
        if float(m) > 20.0
    ]
    if overfull:
        result.add_warning(
            f"Overfull hboxes >20pt: {len(overfull)} instances, "
            f"max {max(overfull):.1f}pt"
        )

    # Check for missing citations
    citation_pattern = r"LaTeX Warning: Citation .* undefined"
    missing_cites = re.findall(citation_pattern, content)
    if missing_cites:
        result.add_error(f"Missing citations: {len(missing_cites)}")


def check_figures_exist(workspace_dir: Path, tex_file: str, result: ValidationResult):
    """Check that all referenced figures exist."""
    tex_path = workspace_dir / tex_file
    if not tex_path.exists():
        result.add_warning(f"TeX file not found: {tex_path}")
        return

    content = tex_path.read_text()

    # Find \input{figs/...} or \includegraphics{figs/...}
    input_pattern = r"\\(?:input|includegraphics)(?:\[[^\]]*\])?\{(figs/[^}]+)\}"
    refs = re.findall(input_pattern, content)

    missing = []
    for ref in refs:
        # Handle both with and without .tex extension
        fig_path = workspace_dir / ref
        if not fig_path.exists():
            # Try with .tex extension if not specified
            if not ref.endswith(".tex"):
                fig_path = workspace_dir / f"{ref}.tex"
        if not fig_path.exists():
            missing.append(ref)

    if missing:
        result.add_error(f"Missing figure files: {len(missing)}")
        for m in missing[:5]:  # Show first 5
            result.add_error(f"  {m}")
    else:
        result.add_info(f"All {len(refs)} referenced figures exist")


def check_stale_outputs(base_dir: Path, result: ValidationResult):
    """Check for stale outputs that should be in archive/."""
    # Common stale patterns
    stale_patterns = [
        "*_els.tex",  # Elsevier versions (now in archive)
        "*_els.pdf",
        "submission_default/",  # Replaced by venue-specific submissions
        "backup/*.tex",  # Should be in git history
        "review/*.png",  # Stale raster crops
    ]

    found_stale = []
    for pattern in stale_patterns:
        if "/" in pattern:
            # Directory check
            stale_dir = base_dir / pattern.rstrip("/")
            if stale_dir.exists() and stale_dir.name != "archive":
                found_stale.append(str(stale_dir.relative_to(base_dir)))
        else:
            # File pattern check
            matches = list(base_dir.glob(pattern))
            # Exclude archive
            matches = [m for m in matches if "archive" not in m.parts]
            found_stale.extend(str(m.relative_to(base_dir)) for m in matches)

    if found_stale:
        result.add_warning(
            f"Stale outputs found (consider moving to archive/): "
            f"{', '.join(found_stale[:5])}"
        )


def validate_workspace(name: str, config: Dict, base_dir: Path) -> ValidationResult:
    """Validate a single workspace."""
    result = ValidationResult()
    result.add_info(f"Validating {name} workspace")

    workspace_dir = base_dir / config["dir"]
    pdf_path = workspace_dir / config["pdf"]
    log_path = workspace_dir / f"{config['tex'].replace('.tex', '.log')}"

    # Core checks
    if check_pdf_exists(pdf_path, result):
        # Page count
        expected_key = name if name != "ieee" else "ieee_main"
        if expected_key in EXPECTED_PAGES:
            check_page_count(pdf_path, EXPECTED_PAGES[expected_key], result)

        # Font embedding
        check_fonts_embedded(pdf_path, config["expected_fonts"], result)

    # Log checks
    check_log_for_issues(log_path, result)

    # Figure existence
    check_figures_exist(workspace_dir, config["tex"], result)

    return result


def main():
    parser = argparse.ArgumentParser(description="Validate paper builds")
    parser.add_argument(
        "--workspace",
        choices=["main", "peerj", "ieee", "all"],
        default="all",
        help="Which workspace to validate"
    )
    parser.add_argument(
        "--check-stale",
        action="store_true",
        help="Check for stale outputs"
    )
    args = parser.parse_args()

    base_dir = Path(__file__).parent

    workspaces = (
        BUILD_CONFIGS.keys() if args.workspace == "all"
        else [args.workspace]
    )

    all_ok = True
    for ws in workspaces:
        result = validate_workspace(ws, BUILD_CONFIGS[ws], base_dir)
        print(f"\n{'='*60}")
        print(f"{ws.upper()} WORKSPACE")
        print('='*60)
        print(result.report())
        all_ok = all_ok and result.is_ok()

    if args.check_stale:
        print(f"\n{'='*60}")
        print("STALE OUTPUT CHECK")
        print('='*60)
        stale_result = ValidationResult()
        check_stale_outputs(base_dir, stale_result)
        print(stale_result.report())

    print(f"\n{'='*60}")
    if all_ok:
        print("✓ All workspaces validated successfully")
        return 0
    else:
        print("✗ Validation failed (see errors above)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
