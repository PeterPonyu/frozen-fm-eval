#!/usr/bin/env Rscript
# Layout validation helpers for sc-fm-benchmark figures
#
# Checks for common layout issues:
# - Y-axis label extent violations
# - Panel size inconsistencies
# - Multiple legend owners in composites
# - Title/annotation collisions
#
# Usage:
#   Rscript check_figure_layout.R [figure_name.tex]
#
# If no argument, checks all .tex files in figs/

library(tools)

# Configuration
MAX_YTICK_LABEL_WIDTH_PT <- 40  # pt
MAX_PANEL_HEIGHT_VARIATION <- 0.15  # 15% tolerance for panel heights
KNOWN_COMPOSITE_LEGEND_OWNERS <- list(
  "fig8_integration" = 1,  # Single collected legend
  "figfair_ieee" = 1,
  "figspatial_ieee" = 1
)

check_ytick_labels <- function(tex_content, filename) {
  # Check for long y-axis tick labels that might overflow
  issues <- c()

  # Pattern: \node[...,anchor=east] at (axis cs:...) {long text};
  ytick_pattern <- "\\\\node\\[.*anchor=east.*\\].*at.*\\{([^}]+)\\}"
  matches <- gregexpr(ytick_pattern, tex_content, perl = TRUE)

  if (matches[[1]][1] != -1) {
    labels <- regmatches(tex_content, matches)[[1]]

    for (label in labels) {
      # Extract text content
      text_match <- regexpr("\\{([^}]+)\\}$", label, perl = TRUE)
      if (text_match != -1) {
        text <- sub(".*\\{([^}]+)\\}$", "\\1", label)

        # Rough heuristic: >30 chars or contains long words
        if (nchar(text) > 30 || any(nchar(strsplit(text, " ")[[1]]) > 15)) {
          issues <- c(issues, sprintf(
            "  WARNING: Long y-tick label may overflow: \"%s\"",
            substr(text, 1, 50)
          ))
        }
      }
    }
  }

  issues
}

check_panel_sizing <- function(tex_content, filename) {
  # Check for inconsistent panel heights in composites
  issues <- c()

  # Pattern: \begin{subfigure}[t]{width}
  # Often followed by \includegraphics or tikzpicture with height
  height_pattern <- "height\\s*=\\s*([0-9.]+)(in|cm|pt)"
  matches <- gregexpr(height_pattern, tex_content, perl = TRUE)

  if (matches[[1]][1] != -1) {
    heights_raw <- regmatches(tex_content, matches)[[1]]

    # Extract numeric values and units
    heights <- lapply(heights_raw, function(h) {
      num <- as.numeric(sub("height\\s*=\\s*([0-9.]+).*", "\\1", h))
      unit <- sub(".*height\\s*=\\s*[0-9.]+([a-z]+).*", "\\1", h, perl = TRUE)
      list(value = num, unit = unit, raw = h)
    })

    if (length(heights) > 1) {
      # Convert all to inches for comparison
      to_inches <- function(h) {
        if (h$unit == "cm") return(h$value / 2.54)
        if (h$unit == "pt") return(h$value / 72.27)
        return(h$value)
      }

      heights_in <- sapply(heights, to_inches)

      if (max(heights_in) - min(heights_in) > max(heights_in) * MAX_PANEL_HEIGHT_VARIATION) {
        issues <- c(issues, sprintf(
          "  WARNING: Panel height variation >%.0f%%: %.2f-%.2fin",
          MAX_PANEL_HEIGHT_VARIATION * 100,
          min(heights_in),
          max(heights_in)
        ))
      }
    }
  }

  issues
}

check_legend_ownership <- function(tex_content, filename) {
  # Check for multiple legends in composites (potential confusion)
  issues <- c()

  count_matches <- function(pattern, text, fixed = TRUE) {
    matches <- gregexpr(pattern, text, fixed = fixed)[[1]]
    if (length(matches) == 1 && matches[[1]] == -1) return(0L)
    length(matches)
  }

  # Count legend directives. A missing pattern must count as zero, not one.
  legend_patterns <- c(
    "legend pos=",
    "legend style=",
    "\\\\legend\\{",
    "legend entries="
  )
  legend_count <- sum(vapply(
    legend_patterns,
    function(pattern) count_matches(pattern, tex_content),
    integer(1)
  ))
  axis_count <- count_matches("\\\\begin\\{axis\\}", tex_content)

  if (legend_count > axis_count && legend_count > 2) {
    # Check if this is a known composite with collected legend
    basename <- sub("\\.tex$", "", basename(filename))
    expected <- KNOWN_COMPOSITE_LEGEND_OWNERS[[basename]]

    if (is.null(expected) || legend_count > expected) {
      issues <- c(issues, sprintf(
        "  INFO: Multiple legend directives (%d) across %d axes - verify only one legend is visible",
        legend_count, axis_count
      ))
    }
  }

  issues
}

check_annotation_collisions <- function(tex_content, filename) {
  # Check for overlapping annotations (panel tags, text boxes)
  issues <- c()

  # Count panel tags
  tag_pattern <- "\\\\node\\[.*panel.?tag.*\\]"
  tag_matches <- gregexpr(tag_pattern, tex_content, perl = TRUE)
  tag_count <- ifelse(tag_matches[[1]][1] == -1, 0, length(tag_matches[[1]]))

  # Count subfigures
  subfig_pattern <- "\\\\begin\\{subfigure\\}"
  subfig_matches <- gregexpr(subfig_pattern, tex_content, fixed = TRUE)
  subfig_count <- ifelse(subfig_matches[[1]][1] == -1, 0, length(subfig_matches[[1]]))

  if (tag_count > 0 && tag_count != subfig_count && subfig_count > 0) {
    issues <- c(issues, sprintf(
      "  WARNING: Panel tag count (%d) != subfigure count (%d)",
      tag_count, subfig_count
    ))
  }

  issues
}

check_figure <- function(filepath) {
  cat(sprintf("\n%s\n", filepath))
  cat(strrep("=", nchar(filepath)), "\n")

  if (!file.exists(filepath)) {
    cat("  ERROR: File not found\n")
    return(FALSE)
  }

  content <- paste(readLines(filepath, warn = FALSE), collapse = "\n")

  all_issues <- c(
    check_ytick_labels(content, filepath),
    check_panel_sizing(content, filepath),
    check_legend_ownership(content, filepath),
    check_annotation_collisions(content, filepath)
  )

  if (length(all_issues) == 0) {
    cat("  ✓ No layout issues detected\n")
    return(TRUE)
  } else {
    cat(paste(all_issues, collapse = "\n"), "\n")
    return(FALSE)
  }
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  file_arg <- grep("^--file=", commandArgs(), value = TRUE)
  script_dir <- if (length(file_arg) == 1) {
    normalizePath(dirname(sub("^--file=", "", file_arg)), mustWork = TRUE)
  } else {
    getwd()
  }

  resolve_figure_path <- function(path) {
    if (file.exists(path)) return(path)
    candidate <- file.path(script_dir, path)
    if (file.exists(candidate)) return(candidate)
    path
  }

  if (length(args) == 0) {
    # Resolve relative to this script so the checker is runnable from repo root.
    figs_dir <- file.path(script_dir, "figs")
    if (!dir.exists(figs_dir)) {
      cat(sprintf("ERROR: figs/ directory not found: %s\n", figs_dir))
      quit(status = 1)
    }

    tex_files <- list.files(figs_dir, pattern = "\\.tex$", full.names = TRUE)

    if (length(tex_files) == 0) {
      cat(sprintf("No .tex files found in %s\n", figs_dir))
      quit(status = 0)
    }

    cat(sprintf("Checking %d figure files...\n", length(tex_files)))

    all_ok <- TRUE
    for (f in tex_files) {
      ok <- check_figure(f)
      all_ok <- all_ok && ok
    }

    cat("\n")
    cat(strrep("=", 60), "\n")
    if (all_ok) {
      cat("✓ All figures passed layout checks\n")
      quit(status = 0)
    } else {
      cat("✗ Some figures have layout warnings (see above)\n")
      quit(status = 0)  # Warnings, not errors
    }

  } else {
    # Respect paths from the caller's working directory, then script directory.
    all_ok <- TRUE
    for (f in args) {
      ok <- check_figure(resolve_figure_path(f))
      all_ok <- all_ok && ok
    }

    quit(status = ifelse(all_ok, 0, 1))
  }
}

if (!interactive() && length(grep("^--file=", commandArgs())) > 0) {
  main()
}
