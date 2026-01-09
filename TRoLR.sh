#!/bin/bash
#
#   TANDEM REPEAT OUTLIERS identified via LONG READS (TRoLR)
# TRoLR.sh - per-sample pipeline for two haplotype-aware BAMs
#
# Usage: TRoLR.sh <BAM_HP1> <BAM_HP2> <KARYOTYPE: XX|XY> [OUTPUT_DIR]
#
set -euo pipefail


if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <BAM_HP1> <BAM_HP2> <KARYOTYPE: XX|XY> [OUTPUT_DIR]"
  exit 1
fi

BAM1=$1
BAM2=$2
KARYOTYPE=${3:-XX}
OUTPUT_DIR=${4:-$(pwd)}

if [[ ! -f "$BAM1" ]]; then echo "BAM not found: $BAM1"; exit 1; fi
if [[ ! -f "$BAM2" ]]; then echo "BAM not found: $BAM2"; exit 1; fi

# Allow overriding repo location with REPO_ROOT; otherwise infer from this script's directory
# This script resides at the repo root, so no parent traversal
REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

# Resources - relative to repo root
MOTIFS="$REPO_ROOT/reference_data/vamos.motif.hg38.v2.1.e0.1.noSTRCHIVE.nohp.bed"
STRCHIVE="$REPO_ROOT/reference_data/vamos_strchive.B2FLLAIV.20250520.bed"
LPS="$REPO_ROOT/scripts/vamos_lps.py"
ANNO="$REPO_ROOT/reference_data/GENCODE_v.45_CANONICAL.bed"
STRCHIVE_INFO="$REPO_ROOT/reference_data/STRchive-disease-loci-v2.4.3.hg38.CE2vK2zA.tsv"
PATHOGENIC_DETECTOR="$REPO_ROOT/scripts/strchive_pathogenic_detector.py"
TEST_OUTLIERS_R="$REPO_ROOT/scripts/outliers.R"
CONTROL_FILE="$REPO_ROOT/reference_data/vamos_asm_lps_e0.1_247_catalog_control_length_counts.tsv.gz"
REF_FILE="$REPO_ROOT/reference_data/vamos_asm_lps_control_summary.tsv.gz"

#module load bedtools/2.31.1
#module load samtools/1.22
#module load R
#source activate vamos-env || source activate vamos-env || true

#mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# derive sample base name (strip _hp1/_hp2 or _hap1/_hap2)
base1=$(basename "$BAM1" .hg38.bam)
base2=$(basename "$BAM2" .hg38.bam)
sample1=$(echo "$base1" | cut -d"_" -f1)
sample2=$(echo "$base2" | cut -d"_" -f1)

if [[ "$sample1" != "$sample2" ]]; then
  echo "Sample base names do not match: $sample1 vs $sample2"
  exit 1
fi
sample="$sample1"
sample_dir="${OUTPUT_DIR}/${sample}"
mkdir -p "$sample_dir"

echo "Running TRoLR for sample: $sample (karyotype=$KARYOTYPE)"
echo "BAMs: $BAM1 , $BAM2"

# Map bam -> hap number by checking filename suffixes; fallback to 1/2 ordering
declare -A BAM_MAP
for b in "$BAM1" "$BAM2"; do
  bn=$(basename "$b" .bam)
  if echo "$bn" | grep -qE 'hp1|hap1'; then
    BAM_MAP[1]="$b"
  elif echo "$bn" | grep -qE 'hp2|hap2'; then
    BAM_MAP[2]="$b"
  fi
done
# if not detected, assign in order
if [[ -z "${BAM_MAP[1]:-}" ]]; then BAM_MAP[1]="$BAM1"; fi
if [[ -z "${BAM_MAP[2]:-}" ]]; then BAM_MAP[2]="$BAM2"; fi

# run vamos for each hap (skip if outputs exist)
for hap in 1 2; do
  bam="${BAM_MAP[$hap]}"
  vcf_plain="${sample_dir}/${sample}_hp${hap}.vcf"
  vcf_gz="${vcf_plain}.gz"
  strchive_plain="${sample_dir}/${sample}_hp${hap}.strchive.vcf"
  #strchive_gz="${strchive_plain}.gz"

  if [[ ! -f "$vcf_gz" && ! -f "$vcf_plain" ]]; then
    echo "Running vamos motifs on hap${hap}"
    vamos --contig -b "$bam" -r "$MOTIFS" -s "${sample}_hp${hap}" -S -o "$vcf_plain" -t 20 || true
    if [[ -f "$vcf_plain" ]]; then bgzip -f "$vcf_plain"; fi
  else
    echo "VCF exists for ${sample}_hp${hap}, skipping motifs call"
  fi

  if [[ ! -f "$strchive_plain" ]]; then
    echo "Running vamos STRchive on hap${hap}"
    vamos --contig -b "$bam" -r "$STRCHIVE" -s "${sample}_hp${hap}" -S -o "$strchive_plain" -t 20 || true
  else
    echo "STRchive VCF exists for ${sample}_hp${hap}, skipping"
  fi
done

#source deactivate vamos-env || true

# Run vamos_lps.py to produce LPS bed using non-strchive vcfs (*.vcf.gz)
lps_bed="${sample_dir}/${sample}_lps.bed"
if [[ ! -f "$lps_bed" ]]; then
  python "$LPS" -d "$sample_dir" -o "$lps_bed" -p "*.vcf.gz"
fi

# Annotate with gene annotations
annotated_bed="${sample_dir}/${sample}_lps_annotated.bed"
if [[ ! -f "$annotated_bed" ]]; then
  bedtools intersect -wa -loj -a "$lps_bed" -b "$ANNO" > "$annotated_bed"
fi

# Run test_outliers.R to produce outliers table (passes karyotype)
outliers_file="${sample_dir}/${sample}_lps_annotated_outliers.bed"
if [[ ! -f "$outliers_file" ]]; then
  Rscript "$TEST_OUTLIERS_R" "$annotated_bed" "$outliers_file" "$KARYOTYPE" "$REF_FILE"
fi

# Run pathogenic detector on strchive vcfs if available
pathogenic_file="${sample_dir}/${sample}_pathogenic_results.tsv"
hp1_vcf="${sample_dir}/${sample}_hp1.strchive.vcf.gz"
hp2_vcf="${sample_dir}/${sample}_hp2.strchive.vcf.gz"
# fall back to plain strchive names if different extension
hp1_vcf_alt="${sample_dir}/${sample}_hp1.strchive.vcf"
hp2_vcf_alt="${sample_dir}/${sample}_hp2.strchive.vcf"

vcf1=""
vcf2=""
if [[ -f "$hp1_vcf" ]]; then vcf1="$hp1_vcf"; elif [[ -f "$hp1_vcf_alt" ]]; then vcf1="$hp1_vcf_alt"; fi
if [[ -f "$hp2_vcf" ]]; then vcf2="$hp2_vcf"; elif [[ -f "$hp2_vcf_alt" ]]; then vcf2="$hp2_vcf_alt"; fi

if [[ -n "$vcf1" && -n "$vcf2" && ! -f "$pathogenic_file" ]]; then
  python "$PATHOGENIC_DETECTOR" "$STRCHIVE_INFO" "$vcf1" "$vcf2" "$pathogenic_file"
fi

# Generate histogram plots for exon/UTR/intron outliers (motif_len >= 3)
# Skip if plots directory already exists with plots
plots_dir="${sample_dir}/plots"

if [[ -d "$plots_dir" ]] && [[ $(find "$plots_dir" -name "*.png" 2>/dev/null | wc -l) -gt 0 ]]; then
  echo "Plots directory exists with PNG files, skipping plot generation"
else
  echo "Generating histogram plots..."
  mkdir -p "$plots_dir"
  
  temp_r_plot="${sample_dir}/make_plots_from_controls.R"
  cat > "$temp_r_plot" <<'R_PLOT'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
outliers_path <- args[1]
control_file  <- args[2]
sample_name   <- args[3]
plots_dir     <- args[4]

library(data.table)
library(dplyr)
library(ggplot2)
library(digest)

if(!file.exists(outliers_path)) {
  stop("Outliers file not found: ", outliers_path)
}
out <- tryCatch(fread(outliers_path, header=TRUE, sep="\t", data.table=FALSE), error=function(e) NULL)
if(is.null(out) || nrow(out)==0) {
  message("No outliers to plot.")
  quit(status=0)
}

# ensure locus column exists
if(!("locus" %in% colnames(out))) {
  if(all(c("chr","start","end") %in% colnames(out))) {
    out$locus <- paste0(out$chr, ":", out$start, "-", out$end)
  } else {
    stop("Cannot construct locus: missing locus or chr/start/end in outliers file")
  }
}

out$motif <- as.character(out$motif)
out$count <- suppressWarnings(as.numeric(out$count))
out$motif_len <- nchar(out$motif)

# Build type classification if not present
if(!("type" %in% colnames(out))) {
  # Try to infer from other columns
  char_cols <- names(out)[sapply(out, is.character)]
  if(length(char_cols) > 0) {
    type_text <- apply(out[, char_cols, drop=FALSE], 1, function(r) paste(na.omit(r), collapse=" "))
    out$type <- case_when(
      grepl("exon", type_text, ignore.case = TRUE) ~ "exon",
      grepl("UTR", type_text, ignore.case = TRUE) ~ "UTR",
      grepl("intron", type_text, ignore.case = TRUE) ~ "intron",
      TRUE ~ "other"
    )
  } else {
    message("Warning: Cannot determine type column, skipping plot generation")
    quit(status=0)
  }
}

# Filter outliers: for introns, only plot if count >= 100
out_filt <- out %>%
  filter(!is.na(motif_len) & motif_len >= 3) %>%
  filter(
    (tolower(type) %in% c("exon", "utr")) |
    (tolower(type) == "intron" & !is.na(count) & count >= 100)
  )

if(nrow(out_filt) == 0) {
  message("No outliers meeting criteria to plot.")
  quit(status=0)
}

if(!file.exists(control_file)) {
  message("Warning: Control file not found, skipping plot generation: ", control_file)
  quit(status=0)
}
ctrl <- tryCatch(fread(control_file, header=TRUE, sep="\t", data.table=FALSE), error=function(e) NULL)
if(is.null(ctrl) || nrow(ctrl) == 0) {
  message("Control file empty or unreadable: ", control_file)
  quit(status=0)
}
if(!all(c("locus","motif","count") %in% colnames(ctrl))) {
  message("Control file must contain columns: locus, motif, count")
  quit(status=0)
}
ctrl$motif <- as.character(ctrl$motif)
ctrl$count <- suppressWarnings(as.numeric(ctrl$count))

if(!dir.exists(plots_dir)) dir.create(plots_dir, recursive=TRUE)

out_pairs <- out_filt %>% select(locus, motif, type, gene) %>% distinct()

for(i in seq_len(nrow(out_pairs))) {
  locus_i <- out_pairs$locus[i]
  motif_i <- out_pairs$motif[i]
  type_i  <- out_pairs$type[i]
  gene_i <- out_pairs$gene[i]
  sample_count <- max(out_filt$count[out_filt$locus == locus_i & out_filt$motif == motif_i], na.rm=TRUE)
  if(is.infinite(sample_count)) sample_count <- NA_real_

  ctrl_sub <- ctrl %>% filter(locus == locus_i & motif == motif_i & type == type_i, gene == gene_i)
  if(nrow(ctrl_sub) == 0) next

  p <- ggplot(ctrl_sub, aes(x = count, y=n_alleles)) +
    geom_bar(stat="identity",fill = "#48379E", color = "black", alpha = 0.8) +
    theme_classic() +
    labs(title = paste0(type_i, " outlier: ", locus_i, " (", motif_i, ")"),
         subtitle = paste0("Sample: ", sample_name, " — sample count: ", ifelse(is.na(sample_count), "NA", sample_count)),
         x = "Copy number (count)", y = "Frequency")
  if(!is.na(sample_count)) p <- p + geom_vline(xintercept = sample_count, color = "red", linetype = "dashed", linewidth = 0.8)

  # create short, safe filename
  locus_hash <- substr(digest(locus_i, algo="sha1"), 1, 8)
  motif_hash <- substr(digest(motif_i, algo="sha1"), 1, 8)
  short_locus <- gsub("[^A-Za-z0-9]", "_", substr(locus_i, 1, 30))
  short_motif <- gsub("[^A-Za-z0-9]", "_", substr(motif_i, 1, 20))
  out_fn <- file.path(plots_dir, paste0(sample_name, "_", type_i, "_", short_locus, "_", locus_hash, "_", short_motif, "_", motif_hash, "_histogram.png"))

  # Skip regenerating if the file already exists
  if (file.exists(out_fn)) {
    message("Skipping existing plot: ", out_fn)
    next
  }

  ggsave(filename = out_fn, plot = p, width = 8, height = 5, dpi = 150)
  message("Wrote plot: ", out_fn)
}
message("Plot generation complete")
R_PLOT

  # run plot maker
  Rscript "$temp_r_plot" "$outliers_file" "$CONTROL_FILE" "$sample" "$plots_dir" || echo "Plot generation failed or skipped"
  rm -f "$temp_r_plot"
fi

# R Markdown report generation with integrated plot display
report_rmd="${sample_dir}/${sample}_outlier_report.Rmd"
report_html="${sample_dir}/${sample}_outlier_report.html"

cat > "$report_rmd" <<'RMD'
---
title: "TRoLR Outlier Report - SAMPLE_PLACEHOLDER"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
date: "`r Sys.Date()`"
---
*Outliers identified in TRoLR analysis pipeline. May not reflect all repeat expansions in the sample. Please refer to the full TRoLR pipeline outputs for comprehensive variant data. Not indented for clinical use.*

```{r setup, include=FALSE}
library(knitr)
library(dplyr)
library(DT)
library(htmltools)
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
outliers_path <- Sys.getenv("OUTLIERS_PATH")
pathogenic_path <- Sys.getenv("PATHOGENIC_PATH")
sample_name <- Sys.getenv("SAMPLE_NAME")
plots_dir_path <- Sys.getenv("PLOTS_DIR")

# Function to find and display plots for a specific type
display_plots_for_type <- function(type_name, plot_files) {
  # Find plots matching this type
  type_plots <- plot_files[grepl(paste0("_", type_name, "_"), plot_files, ignore.case = TRUE)]
  
  if(length(type_plots) > 0) {
    cat("\n\n#### Distribution plots\n\n")
    cat("Control population distributions with sample value (red dashed line):\n\n")
    for(plot_file in type_plots) {
      # Display plot without filename
      cat(sprintf("![](%s){width=85%%}\n\n", plot_file))
    }
  }
}
```

## Summary

```{r summary-table, results='asis'}
if(file.exists(outliers_path)){
  out <- read.delim(outliers_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if(nrow(out) > 0){
    # motif length
    out$motif_len <- ifelse(is.na(out$motif), 0L, nchar(as.character(out$motif)))
    out$count <- suppressWarnings(as.numeric(out$count))
    
    # Build type classification
    char_cols <- names(out)[sapply(out, is.character)]
    if(length(char_cols) > 0){
      type_text <- apply(out[, char_cols, drop=FALSE], 1, function(r) paste(na.omit(r), collapse=" "))
      is_exon   <- grepl("exon", type_text, ignore.case = TRUE, perl = TRUE)
      is_utr    <- grepl("UTR", type_text, ignore.case = TRUE, perl = TRUE)
      is_intron <- grepl("intron", type_text, ignore.case = TRUE, perl = TRUE)
    } else {
      is_exon <- is_utr <- is_intron <- rep(FALSE, nrow(out))
    }
    
    # Count each category
    short_count <- sum(out$motif_len < 3)
    long_mask <- out$motif_len >= 3
    
    exon_count <- sum(is_exon & long_mask)
    utr_count <- sum(!is_exon & is_utr & long_mask)
    intron_high <- sum(!is_exon & !is_utr & is_intron & !is.na(out$count) & out$count >= 100 & long_mask)
    intron_low <- sum(!is_exon & !is_utr & is_intron & (is.na(out$count) | out$count < 100) & long_mask)
    other_count <- sum(!is_exon & !is_utr & !is_intron & long_mask)
    
    summary_df <- data.frame(
      Type = c("Exon", "UTR", "Intron (count ≥ 100)", "Intron (count < 100)", "Other/Unclassified", "Short motifs (<3 bp)"),
      Count = c(exon_count, utr_count, intron_high, intron_low, other_count, short_count)
    )
    
    cat("**Outlier type summary:**\n\n")
    DT::datatable(summary_df, 
                  options = list(paging = FALSE, searching = FALSE, info = FALSE, 
                                dom = 't', columnDefs = list(list(className = 'dt-center', targets = 1))),
                  rownames = FALSE)
  } else {
    cat("No outliers found in analysis.\n")
  }
} else {
  cat("Outliers data not available.\n")
}
```

## Pathogenic calls

```{r pathogenic, results='asis'}
if(nzchar(pathogenic_path) && file.exists(pathogenic_path)){
  pat <- read.delim(pathogenic_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if(nrow(pat)>0){
    DT::datatable(pat, options=list(pageLength=25, scrollX=TRUE))
  } else {
    cat("No pathogenic calls found.\n")
  }
} else {
  cat("Pathogenic results not available.\n")
}
```

## Outliers

```{r outliers-preprocessing, results='asis'}
# Load and process outliers data
plot_files <- character(0)
if(nzchar(plots_dir_path) && dir.exists(plots_dir_path)){
  plot_files <- list.files(plots_dir_path, pattern="\\.png$", full.names=TRUE)
}

if(file.exists(outliers_path)){
  out <- read.delim(outliers_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if(nrow(out)>0){
    # motif length (robust)
    out$motif_len <- ifelse(is.na(out$motif), 0L, nchar(as.character(out$motif)))
    # Convert count to numeric
    out$count <- suppressWarnings(as.numeric(out$count))
    
    short_motifs <- out %>% filter(motif_len < 3)
    long_motifs  <- out %>% filter(motif_len >= 3)

    # Build a single searchable text per row from all character columns (use long_motifs for classification)
    char_cols <- names(long_motifs)[sapply(long_motifs, is.character)]
    if(length(char_cols) == 0){
      type_text <- rep("", nrow(long_motifs))
    } else {
      type_text <- apply(long_motifs[, char_cols, drop=FALSE], 1, function(r) paste(na.omit(r), collapse=" "))
    }

    # classify by searching the combined annotation text (robust to different column names)
    is_exon   <- grepl("exon", type_text, ignore.case = TRUE, perl = TRUE)
    is_utr    <- grepl("UTR", type_text, ignore.case = TRUE, perl = TRUE)
    is_intron <- grepl("intron", type_text, ignore.case = TRUE, perl = TRUE)

    # Select only loci with motif_len >= 3 for these tables
    exon_tab   <- long_motifs[is_exon, , drop=FALSE]
    utr_tab    <- long_motifs[!is_exon & is_utr, , drop=FALSE]
    
    # For introns: only include if count >= 100
    intron_candidates <- long_motifs[!is_exon & !is_utr & is_intron, , drop=FALSE]
    intron_tab <- intron_candidates %>% filter(!is.na(count) & count >= 100)
    intron_low_count <- intron_candidates %>% filter(is.na(count) | count < 100)
    
    # Other table includes non-annotated plus low-count introns
    other_tab_base <- long_motifs[!is_exon & !is_utr & !is_intron, , drop=FALSE]
    other_tab <- rbind(other_tab_base, intron_low_count)
  } else {
    cat("No outliers found.\n")
  }
} else {
  cat("Outliers file not found:", outliers_path, "\n")
}
```

### Exon outliers

```{r exon-outliers-table, results='asis'}
if(exists("exon_tab") && nrow(exon_tab) > 0){
  DT::datatable(exon_tab, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No exon outliers found.\n\n")
}
```

```{r exon-outliers-plots, results='asis'}
if(exists("exon_tab") && nrow(exon_tab) > 0 && length(plot_files) > 0) {
  display_plots_for_type("exon", plot_files)
}
```

### UTR outliers

```{r utr-outliers-table, results='asis'}
if(exists("utr_tab") && nrow(utr_tab) > 0){
  DT::datatable(utr_tab, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No UTR outliers found.\n\n")
}
```

```{r utr-outliers-plots, results='asis'}
if(exists("utr_tab") && nrow(utr_tab) > 0 && length(plot_files) > 0) {
  display_plots_for_type("UTR", plot_files)
}
```

### Intron outliers (count ≥ 100)

```{r intron-outliers-table, results='asis'}
if(exists("intron_tab") && nrow(intron_tab) > 0){
  cat("*Note: Only showing intron outliers with count ≥ 100. Lower count intron outliers are in the Other/Unclassified section.*\n\n")
  DT::datatable(intron_tab, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No intron outliers with count ≥ 100 found.\n\n")
}
```

```{r intron-outliers-plots, results='asis'}
if(exists("intron_tab") && nrow(intron_tab) > 0 && length(plot_files) > 0) {
  display_plots_for_type("intron", plot_files)
}
```

### Other/Unclassified outliers

```{r other-outliers-table, results='asis'}
if(exists("other_tab") && nrow(other_tab) > 0){
  cat("*Includes non-annotated outliers and intron outliers with count < 100.*\n\n")
  DT::datatable(other_tab, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No other/unclassified outliers found.\n\n")
}
```

```{r other-outliers-plots, results='asis'}
if(exists("other_tab") && nrow(other_tab) > 0 && length(plot_files) > 0) {
  # Look for plots that don't match exon/UTR/intron
  other_plots <- plot_files[!grepl("_exon_|_UTR_|_utr_|_intron_", plot_files)]
  if(length(other_plots) > 0) {
    cat("\n\n#### Distribution plots\n\n")
    cat("Control population distributions with sample value (red dashed line):\n\n")
    for(plot_file in other_plots) {
      # Display plot without filename
      cat(sprintf("![](%s){width=85%%}\n\n", plot_file))
    }
  }
}
```

### Short motifs (<3 bp)

```{r short-motifs-table, results='asis'}
if(exists("short_motifs") && nrow(short_motifs) > 0){
  cat("These short motif outliers are reported separately and do not have associated distribution plots.\n\n")
  DT::datatable(short_motifs, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No short-motif outliers found.\n\n")
}
```

---

*Report generated on `r Sys.Date()` for sample `r sample_name`*
RMD

sed -i "s/SAMPLE_PLACEHOLDER/${sample}/g" "$report_rmd"
export OUTLIERS_PATH="$outliers_file"
export PATHOGENIC_PATH="${pathogenic_file:-}"
export SAMPLE_NAME="$sample"
export PLOTS_DIR="$plots_dir"

# render the report
echo "Rendering HTML report..."
Rscript -e "rmarkdown::render('$report_rmd', output_file='$report_html', quiet=TRUE)" || echo "Rmd render failed"

# cleanup
#conda deactivate || true

echo "Pipeline complete for sample: $sample"
echo "Outputs:"
echo " - LPS: $lps_bed"
echo " - Annotated LPS: $annotated_bed"
echo " - Outliers: $outliers_file"
echo " - Pathogenic calls: $pathogenic_file"
echo " - Histogram plots: $plots_dir"
echo " - HTML report: $report_html"