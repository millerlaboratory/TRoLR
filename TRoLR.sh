#!/bin/bash
#
#   TANDEM REPEAT OUTLIERS identified via LONG READS (TRoLR)
# TRoLR.sh - per-sample pipeline for two haplotype-aware BAMs
#
# Usage: TRoLR.sh <BAM_HP1> <BAM_HP2> <KARYOTYPE: XX|XY> [OUTPUT_DIR]
#
set -euo pipefail

# Performance tracking
SCRIPT_START_TIME=$(date +%s)
SCRIPT_START_DATE=$(date)
PERF_LOG_FILE=""
MAX_PEAK_MEMORY_KB=0  # Track maximum peak memory across all runs

# Function to get current memory usage in MB
get_memory_mb() {
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux: use /proc/self/status
    awk '/VmRSS/ {print int($2/1024)}' /proc/self/status 2>/dev/null || echo "0"
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS: use ps
    ps -o rss= -p $$ 2>/dev/null | awk '{print int($1/1024)}' || echo "0"
  else
    echo "0"
  fi
}

# Function to log timing checkpoint
log_checkpoint() {
  local stage_name="$1"
  local peak_mem_kb="${2:-}"  # Optional peak memory in KB from /usr/bin/time
  local current_time=$(date +%s)
  local elapsed=$((current_time - SCRIPT_START_TIME))
  
  if [[ -n "$peak_mem_kb" ]]; then
    # Update global max peak memory if this is higher
    if (( peak_mem_kb > MAX_PEAK_MEMORY_KB )); then
      MAX_PEAK_MEMORY_KB=$peak_mem_kb
    fi
    # Convert KB to MB
    local mem_mb=$((peak_mem_kb / 1024))
    echo "CHECKPOINT: $stage_name | Elapsed: ${elapsed}s | Memory: ${mem_mb}MB" | tee -a "$PERF_LOG_FILE"
  else
    # Fallback to current memory
    local mem_mb=$(get_memory_mb)
    echo "CHECKPOINT: $stage_name | Elapsed: ${elapsed}s | Memory: ${mem_mb}MB" | tee -a "$PERF_LOG_FILE"
  fi
}

# Function to run command with /usr/bin/time wrapper for resource tracking
run_with_timing() {
  local stage_name="$1"
  shift
  local cmd=("$@")

  # Write status to log file and stderr (not stdout)
  {
    echo ">>> Running: $stage_name"
    echo "CHECKPOINT: $stage_name (started) | $(date)"
  } | tee -a "$PERF_LOG_FILE" >&2

  # Run command: only stdout from the command itself (not the status messages)
  local peak_mem_kb=""
  if command -v /usr/bin/time &> /dev/null; then
    local time_file=$(mktemp)
    /usr/bin/time -f "Time: %Es | Peak Memory: %MKB | CPU: %P" "${cmd[@]}" 2>"$time_file"
    local time_output=$(cat "$time_file")
    echo "  $time_output" >> "$PERF_LOG_FILE"
    # Parse peak memory from time output (extract number after "Peak Memory: ")
    peak_mem_kb=$(echo "$time_output" | grep -oP 'Peak Memory: \K\d+' || echo "")
    rm -f "$time_file"
  else
    # Fallback: just run the command
    "${cmd[@]}" 2>>"$PERF_LOG_FILE"
  fi

  # Write completion status to log and stderr with parsed peak memory
  log_checkpoint "$stage_name (completed)" "$peak_mem_kb" >&2
}


if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <BAM_HP1> <BAM_HP2> <KARYOTYPE: XX|XY> [EXISTING_OUTPUT_DIR]"
  exit 1
fi

SAMPLE_ID=$1
BAM1=$2
BAM2=$3
KARYOTYPE=${4:-XX}
OUTPUT_DIR=${5:-$(pwd)}

if [[ ! -f "$BAM1" ]]; then echo "BAM not found: $BAM1"; exit 1; fi
if [[ ! -f "$BAM2" ]]; then echo "BAM not found: $BAM2"; exit 1; fi

# Allow overriding repo location with REPO_ROOT; otherwise infer from this script's directory
# This script resides at the repo root, so no parent traversal
REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}

# Resources - relative to repo root
MOTIFS="$REPO_ROOT/TRoLR_reference_data/DATA_FILE_1_Curated_TRCompV2_eff0.1_noHp_noSTRchive.tsv"
STRCHIVE="$REPO_ROOT/TRoLR_reference_data/DATA_FILE_2_Curated_STRchive_loci_v1.0_20260330.bed"
LPS="$REPO_ROOT/scripts/vamos_lps.py"
ANNO="$REPO_ROOT/TRoLR_reference_data/TR_GENCODE_v.45_ANNOTATION_PHENOTYPE_SORTED.bed"
STRCHIVE_INFO="$REPO_ROOT/TRoLR_reference_data/STRCHIVE_locus_information_annotation.tsv"
PATHOGENIC_DETECTOR="$REPO_ROOT/scripts/strchive_pathogenic_detector.py"
TEST_OUTLIERS_R="$REPO_ROOT/scripts/outliers.R"

#Update with final file paths later
CONTROL_FILE="$REPO_ROOT/TRoLR_reference_data/Supp_LPS_MOTIF_ALLELE_COUNTS.tsv.gz"
REF_FILE="$REPO_ROOT/TRoLR_reference_data/DATA_FILE_5_ALL_LOCI_SUMMARY_STATS_PER_MOTIF.tsv.gz"

cd "$OUTPUT_DIR"


sample_dir="${OUTPUT_DIR}/${SAMPLE_ID}"
mkdir -p "$sample_dir"

# Initialize performance log file
PERF_LOG_FILE="${sample_dir}/${SAMPLE_ID}_performance.log"
{
  echo "======================================"
  echo "TRoLR Performance Monitoring Log"
  echo "======================================"
  echo "Sample: $SAMPLE_ID"
  echo "Karyotype: $KARYOTYPE"
  echo "Start Time: $SCRIPT_START_DATE"
  echo "Start Timestamp: $SCRIPT_START_TIME"
  echo "======================================" 
} > "$PERF_LOG_FILE"

log_checkpoint "Script initialization"

echo "Running TRoLR for sample: $SAMPLE_ID (karyotype=$KARYOTYPE)"
echo "BAMs: $BAM1 , $BAM2"

# Map bam -> hap number by checking filename suffixes; fallback to 1/2 ordering
declare -A BAM_MAP
for b in "$BAM1" "$BAM2"; do
  bn=$(basename "$b" .bam)
  if echo "$bn" | grep -qE 'hp1|hap1|mat'; then
    BAM_MAP[1]="$b"
  elif echo "$bn" | grep -qE 'hp2|hap2|pat'; then
    BAM_MAP[2]="$b"
  fi
done
# if not detected, assign in order
if [[ -z "${BAM_MAP[1]:-}" ]]; then BAM_MAP[1]="$BAM1"; fi
if [[ -z "${BAM_MAP[2]:-}" ]]; then BAM_MAP[2]="$BAM2"; fi




# run vamos for each hap (skip if outputs exist)
for hap in 1 2; do
  bam="${BAM_MAP[$hap]}"
  vcf_plain="${sample_dir}/${SAMPLE_ID}_hp${hap}.vcf"
  vcf_gz="${vcf_plain}.gz"
  strchive_plain="${sample_dir}/${SAMPLE_ID}_hp${hap}.strchive.vcf"
  #strchive_gz="${strchive_plain}.gz"

  if [[ ! -f "$vcf_gz" && ! -f "$vcf_plain" ]]; then
    echo "Running vamos motifs on hap${hap}"
    run_with_timing "VAMOS motifs hap${hap}" vamos --contig -b "$bam" -r "$MOTIFS" -s "${SAMPLE_ID}_hp${hap}" -E -S -o "$vcf_plain" -t 20 || true
    if [[ -f "$vcf_plain" ]]; then bgzip -f "$vcf_plain"; fi
  else
    echo "VCF exists for ${SAMPLE_ID}_hp${hap}, skipping motifs call"
    log_checkpoint "VCF skipped for hap${hap} (already exists)"
  fi

  if [[ ! -f "$strchive_plain" ]]; then
    echo "Running vamos STRchive on hap${hap}"
    run_with_timing "VAMOS STRchive hap${hap}" vamos --contig -b "$bam" -r "$STRCHIVE" -s "${SAMPLE_ID}_hp${hap}" -Z -E -S -o "$strchive_plain" -t 20 || true
  else
    echo "STRchive VCF exists for ${SAMPLE_ID}_hp${hap}, skipping"
    log_checkpoint "STRchive VCF skipped for hap${hap} (already exists)"
  fi
done



# Run vamos_lps.py to produce LPS bed using non-strchive vcfs (*.vcf.gz)
lps_bed="${sample_dir}/${SAMPLE_ID}_lps.bed"
if [[ ! -f "$lps_bed" ]]; then
  run_with_timing "LPS bed generation" python "$LPS" -d "$sample_dir" -o "$lps_bed" -p "*.vcf.gz"
else
  log_checkpoint "LPS bed skipped (already exists)"
fi

# Annotate with gene annotations
annotated_bed="${sample_dir}/${SAMPLE_ID}_lps_annotated.bed"
if [[ ! -f "$annotated_bed" ]]; then
  run_with_timing "Gene annotation" bedtools intersect -wa -loj -a "$lps_bed" -b "$ANNO" > "$annotated_bed"
else
  log_checkpoint "Gene annotation skipped (already exists)"
fi

# Run test_outliers.R to produce outliers table (passes karyotype)
outliers_file="${sample_dir}/${SAMPLE_ID}_lps_annotated_outliers.bed"
if [[ ! -f "$outliers_file" ]]; then
  run_with_timing "Outlier detection (R)" Rscript "$TEST_OUTLIERS_R" "$annotated_bed" "$outliers_file" "$KARYOTYPE" "$REF_FILE"
else
  log_checkpoint "Outlier detection skipped (already exists)"
fi

# Run pathogenic detector on strchive vcfs if available
pathogenic_file="${sample_dir}/${SAMPLE_ID}_pathogenic_results.tsv"
hp1_vcf="${sample_dir}/${SAMPLE_ID}_hp1.strchive.vcf.gz"
hp2_vcf="${sample_dir}/${SAMPLE_ID}_hp2.strchive.vcf.gz"
# fall back to plain strchive names if different extension
hp1_vcf_alt="${sample_dir}/${SAMPLE_ID}_hp1.strchive.vcf"
hp2_vcf_alt="${sample_dir}/${SAMPLE_ID}_hp2.strchive.vcf"

vcf1=""
vcf2=""
if [[ -f "$hp1_vcf" ]]; then vcf1="$hp1_vcf"; elif [[ -f "$hp1_vcf_alt" ]]; then vcf1="$hp1_vcf_alt"; fi
if [[ -f "$hp2_vcf" ]]; then vcf2="$hp2_vcf"; elif [[ -f "$hp2_vcf_alt" ]]; then vcf2="$hp2_vcf_alt"; fi

if [[ -n "$vcf1" && -n "$vcf2" && ! -f "$pathogenic_file" ]]; then
  run_with_timing "Pathogenic detection" python "$PATHOGENIC_DETECTOR" "$STRCHIVE_INFO" "$vcf1" "$vcf2" "$pathogenic_file"
elif [[ ! -n "$vcf1" || ! -n "$vcf2" ]]; then
  echo "STRchive VCFs not available, skipping pathogenic detection" | tee -a "$PERF_LOG_FILE"
else
  log_checkpoint "Pathogenic detection skipped (already exists)"
fi

# Generate histogram plots for exon/UTR/intron outliers (motif_len >= 3)
# Skip if plots directory already exists with plots
plots_dir="${sample_dir}/plots"

if [[ -d "$plots_dir" ]] && [[ $(find "$plots_dir" -name "*.png" 2>/dev/null | wc -l) -gt 0 ]]; then
  echo "Plots directory exists with PNG files, skipping plot generation"
  log_checkpoint "Plot generation skipped (already exists)"
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
if(!("feature" %in% colnames(out))) {
  # Try to infer from other columns
  char_cols <- names(out)[sapply(out, is.character)]
  if(length(char_cols) > 0) {
    feature_text <- apply(out[, char_cols, drop=FALSE], 1, function(r) paste(na.omit(r), collapse=" "))
    out$feature <- case_when(
      grepl("exon", feature_text, ignore.case = TRUE) ~ "exon",
      grepl("5'UTR|5UTR", feature_text, ignore.case = TRUE) ~ "5UTR",
      grepl("3'UTR|3UTR", feature_text, ignore.case = TRUE) ~ "3UTR",
      grepl("UTR", feature_text, ignore.case = TRUE) ~ "UTR",
      grepl("intron", feature_text, ignore.case = TRUE) ~ "intron",
      TRUE ~ "other"
    )
  } else {
    message("Warning: Cannot determine feature column, skipping plot generation")
    quit(status=0)
  }
}

# Filter outliers: for introns, only plot if count >= 100
out_filt <- out %>%
  filter(!is.na(motif_len) & motif_len >= 3) %>%
  filter(
    (tolower(feature) %in% c("exon", "utr", "5utr", "3utr")) |
    (tolower(feature) == "intron" & !is.na(count) & count >= 100)
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

out_pairs <- out_filt %>% select(locus, motif, feature, gene) %>% distinct()

plot_order_file <- file.path(plots_dir, ".plot_order.txt")
plot_connections <- file(plot_order_file, "w")

for(i in seq_len(nrow(out_pairs))) {
  locus_i <- out_pairs$locus[i]
  motif_i <- out_pairs$motif[i]
  feature_i  <- out_pairs$feature[i]
  gene_i <- out_pairs$gene[i]
  sample_count <- max(out_filt$count[out_filt$locus == locus_i & out_filt$motif == motif_i], na.rm=TRUE)
  if(is.infinite(sample_count)) sample_count <- NA_real_

  ctrl_sub <- ctrl %>% filter(locus == locus_i & motif == motif_i)
  if(nrow(ctrl_sub) == 0) next

  p <- ggplot(ctrl_sub, aes(x = count, y=n_alleles)) +
    geom_bar(stat="identity",fill = "#48379E", color = "black", alpha = 0.8) +
    theme_classic() +
    labs(title = paste0(feature_i, " outlier: ", locus_i, " (", motif_i, ")"),
         subtitle = paste0("Sample: ", sample_name, " — sample count: ", ifelse(is.na(sample_count), "NA", sample_count)),
         x = "Copy number (count)", y = "Frequency")
  if(!is.na(sample_count)) p <- p + geom_vline(xintercept = sample_count, color = "red", linetype = "dashed", linewidth = 0.8)

  # create short, safe filename
  locus_hash <- substr(digest(locus_i, algo="sha1"), 1, 8)
  motif_hash <- substr(digest(motif_i, algo="sha1"), 1, 8)
  short_locus <- gsub("[^A-Za-z0-9]", "_", substr(locus_i, 1, 30))
  short_motif <- gsub("[^A-Za-z0-9]", "_", substr(motif_i, 1, 20))
  out_fn <- file.path(plots_dir, paste0(sprintf("%04d", i), "_", sample_name, "_", feature_i, "_", short_locus, "_", locus_hash, "_", short_motif, "_", motif_hash, "_histogram.png"))

  # Record plot metadata in order (locus, motif, feature, plot_file)
  writeLines(paste(locus_i, motif_i, feature_i, basename(out_fn), sep="\t"), plot_connections)

  # Skip regenerating if the file already exists
  if (file.exists(out_fn)) {
    message("Skipping existing plot: ", out_fn)
    next
  }

  ggsave(filename = out_fn, plot = p, width = 8, height = 5, dpi = 150)
  message("Wrote plot: ", out_fn)
}
close(plot_connections)
message("Plot generation complete")
R_PLOT

  # run plot maker
  run_with_timing "Histogram plot generation" Rscript "$temp_r_plot" "$outliers_file" "$CONTROL_FILE" "$SAMPLE_ID" "$plots_dir" || echo "Plot generation failed or skipped"
  rm -f "$temp_r_plot"
fi

# R Markdown report generation with integrated plot display
report_rmd="${sample_dir}/${SAMPLE_ID}_outlier_report.Rmd"
report_html="${sample_dir}/${SAMPLE_ID}_outlier_report.html"

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
library(ggplot2)
library(digest)
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
outliers_path <- Sys.getenv("OUTLIERS_PATH")
pathogenic_path <- Sys.getenv("PATHOGENIC_PATH")
sample_name <- Sys.getenv("SAMPLE_NAME")
plots_dir_path <- Sys.getenv("PLOTS_DIR")
control_file <- Sys.getenv("CONTROL_FILE")

# Function to find and display plots for a specific feature
display_plots_for_feature <- function(feature_name, plot_connections_list, plot_files_legacy) {
  # First try using the ordered list
  if(!is.null(plot_connections_list[[feature_name]]) && length(plot_connections_list[[feature_name]]) > 0){
    feature_plots <- plot_connections_list[[feature_name]]
  } else {
    # Fallback to legacy file search
    feature_plots <- plot_files_legacy[grepl(paste0("_", feature_name, "_"), plot_files_legacy, ignore.case = TRUE)]
  }
  
  if(length(feature_plots) > 0) {
    cat("\n\n#### Distribution plots\n\n")
    cat("Control population distributions with sample value (red dashed line):\n\n")
    for(plot_file in feature_plots) {
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
      feature_text <- apply(out[, char_cols, drop=FALSE], 1, function(r) paste(na.omit(r), collapse=" "))
      is_exon   <- grepl("exon", feature_text, ignore.case = TRUE, perl = TRUE)
      is_utr    <- grepl("UTR", feature_text, ignore.case = TRUE, perl = TRUE)
      is_intron <- grepl("intron", feature_text, ignore.case = TRUE, perl = TRUE)
    } else {
      is_exon <- is_utr <- is_intron <- rep(FALSE, nrow(out))
    }
    
    # Count each category (distinguish 5UTR and 3UTR)
    short_count <- sum(out$motif_len < 3)
    long_mask <- out$motif_len >= 3
    
    exon_count <- sum(is_exon & long_mask)
    is_5utr <- grepl("5'UTR|5UTR", feature_text, ignore.case = TRUE)
    is_3utr <- grepl("3'UTR|3UTR", feature_text, ignore.case = TRUE)
    utr_5_count <- sum(!is_exon & is_5utr & long_mask)
    utr_3_count <- sum(!is_exon & is_3utr & long_mask)
    intron_high <- sum(!is_exon & !is_utr & is_intron & !is.na(out$count) & out$count >= 100 & long_mask)
    intron_low <- sum(!is_exon & !is_utr & is_intron & (is.na(out$count) | out$count < 100) & long_mask)
    other_count <- sum(!is_exon & !is_utr & !is_intron & long_mask)
    
    # Count pathogenic calls if available
    pathogenic_count <- 0
    if(nzchar(pathogenic_path) && file.exists(pathogenic_path)){
      pat_temp <- tryCatch(read.delim(pathogenic_path, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
      if(!is.null(pat_temp) && nrow(pat_temp) > 0){
        pathogenic_count <- nrow(pat_temp)
      }
    }
    
    summary_df <- data.frame(
      Type = c("Exon", "5'UTR", "3'UTR", "Intron (count ≥ 100)", "Intron (count < 100)", "Other/Unclassified", "Short motifs (<3 bp)", "Pathogenic calls"),
      Count = c(exon_count, utr_5_count, utr_3_count, intron_high, intron_low, other_count, short_count, pathogenic_count)
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
*Pathogenic calls count the number of a pathogenic-associated motif regardless of interruption, resulting in potential false positives. If a count is at or near the pathogenic minimum, please consult the genotype files.*

```{r pathogenic-table, results='asis'}
if(nzchar(pathogenic_path) && file.exists(pathogenic_path)){
  pat <- read.delim(pathogenic_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if(nrow(pat)>0){
    cat("**Number of pathogenic calls:** ", nrow(pat), "\n\n")
    DT::datatable(pat, options=list(pageLength=25, scrollX=TRUE))
  } else {
    cat("No pathogenic calls found.\n")
  }
} else {
  cat("Pathogenic results not available.\n")
}
```

```{r pathogenic-plots, results='asis'}
# Generate plots for pathogenic calls using control data
if(nzchar(pathogenic_path) && file.exists(pathogenic_path) && file.exists(control_file)){
  pat <- read.delim(pathogenic_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  ctrl_all <- tryCatch(read.delim(control_file, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
  
  if(!is.null(ctrl_all) && nrow(ctrl_all) > 0 && nrow(pat) > 0){
    # Ensure necessary columns
    if(all(c("locus","motif","count") %in% colnames(ctrl_all))){
      ctrl_all$motif <- as.character(ctrl_all$motif)
      ctrl_all$count <- suppressWarnings(as.numeric(ctrl_all$count))
      
      # Check if n_alleles exists; if not, create it from row frequencies
      if(!("n_alleles" %in% colnames(ctrl_all))){
        ctrl_all <- ctrl_all %>% group_by(locus, motif, count) %>% mutate(n_alleles = n()) %>% ungroup()
      } else {
        ctrl_all$n_alleles <- suppressWarnings(as.numeric(ctrl_all$n_alleles))
      }
      
      cat("\n### Population distributions for pathogenic calls\n\n")
      cat("Control population distributions with sample value (red dashed line):\n\n")
      
      for(i in seq_len(nrow(pat))){
        locus_i <- pat$locus[i]
        motif_i <- as.character(pat$motif[i])
        sample_count <- suppressWarnings(as.numeric(pat$count[i]))
        gene_i <- pat$gene[i]
        
        # Match control data for this locus and motif
        ctrl_sub <- ctrl_all %>% 
          filter(locus == locus_i & motif == motif_i) %>%
          distinct()
        
        if(nrow(ctrl_sub) > 0){
          p <- ggplot(ctrl_sub, aes(x = count, y = n_alleles)) +
            geom_bar(stat="identity", fill = "#48379E", color = "black", alpha = 0.8) +
            theme_classic() +
            labs(title = paste0("Pathogenic: ", locus_i),
                 subtitle = paste0("Gene: ", gene_i, " | Motif: ", motif_i, " | Sample: ", sample_name),
                 x = "Copy number (count)", y = "Frequency")
          
          if(!is.na(sample_count)){
            p <- p + geom_vline(xintercept = sample_count, color = "red", linetype = "dashed", linewidth = 0.8)
          }
          
          # Create safe filename with hash
          locus_hash <- substr(digest(locus_i, algo="sha1"), 1, 8)
          motif_hash <- substr(digest(motif_i, algo="sha1"), 1, 8)
          plot_file <- tempfile(paste0("pathogenic_", locus_hash, "_", motif_hash, "_"), fileext=".png")
          
          ggsave(filename = plot_file, plot = p, width = 8, height = 5, dpi = 150)
          cat(sprintf("![](%s){width=85%%}\n\n", plot_file))
        }
      }
    } else {
      cat("Control file missing required columns (locus, motif, count).\n")
    }
  } else if(nzchar(pathogenic_path) && file.exists(pathogenic_path)){
    pat <- read.delim(pathogenic_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    if(nrow(pat) > 0){
      cat("No control data available for pathogenic call plots.\n")
    }
  }
}
```

## Outliers

```{r outliers-preprocessing, results='asis'}
# Load and process outliers data
plot_files <- character(0)
plot_order_list <- list()

# Load plot order from file if it exists
plot_order_file <- file.path(plots_dir_path, ".plot_order.txt")
if(nzchar(plots_dir_path) && dir.exists(plots_dir_path) && file.exists(plot_order_file)){
  plot_order <- read.delim(plot_order_file, header=FALSE, col.names=c("locus","motif","feature","filename"), stringsAsFactors=FALSE)
  # Create lookup lists by feature to maintain order
  for(feat in unique(plot_order$feature)){
    plot_order_list[[feat]] <- file.path(plots_dir_path, plot_order$filename[plot_order$feature == feat])
  }
} else if(nzchar(plots_dir_path) && dir.exists(plots_dir_path)){
  # Fallback to listing files if order file doesn't exist
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
    is_5utr   <- grepl("5'UTR|5UTR", type_text, ignore.case = TRUE, perl = TRUE)
    is_3utr   <- grepl("3'UTR|3UTR", type_text, ignore.case = TRUE, perl = TRUE)
    is_utr    <- grepl("UTR", type_text, ignore.case = TRUE, perl = TRUE)
    is_intron <- grepl("intron", type_text, ignore.case = TRUE, perl = TRUE)

    # Select only loci with motif_len >= 3 for these tables
    exon_tab   <- long_motifs[is_exon, , drop=FALSE]
    utr_5_tab  <- long_motifs[!is_exon & is_5utr, , drop=FALSE]
    utr_3_tab  <- long_motifs[!is_exon & is_3utr, , drop=FALSE]
    
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
if(exists("exon_tab") && nrow(exon_tab) > 0 && (length(plot_files) > 0 || length(plot_order_list) > 0)) {
  display_plots_for_feature("exon", plot_order_list, plot_files)
}
```

### 5'UTR outliers

```{r utr-5-outliers-table, results='asis'}
if(exists("utr_5_tab") && nrow(utr_5_tab) > 0){
  DT::datatable(utr_5_tab, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No 5'UTR outliers found.\n\n")
}
```

```{r utr-5-outliers-plots, results='asis'}
if(exists("utr_5_tab") && nrow(utr_5_tab) > 0 && (length(plot_files) > 0 || length(plot_order_list) > 0)) {
  display_plots_for_feature("5UTR", plot_order_list, plot_files)
}
```

### 3'UTR outliers

```{r utr-3-outliers-table, results='asis'}
if(exists("utr_3_tab") && nrow(utr_3_tab) > 0){
  DT::datatable(utr_3_tab, options=list(pageLength=25, scrollX=TRUE))
} else {
  cat("No 3'UTR outliers found.\n\n")
}
```

```{r utr-3-outliers-plots, results='asis'}
if(exists("utr_3_tab") && nrow(utr_3_tab) > 0 && (length(plot_files) > 0 || length(plot_order_list) > 0)) {
  display_plots_for_feature("3UTR", plot_order_list, plot_files)
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
if(exists("intron_tab") && nrow(intron_tab) > 0 && (length(plot_files) > 0 || length(plot_order_list) > 0)) {
  display_plots_for_feature("intron", plot_order_list, plot_files)
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
if(exists("other_tab") && nrow(other_tab) > 0) {
  # Look for plots that don't match exon/5UTR/3UTR/intron patterns
  if(length(plot_order_list) > 0){
    # Get all plots not in the known feature categories
    all_features <- names(plot_order_list)
    known_features <- c("exon", "5UTR", "3UTR", "intron")
    other_features <- setdiff(all_features, known_features)
    
    for(feat in other_features){
      if(!is.null(plot_order_list[[feat]]) && length(plot_order_list[[feat]]) > 0){
        if(!exists("other_plots_shown")){
          cat("\n\n#### Distribution plots\n\n")
          cat("Control population distributions with sample value (red dashed line):\n\n")
          other_plots_shown <- TRUE
        }
        for(plot_file in plot_order_list[[feat]]){
          cat(sprintf("![](%s){width=85%%}\n\n", plot_file))
        }
      }
    }
  } else if(length(plot_files) > 0){
    # Fallback to legacy filtering
    other_plots <- plot_files[!grepl("_exon_|_5UTR_|_3UTR_|_intron_", plot_files)]
    if(length(other_plots) > 0) {
      cat("\n\n#### Distribution plots\n\n")
      cat("Control population distributions with sample value (red dashed line):\n\n")
      for(plot_file in other_plots) {
        cat(sprintf("![](%s){width=85%%}\n\n", plot_file))
      }
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

sed -i "s/SAMPLE_PLACEHOLDER/${SAMPLE_ID}/g" "$report_rmd"
export OUTLIERS_PATH="$outliers_file"
export PATHOGENIC_PATH="${pathogenic_file:-}"
export SAMPLE_NAME="$SAMPLE_ID"
export PLOTS_DIR="$plots_dir"
export CONTROL_FILE="$CONTROL_FILE"

# render the report
echo "Rendering HTML report..."
run_with_timing "HTML report generation" Rscript -e "rmarkdown::render('$report_rmd', output_file='$report_html', quiet=TRUE)" || echo "Rmd render failed"

# cleanup
#conda deactivate || true

# Final summary
SCRIPT_END_TIME=$(date +%s)
SCRIPT_END_DATE=$(date)
TOTAL_ELAPSED=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
TOTAL_MINUTES=$((TOTAL_ELAPSED / 60))
TOTAL_SECONDS=$((TOTAL_ELAPSED % 60))
# Use tracked peak memory instead of current memory
FINAL_MEMORY=$((MAX_PEAK_MEMORY_KB / 1024))
if (( FINAL_MEMORY == 0 )); then
  # Fallback to current memory if no peak was tracked
  FINAL_MEMORY=$(get_memory_mb)
fi

{
  echo ""
  echo "======================================"
  echo "TRoLR Pipeline Complete - Summary"
  echo "======================================"
  echo "Sample: $SAMPLE_ID"
  echo "End Time: $SCRIPT_END_DATE"
  echo "End Timestamp: $SCRIPT_END_TIME"
  echo "Total Runtime: ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s (${TOTAL_ELAPSED}s)"
  echo "Peak Memory Usage: ${FINAL_MEMORY}MB"
  echo "======================================"
} | tee -a "$PERF_LOG_FILE"

echo "Pipeline complete for sample: $SAMPLE_ID"
echo "Outputs:"
echo " - LPS: $lps_bed"
echo " - Annotated LPS: $annotated_bed"
echo " - Outliers: $outliers_file"
echo " - Pathogenic calls: $pathogenic_file"
echo " - Histogram plots: $plots_dir"
echo " - HTML report: $report_html"
echo " - Performance log: $PERF_LOG_FILE"