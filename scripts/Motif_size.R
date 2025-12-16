
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript find_motif_size_conflicts.R <input.tsv> [output.tsv]")
}
input_path <- args[1]
output_path <- if (length(args) >= 2) args[2] else "motif_size_conflicts.tsv"

dt <- data.table::fread(input_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))

# normalize expected column names
colnames(dt) <- make.names(colnames(dt))
if (!all(c("locus", "motif", "count") %in% tolower(colnames(dt)))) {
  # try to map by name ignoring case
  names(dt) <- tolower(names(dt))
}
if (!all(c("locus", "motif", "count") %in% names(dt))) {
  stop("Input must contain columns: locus, motif, count (case-insensitive)")
}

# ensure count is numeric
dt$count <- as.numeric(dt$count)

# Summary per locus (optionally per sample/hp if present)
group_cols <- c("locus")

summary_tbl <- dt %>%
  group_by(across(all_of(group_cols))) %>%
  summarize(
    n_alleles = sum(!is.na(count)),                            # n = number of alleles with size
    n_motifs = n_distinct(motif),
    motifs = paste(sort(unique(na.omit(motif))), collapse = ";"),
    n_distinct_counts = n_distinct(na.omit(count)),
    counts = paste(sort(unique(na.omit(as.character(count)))), collapse = ";"),
    differing = (n_motifs > 1 & n_distinct_counts > 1),
    .groups = "drop"
  ) %>%
  arrange(desc(differing), locus)

# also produce a locus-level rollup (across samples/hp)
locus_rollup <- summary_tbl %>%
  group_by(locus) %>%
  summarize(
    total_samples = n(),
    total_alleles = sum(n_alleles),
    n_motif_variants = n_distinct(unlist(strsplit(motifs, ";"))),
    n_distinct_counts_any = n_distinct(unlist(strsplit(counts, ";"))),
    any_differing = any(differing),
    motifs_all = paste(sort(unique(unlist(strsplit(motifs, ";")))), collapse = ";"),
    counts_all = paste(sort(unique(unlist(strsplit(counts, ";")))), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(flag = (n_motif_variants > 1 & n_distinct_counts_any > 1)) %>%
  arrange(desc(flag), locus)

# write outputs
fwrite(summary_tbl, file = sub("\\.tsv$","_per_sample.tsv", output_path), sep = "\t")
fwrite(locus_rollup, file = output_path, sep = "\t")

cat("Wrote:\n -", sub("\\.tsv$","_per_sample.tsv", output_path), "\n -", output_path, "\n")
