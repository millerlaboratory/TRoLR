library(dplyr)
library(data.table)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
karyotype <- if (length(args) >= 3) {
  val <- args[3]
  val <- gsub("\\r", "", val)   # remove CR if present
  val <- trimws(val)           # trim whitespace
  val <- toupper(val)          # uppercase
  if (nchar(val) == 0) "XX" else val
} else {
  "XX"
}
if (!karyotype %in% c("XX", "XY")) stop("karyotype must be 'XX' or 'XY'")

if (length(args) < 2) {
  stop("Usage: Rscript test_outliers.R <input_data.tsv> <output_outliers.tsv> [karyotype: XX|XY]")
}
input_path <- args[1]
output_path <- args[2]
ref_data <- data.table::fread(args[4])

data <- data.table::fread(input_path)

formatted_data <- data %>%
  dplyr::mutate(locus = paste0(V1, ":", V2, "-", V3)) %>%
  tidyr::separate(V4, into = c("sample_hp", "motif", "count"), sep = "\\|") %>%
  tidyr::separate(sample_hp, into = c("sample", "hp"), sep = "_") %>%
  dplyr::mutate(count = as.numeric(count)) %>%
  dplyr::select(locus, sample, hp, motif, count)

#ref_data <- data.table::fread("/n/users/sgibson/Projects/TANDEM_REPEATS/TRoLR/reference_data/vamos_asm_lps_control_summary.tsv") #need to change this path

joined <- formatted_data %>%
  dplyr::left_join(ref_data, by = c("locus", "motif"), relationship = "many-to-many")

if (karyotype == "XY") {
  # For XY samples, keep at most one allele per sample at chrX loci by taking the larger count
  chrX_XY <- joined %>%
    dplyr::filter(grepl("^chrX", locus)) %>%
    dplyr::group_by(locus, sample) %>%
    dplyr::slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(hp = "hp1")

  # drop those sample+locus pairs from full set (only remove matching locus+sample pairs)
  data_dropped <- joined %>%
    dplyr::anti_join(chrX_XY %>% dplyr::select(locus, sample), by = c("locus", "sample"))

  # combine
  data_final <- dplyr::bind_rows(data_dropped, chrX_XY)
} else {
  data_final <- joined
}

outliers <- data_final %>%
  dplyr::mutate(diff = as.numeric(count) - per99) %>%
  dplyr::filter((as.numeric(count) > per99) & diff >= 1) %>%
  dplyr::arrange(dplyr::desc(diff)) %>%
  dplyr::group_by(locus, sample, hp, motif, count, gene, feature, n_alleles, mean_lps, sd_lps, range_lps, per99, OMIM, phenotype, diff) %>%
  dplyr::distinct() %>%
  dplyr::ungroup()

# write to output path provided by caller
write.table(outliers, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
