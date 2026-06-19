# Step 35: beta-cell identity audit and optional CELLxGENE/author-label comparison
#
# Purpose:
# 1. Check whether the selected beta-cell pseudo-bulk profiles show high beta markers.
# 2. Check whether non-beta markers are low.
# 3. Optionally compare your beta-cell counts with an external CELLxGENE/author metadata file.
#
# Required input:
#   results/06_beta_pseudobulk/beta_pseudobulk_counts_matrix.csv
#   results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv
#
# Optional input:
#   data/reference/CELLxGENE/cell_metadata.csv
#
# Outputs:
#   results/35_beta_identity_audit/beta_identity_marker_logCPM_by_sample.csv
#   results/35_beta_identity_audit/beta_identity_marker_logCPM_by_group.csv
#   results/35_beta_identity_audit/beta_identity_score_summary_by_sample.csv
#   results/35_beta_identity_audit/non_beta_marker_contamination_summary.csv
#   results/35_beta_identity_audit/optional_cellxgene_beta_count_comparison.csv

out_dir <- "results/35_beta_identity_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- "results/06_beta_pseudobulk/beta_pseudobulk_counts_matrix.csv"
sample_file <- "results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv"

if (!file.exists(counts_file)) stop("Missing file: ", counts_file)
if (!file.exists(sample_file)) stop("Missing file: ", sample_file)

find_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) return(hit[1])
  lower_names <- tolower(colnames(df))
  lower_candidates <- tolower(candidates)
  idx <- match(lower_candidates, lower_names)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(colnames(df)[idx[1]])
  if (required) stop("Cannot find required column. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

safe_genes <- function(gene_set, available_genes) {
  gene_set <- toupper(gene_set)
  available_upper <- toupper(available_genes)
  available_genes[available_upper %in% gene_set]
}

mean_marker <- function(mat, genes) {
  genes2 <- safe_genes(genes, rownames(mat))
  if (length(genes2) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[genes2, , drop = FALSE], na.rm = TRUE)
}

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
counts <- round(counts)

sample_info <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)

sample_col <- find_col(sample_info, c("sample_prefix", "sample", "Sample", "sample_id", "Sample_ID", "library", "LibraryID"))
group_col <- find_col(sample_info, c("disease_group", "group", "Group", "disease", "Disease", "condition", "DiseaseState"), required = FALSE)
beta_col <- find_col(sample_info, c("n_beta_cells", "beta_cells", "Beta cells", "n_beta", "beta_cell_count", "selected_beta_cells"), required = FALSE)

sample_info <- sample_info[sample_info[[sample_col]] %in% colnames(counts), ]
sample_info <- sample_info[match(colnames(counts), sample_info[[sample_col]]), ]

if (!identical(as.character(sample_info[[sample_col]]), colnames(counts))) {
  stop("Sample metadata could not be aligned with count matrix columns.")
}

lib_size <- colSums(counts)
logcpm <- log2(t(t(counts) / lib_size * 1e6) + 1)

marker_sets <- list(
  beta = c("INS", "IAPP", "PCSK1", "PCSK2", "MAFA", "PDX1", "NKX6-1", "NKX6_1", "SLC30A8"),
  alpha = c("GCG", "ARX"),
  delta = c("SST"),
  pp = c("PPY"),
  epsilon = c("GHRL"),
  acinar = c("PRSS1", "CPA1", "CTRB2", "CELA3A", "CLPS"),
  ductal = c("KRT19", "KRT8", "KRT18", "CFTR"),
  endothelial = c("PECAM1", "VWF", "KDR"),
  immune = c("PTPRC", "LST1", "C1QA"),
  stellate = c("COL1A1", "COL1A2", "DCN", "LUM")
)

marker_by_sample <- data.frame(sample = colnames(logcpm), stringsAsFactors = FALSE)

if (!is.na(group_col)) marker_by_sample$group <- as.character(sample_info[[group_col]])
if (!is.na(beta_col)) marker_by_sample$n_beta_cells <- as.numeric(sample_info[[beta_col]])

for (set_name in names(marker_sets)) {
  marker_by_sample[[paste0(set_name, "_marker_mean_logCPM")]] <- mean_marker(logcpm, marker_sets[[set_name]])
}

marker_by_sample$beta_to_alpha_delta_ratio <- marker_by_sample$beta_marker_mean_logCPM - pmax(
  marker_by_sample$alpha_marker_mean_logCPM,
  marker_by_sample$delta_marker_mean_logCPM,
  na.rm = TRUE
)

marker_by_sample$beta_to_exocrine_ratio <- marker_by_sample$beta_marker_mean_logCPM - pmax(
  marker_by_sample$acinar_marker_mean_logCPM,
  marker_by_sample$ductal_marker_mean_logCPM,
  na.rm = TRUE
)

marker_by_sample$beta_to_non_endocrine_ratio <- marker_by_sample$beta_marker_mean_logCPM - pmax(
  marker_by_sample$immune_marker_mean_logCPM,
  marker_by_sample$endothelial_marker_mean_logCPM,
  marker_by_sample$stellate_marker_mean_logCPM,
  na.rm = TRUE
)

write.csv(marker_by_sample, file.path(out_dir, "beta_identity_marker_logCPM_by_sample.csv"), row.names = FALSE)

if (!is.na(group_col)) {
  numeric_cols <- sapply(marker_by_sample, is.numeric)
  marker_by_group <- aggregate(
    marker_by_sample[, numeric_cols, drop = FALSE],
    by = list(group = marker_by_sample$group),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  write.csv(marker_by_group, file.path(out_dir, "beta_identity_marker_logCPM_by_group.csv"), row.names = FALSE)
}

score_summary <- data.frame(
  metric = c(
    "median_beta_marker_mean_logCPM",
    "median_alpha_marker_mean_logCPM",
    "median_delta_marker_mean_logCPM",
    "median_acinar_marker_mean_logCPM",
    "median_ductal_marker_mean_logCPM",
    "median_immune_marker_mean_logCPM",
    "median_endothelial_marker_mean_logCPM",
    "median_beta_to_alpha_delta_ratio",
    "median_beta_to_exocrine_ratio",
    "median_beta_to_non_endocrine_ratio"
  ),
  value = c(
    median(marker_by_sample$beta_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$alpha_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$delta_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$acinar_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$ductal_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$immune_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$endothelial_marker_mean_logCPM, na.rm = TRUE),
    median(marker_by_sample$beta_to_alpha_delta_ratio, na.rm = TRUE),
    median(marker_by_sample$beta_to_exocrine_ratio, na.rm = TRUE),
    median(marker_by_sample$beta_to_non_endocrine_ratio, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(score_summary, file.path(out_dir, "beta_identity_score_summary_by_sample.csv"), row.names = FALSE)

contam <- marker_by_sample
contam$possible_alpha_delta_signal <- contam$beta_to_alpha_delta_ratio < 1
contam$possible_exocrine_signal <- contam$beta_to_exocrine_ratio < 1
contam$possible_non_endocrine_signal <- contam$beta_to_non_endocrine_ratio < 1

contam_summary <- data.frame(
  metric = c(
    "n_samples",
    "n_possible_alpha_delta_signal",
    "n_possible_exocrine_signal",
    "n_possible_non_endocrine_signal"
  ),
  value = c(
    nrow(contam),
    sum(contam$possible_alpha_delta_signal, na.rm = TRUE),
    sum(contam$possible_exocrine_signal, na.rm = TRUE),
    sum(contam$possible_non_endocrine_signal, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(contam_summary, file.path(out_dir, "non_beta_marker_contamination_summary.csv"), row.names = FALSE)

external_file <- "data/reference/CELLxGENE/cell_metadata.csv"

if (file.exists(external_file)) {
  ext <- read.csv(external_file, stringsAsFactors = FALSE, check.names = FALSE)
  ext_sample_col <- find_col(ext, c("sample", "sample_id", "Sample", "Sample_ID", "library", "LibraryID", "sample_prefix"))
  ext_celltype_col <- find_col(ext, c("cell_type", "celltype", "CellType", "annotation", "author_cell_type", "cell_type_original", "cell_type_ontology_term_label"))

  ext$sample_for_compare <- as.character(ext[[ext_sample_col]])
  ext$celltype_for_compare <- tolower(as.character(ext[[ext_celltype_col]]))
  ext$is_beta_external <- grepl("beta|β", ext$celltype_for_compare, ignore.case = TRUE)

  ext_counts <- aggregate(
    ext$is_beta_external,
    by = list(sample = ext$sample_for_compare),
    FUN = function(x) sum(x, na.rm = TRUE)
  )
  colnames(ext_counts)[2] <- "external_beta_cells"

  our_counts <- data.frame(
    sample = marker_by_sample$sample,
    our_beta_cells = if (!is.na(beta_col)) marker_by_sample$n_beta_cells else NA_real_,
    stringsAsFactors = FALSE
  )

  cmp <- merge(our_counts, ext_counts, by = "sample", all = TRUE)
  cmp$difference_external_minus_ours <- cmp$external_beta_cells - cmp$our_beta_cells
  cmp$ratio_external_to_ours <- cmp$external_beta_cells / cmp$our_beta_cells

  write.csv(cmp, file.path(out_dir, "optional_cellxgene_beta_count_comparison.csv"), row.names = FALSE)
  cat("\nOptional external metadata comparison completed.\n")
} else {
  cat("\nNo external CELLxGENE/author metadata file found at: ", external_file, "\n", sep = "")
  cat("Place a CSV there if you want sample-level beta-cell count comparison.\n")
}

cat("\nMarker identity audit summary:\n")
print(score_summary)

cat("\nPossible contamination summary:\n")
print(contam_summary)

cat("\nDone: beta-cell identity audit completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
