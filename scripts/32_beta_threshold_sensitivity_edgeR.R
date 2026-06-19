# Step 32: beta-cell threshold sensitivity analysis
#
# Purpose:
# Test whether the main beta-cell pseudo-bulk DE result depends on the beta-cell count threshold.
#
# Thresholds tested:
#   n_beta_cells >= 50
#   n_beta_cells >= 100
#   n_beta_cells >= 200
#
# Main outputs:
#   results/32_beta_threshold_sensitivity/threshold_DE_summary.csv
#   results/32_beta_threshold_sensitivity/threshold_sample_summary.csv
#   results/32_beta_threshold_sensitivity/threshold_DEG_overlap_vs_100.csv
#   results/32_beta_threshold_sensitivity/<threshold>/<contrast>_edgeR_all.csv
#   results/32_beta_threshold_sensitivity/<threshold>/<contrast>_edgeR_FDR005.csv

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    BiocManager::install("edgeR", ask = FALSE, update = FALSE)
  }
  library(edgeR)
})

out_dir <- "results/32_beta_threshold_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- "results/06_beta_pseudobulk/beta_pseudobulk_counts_matrix.csv"
sample_file <- "results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv"

if (!file.exists(counts_file)) {
  stop("Missing count matrix: ", counts_file)
}

if (!file.exists(sample_file)) {
  stop("Missing sample summary: ", sample_file)
}

counts <- read.csv(
  counts_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
counts <- round(counts)

sample_info <- read.csv(
  sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

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

sample_col <- find_col(
  sample_info,
  c("sample_prefix", "sample", "Sample", "sample_id", "Sample_ID", "library", "LibraryID")
)

group_col <- find_col(
  sample_info,
  c("disease_group", "group", "Group", "disease", "Disease", "condition", "DiseaseState")
)

beta_col <- find_col(
  sample_info,
  c("n_beta_cells", "beta_cells", "Beta cells", "n_beta", "beta_cell_count", "selected_beta_cells")
)

sample_info <- sample_info[sample_info[[sample_col]] %in% colnames(counts), ]
sample_info <- sample_info[match(colnames(counts), sample_info[[sample_col]]), ]

if (!identical(as.character(sample_info[[sample_col]]), colnames(counts))) {
  stop("Sample metadata could not be aligned to count matrix columns.")
}

sample_info$sample_id_for_model <- as.character(sample_info[[sample_col]])
sample_info$disease_group_for_model <- as.character(sample_info[[group_col]])
sample_info$n_beta_cells_for_model <- as.numeric(sample_info[[beta_col]])

sample_info$disease_group_for_model <- toupper(sample_info$disease_group_for_model)
sample_info$disease_group_for_model <- gsub("NON.*DIABETIC|NON-DIABETIC|NODIABETIC", "ND", sample_info$disease_group_for_model)
sample_info$disease_group_for_model <- gsub("PRE.*DIABETIC|PRE-DIABETIC|PREDIABETIC", "PD", sample_info$disease_group_for_model)
sample_info$disease_group_for_model <- gsub("TYPE.*2.*DIABETES|TYPE 2 DIABETES", "T2D", sample_info$disease_group_for_model)

valid_groups <- c("ND", "PD", "T2D")
thresholds <- c(50, 100, 200)
contrasts_to_test <- list(
  T2D_vs_ND = c("T2D", "ND"),
  PD_vs_ND = c("PD", "ND"),
  T2D_vs_PD = c("T2D", "PD")
)

run_edgeR_threshold <- function(threshold) {
  label <- paste0("beta_min", threshold)
  threshold_dir <- file.path(out_dir, label)
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)

  keep_samples <- sample_info$disease_group_for_model %in% valid_groups &
    !is.na(sample_info$n_beta_cells_for_model) &
    sample_info$n_beta_cells_for_model >= threshold

  meta <- sample_info[keep_samples, ]
  count_mat <- counts[, meta$sample_id_for_model, drop = FALSE]
  group <- factor(meta$disease_group_for_model, levels = c("ND", "PD", "T2D"))

  sample_summary <- data.frame(
    threshold = threshold,
    n_samples_total = nrow(meta),
    n_ND = sum(group == "ND"),
    n_PD = sum(group == "PD"),
    n_T2D = sum(group == "T2D"),
    min_beta_cells = min(meta$n_beta_cells_for_model, na.rm = TRUE),
    median_beta_cells = median(meta$n_beta_cells_for_model, na.rm = TRUE),
    max_beta_cells = max(meta$n_beta_cells_for_model, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  write.csv(
    meta,
    file.path(threshold_dir, paste0(label, "_samples_used.csv")),
    row.names = FALSE
  )

  y <- DGEList(counts = count_mat, group = group)
  keep_genes <- filterByExpr(y, group = group)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)

  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  write.csv(
    logcpm,
    file.path(threshold_dir, paste0(label, "_logCPM_TMM.csv")),
    row.names = TRUE
  )

  contrast_summary_list <- list()
  deg_sets <- list()

  for (contrast_name in names(contrasts_to_test)) {
    first_group <- contrasts_to_test[[contrast_name]][1]
    second_group <- contrasts_to_test[[contrast_name]][2]

    contrast_vec <- rep(0, ncol(design))
    names(contrast_vec) <- colnames(design)
    contrast_vec[first_group] <- 1
    contrast_vec[second_group] <- -1

    qlf <- glmQLFTest(fit, contrast = contrast_vec)
    tab <- topTags(qlf, n = Inf)$table
    tab$gene <- rownames(tab)
    tab <- tab[, c("gene", setdiff(colnames(tab), "gene"))]
    tab$direction <- ifelse(
      tab$logFC > 0,
      paste0("higher_in_", first_group),
      paste0("lower_in_", first_group)
    )
    tab <- tab[order(tab$FDR, tab$PValue), ]

    write.csv(tab, file.path(threshold_dir, paste0(contrast_name, "_edgeR_all.csv")), row.names = FALSE)

    sig005 <- tab[tab$FDR < 0.05, ]
    sig010 <- tab[tab$FDR < 0.10, ]

    write.csv(sig005, file.path(threshold_dir, paste0(contrast_name, "_edgeR_FDR005.csv")), row.names = FALSE)
    write.csv(sig010, file.path(threshold_dir, paste0(contrast_name, "_edgeR_FDR010.csv")), row.names = FALSE)

    deg_sets[[contrast_name]] <- sig005$gene

    contrast_summary_list[[contrast_name]] <- data.frame(
      threshold = threshold,
      contrast = contrast_name,
      n_samples = nrow(meta),
      n_ND = sum(group == "ND"),
      n_PD = sum(group == "PD"),
      n_T2D = sum(group == "T2D"),
      n_tested_genes = nrow(tab),
      n_FDR005 = nrow(sig005),
      n_up_FDR005 = sum(sig005$logFC > 0),
      n_down_FDR005 = sum(sig005$logFC < 0),
      n_FDR010 = nrow(sig010),
      top_gene = tab$gene[1],
      top_logFC = tab$logFC[1],
      top_P_value = tab$PValue[1],
      top_FDR = tab$FDR[1],
      top_direction = tab$direction[1],
      stringsAsFactors = FALSE
    )
  }

  contrast_summary <- do.call(rbind, contrast_summary_list)
  write.csv(
    contrast_summary,
    file.path(threshold_dir, paste0(label, "_DE_summary.csv")),
    row.names = FALSE
  )

  list(
    threshold = threshold,
    sample_summary = sample_summary,
    contrast_summary = contrast_summary,
    deg_sets = deg_sets
  )
}

all_results <- lapply(thresholds, run_edgeR_threshold)
names(all_results) <- paste0("beta_min", thresholds)

sample_summary_all <- do.call(rbind, lapply(all_results, function(x) x$sample_summary))
de_summary_all <- do.call(rbind, lapply(all_results, function(x) x$contrast_summary))

write.csv(sample_summary_all, file.path(out_dir, "threshold_sample_summary.csv"), row.names = FALSE)
write.csv(de_summary_all, file.path(out_dir, "threshold_DE_summary.csv"), row.names = FALSE)

overlap_rows <- list()
ref <- all_results[["beta_min100"]]$deg_sets

for (threshold in thresholds) {
  cur <- all_results[[paste0("beta_min", threshold)]]$deg_sets

  for (contrast_name in names(contrasts_to_test)) {
    a <- unique(ref[[contrast_name]])
    b <- unique(cur[[contrast_name]])
    if (is.null(a)) a <- character()
    if (is.null(b)) b <- character()

    inter <- intersect(a, b)
    uni <- union(a, b)

    overlap_rows[[paste(threshold, contrast_name, sep = "_")]] <- data.frame(
      comparison = paste0("beta_min", threshold, "_vs_beta_min100"),
      threshold = threshold,
      contrast = contrast_name,
      n_DEG_threshold100 = length(a),
      n_DEG_current_threshold = length(b),
      n_overlap = length(inter),
      jaccard = ifelse(length(uni) > 0, length(inter) / length(uni), NA),
      fraction_of_threshold100_recovered = ifelse(length(a) > 0, length(inter) / length(a), NA),
      fraction_of_current_recovered_in_threshold100 = ifelse(length(b) > 0, length(inter) / length(b), NA),
      stringsAsFactors = FALSE
    )
  }
}

overlap_df <- do.call(rbind, overlap_rows)
write.csv(overlap_df, file.path(out_dir, "threshold_DEG_overlap_vs_100.csv"), row.names = FALSE)

stable_list <- list()

for (contrast_name in names(contrasts_to_test)) {
  sets <- lapply(all_results, function(x) unique(x$deg_sets[[contrast_name]]))
  sets <- lapply(sets, function(x) if (is.null(x)) character() else x)
  stable <- Reduce(intersect, sets)
  union_all <- Reduce(union, sets)

  stable_list[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    n_stable_DEGs_all_thresholds = length(stable),
    n_union_DEGs_any_threshold = length(union_all),
    stable_fraction_of_union = ifelse(length(union_all) > 0, length(stable) / length(union_all), NA),
    stable_genes = paste(head(stable, 50), collapse = ";"),
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(gene = stable),
    file.path(out_dir, paste0(contrast_name, "_stable_DEGs_all_thresholds.csv")),
    row.names = FALSE
  )
}

stable_summary <- do.call(rbind, stable_list)
write.csv(stable_summary, file.path(out_dir, "stable_DEGs_all_thresholds_summary.csv"), row.names = FALSE)

print(sample_summary_all)
print(de_summary_all)
print(overlap_df)
print(stable_summary)

cat("\nDone: beta-cell threshold sensitivity analysis completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
