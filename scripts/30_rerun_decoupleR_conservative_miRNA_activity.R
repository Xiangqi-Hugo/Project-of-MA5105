# Conservative decoupleR re-run for miRNA activity.
# This uses the same beta-cell pseudo-bulk expression matrix, but filters the miRNA target network
# by target-set size before activity inference.
#
# Filters:
#   min_targets = 15
#   max_targets = 200, 300, 500
#
# Method:
#   decoupleR::run_wmean with mor = -1 for miRNA repression.
#   limma tests predicted miRNA activity across ND, PD, and T2D groups.
#   BH global FDR is recalculated across all miRNA activity tests and all three contrasts
#   within each target-size filter.

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    BiocManager::install("decoupleR", ask = FALSE, update = FALSE)
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma", ask = FALSE, update = FALSE)
  }
  library(decoupleR)
  library(limma)
})

strict_dir <- "results/14_strict_miRNA_target_GSEA"
out_dir <- "results/30_decoupleR_conservative_miRNA_activity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"
sample_file <- "results/07_beta_DE_edgeR/samples_used_for_DE.csv"
term2gene_file <- file.path(strict_dir, "strict_miRNA_TERM2GENE.csv")

if (!file.exists(expr_file)) stop("Missing expression file: ", expr_file)
if (!file.exists(sample_file)) stop("Missing sample file: ", sample_file)
if (!file.exists(term2gene_file)) stop("Missing TERM2GENE file: ", term2gene_file)

# --------------------------
# Read expression and metadata
# --------------------------
expr <- read.csv(
  expr_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

sample_info <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)

# Robust metadata columns.
sample_col <- intersect(
  c("sample_prefix", "sample", "Sample", "sample_id", "LibraryID"),
  colnames(sample_info)
)[1]

group_col <- intersect(
  c("disease_group", "group", "disease", "Disease", "condition"),
  colnames(sample_info)
)[1]

if (is.na(sample_col)) stop("Cannot find sample ID column in samples_used_for_DE.csv")
if (is.na(group_col)) stop("Cannot find disease group column in samples_used_for_DE.csv")

sample_info <- sample_info[match(colnames(expr), sample_info[[sample_col]]), ]

if (any(is.na(sample_info[[sample_col]]))) {
  stop("Some expression columns are missing from sample metadata.")
}

sample_info$disease_group_for_model <- as.character(sample_info[[group_col]])
sample_info$disease_group_for_model <- factor(
  sample_info$disease_group_for_model,
  levels = c("ND", "PD", "T2D")
)

if (any(is.na(sample_info$disease_group_for_model))) {
  stop("Disease groups must include only ND, PD, and T2D for this model.")
}

# --------------------------
# Read and clean TERM2GENE
# --------------------------
term2gene <- read.csv(term2gene_file, stringsAsFactors = FALSE, check.names = FALSE)

cn <- colnames(term2gene)
term_col <- cn[tolower(cn) %in% c("term", "id", "mirna", "mirna_id", "source", "gs_name", "pathway")][1]
gene_col <- cn[tolower(cn) %in% c("gene", "target", "target_gene", "gene_symbol", "symbol")][1]

if (is.na(term_col)) term_col <- cn[1]
if (is.na(gene_col)) gene_col <- cn[2]

term2gene <- term2gene[, c(term_col, gene_col)]
colnames(term2gene) <- c("source", "target")

term2gene <- term2gene[!is.na(term2gene$source) & !is.na(term2gene$target), ]
term2gene <- term2gene[term2gene$source != "" & term2gene$target != "", ]
term2gene <- unique(term2gene)

# Keep only genes present in expression matrix.
term2gene <- term2gene[term2gene$target %in% rownames(expr), ]

target_sizes <- aggregate(
  target ~ source,
  data = term2gene,
  FUN = function(x) length(unique(x))
)
colnames(target_sizes) <- c("source", "n_targets_in_expr")

write.csv(
  target_sizes[order(target_sizes$n_targets_in_expr, decreasing = TRUE), ],
  file.path(out_dir, "miRNA_target_set_sizes_after_gene_matching.csv"),
  row.names = FALSE
)

primary_terms <- c(
  "MIR195_5P",
  "MIR16_5P",
  "MIR15A_5P",
  "MIR15B_5P",
  "MIR649",
  "MIR6838_5P"
)

# --------------------------
# decoupleR output parser
# --------------------------
activity_from_decoupleR <- function(act_long) {
  act_long <- as.data.frame(act_long)

  # Some decoupleR versions return columns: statistic/source/condition/score.
  if ("statistic" %in% colnames(act_long)) {
    stat_values <- unique(as.character(act_long$statistic))
    preferred <- c("norm_wmean", "wmean", "corr_wmean", "score")
    chosen <- intersect(preferred, stat_values)[1]
    if (is.na(chosen)) chosen <- stat_values[1]
    act_long <- act_long[as.character(act_long$statistic) == chosen, , drop = FALSE]
    score_col <- "score"
  } else {
    candidates <- c("norm_wmean", "score", "wmean", "corr_wmean")
    score_col <- intersect(candidates, colnames(act_long))[1]
  }

  if (is.na(score_col) || !(score_col %in% colnames(act_long))) {
    stop("Cannot detect activity score column from decoupleR output.")
  }

  source_col <- intersect(c("source", "regulator", "miRNA", "ID"), colnames(act_long))[1]
  condition_col <- intersect(c("condition", "sample", "Sample", "sample_id"), colnames(act_long))[1]

  if (is.na(source_col)) stop("Cannot detect source column from decoupleR output.")
  if (is.na(condition_col)) stop("Cannot detect condition/sample column from decoupleR output.")

  sources <- unique(as.character(act_long[[source_col]]))
  conditions <- unique(as.character(act_long[[condition_col]]))

  mat <- matrix(
    NA_real_,
    nrow = length(sources),
    ncol = length(conditions),
    dimnames = list(sources, conditions)
  )

  idx <- cbind(
    match(as.character(act_long[[source_col]]), sources),
    match(as.character(act_long[[condition_col]]), conditions)
  )

  mat[idx] <- as.numeric(act_long[[score_col]])

  mat
}

# --------------------------
# Limma testing function
# --------------------------
run_limma_contrasts <- function(activity_mat, filter_label, filter_out_dir) {
  activity_mat <- activity_mat[, colnames(expr), drop = FALSE]

  design <- model.matrix(~ 0 + sample_info$disease_group_for_model)
  colnames(design) <- levels(sample_info$disease_group_for_model)

  fit <- lmFit(activity_mat, design)

  contrast_matrix <- makeContrasts(
    T2D_vs_ND = T2D - ND,
    PD_vs_ND = PD - ND,
    T2D_vs_PD = T2D - PD,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  contrast_names <- colnames(contrast_matrix)
  all_list <- list()
  within_summary_list <- list()

  for (contrast in contrast_names) {
    tab <- topTable(fit2, coef = contrast, number = Inf, sort.by = "P")
    tab$miRNA <- rownames(tab)
    tab <- tab[, c("miRNA", setdiff(colnames(tab), "miRNA"))]

    first_group <- sub("_vs_.*$", "", contrast)
    tab$activity_direction <- ifelse(
      tab$logFC > 0,
      paste0("higher_in_", first_group),
      paste0("lower_in_", first_group)
    )

    tab$contrast <- contrast
    tab$filter <- filter_label

    write.csv(
      tab,
      file.path(filter_out_dir, paste0(contrast, "_decoupleR_activity_all.csv")),
      row.names = FALSE
    )

    write.csv(
      tab[tab$adj.P.Val < 0.05, ],
      file.path(filter_out_dir, paste0(contrast, "_decoupleR_activity_FDR005.csv")),
      row.names = FALSE
    )

    write.csv(
      tab[tab$adj.P.Val < 0.10, ],
      file.path(filter_out_dir, paste0(contrast, "_decoupleR_activity_FDR010.csv")),
      row.names = FALSE
    )

    within_summary_list[[contrast]] <- data.frame(
      filter = filter_label,
      contrast = contrast,
      n_miRNAs_tested = nrow(tab),
      n_within_FDR005 = sum(tab$adj.P.Val < 0.05),
      n_within_FDR010 = sum(tab$adj.P.Val < 0.10),
      top_miRNA = tab$miRNA[1],
      top_logFC = tab$logFC[1],
      top_P_value = tab$P.Value[1],
      top_within_FDR = tab$adj.P.Val[1],
      top_direction = tab$activity_direction[1],
      stringsAsFactors = FALSE
    )

    all_list[[contrast]] <- tab
  }

  all_tests <- do.call(rbind, all_list)
  all_tests$global_FDR_across_all_contrasts <- p.adjust(all_tests$P.Value, method = "BH")
  all_tests <- all_tests[order(all_tests$global_FDR_across_all_contrasts, all_tests$P.Value), ]

  write.csv(
    all_tests,
    file.path(filter_out_dir, "decoupleR_all_contrasts_global_FDR.csv"),
    row.names = FALSE
  )

  global_summary <- do.call(
    rbind,
    lapply(contrast_names, function(contrast) {
      x <- all_tests[all_tests$contrast == contrast, ]
      x <- x[order(x$global_FDR_across_all_contrasts, x$P.Value), ]

      data.frame(
        filter = filter_label,
        contrast = contrast,
        n_tests = nrow(x),
        n_global_FDR005 = sum(x$global_FDR_across_all_contrasts < 0.05),
        n_global_FDR010 = sum(x$global_FDR_across_all_contrasts < 0.10),
        top_miRNA = x$miRNA[1],
        top_logFC = x$logFC[1],
        top_P_value = x$P.Value[1],
        top_within_FDR = x$adj.P.Val[1],
        top_global_FDR = x$global_FDR_across_all_contrasts[1],
        top_direction = x$activity_direction[1],
        stringsAsFactors = FALSE
      )
    })
  )

  write.csv(
    global_summary,
    file.path(filter_out_dir, "decoupleR_global_FDR_summary.csv"),
    row.names = FALSE
  )

  primary <- all_tests[all_tests$miRNA %in% primary_terms, ]

  primary <- merge(
    primary,
    target_sizes,
    by.x = "miRNA",
    by.y = "source",
    all.x = TRUE
  )

  primary <- primary[order(primary$contrast, primary$global_FDR_across_all_contrasts, primary$P.Value), ]

  write.csv(
    primary,
    file.path(filter_out_dir, "primary_candidates_decoupleR_conservative_network.csv"),
    row.names = FALSE
  )

  list(
    within_summary = do.call(rbind, within_summary_list),
    global_summary = global_summary,
    primary = primary,
    all_tests = all_tests
  )
}

# --------------------------
# Main loop
# --------------------------
min_targets <- 15
max_targets_vec <- c(200, 300, 500)

overall_rows <- list()
global_summary_list <- list()
within_summary_list <- list()
primary_list <- list()

for (max_targets in max_targets_vec) {
  filter_label <- paste0("min", min_targets, "_max", max_targets)
  filter_out_dir <- file.path(out_dir, filter_label)
  dir.create(filter_out_dir, recursive = TRUE, showWarnings = FALSE)

  keep_sources <- target_sizes$source[
    target_sizes$n_targets_in_expr >= min_targets &
      target_sizes$n_targets_in_expr <= max_targets
  ]

  net <- term2gene[term2gene$source %in% keep_sources, ]
  net$mor <- -1

  write.csv(
    net,
    file.path(filter_out_dir, paste0("network_", filter_label, ".csv")),
    row.names = FALSE
  )

  network_summary <- data.frame(
    filter = filter_label,
    min_targets = min_targets,
    max_targets = max_targets,
    n_miRNA_sources = length(unique(net$source)),
    n_network_edges = nrow(net),
    n_unique_target_genes = length(unique(net$target)),
    n_primary_terms_in_network = sum(primary_terms %in% unique(net$source)),
    primary_terms_in_network = paste(primary_terms[primary_terms %in% unique(net$source)], collapse = ";"),
    stringsAsFactors = FALSE
  )

  write.csv(
    network_summary,
    file.path(filter_out_dir, "network_summary.csv"),
    row.names = FALSE
  )

  print(network_summary)

  message("Running decoupleR::run_wmean for ", filter_label)

  set.seed(123)

  act_long <- decoupleR::run_wmean(
    mat = expr,
    net = net,
    .source = "source",
    .target = "target",
    .mor = "mor",
    times = 100,
    minsize = min_targets
  )

  activity_mat <- activity_from_decoupleR(act_long)

  # Keep only samples used in DE and in the correct order.
  activity_mat <- activity_mat[, colnames(expr), drop = FALSE]

  write.csv(
    activity_mat,
    file.path(filter_out_dir, "decoupleR_activity_matrix.csv"),
    row.names = TRUE
  )

  saveRDS(
    activity_mat,
    file.path(filter_out_dir, "decoupleR_activity_matrix.rds")
  )

  res <- run_limma_contrasts(activity_mat, filter_label, filter_out_dir)

  global_summary_list[[filter_label]] <- res$global_summary
  within_summary_list[[filter_label]] <- res$within_summary
  primary_list[[filter_label]] <- res$primary

  overall_rows[[filter_label]] <- network_summary
}

overall_network_summary <- do.call(rbind, overall_rows)
overall_global_summary <- do.call(rbind, global_summary_list)
overall_within_summary <- do.call(rbind, within_summary_list)
overall_primary <- do.call(rbind, primary_list)

write.csv(
  overall_network_summary,
  file.path(out_dir, "conservative_decoupleR_network_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_within_summary,
  file.path(out_dir, "conservative_decoupleR_within_contrast_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_global_summary,
  file.path(out_dir, "conservative_decoupleR_global_FDR_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_primary,
  file.path(out_dir, "conservative_decoupleR_primary_candidates_all_filters.csv"),
  row.names = FALSE
)

print(overall_global_summary)
print(overall_primary)

cat("\nDone: conservative decoupleR miRNA activity re-run completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
