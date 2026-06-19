# Conservative re-filtering of strict miRNA target GSEA results.
# This script does not claim that the original GSEA is wrong.
# It tests whether the significant signal is driven by very large and/or redundant miRNA target sets.

strict_dir <- "results/14_strict_miRNA_target_GSEA"
audit_dir <- "results/28_miRNA_GSEA_audit"
out_dir <- "results/29_conservative_miRNA_GSEA_filtering"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

term2gene_file <- file.path(strict_dir, "strict_miRNA_TERM2GENE.csv")

if (!file.exists(term2gene_file)) {
  stop("Cannot find strict_miRNA_TERM2GENE.csv")
}

term2gene <- read.csv(term2gene_file, stringsAsFactors = FALSE, check.names = FALSE)

detect_term_gene_cols <- function(x) {
  cn <- colnames(x)
  term_col <- cn[tolower(cn) %in% c("term", "id", "mirna", "mirna_id", "gs_name", "pathway")][1]
  gene_col <- cn[tolower(cn) %in% c("gene", "target", "target_gene", "gene_symbol", "symbol")][1]

  if (is.na(term_col)) term_col <- cn[1]
  if (is.na(gene_col)) gene_col <- cn[2]

  c(term_col, gene_col)
}

tg_cols <- detect_term_gene_cols(term2gene)

term2gene <- term2gene[, tg_cols]
colnames(term2gene) <- c("ID", "gene")

term2gene <- term2gene[!is.na(term2gene$ID) & !is.na(term2gene$gene), ]
term2gene <- unique(term2gene)

target_list <- split(term2gene$gene, term2gene$ID)
target_list <- lapply(target_list, unique)

target_size <- data.frame(
  ID = names(target_list),
  n_targets = vapply(target_list, length, numeric(1)),
  stringsAsFactors = FALSE
)

primary_terms <- c(
  "MIR195_5P",
  "MIR16_5P",
  "MIR15A_5P",
  "MIR15B_5P",
  "MIR649",
  "MIR6838_5P"
)

primary_size <- target_size[target_size$ID %in% primary_terms, ]
primary_size$passes_min15_max200 <- primary_size$n_targets >= 15 & primary_size$n_targets <= 200
primary_size$passes_min15_max300 <- primary_size$n_targets >= 15 & primary_size$n_targets <= 300
primary_size$passes_min15_max500 <- primary_size$n_targets >= 15 & primary_size$n_targets <= 500

write.csv(
  primary_size[order(primary_size$n_targets, decreasing = TRUE), ],
  file.path(out_dir, "primary_candidate_target_size_filter_status.csv"),
  row.names = FALSE
)

find_all_gsea_file <- function(contrast) {
  files <- list.files(strict_dir, pattern = "\\.csv$", full.names = TRUE)
  base <- basename(files)

  keep <- grepl(contrast, base, fixed = TRUE)
  keep <- keep & grepl("GSEA", base, ignore.case = TRUE)
  keep <- keep & !grepl("FDR005|FDR010|summary|TERM2GENE|target_set|pairwise|overlap", base, ignore.case = TRUE)

  candidates <- files[keep]

  # Prefer names containing all or full.
  preferred <- candidates[grepl("all|full", basename(candidates), ignore.case = TRUE)]
  if (length(preferred) > 0) return(preferred[1])

  if (length(candidates) > 0) return(candidates[1])

  return(NA_character_)
}

detect_id_col <- function(x) {
  cn <- colnames(x)
  id_col <- cn[tolower(cn) %in% c("id", "term", "pathway", "gs_name", "mirna")][1]
  if (is.na(id_col)) id_col <- cn[1]
  id_col
}

detect_p_col <- function(x) {
  cn <- colnames(x)
  p_col <- cn[tolower(cn) %in% c("pval", "pvalue", "p.value", "p_value", "p")][1]
  p_col
}

detect_fdr_col <- function(x) {
  cn <- colnames(x)
  fdr_col <- cn[tolower(cn) %in% c("padj", "p.adjust", "p_adj", "fdr", "qvalue", "q_value")][1]
  fdr_col
}

detect_nes_col <- function(x) {
  cn <- colnames(x)
  nes_col <- cn[tolower(cn) %in% c("nes", "normalized_enrichment_score")][1]
  nes_col
}

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

collapse_redundant <- function(df, fdr_col = "conservative_FDR", threshold = 0.25) {
  if (nrow(df) == 0) return(df)

  df <- df[order(df[[fdr_col]], abs(df$NES), decreasing = c(FALSE, TRUE)), ]

  kept <- character()
  keep_idx <- logical(nrow(df))

  for (i in seq_len(nrow(df))) {
    term <- df$ID[i]

    if (length(kept) == 0) {
      keep_idx[i] <- TRUE
      kept <- c(kept, term)
    } else {
      overlaps <- vapply(
        kept,
        function(k) jaccard(target_list[[term]], target_list[[k]]),
        numeric(1)
      )

      if (all(is.na(overlaps)) || max(overlaps, na.rm = TRUE) <= threshold) {
        keep_idx[i] <- TRUE
        kept <- c(kept, term)
      }
    }
  }

  df[keep_idx, , drop = FALSE]
}

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")
max_sizes <- c(200, 300, 500)

summary_rows <- list()
missing_files <- character()

for (contrast in contrasts) {
  gsea_file <- find_all_gsea_file(contrast)

  if (is.na(gsea_file) || !file.exists(gsea_file)) {
    missing_files <- c(missing_files, contrast)
    next
  }

  message("Reading ", contrast, " GSEA file: ", gsea_file)

  res <- read.csv(gsea_file, stringsAsFactors = FALSE, check.names = FALSE)

  id_col <- detect_id_col(res)
  p_col <- detect_p_col(res)
  fdr_col <- detect_fdr_col(res)
  nes_col <- detect_nes_col(res)

  if (is.na(nes_col)) {
    stop("Cannot find NES column in ", gsea_file)
  }

  res$ID <- res[[id_col]]
  res$NES <- res[[nes_col]]

  if (!is.na(p_col)) {
    res$P_value_for_refilter <- res[[p_col]]
  } else if (!is.na(fdr_col)) {
    warning("No raw p-value column found for ", contrast, ". Using existing FDR only. This is less ideal.")
    res$P_value_for_refilter <- NA_real_
  } else {
    stop("Cannot find p-value or FDR column in ", gsea_file)
  }

  if (!is.na(fdr_col)) {
    res$original_FDR <- res[[fdr_col]]
  } else {
    res$original_FDR <- NA_real_
  }

  res <- merge(res, target_size, by = "ID", all.x = TRUE)

  write.csv(
    res,
    file.path(out_dir, paste0(contrast, "_GSEA_all_with_target_sizes.csv")),
    row.names = FALSE
  )

  for (max_size in max_sizes) {
    filtered <- res[
      !is.na(res$n_targets) &
        res$n_targets >= 15 &
        res$n_targets <= max_size,
    ]

    if (!all(is.na(filtered$P_value_for_refilter))) {
      filtered$conservative_FDR <- p.adjust(filtered$P_value_for_refilter, method = "BH")
    } else {
      filtered$conservative_FDR <- filtered$original_FDR
    }

    filtered <- filtered[order(filtered$conservative_FDR, filtered$P_value_for_refilter), ]

    fdr005 <- filtered[filtered$conservative_FDR < 0.05, ]
    fdr010 <- filtered[filtered$conservative_FDR < 0.10, ]

    collapsed_025 <- collapse_redundant(fdr005, threshold = 0.25)
    collapsed_050 <- collapse_redundant(fdr005, threshold = 0.50)

    tag <- paste0("min15_max", max_size)

    write.csv(
      filtered,
      file.path(out_dir, paste0(contrast, "_", tag, "_all_refiltered.csv")),
      row.names = FALSE
    )

    write.csv(
      fdr005,
      file.path(out_dir, paste0(contrast, "_", tag, "_FDR005.csv")),
      row.names = FALSE
    )

    write.csv(
      fdr010,
      file.path(out_dir, paste0(contrast, "_", tag, "_FDR010.csv")),
      row.names = FALSE
    )

    write.csv(
      collapsed_025,
      file.path(out_dir, paste0(contrast, "_", tag, "_FDR005_redundancy_collapsed_jaccard025.csv")),
      row.names = FALSE
    )

    write.csv(
      collapsed_050,
      file.path(out_dir, paste0(contrast, "_", tag, "_FDR005_redundancy_collapsed_jaccard050.csv")),
      row.names = FALSE
    )

    top_id <- if (nrow(filtered) > 0) filtered$ID[1] else NA
    top_nes <- if (nrow(filtered) > 0) filtered$NES[1] else NA
    top_fdr <- if (nrow(filtered) > 0) filtered$conservative_FDR[1] else NA
    top_size <- if (nrow(filtered) > 0) filtered$n_targets[1] else NA

    summary_rows[[paste(contrast, max_size, sep = "_")]] <- data.frame(
      contrast = contrast,
      filter = tag,
      n_terms_tested_after_size_filter = nrow(filtered),
      n_FDR005_after_refilter = nrow(fdr005),
      n_FDR010_after_refilter = nrow(fdr010),
      n_FDR005_after_jaccard025_collapse = nrow(collapsed_025),
      n_FDR005_after_jaccard050_collapse = nrow(collapsed_050),
      top_ID = top_id,
      top_NES = top_nes,
      top_conservative_FDR = top_fdr,
      top_n_targets = top_size,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)

write.csv(
  summary_df,
  file.path(out_dir, "conservative_miRNA_GSEA_filtering_summary.csv"),
  row.names = FALSE
)

if (length(missing_files) > 0) {
  writeLines(
    missing_files,
    con = file.path(out_dir, "missing_GSEA_all_files.txt")
  )

  message("Missing all-GSEA files for contrasts: ", paste(missing_files, collapse = ", "))
  message("Please run: ls results/14_strict_miRNA_target_GSEA")
}

print(summary_df)
cat("\nDone: conservative strict miRNA GSEA filtering completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
