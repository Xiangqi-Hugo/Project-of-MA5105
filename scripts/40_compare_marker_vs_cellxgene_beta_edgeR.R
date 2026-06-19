# Step 40: Compare marker-based beta analysis with CELLxGENE beta-only analysis
#
# Purpose:
#   Quantify how similar the original marker-based pseudo-bulk edgeR result is to
#   the CELLxGENE beta-only pseudo-bulk edgeR result.
#
# Inputs:
#   Original marker-based edgeR results:
#     results/07_beta_DE_edgeR/
#   CELLxGENE beta-only edgeR results:
#     results/39_cellxgene_beta_matched41_edgeR_GO/
#
# Outputs:
#   results/40_marker_vs_cellxgene_beta_comparison/
#     marker_vs_cellxgene_edgeR_summary.csv
#     marker_vs_cellxgene_DEG_overlap_summary.csv
#     marker_vs_cellxgene_rank_correlation_summary.csv
#     <contrast>_merged_gene_level_comparison.csv

out_dir <- "results/40_marker_vs_cellxgene_beta_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

orig_dir <- "results/07_beta_DE_edgeR"
cxg_dir <- "results/39_cellxgene_beta_matched41_edgeR_GO"

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

find_result_file <- function(dir_path, contrast) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  files <- files[grepl(contrast, basename(files), fixed = TRUE)]

  # Prefer all-gene files.
  preferred <- files[grepl("all", basename(files), ignore.case = TRUE) |
                       grepl("edgeR", basename(files), ignore.case = TRUE)]
  if (length(preferred) > 0) {
    # Avoid DEG-only files if possible.
    preferred2 <- preferred[!grepl("FDR|DEG", basename(preferred), ignore.case = TRUE)]
    if (length(preferred2) > 0) return(preferred2[1])
    return(preferred[1])
  }
  if (length(files) > 0) return(files[1])
  NA_character_
}

standardize <- function(tab) {
  # Standardize gene column.
  if (!"gene" %in% colnames(tab)) {
    possible_gene_cols <- c("Gene", "genes", "symbol", "SYMBOL", "rowname")
    hit <- intersect(possible_gene_cols, colnames(tab))
    if (length(hit) > 0) {
      colnames(tab)[match(hit[1], colnames(tab))] <- "gene"
    } else {
      tab$gene <- rownames(tab)
    }
  }

  # Standardize FDR column.
  if (!"FDR" %in% colnames(tab)) {
    hit <- intersect(c("adj.P.Val", "padj", "PAdj", "qvalue", "p.adjust"), colnames(tab))
    if (length(hit) > 0) colnames(tab)[match(hit[1], colnames(tab))] <- "FDR"
  }

  # Standardize PValue column.
  if (!"PValue" %in% colnames(tab)) {
    hit <- intersect(c("P.Value", "Pvalue", "pvalue", "p.value", "P"), colnames(tab))
    if (length(hit) > 0) colnames(tab)[match(hit[1], colnames(tab))] <- "PValue"
  }

  required <- c("gene", "logFC", "FDR")
  missing <- setdiff(required, colnames(tab))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  tab <- tab[!is.na(tab$gene) & tab$gene != "", ]
  tab <- tab[!duplicated(tab$gene), ]
  tab
}

overlap_rows <- list()
cor_rows <- list()

for (contrast in contrasts) {
  orig_file <- find_result_file(orig_dir, contrast)
  cxg_file <- file.path(cxg_dir, paste0(contrast, "_edgeR_all_genes.csv"))

  if (is.na(orig_file) || !file.exists(orig_file)) {
    warning("Could not find original file for ", contrast)
    next
  }
  if (!file.exists(cxg_file)) {
    warning("Could not find CELLxGENE file for ", contrast)
    next
  }

  orig <- standardize(read.csv(orig_file, stringsAsFactors = FALSE, check.names = FALSE))
  cxg <- standardize(read.csv(cxg_file, stringsAsFactors = FALSE, check.names = FALSE))

  merged <- merge(
    orig[, c("gene", "logFC", "FDR", if ("PValue" %in% colnames(orig)) "PValue" else NULL)],
    cxg[, c("gene", "logFC", "FDR", if ("PValue" %in% colnames(cxg)) "PValue" else NULL)],
    by = "gene",
    suffixes = c("_marker", "_cellxgene")
  )

  if ("PValue_marker" %in% colnames(merged)) {
    merged$rank_marker <- sign(merged$logFC_marker) * (-log10(pmax(merged$PValue_marker, .Machine$double.xmin)))
  } else {
    merged$rank_marker <- sign(merged$logFC_marker) * (-log10(pmax(merged$FDR_marker, .Machine$double.xmin)))
  }

  if ("PValue_cellxgene" %in% colnames(merged)) {
    merged$rank_cellxgene <- sign(merged$logFC_cellxgene) * (-log10(pmax(merged$PValue_cellxgene, .Machine$double.xmin)))
  } else {
    merged$rank_cellxgene <- sign(merged$logFC_cellxgene) * (-log10(pmax(merged$FDR_cellxgene, .Machine$double.xmin)))
  }

  merged$marker_sig_FDR005 <- merged$FDR_marker < 0.05
  merged$cellxgene_sig_FDR005 <- merged$FDR_cellxgene < 0.05
  merged$marker_sig_FDR010 <- merged$FDR_marker < 0.10
  merged$cellxgene_sig_FDR010 <- merged$FDR_cellxgene < 0.10
  merged$direction_match <- sign(merged$logFC_marker) == sign(merged$logFC_cellxgene)

  merged <- merged[order(merged$FDR_marker + merged$FDR_cellxgene), ]
  write.csv(
    merged,
    file.path(out_dir, paste0(contrast, "_merged_gene_level_comparison.csv")),
    row.names = FALSE
  )

  for (thr in c("FDR005", "FDR010")) {
    marker_col <- paste0("marker_sig_", thr)
    cxg_col <- paste0("cellxgene_sig_", thr)
    marker_genes <- merged$gene[merged[[marker_col]]]
    cxg_genes <- merged$gene[merged[[cxg_col]]]
    ov <- intersect(marker_genes, cxg_genes)
    union <- union(marker_genes, cxg_genes)

    overlap_rows[[paste(contrast, thr, sep = "_")]] <- data.frame(
      contrast = contrast,
      threshold = thr,
      n_marker_DEGs = length(marker_genes),
      n_cellxgene_DEGs = length(cxg_genes),
      n_overlap = length(ov),
      jaccard = ifelse(length(union) > 0, length(ov) / length(union), NA_real_),
      fraction_marker_recovered_in_cellxgene = ifelse(length(marker_genes) > 0, length(ov) / length(marker_genes), NA_real_),
      fraction_cellxgene_recovered_in_marker = ifelse(length(cxg_genes) > 0, length(ov) / length(cxg_genes), NA_real_),
      same_direction_overlap_fraction = ifelse(
        length(ov) > 0,
        mean(merged$direction_match[merged$gene %in% ov], na.rm = TRUE),
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
  }

  # Top 100 ranked overlap.
  top_marker <- merged$gene[order(-abs(merged$rank_marker))][1:min(100, nrow(merged))]
  top_cxg <- merged$gene[order(-abs(merged$rank_cellxgene))][1:min(100, nrow(merged))]

  cor_rows[[contrast]] <- data.frame(
    contrast = contrast,
    n_common_genes = nrow(merged),
    logFC_pearson = suppressWarnings(cor(merged$logFC_marker, merged$logFC_cellxgene, method = "pearson", use = "complete.obs")),
    logFC_spearman = suppressWarnings(cor(merged$logFC_marker, merged$logFC_cellxgene, method = "spearman", use = "complete.obs")),
    rank_spearman = suppressWarnings(cor(merged$rank_marker, merged$rank_cellxgene, method = "spearman", use = "complete.obs")),
    global_direction_match_fraction = mean(merged$direction_match, na.rm = TRUE),
    top100_overlap = length(intersect(top_marker, top_cxg)),
    top100_jaccard = length(intersect(top_marker, top_cxg)) / length(union(top_marker, top_cxg)),
    stringsAsFactors = FALSE
  )
}

overlap_summary <- do.call(rbind, overlap_rows)
cor_summary <- do.call(rbind, cor_rows)

write.csv(overlap_summary, file.path(out_dir, "marker_vs_cellxgene_DEG_overlap_summary.csv"), row.names = FALSE)
write.csv(cor_summary, file.path(out_dir, "marker_vs_cellxgene_rank_correlation_summary.csv"), row.names = FALSE)

# Combine edgeR summaries.
orig_sum_file <- file.path(orig_dir, "DE_contrast_summary.csv")
cxg_sum_file <- file.path(cxg_dir, "CELLxGENE_beta_matched41_edgeR_summary.csv")

if (file.exists(orig_sum_file) && file.exists(cxg_sum_file)) {
  orig_sum <- read.csv(orig_sum_file, stringsAsFactors = FALSE, check.names = FALSE)
  cxg_sum <- read.csv(cxg_sum_file, stringsAsFactors = FALSE, check.names = FALSE)

  # Try to standardize original summary column names.
  names(orig_sum) <- gsub(" ", "_", names(orig_sum))
  names(cxg_sum) <- gsub(" ", "_", names(cxg_sum))

  orig_sum$source <- "marker_based_beta"
  cxg_sum$source <- "CELLxGENE_beta_only"

  common_cols <- intersect(colnames(orig_sum), colnames(cxg_sum))
  combined <- rbind(orig_sum[, common_cols, drop = FALSE], cxg_sum[, common_cols, drop = FALSE])
  write.csv(combined, file.path(out_dir, "marker_vs_cellxgene_edgeR_summary.csv"), row.names = FALSE)
}

cat("\nDEG overlap summary:\n")
print(overlap_summary)

cat("\nRank/logFC correlation summary:\n")
print(cor_summary)

cat("\nDone.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
