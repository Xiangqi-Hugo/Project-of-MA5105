# Audit script for strict miRNA target GSEA.
# Purpose:
# 1. Check whether strict miRNA GSEA is overcalling because of large/redundant target sets.
# 2. Check target-set sizes for significant terms and primary candidates.
# 3. Check overlap among significant miRNA target sets.
# 4. Create simple diagnostic tables for the report/method discussion.

out_dir <- "results/28_miRNA_GSEA_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strict_dir <- "results/14_strict_miRNA_target_GSEA"
term2gene_file <- file.path(strict_dir, "strict_miRNA_TERM2GENE.csv")

if (!file.exists(term2gene_file)) {
  stop("Cannot find strict_miRNA_TERM2GENE.csv in results/14_strict_miRNA_target_GSEA")
}

term2gene <- read.csv(term2gene_file, stringsAsFactors = FALSE, check.names = FALSE)

# Robustly detect TERM and GENE columns.
cn <- colnames(term2gene)
term_col <- cn[tolower(cn) %in% c("term", "id", "miRNA", "mirna", "gs_name", "pathway")][1]
gene_col <- cn[tolower(cn) %in% c("gene", "target", "target_gene", "gene_symbol")][1]

if (is.na(term_col) || is.na(gene_col)) {
  term_col <- cn[1]
  gene_col <- cn[2]
}

term2gene <- term2gene[, c(term_col, gene_col)]
colnames(term2gene) <- c("term", "gene")

term2gene <- term2gene[!is.na(term2gene$term) & !is.na(term2gene$gene), ]
term2gene <- unique(term2gene)

target_size <- aggregate(
  gene ~ term,
  data = term2gene,
  FUN = function(x) length(unique(x))
)

colnames(target_size) <- c("term", "n_targets")

target_size <- target_size[order(target_size$n_targets, decreasing = TRUE), ]

write.csv(
  target_size,
  file.path(out_dir, "strict_miRNA_target_set_sizes_all.csv"),
  row.names = FALSE
)

size_summary <- data.frame(
  metric = c(
    "n_miRNA_terms",
    "min_targets",
    "median_targets",
    "mean_targets",
    "max_targets",
    "n_terms_gt_200_targets",
    "n_terms_gt_500_targets",
    "n_terms_gt_1000_targets"
  ),
  value = c(
    nrow(target_size),
    min(target_size$n_targets),
    median(target_size$n_targets),
    mean(target_size$n_targets),
    max(target_size$n_targets),
    sum(target_size$n_targets > 200),
    sum(target_size$n_targets > 500),
    sum(target_size$n_targets > 1000)
  )
)

write.csv(
  size_summary,
  file.path(out_dir, "strict_miRNA_target_set_size_summary.csv"),
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

primary_sizes <- target_size[target_size$term %in% primary_terms, ]
write.csv(
  primary_sizes,
  file.path(out_dir, "primary_candidate_target_set_sizes.csv"),
  row.names = FALSE
)

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

sig_size_list <- list()
overlap_summary_list <- list()

target_list <- split(term2gene$gene, term2gene$term)
target_list <- lapply(target_list, unique)

for (contrast in contrasts) {
  fdr005_file <- file.path(strict_dir, paste0(contrast, "_strict_miRNA_GSEA_FDR005.csv"))

  if (!file.exists(fdr005_file)) {
    warning("Missing file: ", fdr005_file)
    next
  }

  sig <- read.csv(fdr005_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (!"ID" %in% colnames(sig)) {
    id_col <- colnames(sig)[1]
    sig$ID <- sig[[id_col]]
  }

  sig <- sig[!is.na(sig$ID), ]

  sig_sizes <- merge(
    sig,
    target_size,
    by.x = "ID",
    by.y = "term",
    all.x = TRUE
  )

  sig_sizes$contrast <- contrast

  write.csv(
    sig_sizes,
    file.path(out_dir, paste0(contrast, "_significant_miRNA_target_set_sizes.csv")),
    row.names = FALSE
  )

  sig_size_list[[contrast]] <- sig_sizes

  sig_terms <- unique(sig$ID)
  sig_terms <- sig_terms[sig_terms %in% names(target_list)]

  if (length(sig_terms) >= 2) {
    pair_rows <- list()
    k <- 1

    for (i in seq_len(length(sig_terms) - 1)) {
      for (j in (i + 1):length(sig_terms)) {
        a <- target_list[[sig_terms[i]]]
        b <- target_list[[sig_terms[j]]]
        inter <- length(intersect(a, b))
        union <- length(union(a, b))
        jaccard <- ifelse(union > 0, inter / union, NA)

        pair_rows[[k]] <- data.frame(
          contrast = contrast,
          term1 = sig_terms[i],
          term2 = sig_terms[j],
          n_targets_term1 = length(a),
          n_targets_term2 = length(b),
          n_overlap = inter,
          n_union = union,
          jaccard = jaccard,
          stringsAsFactors = FALSE
        )
        k <- k + 1
      }
    }

    overlap_df <- do.call(rbind, pair_rows)
    overlap_df <- overlap_df[order(overlap_df$jaccard, decreasing = TRUE), ]

    write.csv(
      overlap_df,
      file.path(out_dir, paste0(contrast, "_significant_miRNA_target_set_pairwise_overlap.csv")),
      row.names = FALSE
    )

    overlap_summary_list[[contrast]] <- data.frame(
      contrast = contrast,
      n_significant_terms = length(sig_terms),
      median_pairwise_jaccard = median(overlap_df$jaccard, na.rm = TRUE),
      mean_pairwise_jaccard = mean(overlap_df$jaccard, na.rm = TRUE),
      max_pairwise_jaccard = max(overlap_df$jaccard, na.rm = TRUE),
      n_pairs_jaccard_gt_0_25 = sum(overlap_df$jaccard > 0.25, na.rm = TRUE),
      n_pairs_jaccard_gt_0_50 = sum(overlap_df$jaccard > 0.50, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  } else {
    overlap_summary_list[[contrast]] <- data.frame(
      contrast = contrast,
      n_significant_terms = length(sig_terms),
      median_pairwise_jaccard = NA,
      mean_pairwise_jaccard = NA,
      max_pairwise_jaccard = NA,
      n_pairs_jaccard_gt_0_25 = NA,
      n_pairs_jaccard_gt_0_50 = NA,
      stringsAsFactors = FALSE
    )
  }
}

if (length(sig_size_list) > 0) {
  sig_sizes_all <- do.call(rbind, sig_size_list)

  write.csv(
    sig_sizes_all,
    file.path(out_dir, "all_contrasts_significant_miRNA_target_set_sizes.csv"),
    row.names = FALSE
  )

  sig_size_summary <- do.call(
    rbind,
    lapply(split(sig_sizes_all, sig_sizes_all$contrast), function(x) {
      data.frame(
        contrast = unique(x$contrast),
        n_significant_terms = nrow(x),
        min_targets = min(x$n_targets, na.rm = TRUE),
        median_targets = median(x$n_targets, na.rm = TRUE),
        mean_targets = mean(x$n_targets, na.rm = TRUE),
        max_targets = max(x$n_targets, na.rm = TRUE),
        n_sig_terms_gt_200_targets = sum(x$n_targets > 200, na.rm = TRUE),
        n_sig_terms_gt_500_targets = sum(x$n_targets > 500, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  write.csv(
    sig_size_summary,
    file.path(out_dir, "significant_miRNA_target_set_size_summary_by_contrast.csv"),
    row.names = FALSE
  )
}

if (length(overlap_summary_list) > 0) {
  overlap_summary <- do.call(rbind, overlap_summary_list)

  write.csv(
    overlap_summary,
    file.path(out_dir, "significant_miRNA_target_set_overlap_summary.csv"),
    row.names = FALSE
  )
}

cat("\nDone: strict miRNA GSEA audit completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
