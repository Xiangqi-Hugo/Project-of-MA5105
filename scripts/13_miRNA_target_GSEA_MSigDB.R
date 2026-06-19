suppressPackageStartupMessages({
  library(msigdbr)
  library(clusterProfiler)
})

in_dir <- "results/07_beta_DE_edgeR"
out_dir <- "results/13_miRNA_target_GSEA"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading MSigDB gene sets...")

msig_all <- msigdbr(species = "Homo sapiens")

message("Available MSigDB columns:")
print(colnames(msig_all))

# Robustly detect miRNA target gene sets across msigdbr versions.
possible_cols <- intersect(
  c("gs_collection", "gs_subcollection", "gs_cat", "gs_subcat", "gs_name"),
  colnames(msig_all)
)

mir_idx <- rep(FALSE, nrow(msig_all))

for (cc in possible_cols) {
  mir_idx <- mir_idx | grepl("MIR|miR|microRNA|MIRNA", msig_all[[cc]])
}

mir_sets <- msig_all[mir_idx, ]

if (nrow(mir_sets) == 0) {
  message("No miRNA-related gene sets detected automatically.")
  message("MSigDB collection table:")
  print(unique(msig_all[, possible_cols]))
  stop("No miRNA gene sets found. Please inspect msigdbr collections.")
}

message("Number of miRNA target gene-set rows:")
print(nrow(mir_sets))

message("Example miRNA gene-set names:")
print(head(unique(mir_sets$gs_name), 20))

TERM2GENE <- unique(
  data.frame(
    term = mir_sets$gs_name,
    gene = mir_sets$gene_symbol,
    stringsAsFactors = FALSE
  )
)

write.csv(
  unique(mir_sets[, intersect(c("gs_name", "gs_collection", "gs_subcollection", "gs_cat", "gs_subcat"), colnames(mir_sets))]),
  file.path(out_dir, "miRNA_gene_sets_detected.csv"),
  row.names = FALSE
)

run_mirna_gsea <- function(contrast_name) {
  message("Running miRNA target GSEA for: ", contrast_name)
  
  file <- file.path(in_dir, paste0(contrast_name, "_edgeR_all_genes.csv"))
  res <- read.csv(file, stringsAsFactors = FALSE)
  
  res <- res[!is.na(res$gene) & !is.na(res$logFC) & !is.na(res$PValue), ]
  res <- res[res$gene != "", ]
  
  # Remove duplicate gene symbols by smallest p-value.
  res <- res[order(res$gene, res$PValue), ]
  res <- res[!duplicated(res$gene), ]
  
  # Ranked metric.
  # Positive = higher in first group of contrast.
  # Example: T2D_vs_ND positive means higher in T2D.
  res$rank_metric <- sign(res$logFC) * (-log10(res$PValue + 1e-300))
  res <- res[is.finite(res$rank_metric), ]
  
  gene_list <- res$rank_metric
  names(gene_list) <- res$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea <- GSEA(
    geneList = gene_list,
    TERM2GENE = TERM2GENE,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  gsea_df <- as.data.frame(gsea)
  
  # Infer miRNA activity direction from target direction.
  # NES < 0 means targets lower in the first group.
  # That is compatible with higher miRNA activity in the first group.
  gsea_df$inferred_miRNA_activity <- ifelse(
    gsea_df$NES < 0,
    paste0("higher_in_", sub("_vs_.*$", "", contrast_name)),
    paste0("lower_in_", sub("_vs_.*$", "", contrast_name))
  )
  
  write.csv(
    gsea_df,
    file.path(out_dir, paste0(contrast_name, "_miRNA_target_GSEA_all.csv")),
    row.names = FALSE
  )
  
  sig005 <- gsea_df[gsea_df$p.adjust < 0.05, ]
  sig010 <- gsea_df[gsea_df$p.adjust < 0.10, ]
  
  write.csv(
    sig005,
    file.path(out_dir, paste0(contrast_name, "_miRNA_target_GSEA_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig010,
    file.path(out_dir, paste0(contrast_name, "_miRNA_target_GSEA_FDR010.csv")),
    row.names = FALSE
  )
  
  positive <- sig005[sig005$NES > 0, ]
  negative <- sig005[sig005$NES < 0, ]
  
  positive <- positive[order(positive$p.adjust, -positive$NES), ]
  negative <- negative[order(negative$p.adjust, negative$NES), ]
  
  write.csv(
    head(positive, 50),
    file.path(out_dir, paste0(contrast_name, "_top50_targets_up_in_first_group.csv")),
    row.names = FALSE
  )
  
  write.csv(
    head(negative, 50),
    file.path(out_dir, paste0(contrast_name, "_top50_targets_down_in_first_group.csv")),
    row.names = FALSE
  )
  
  summary_one <- data.frame(
    contrast = contrast_name,
    n_ranked_genes = length(gene_list),
    n_miRNA_terms_all = nrow(gsea_df),
    n_miRNA_FDR005 = nrow(sig005),
    n_targets_up_in_first_group = nrow(positive),
    n_targets_down_in_first_group = nrow(negative),
    top_targets_up_term = ifelse(nrow(positive) > 0, positive$Description[1], NA),
    top_targets_up_NES = ifelse(nrow(positive) > 0, positive$NES[1], NA),
    top_targets_up_FDR = ifelse(nrow(positive) > 0, positive$p.adjust[1], NA),
    top_targets_down_term = ifelse(nrow(negative) > 0, negative$Description[1], NA),
    top_targets_down_NES = ifelse(nrow(negative) > 0, negative$NES[1], NA),
    top_targets_down_FDR = ifelse(nrow(negative) > 0, negative$p.adjust[1], NA),
    stringsAsFactors = FALSE
  )
  
  summary_one
}

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

summary_list <- lapply(contrasts, run_mirna_gsea)
summary_df <- do.call(rbind, summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "miRNA_target_GSEA_summary.csv"),
  row.names = FALSE
)

print(summary_df)

message("Done: miRNA target GSEA completed.")
