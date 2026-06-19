suppressPackageStartupMessages({
  library(msigdbr)
  library(clusterProfiler)
})

in_dir <- "results/07_beta_DE_edgeR"
out_dir <- "results/14_strict_miRNA_target_GSEA"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msig_all <- msigdbr(species = "Homo sapiens")

# Keep only gene sets whose names look like direct miRNA IDs.
# Examples: MIR195_5P, MIR16_5P, MIR15A_5P, MIR144_3P, LET7A_5P
strict_idx <- grepl("^MIR[0-9A-Z]+(_[35]P)?$", msig_all$gs_name) |
  grepl("^MIRLET[0-9A-Z]+(_[35]P)?$", msig_all$gs_name) |
  grepl("^LET7[0-9A-Z]+(_[35]P)?$", msig_all$gs_name)

mir_sets <- msig_all[strict_idx, ]

message("Strict miRNA gene-set rows:")
print(nrow(mir_sets))

message("Strict miRNA gene-set number:")
print(length(unique(mir_sets$gs_name)))

message("Example strict miRNA sets:")
print(head(unique(mir_sets$gs_name), 30))

if (nrow(mir_sets) == 0) {
  stop("No strict miRNA gene sets found. Need to inspect MSigDB names.")
}

TERM2GENE <- unique(
  data.frame(
    term = mir_sets$gs_name,
    gene = mir_sets$gene_symbol,
    stringsAsFactors = FALSE
  )
)

write.csv(
  unique(mir_sets[, c("gs_name", "gene_symbol")]),
  file.path(out_dir, "strict_miRNA_TERM2GENE.csv"),
  row.names = FALSE
)

run_gsea <- function(contrast_name) {
  message("Running strict miRNA target GSEA for: ", contrast_name)
  
  file <- file.path(in_dir, paste0(contrast_name, "_edgeR_all_genes.csv"))
  res <- read.csv(file, stringsAsFactors = FALSE)
  
  res <- res[!is.na(res$gene) & !is.na(res$logFC) & !is.na(res$PValue), ]
  res <- res[res$gene != "", ]
  
  res <- res[order(res$gene, res$PValue), ]
  res <- res[!duplicated(res$gene), ]
  
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
  
  first_group <- sub("_vs_.*$", "", contrast_name)
  
  gsea_df$target_direction <- ifelse(
    gsea_df$NES < 0,
    paste0("targets_down_in_", first_group),
    paste0("targets_up_in_", first_group)
  )
  
  gsea_df$inferred_miRNA_activity <- ifelse(
    gsea_df$NES < 0,
    paste0("higher_in_", first_group),
    paste0("lower_in_", first_group)
  )
  
  write.csv(
    gsea_df,
    file.path(out_dir, paste0(contrast_name, "_strict_miRNA_GSEA_all.csv")),
    row.names = FALSE
  )
  
  sig005 <- gsea_df[gsea_df$p.adjust < 0.05, ]
  sig010 <- gsea_df[gsea_df$p.adjust < 0.10, ]
  
  write.csv(
    sig005,
    file.path(out_dir, paste0(contrast_name, "_strict_miRNA_GSEA_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig010,
    file.path(out_dir, paste0(contrast_name, "_strict_miRNA_GSEA_FDR010.csv")),
    row.names = FALSE
  )
  
  negative <- sig005[sig005$NES < 0, ]
  positive <- sig005[sig005$NES > 0, ]
  
  negative <- negative[order(negative$p.adjust, negative$NES), ]
  positive <- positive[order(positive$p.adjust, -positive$NES), ]
  
  data.frame(
    contrast = contrast_name,
    n_ranked_genes = length(gene_list),
    n_strict_miRNA_terms_all = nrow(gsea_df),
    n_strict_miRNA_FDR005 = nrow(sig005),
    n_targets_down_in_first_group = nrow(negative),
    n_targets_up_in_first_group = nrow(positive),
    top_targets_down_term = ifelse(nrow(negative) > 0, negative$Description[1], NA),
    top_targets_down_NES = ifelse(nrow(negative) > 0, negative$NES[1], NA),
    top_targets_down_FDR = ifelse(nrow(negative) > 0, negative$p.adjust[1], NA),
    top_targets_up_term = ifelse(nrow(positive) > 0, positive$Description[1], NA),
    top_targets_up_NES = ifelse(nrow(positive) > 0, positive$NES[1], NA),
    top_targets_up_FDR = ifelse(nrow(positive) > 0, positive$p.adjust[1], NA),
    stringsAsFactors = FALSE
  )
}

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

summary_list <- lapply(contrasts, run_gsea)
summary_df <- do.call(rbind, summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "strict_miRNA_GSEA_summary.csv"),
  row.names = FALSE
)

print(summary_df)

message("Done: strict miRNA target GSEA completed.")
