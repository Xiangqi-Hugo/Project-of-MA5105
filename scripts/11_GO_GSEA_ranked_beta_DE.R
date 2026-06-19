suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

in_dir <- "results/07_beta_DE_edgeR"
out_dir <- "results/11_GO_GSEA"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_gsea <- function(contrast_name) {
  message("Running ranked GO GSEA for: ", contrast_name)
  
  file <- file.path(in_dir, paste0(contrast_name, "_edgeR_all_genes.csv"))
  res <- read.csv(file, stringsAsFactors = FALSE)
  
  res <- res[!is.na(res$gene) & !is.na(res$logFC) & !is.na(res$PValue), ]
  res <- res[res$gene != "", ]
  
  # Remove duplicate gene symbols.
  # Keep the row with the smallest p-value.
  res <- res[order(res$gene, res$PValue), ]
  res <- res[!duplicated(res$gene), ]
  
  # Ranked metric.
  # Positive value means higher in the first group of the contrast.
  # Example: T2D_vs_ND positive means higher in T2D.
  res$rank_metric <- sign(res$logFC) * (-log10(res$PValue + 1e-300))
  
  res <- res[is.finite(res$rank_metric), ]
  res <- res[order(res$rank_metric, decreasing = TRUE), ]
  
  gene_list <- res$rank_metric
  names(gene_list) <- res$gene
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  gsea_df <- as.data.frame(gsea)
  
  write.csv(
    gsea_df,
    file.path(out_dir, paste0(contrast_name, "_GO_BP_GSEA_all.csv")),
    row.names = FALSE
  )
  
  sig005 <- gsea_df[gsea_df$p.adjust < 0.05, ]
  sig010 <- gsea_df[gsea_df$p.adjust < 0.10, ]
  
  write.csv(
    sig005,
    file.path(out_dir, paste0(contrast_name, "_GO_BP_GSEA_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig010,
    file.path(out_dir, paste0(contrast_name, "_GO_BP_GSEA_FDR010.csv")),
    row.names = FALSE
  )
  
  summary_one <- data.frame(
    contrast = contrast_name,
    n_ranked_genes = length(gene_list),
    n_GSEA_terms_all = nrow(gsea_df),
    n_GSEA_FDR005 = nrow(sig005),
    n_GSEA_FDR010 = nrow(sig010),
    top_positive_term = ifelse(
      any(gsea_df$NES > 0),
      gsea_df$Description[which.min(ifelse(gsea_df$NES > 0, gsea_df$p.adjust, Inf))],
      NA
    ),
    top_negative_term = ifelse(
      any(gsea_df$NES < 0),
      gsea_df$Description[which.min(ifelse(gsea_df$NES < 0, gsea_df$p.adjust, Inf))],
      NA
    ),
    stringsAsFactors = FALSE
  )
  
  # Optional dotplot
  if (requireNamespace("enrichplot", quietly = TRUE) &&
      requireNamespace("ggplot2", quietly = TRUE) &&
      nrow(gsea_df) > 0) {
    
    png(
      filename = file.path(out_dir, paste0(contrast_name, "_GO_BP_GSEA_dotplot.png")),
      width = 1400,
      height = 1000,
      res = 150
    )
    
    print(
      enrichplot::dotplot(gsea, showCategory = 20) +
        ggplot2::ggtitle(paste0(contrast_name, " GO BP GSEA"))
    )
    
    dev.off()
  }
  
  summary_one
}

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

summary_list <- lapply(contrasts, run_gsea)
summary_df <- do.call(rbind, summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "GO_BP_GSEA_summary.csv"),
  row.names = FALSE
)

print(summary_df)

message("Done: ranked GO GSEA completed.")
