out_dir <- "results/08_DE_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

for (contrast in contrasts) {
  file <- paste0("results/07_beta_DE_edgeR/", contrast, "_edgeR_all_genes.csv")
  res <- read.csv(file, stringsAsFactors = FALSE)
  
  res <- res[order(res$FDR, -abs(res$logFC)), ]
  
  sig005 <- res[res$FDR < 0.05, ]
  sig010 <- res[res$FDR < 0.10, ]
  
  top_up005 <- sig005[sig005$logFC > 0, ]
  top_down005 <- sig005[sig005$logFC < 0, ]
  
  top_up005 <- top_up005[order(top_up005$FDR, -top_up005$logFC), ]
  top_down005 <- top_down005[order(top_down005$FDR, top_down005$logFC), ]
  
  write.csv(
    head(top_up005, 50),
    file.path(out_dir, paste0(contrast, "_top50_up_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    head(top_down005, 50),
    file.path(out_dir, paste0(contrast, "_top50_down_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig005,
    file.path(out_dir, paste0(contrast, "_all_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig010,
    file.path(out_dir, paste0(contrast, "_all_FDR010.csv")),
    row.names = FALSE
  )
  
  summary <- data.frame(
    contrast = contrast,
    n_tested = nrow(res),
    n_FDR005 = nrow(sig005),
    n_up_FDR005 = sum(sig005$logFC > 0),
    n_down_FDR005 = sum(sig005$logFC < 0),
    n_FDR010 = nrow(sig010),
    n_up_FDR010 = sum(sig010$logFC > 0),
    n_down_FDR010 = sum(sig010$logFC < 0),
    top_up_gene = ifelse(nrow(top_up005) > 0, top_up005$gene[1], NA),
    top_down_gene = ifelse(nrow(top_down005) > 0, top_down005$gene[1], NA)
  )
  
  write.csv(
    summary,
    file.path(out_dir, paste0(contrast, "_summary.csv")),
    row.names = FALSE
  )
}

all_summaries <- do.call(
  rbind,
  lapply(
    contrasts,
    function(x) read.csv(file.path(out_dir, paste0(x, "_summary.csv")))
  )
)

write.csv(
  all_summaries,
  file.path(out_dir, "all_contrast_DE_summary.csv"),
  row.names = FALSE
)

print(all_summaries)
