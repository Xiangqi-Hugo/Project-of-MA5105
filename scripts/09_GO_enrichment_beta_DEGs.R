suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

in_dir <- "results/07_beta_DE_edgeR"
out_dir <- "results/10_GO_enrichment"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_go <- function(contrast_name) {
  message("Running GO enrichment for: ", contrast_name)
  
  file <- file.path(in_dir, paste0(contrast_name, "_edgeR_all_genes.csv"))
  res <- read.csv(file, stringsAsFactors = FALSE)
  
  sig <- res[res$FDR < 0.05, ]
  
  up_genes <- sig$gene[sig$logFC > 0]
  down_genes <- sig$gene[sig$logFC < 0]
  
  gene_sets <- list(
    up = up_genes,
    down = down_genes
  )
  
  for (direction in names(gene_sets)) {
    genes <- unique(gene_sets[[direction]])
    
    genes <- genes[!grepl("^AC[0-9]|^AL[0-9]|^LINC", genes)]
    
    if (length(genes) < 10) {
      message("Skipping ", contrast_name, " ", direction, ": fewer than 10 genes.")
      next
    }
    
    mapped <- bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    mapped <- unique(mapped)
    
    if (nrow(mapped) < 10) {
      message("Skipping ", contrast_name, " ", direction, ": fewer than 10 mapped genes.")
      next
    }
    
    ego <- enrichGO(
      gene = mapped$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20,
      readable = TRUE
    )
    
    ego_df <- as.data.frame(ego)
    
    out_csv <- file.path(
      out_dir,
      paste0(contrast_name, "_", direction, "_GO_BP.csv")
    )
    
    write.csv(
      ego_df,
      out_csv,
      row.names = FALSE
    )
    
    if (nrow(ego_df) > 0) {
      png(
        filename = file.path(out_dir, paste0(contrast_name, "_", direction, "_GO_BP_dotplot.png")),
        width = 1400,
        height = 1000,
        res = 150
      )
      
      print(dotplot(ego, showCategory = 20) + ggplot2::ggtitle(paste0(contrast_name, " ", direction, " GO BP")))
      
      dev.off()
    }
    
    message("Saved: ", out_csv)
  }
}

run_go("T2D_vs_ND")
run_go("T2D_vs_PD")

message("Done: GO enrichment completed.")
