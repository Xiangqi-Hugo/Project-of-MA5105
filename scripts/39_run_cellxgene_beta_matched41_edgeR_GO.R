# Step 39: Run edgeR and optional GO analysis on CELLxGENE beta-cell matched41 pseudo-bulk
#
# Inputs from Step 38:
#   results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_counts_matrix_matched41.csv
#   results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_sample_summary_matched41.csv
#
# Outputs:
#   results/39_cellxgene_beta_matched41_edgeR_GO/
#
# Main purpose:
#   Repeat the same sample-level beta-cell pseudo-bulk DE workflow using CELLxGENE beta-cell annotation.
#   This is an annotation sensitivity analysis against the original marker-based beta-cell workflow.

out_dir <- "results/39_cellxgene_beta_matched41_edgeR_GO"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- "results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_counts_matrix_matched41.csv"
sample_file <- "results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_sample_summary_matched41.csv"

if (!file.exists(counts_file)) stop("Missing counts file: ", counts_file)
if (!file.exists(sample_file)) stop("Missing sample file: ", sample_file)

suppressPackageStartupMessages({
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR is not installed. Install edgeR first.")
  }
  library(edgeR)
})

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
counts <- round(counts)

samples <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!"LibraryID" %in% colnames(samples)) stop("sample summary must contain LibraryID.")
if (!"disease_group_for_DE" %in% colnames(samples)) stop("sample summary must contain disease_group_for_DE.")

samples <- samples[match(colnames(counts), samples$LibraryID), ]
if (!identical(samples$LibraryID, colnames(counts))) {
  stop("Sample order could not be aligned between counts matrix and sample table.")
}

samples$group <- factor(samples$disease_group_for_DE, levels = c("ND", "PD", "T2D"))

# Remove genes with missing or zero library sizes.
keep_samples <- !is.na(samples$group)
counts <- counts[, keep_samples, drop = FALSE]
samples <- samples[keep_samples, ]

dge <- DGEList(counts = counts, samples = samples, group = samples$group)
keep_genes <- filterByExpr(dge, group = samples$group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0 + group, data = samples)
colnames(design) <- levels(samples$group)

fit <- glmQLFit(dge, design, robust = TRUE)

contrast_list <- list(
  T2D_vs_ND = makeContrasts(T2D - ND, levels = design),
  PD_vs_ND = makeContrasts(PD - ND, levels = design),
  T2D_vs_PD = makeContrasts(T2D - PD, levels = design)
)

summary_rows <- list()

plot_volcano <- function(tab, contrast_name, out_png) {
  fdr <- tab$FDR
  logfc <- tab$logFC
  neglog <- -log10(pmax(tab$PValue, .Machine$double.xmin))
  sig <- fdr < 0.05

  png(out_png, width = 1800, height = 1600, res = 220)
  plot(
    logfc, neglog,
    pch = 16,
    cex = 0.45,
    xlab = "logFC",
    ylab = "-log10(P-value)",
    main = paste0("CELLxGENE beta matched41: ", contrast_name)
  )
  points(logfc[sig], neglog[sig], pch = 16, cex = 0.45)
  abline(h = -log10(0.05), lty = 2)
  abline(v = c(-1, 1), lty = 3)
  dev.off()
}

for (cn in names(contrast_list)) {
  qlf <- glmQLFTest(fit, contrast = contrast_list[[cn]])
  tab <- topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  tab <- tab[, c("gene", setdiff(colnames(tab), "gene"))]
  tab <- tab[order(tab$FDR, tab$PValue), ]

  write.csv(tab, file.path(out_dir, paste0(cn, "_edgeR_all_genes.csv")), row.names = FALSE)
  write.csv(tab[tab$FDR < 0.05, ], file.path(out_dir, paste0(cn, "_edgeR_DEGs_FDR005.csv")), row.names = FALSE)
  write.csv(tab[tab$FDR < 0.10, ], file.path(out_dir, paste0(cn, "_edgeR_DEGs_FDR010.csv")), row.names = FALSE)

  n_sig <- sum(tab$FDR < 0.05, na.rm = TRUE)
  n_up <- sum(tab$FDR < 0.05 & tab$logFC > 0, na.rm = TRUE)
  n_down <- sum(tab$FDR < 0.05 & tab$logFC < 0, na.rm = TRUE)
  n_sig10 <- sum(tab$FDR < 0.10, na.rm = TRUE)

  top_gene <- if (nrow(tab) > 0) tab$gene[1] else NA
  min_fdr <- if (nrow(tab) > 0) min(tab$FDR, na.rm = TRUE) else NA

  summary_rows[[cn]] <- data.frame(
    contrast = cn,
    tested_genes = nrow(tab),
    FDR005 = n_sig,
    up = n_up,
    down = n_down,
    FDR010 = n_sig10,
    top_gene = top_gene,
    min_FDR = min_fdr,
    stringsAsFactors = FALSE
  )

  plot_volcano(tab, cn, file.path(out_dir, paste0(cn, "_volcano.png")))
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "CELLxGENE_beta_matched41_edgeR_summary.csv"), row.names = FALSE)

# Compare against original marker-based summary, if available.
orig_summary_file <- "results/07_beta_DE_edgeR/DE_contrast_summary.csv"
if (file.exists(orig_summary_file)) {
  orig <- read.csv(orig_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  write.csv(orig, file.path(out_dir, "original_marker_based_edgeR_summary_copied.csv"), row.names = FALSE)
}

# Optional GO analysis.
# This section runs only if clusterProfiler and org.Hs.eg.db are installed.
run_go <- requireNamespace("clusterProfiler", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)

if (run_go) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })

  go_dir <- file.path(out_dir, "GO")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

  for (cn in names(contrast_list)) {
    tab <- read.csv(file.path(out_dir, paste0(cn, "_edgeR_all_genes.csv")), stringsAsFactors = FALSE)

    # ORA for FDR<0.05 up/down genes.
    for (direction in c("up", "down")) {
      if (direction == "up") {
        genes <- tab$gene[tab$FDR < 0.05 & tab$logFC > 0]
      } else {
        genes <- tab$gene[tab$FDR < 0.05 & tab$logFC < 0]
      }

      if (length(genes) >= 10) {
        eg <- tryCatch(
          bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
          error = function(e) NULL
        )
        if (!is.null(eg) && nrow(eg) >= 10) {
          ego <- tryCatch(
            enrichGO(
              gene = unique(eg$ENTREZID),
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              pAdjustMethod = "BH",
              readable = TRUE
            ),
            error = function(e) NULL
          )
          if (!is.null(ego)) {
            out <- as.data.frame(ego)
            write.csv(out, file.path(go_dir, paste0(cn, "_", direction, "_GO_BP_ORA.csv")), row.names = FALSE)
          }
        }
      }
    }

    # Ranked GO GSEA.
    tab <- tab[!is.na(tab$PValue) & !is.na(tab$logFC), ]
    tab$rank_metric <- sign(tab$logFC) * (-log10(pmax(tab$PValue, .Machine$double.xmin)))
    tab <- tab[order(-tab$rank_metric), ]
    tab <- tab[!duplicated(tab$gene), ]

    map <- tryCatch(
      bitr(tab$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) NULL
    )

    if (!is.null(map) && nrow(map) > 100) {
      rank_df <- merge(tab[, c("gene", "rank_metric")], map, by.x = "gene", by.y = "SYMBOL")
      rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]
      gene_list <- rank_df$rank_metric
      names(gene_list) <- rank_df$ENTREZID
      gene_list <- sort(gene_list, decreasing = TRUE)

      gsea <- tryCatch(
        gseGO(
          geneList = gene_list,
          OrgDb = org.Hs.eg.db,
          ont = "BP",
          keyType = "ENTREZID",
          pAdjustMethod = "BH",
          minGSSize = 10,
          maxGSSize = 500,
          eps = 0,
          verbose = FALSE
        ),
        error = function(e) NULL
      )

      if (!is.null(gsea)) {
        gout <- as.data.frame(gsea)
        write.csv(gout, file.path(go_dir, paste0(cn, "_GO_BP_GSEA_all.csv")), row.names = FALSE)
        write.csv(gout[gout$p.adjust < 0.05, ], file.path(go_dir, paste0(cn, "_GO_BP_GSEA_FDR005.csv")), row.names = FALSE)
      }
    }
  }
} else {
  cat("\n[WARN] clusterProfiler/org.Hs.eg.db not available. GO analysis skipped.\n")
  cat("[WARN] edgeR results were still completed.\n")
}

cat("\nCELLxGENE beta matched41 edgeR summary:\n")
print(summary_df)

cat("\nDone.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
