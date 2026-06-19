suppressPackageStartupMessages({
  library(edgeR)
})

count_csv <- "results/06_beta_pseudobulk/beta_pseudobulk_counts_matrix.csv"
sample_file <- "results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv"

out_dir <- "results/07_beta_DE_edgeR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading pseudo-bulk count matrix from CSV...")

count_df <- read.csv(
  count_csv,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(count_df)
storage.mode(counts) <- "numeric"

message("Count matrix dimension:")
print(dim(counts))

sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

sample_info <- sample_info[match(colnames(counts), sample_info$sample_prefix), ]

if (any(is.na(sample_info$sample_prefix))) {
  stop("Sample metadata does not match pseudo-bulk matrix columns.")
}

min_beta_cells <- 100

keep_samples <- sample_info$use_for_conservative_DE &
  sample_info$disease_group %in% c("ND", "PD", "T2D") &
  sample_info$n_beta_cells >= min_beta_cells

counts_de <- counts[, keep_samples, drop = FALSE]
sample_info_de <- sample_info[keep_samples, ]

sample_info_de$disease_group <- factor(
  sample_info_de$disease_group,
  levels = c("ND", "PD", "T2D")
)

write.csv(
  sample_info_de,
  file.path(out_dir, "samples_used_for_DE.csv"),
  row.names = FALSE
)

group_counts <- as.data.frame(table(sample_info_de$disease_group))
colnames(group_counts) <- c("disease_group", "n_samples")

write.csv(
  group_counts,
  file.path(out_dir, "DE_group_sample_counts.csv"),
  row.names = FALSE
)

message("Samples used for DE:")
print(group_counts)

message("Samples excluded from DE:")
excluded <- sample_info[!keep_samples, c(
  "sample_prefix",
  "disease_group",
  "donor_or_pair_id",
  "n_beta_cells",
  "beta_fraction",
  "use_for_conservative_DE"
)]

write.csv(
  excluded,
  file.path(out_dir, "samples_excluded_from_DE.csv"),
  row.names = FALSE
)

print(excluded)

if (any(group_counts$n_samples < 2)) {
  stop("At least one group has fewer than 2 samples after filtering.")
}

dge <- DGEList(
  counts = counts_de,
  samples = sample_info_de,
  group = sample_info_de$disease_group
)

design <- model.matrix(~ 0 + disease_group, data = sample_info_de)
colnames(design) <- levels(sample_info_de$disease_group)

keep_genes <- filterByExpr(dge, design = design)

gene_filter_summary <- data.frame(
  n_genes_before_filtering = nrow(dge),
  n_genes_after_filtering = sum(keep_genes),
  stringsAsFactors = FALSE
)

write.csv(
  gene_filter_summary,
  file.path(out_dir, "gene_filtering_summary.csv"),
  row.names = FALSE
)

dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design, robust = TRUE)

contrasts <- makeContrasts(
  T2D_vs_ND = T2D - ND,
  PD_vs_ND = PD - ND,
  T2D_vs_PD = T2D - PD,
  levels = design
)

run_contrast <- function(contrast_name) {
  message("Running contrast: ", contrast_name)
  
  qlf <- glmQLFTest(fit, contrast = contrasts[, contrast_name])
  res <- topTags(qlf, n = Inf)$table
  
  res$gene <- rownames(res)
  res <- res[, c("gene", setdiff(colnames(res), "gene"))]
  
  res$direction_FDR005 <- ifelse(
    res$FDR < 0.05 & res$logFC > 0,
    "up",
    ifelse(
      res$FDR < 0.05 & res$logFC < 0,
      "down",
      "not_significant"
    )
  )
  
  res$direction_FDR01 <- ifelse(
    res$FDR < 0.10 & res$logFC > 0,
    "up",
    ifelse(
      res$FDR < 0.10 & res$logFC < 0,
      "down",
      "not_significant"
    )
  )
  
  write.csv(
    res,
    file.path(out_dir, paste0(contrast_name, "_edgeR_all_genes.csv")),
    row.names = FALSE
  )
  
  sig005 <- res[res$FDR < 0.05, ]
  sig01 <- res[res$FDR < 0.10, ]
  
  write.csv(
    sig005,
    file.path(out_dir, paste0(contrast_name, "_edgeR_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig01,
    file.path(out_dir, paste0(contrast_name, "_edgeR_FDR010.csv")),
    row.names = FALSE
  )
  
  png(
    filename = file.path(out_dir, paste0(contrast_name, "_volcano.png")),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  plot(
    res$logFC,
    -log10(res$PValue),
    pch = 16,
    cex = 0.45,
    xlab = "log2 fold change",
    ylab = "-log10 p-value",
    main = contrast_name
  )
  
  abline(v = 0, lty = 2)
  abline(h = -log10(0.05), lty = 2)
  
  dev.off()
  
  data.frame(
    contrast = contrast_name,
    n_tested_genes = nrow(res),
    n_FDR005 = nrow(sig005),
    n_up_FDR005 = sum(sig005$logFC > 0),
    n_down_FDR005 = sum(sig005$logFC < 0),
    n_FDR010 = nrow(sig01),
    n_up_FDR010 = sum(sig01$logFC > 0),
    n_down_FDR010 = sum(sig01$logFC < 0),
    top_gene_by_FDR = res$gene[which.min(res$FDR)],
    min_FDR = min(res$FDR, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

de_summary <- do.call(
  rbind,
  lapply(colnames(contrasts), run_contrast)
)

write.csv(
  de_summary,
  file.path(out_dir, "DE_contrast_summary.csv"),
  row.names = FALSE
)

png(
  filename = file.path(out_dir, "MDS_plot_beta_pseudobulk.png"),
  width = 1200,
  height = 1000,
  res = 150
)

plotMDS(
  dge,
  labels = sample_info_de$disease_group,
  col = as.numeric(sample_info_de$disease_group),
  main = "Beta-cell pseudo-bulk MDS"
)

legend(
  "topright",
  legend = levels(sample_info_de$disease_group),
  col = seq_along(levels(sample_info_de$disease_group)),
  pch = 16
)

dev.off()

logcpm <- cpm(dge, log = TRUE, prior.count = 1)

write.csv(
  logcpm,
  file.path(out_dir, "beta_pseudobulk_logCPM_TMM.csv"),
  row.names = TRUE
)

saveRDS(
  dge,
  file.path(out_dir, "edgeR_DGEList_beta_pseudobulk.rds")
)

print(gene_filter_summary)
print(de_summary)

message("Done: edgeR pseudo-bulk DE completed.")
