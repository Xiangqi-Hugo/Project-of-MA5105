# Run miRSCAPE with GTEx Pancreas paired mRNA-miRNA training data
# Project root: ~/T2D_beta_54sample_project

suppressPackageStartupMessages({
  library(Seurat)
  library(xgboost)
  library(limma)
  library(ggplot2)
})

source("tools/miRSCAPE/code/miRSCAPE.R")

in_dir <- "results/21_GTEx_miRSCAPE_training_input"
out_dir <- "results/22_miRSCAPE_prediction"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bulk_mrna_file <- file.path(in_dir, "GTEx_pancreas_bulk_mRNA_TPM_common_genes.rds")
bulk_mirna_file <- file.path(in_dir, "GTEx_pancreas_bulk_miRNA_TPM_miRBase_names.rds")
target_mrna_file <- file.path(in_dir, "GSE221156_beta_pseudobulk_logCPM_common_genes_for_miRSCAPE.rds")
sample_file <- "results/07_beta_DE_edgeR/samples_used_for_DE.csv"

message("Reading matrices...")
bulk_mrna <- as.matrix(readRDS(bulk_mrna_file))
bulk_mirna <- as.matrix(readRDS(bulk_mirna_file))
target_mrna <- as.matrix(readRDS(target_mrna_file))
storage.mode(bulk_mrna) <- "numeric"
storage.mode(bulk_mirna) <- "numeric"
storage.mode(target_mrna) <- "numeric"

message("Input dimensions:")
print(dim(bulk_mrna))
print(dim(bulk_mirna))
print(dim(target_mrna))

# Match genes.
common_genes <- intersect(rownames(bulk_mrna), rownames(target_mrna))
bulk_mrna <- bulk_mrna[common_genes, , drop = FALSE]
target_mrna <- target_mrna[common_genes, , drop = FALSE]
target_mrna <- target_mrna[rownames(bulk_mrna), , drop = FALSE]

# Match GTEx pancreas paired samples.
common_samples <- intersect(colnames(bulk_mrna), colnames(bulk_mirna))
bulk_mrna <- bulk_mrna[, common_samples, drop = FALSE]
bulk_mirna <- bulk_mirna[, common_samples, drop = FALSE]
bulk_mirna <- bulk_mirna[, colnames(bulk_mrna), drop = FALSE]

if (!identical(colnames(bulk_mrna), colnames(bulk_mirna))) {
  stop("GTEx mRNA and miRNA columns are not aligned.")
}

# Remove miRNAs with zero variance.
mirna_var <- apply(bulk_mirna, 1, var, na.rm = TRUE)
bulk_mirna <- bulk_mirna[is.finite(mirna_var) & mirna_var > 0, , drop = FALSE]

# Candidate availability.
primary_candidates <- c(
  "hsa-miR-195-5p",
  "hsa-miR-16-5p",
  "hsa-miR-15a-5p",
  "hsa-miR-15b-5p",
  "hsa-miR-649",
  "hsa-miR-6838-5p"
)

candidate_availability <- data.frame(
  miRNA = primary_candidates,
  available_in_GTEx_miRSCAPE = primary_candidates %in% rownames(bulk_mirna),
  stringsAsFactors = FALSE
)
write.csv(candidate_availability,
          file.path(out_dir, "primary_candidate_availability_in_GTEx_miRSCAPE.csv"),
          row.names = FALSE)
print(candidate_availability)

# Optional quick mode. Use: MIRSCAPE_MODE=candidate Rscript scripts/22_run_miRSCAPE_GTEx_pancreas_prediction.R
run_mode <- Sys.getenv("MIRSCAPE_MODE", unset = "full")
if (run_mode == "candidate") {
  keep <- intersect(primary_candidates, rownames(bulk_mirna))
  if (length(keep) == 0) stop("No primary candidates available for candidate mode.")
  bulk_mirna <- bulk_mirna[keep, , drop = FALSE]
  message("Candidate mode: testing only available primary candidates.")
}

message("Final dimensions used by miRSCAPE:")
print(dim(bulk_mrna))
print(dim(bulk_mirna))
print(dim(target_mrna))

message("Transforming bulk mRNA and bulk miRNA...")
bulkmRNA <- bulkTransform(bulk_mrna)
bulkmiRNA <- bulkTransform(bulk_mirna, justNorm = TRUE)
scmRNA <- target_mrna

message("Running miRSCAPE. This may take time...")
set.seed(123)
pred <- miRSCAPE(bulkmRNA = bulkmRNA, bulkmiRNA = bulkmiRNA, scmRNA = scmRNA)

message("miRSCAPE output class:")
print(class(pred))

if (is.data.frame(pred) || is.matrix(pred)) {
  pred_mat <- as.matrix(pred)
} else if (is.list(pred)) {
  if ("pred" %in% names(pred)) {
    pred_mat <- as.matrix(pred$pred)
  } else if ("prediction" %in% names(pred)) {
    pred_mat <- as.matrix(pred$prediction)
  } else {
    message("List names:")
    print(names(pred))
    stop("Unknown miRSCAPE list output format.")
  }
} else {
  stop("Unknown miRSCAPE output format.")
}

storage.mode(pred_mat) <- "numeric"

# Orient as miRNA x beta-sample.
if (all(colnames(target_mrna) %in% colnames(pred_mat))) {
  pred_mirna_by_sample <- pred_mat
} else if (all(colnames(target_mrna) %in% rownames(pred_mat))) {
  pred_mirna_by_sample <- t(pred_mat)
} else {
  message("Prediction dimensions:")
  print(dim(pred_mat))
  message("Prediction rownames example:")
  print(head(rownames(pred_mat)))
  message("Prediction colnames example:")
  print(head(colnames(pred_mat)))
  stop("Could not orient prediction matrix.")
}

pred_mirna_by_sample <- pred_mirna_by_sample[, colnames(target_mrna), drop = FALSE]

write.csv(pred_mirna_by_sample,
          file.path(out_dir, "miRSCAPE_GTEx_pancreas_predicted_miRNA_expression.csv"),
          row.names = TRUE)
saveRDS(pred_mirna_by_sample,
        file.path(out_dir, "miRSCAPE_GTEx_pancreas_predicted_miRNA_expression.rds"))

message("Prediction matrix dimension:")
print(dim(pred_mirna_by_sample))

# Differential analysis of predicted miRNA expression.
sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)
sample_info <- sample_info[match(colnames(pred_mirna_by_sample), sample_info$sample_prefix), ]
if (any(is.na(sample_info$sample_prefix))) stop("Sample metadata does not match prediction columns.")
sample_info$disease_group <- factor(sample_info$disease_group, levels = c("ND", "PD", "T2D"))

design <- model.matrix(~ 0 + disease_group, data = sample_info)
colnames(design) <- levels(sample_info$disease_group)
fit <- lmFit(pred_mirna_by_sample, design)
contrast_mat <- makeContrasts(
  T2D_vs_ND = T2D - ND,
  PD_vs_ND = PD - ND,
  T2D_vs_PD = T2D - PD,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contrast_mat))

summary_list <- list()
candidate_list <- list()

for (contrast in colnames(contrast_mat)) {
  tab <- topTable(fit2, coef = contrast, number = Inf, sort.by = "P")
  tab$miRNA <- rownames(tab)
  tab <- tab[, c("miRNA", setdiff(colnames(tab), "miRNA"))]
  first_group <- sub("_vs_.*$", "", contrast)
  tab$predicted_expression_direction <- ifelse(tab$logFC > 0,
                                                paste0("higher_in_", first_group),
                                                paste0("lower_in_", first_group))
  write.csv(tab, file.path(out_dir, paste0(contrast, "_miRSCAPE_predicted_miRNA_all.csv")), row.names = FALSE)
  write.csv(tab[tab$adj.P.Val < 0.05, ], file.path(out_dir, paste0(contrast, "_miRSCAPE_predicted_miRNA_FDR005.csv")), row.names = FALSE)
  write.csv(tab[tab$adj.P.Val < 0.10, ], file.path(out_dir, paste0(contrast, "_miRSCAPE_predicted_miRNA_FDR010.csv")), row.names = FALSE)
  summary_list[[contrast]] <- data.frame(
    contrast = contrast,
    n_miRNAs_tested = nrow(tab),
    n_FDR005 = sum(tab$adj.P.Val < 0.05),
    n_FDR010 = sum(tab$adj.P.Val < 0.10),
    top_miRNA = tab$miRNA[1],
    top_logFC = tab$logFC[1],
    top_P_value = tab$P.Value[1],
    top_FDR = tab$adj.P.Val[1],
    top_direction = tab$predicted_expression_direction[1],
    stringsAsFactors = FALSE
  )
  subtab <- tab[tab$miRNA %in% primary_candidates, ]
  subtab$contrast <- contrast
  candidate_list[[contrast]] <- subtab
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_dir, "miRSCAPE_predicted_miRNA_DE_summary.csv"), row.names = FALSE)
print(summary_df)

candidate_df <- do.call(rbind, candidate_list)
write.csv(candidate_df, file.path(out_dir, "primary_candidate_miRSCAPE_predicted_expression_summary.csv"), row.names = FALSE)
print(candidate_df)

# Candidate boxplot.
available_candidates <- intersect(primary_candidates, rownames(pred_mirna_by_sample))
if (length(available_candidates) > 0) {
  long_df <- data.frame(
    miRNA = rep(available_candidates, times = ncol(pred_mirna_by_sample)),
    sample_prefix = rep(colnames(pred_mirna_by_sample), each = length(available_candidates)),
    predicted_expression = as.numeric(pred_mirna_by_sample[available_candidates, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  long_df <- merge(long_df, sample_info[, c("sample_prefix", "disease_group")], by = "sample_prefix", all.x = TRUE)
  long_df$disease_group <- factor(long_df$disease_group, levels = c("ND", "PD", "T2D"))
  p <- ggplot(long_df, aes(x = disease_group, y = predicted_expression)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
    facet_wrap(~ miRNA, scales = "free_y") +
    theme_bw(base_size = 12) +
    labs(title = "miRSCAPE predicted miRNA expression", x = "Disease group", y = "Predicted miRNA expression")
  ggsave(file.path(out_dir, "primary_candidate_miRSCAPE_predicted_expression_boxplots.png"), p, width = 8, height = 5, dpi = 300)
}

cat("\nDone: miRSCAPE GTEx pancreas prediction completed.\n")
