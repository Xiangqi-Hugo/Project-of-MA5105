# Step 42: CELLxGENE beta-only miRSCAPE Model1 prediction using GTEx pancreas training
#
# Purpose:
#   Repeat miRSCAPE Model1 on CELLxGENE beta-only matched41 pseudo-bulk.
#
# Inputs:
#   results/21_GTEx_miRSCAPE_training_input/
#   results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_counts_matrix_matched41.csv
#   results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_sample_summary_matched41.csv
#
# Outputs:
#   results/42_cellxgene_miRSCAPE_GTEx_pancreas_prediction/

suppressPackageStartupMessages({
  library(edgeR)
  library(xgboost)
  library(limma)
  library(ggplot2)
})

source("tools/miRSCAPE/code/miRSCAPE.R")

miRSCAPE_xgbtrain <- function(
  bulkmRNA,
  bulkmiRNA,
  scmRNA,
  nrounds = 30,
  max_depth = 4,
  eta = 0.3,
  subsample = 0.9,
  colsample_bytree = 0.8,
  mode_label = "full"
) {
  commGenes <- intersect(rownames(bulkmRNA), rownames(scmRNA))
  if (length(commGenes) < 1000) stop("Too few common genes between bulk mRNA and target mRNA.")

  bulkmRNA <- bulkmRNA[commGenes, , drop = FALSE]
  scmRNA <- scmRNA[commGenes, , drop = FALSE]

  x_train <- t(bulkmRNA)
  x_test <- t(scmRNA)

  if (nrow(x_train) != ncol(bulkmiRNA)) {
    stop("Training sample mismatch: nrow(t(bulkmRNA)) = ", nrow(x_train), ", ncol(bulkmiRNA) = ", ncol(bulkmiRNA))
  }

  result <- matrix(NA_real_, nrow = nrow(bulkmiRNA), ncol = ncol(scmRNA),
                   dimnames = list(rownames(bulkmiRNA), colnames(scmRNA)))

  dtest <- xgb.DMatrix(data = x_test)

  params <- list(
    booster = "gbtree",
    objective = "reg:squarederror",
    max_depth = max_depth,
    eta = eta,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    eval_metric = "rmse"
  )

  for (i in seq_len(nrow(bulkmiRNA))) {
    y <- as.numeric(bulkmiRNA[i, ])
    dtrain <- xgb.DMatrix(data = x_train, label = y)

    bst <- xgb.train(params = params, data = dtrain, nrounds = nrounds, verbose = 0)
    result[i, ] <- predict(bst, dtest)

    if (i %% 10 == 0 || i == nrow(bulkmiRNA)) {
      message("Predicted ", i, " / ", nrow(bulkmiRNA), " miRNAs [", mode_label, "]")
    }
  }

  result
}

in_dir <- "results/21_GTEx_miRSCAPE_training_input"
out_dir <- "results/42_cellxgene_miRSCAPE_GTEx_pancreas_prediction"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bulk_mrna_file <- file.path(in_dir, "GTEx_pancreas_bulk_mRNA_TPM_common_genes.rds")
bulk_mirna_file <- file.path(in_dir, "GTEx_pancreas_bulk_miRNA_TPM_miRBase_names.rds")
counts_file <- "results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_counts_matrix_matched41.csv"
sample_file <- "results/38_cellxgene_beta_pseudobulk/cellxgene_beta_pseudobulk_sample_summary_matched41.csv"

if (!file.exists(bulk_mrna_file)) stop("Missing: ", bulk_mrna_file)
if (!file.exists(bulk_mirna_file)) stop("Missing: ", bulk_mirna_file)
if (!file.exists(counts_file)) stop("Missing: ", counts_file)
if (!file.exists(sample_file)) stop("Missing: ", sample_file)

message("Reading GTEx matrices...")
bulk_mrna <- as.matrix(readRDS(bulk_mrna_file))
bulk_mirna <- as.matrix(readRDS(bulk_mirna_file))

message("Reading CELLxGENE beta counts...")
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
counts <- round(counts)

sample_info <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)
sample_info <- sample_info[match(colnames(counts), sample_info$LibraryID), ]
if (!identical(sample_info$LibraryID, colnames(counts))) stop("Sample metadata does not match count matrix columns.")

sample_info$disease_group <- factor(sample_info$disease_group_for_DE, levels = c("ND", "PD", "T2D"))

dge <- DGEList(counts = counts, group = sample_info$disease_group)
keep <- filterByExpr(dge, group = sample_info$disease_group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
target_mrna <- cpm(dge, log = TRUE, prior.count = 1)

write.csv(target_mrna, file.path(out_dir, "cellxgene_beta_matched41_logCPM_TMM_for_miRSCAPE_GTEx.csv"))

storage.mode(bulk_mrna) <- "numeric"
storage.mode(bulk_mirna) <- "numeric"
storage.mode(target_mrna) <- "numeric"

common_genes <- intersect(rownames(bulk_mrna), rownames(target_mrna))
bulk_mrna <- bulk_mrna[common_genes, , drop = FALSE]
target_mrna <- target_mrna[common_genes, , drop = FALSE]
target_mrna <- target_mrna[rownames(bulk_mrna), , drop = FALSE]

common_samples <- intersect(colnames(bulk_mrna), colnames(bulk_mirna))
bulk_mrna <- bulk_mrna[, common_samples, drop = FALSE]
bulk_mirna <- bulk_mirna[, common_samples, drop = FALSE]
bulk_mirna <- bulk_mirna[, colnames(bulk_mrna), drop = FALSE]

mirna_var <- apply(bulk_mirna, 1, var, na.rm = TRUE)
bulk_mirna <- bulk_mirna[is.finite(mirna_var) & mirna_var > 0, , drop = FALSE]

primary_candidates <- c("hsa-miR-195-5p", "hsa-miR-16-5p", "hsa-miR-15a-5p", "hsa-miR-15b-5p", "hsa-miR-649", "hsa-miR-6838-5p")

candidate_availability <- data.frame(
  miRNA = primary_candidates,
  available_in_GTEx_miRSCAPE = primary_candidates %in% rownames(bulk_mirna),
  stringsAsFactors = FALSE
)
write.csv(candidate_availability, file.path(out_dir, "cellxgene_primary_candidate_availability_in_GTEx_miRSCAPE.csv"), row.names = FALSE)
print(candidate_availability)

run_mode <- Sys.getenv("MIRSCAPE_MODE", unset = "full")

if (run_mode == "candidate") {
  keep_candidates <- c("hsa-miR-15a-5p", "hsa-miR-15b-5p")
  keep_candidates <- keep_candidates[keep_candidates %in% rownames(bulk_mirna)]
  if (length(keep_candidates) == 0) stop("Candidate mode requested, but no available candidates were found.")
  bulk_mirna <- bulk_mirna[keep_candidates, , drop = FALSE]
}

message("Final dimensions used by CELLxGENE GTEx miRSCAPE:")
print(dim(bulk_mrna))
print(dim(bulk_mirna))
print(dim(target_mrna))

message("Transforming bulk mRNA and bulk miRNA...")
bulkmRNA <- bulkTransform(bulk_mrna)
bulkmiRNA <- bulkTransform(bulk_mirna, justNorm = TRUE)
scmRNA <- target_mrna

set.seed(123)
pred_mirna_by_sample <- miRSCAPE_xgbtrain(
  bulkmRNA = bulkmRNA,
  bulkmiRNA = bulkmiRNA,
  scmRNA = scmRNA,
  nrounds = 30,
  mode_label = paste0("CELLxGENE_GTEx_", run_mode)
)

write.csv(pred_mirna_by_sample, file.path(out_dir, "cellxgene_miRSCAPE_GTEx_pancreas_predicted_miRNA_expression.csv"), row.names = TRUE)
saveRDS(pred_mirna_by_sample, file.path(out_dir, "cellxgene_miRSCAPE_GTEx_pancreas_predicted_miRNA_expression.rds"))

sample_info <- sample_info[match(colnames(pred_mirna_by_sample), sample_info$LibraryID), ]
if (any(is.na(sample_info$LibraryID))) stop("Sample metadata does not match predicted miRNA matrix.")

design <- model.matrix(~ 0 + disease_group, data = sample_info)
colnames(design) <- levels(sample_info$disease_group)

fit <- lmFit(pred_mirna_by_sample, design)

contrasts <- makeContrasts(
  T2D_vs_ND = T2D - ND,
  PD_vs_ND = PD - ND,
  T2D_vs_PD = T2D - PD,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

contrast_names <- colnames(contrasts)
summary_list <- list()
all_list <- list()

for (contrast in contrast_names) {
  tab <- topTable(fit2, coef = contrast, number = Inf, sort.by = "P")
  tab$miRNA <- rownames(tab)
  tab <- tab[, c("miRNA", setdiff(colnames(tab), "miRNA"))]
  first_group <- sub("_vs_.*$", "", contrast)
  tab$predicted_expression_direction <- ifelse(tab$logFC > 0, paste0("higher_in_", first_group), paste0("lower_in_", first_group))
  tab$contrast <- contrast

  write.csv(tab, file.path(out_dir, paste0(contrast, "_cellxgene_miRSCAPE_GTEx_predicted_miRNA_all.csv")), row.names = FALSE)
  write.csv(tab[tab$adj.P.Val < 0.05, ], file.path(out_dir, paste0(contrast, "_cellxgene_miRSCAPE_GTEx_predicted_miRNA_FDR005.csv")), row.names = FALSE)
  write.csv(tab[tab$adj.P.Val < 0.10, ], file.path(out_dir, paste0(contrast, "_cellxgene_miRSCAPE_GTEx_predicted_miRNA_FDR010.csv")), row.names = FALSE)

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
  all_list[[contrast]] <- tab
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(out_dir, "cellxgene_miRSCAPE_GTEx_predicted_miRNA_DE_summary.csv"), row.names = FALSE)

all_tests <- do.call(rbind, all_list)
all_tests$global_FDR_across_all_GTEx_miRSCAPE_tests <- p.adjust(all_tests$P.Value, method = "BH")
all_tests <- all_tests[order(all_tests$global_FDR_across_all_GTEx_miRSCAPE_tests, all_tests$P.Value), ]
write.csv(all_tests, file.path(out_dir, "cellxgene_GTEx_miRSCAPE_all_contrasts_global_FDR.csv"), row.names = FALSE)

global_summary <- do.call(rbind, lapply(contrast_names, function(contrast) {
  x <- all_tests[all_tests$contrast == contrast, ]
  x <- x[order(x$global_FDR_across_all_GTEx_miRSCAPE_tests, x$P.Value), ]
  data.frame(
    contrast = contrast,
    n_tests = nrow(x),
    n_global_FDR005 = sum(x$global_FDR_across_all_GTEx_miRSCAPE_tests < 0.05),
    n_global_FDR010 = sum(x$global_FDR_across_all_GTEx_miRSCAPE_tests < 0.10),
    top_miRNA = x$miRNA[1],
    top_logFC = x$logFC[1],
    top_P_value = x$P.Value[1],
    top_within_contrast_FDR = x$adj.P.Val[1],
    top_global_FDR = x$global_FDR_across_all_GTEx_miRSCAPE_tests[1],
    top_direction = x$predicted_expression_direction[1],
    stringsAsFactors = FALSE
  )
}))

write.csv(global_summary, file.path(out_dir, "cellxgene_GTEx_miRSCAPE_global_FDR_summary.csv"), row.names = FALSE)

candidate_rows <- unique(c(primary_candidates, "hsa-miR-15a-5p", "hsa-miR-15b-5p"))
candidate_list <- list()
for (contrast in contrast_names) {
  tab <- read.csv(file.path(out_dir, paste0(contrast, "_cellxgene_miRSCAPE_GTEx_predicted_miRNA_all.csv")), stringsAsFactors = FALSE)
  subtab <- tab[tab$miRNA %in% candidate_rows, ]
  subtab$contrast <- contrast
  candidate_list[[contrast]] <- subtab
}
candidate_df <- do.call(rbind, candidate_list)
write.csv(candidate_df, file.path(out_dir, "cellxgene_GTEx_primary_candidate_miRSCAPE_predicted_expression_summary.csv"), row.names = FALSE)

print(summary_df)
print(global_summary)
print(candidate_df)

cat("\nDone: CELLxGENE xgb.train miRSCAPE GTEx pancreas prediction completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
