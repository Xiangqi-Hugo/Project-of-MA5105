# Step 44: Summarize CELLxGENE beta-only miRNA robustness results
#
# Purpose:
#   Combine CELLxGENE beta-only decoupleR, GTEx-miRSCAPE, and TCGA-PAAD miRSCAPE outputs.
#   Create report-ready tables and a conservative interpretation sentence.
#
# Inputs:
#   results/41_cellxgene_decoupleR_conservative_miRNA_activity/
#   results/42_cellxgene_miRSCAPE_GTEx_pancreas_prediction/
#   results/43_cellxgene_miRSCAPE_Model2_TCGA_PAAD_prediction/
#
# Outputs:
#   results/44_cellxgene_miRNA_robustness_summary/
#     cellxgene_miRNA_method_level_summary.csv
#     cellxgene_primary_candidate_summary.csv
#     cellxgene_miRNA_report_interpretation.txt

out_dir <- "results/44_cellxgene_miRNA_robustness_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dec_dir <- "results/41_cellxgene_decoupleR_conservative_miRNA_activity"
gtex_dir <- "results/42_cellxgene_miRSCAPE_GTEx_pancreas_prediction"
m2_dir <- "results/43_cellxgene_miRSCAPE_Model2_TCGA_PAAD_prediction"

required_files <- c(
  file.path(dec_dir, "cellxgene_conservative_decoupleR_global_FDR_summary.csv"),
  file.path(dec_dir, "cellxgene_conservative_decoupleR_primary_candidates_all_filters.csv"),
  file.path(gtex_dir, "cellxgene_GTEx_miRSCAPE_global_FDR_summary.csv"),
  file.path(gtex_dir, "cellxgene_GTEx_primary_candidate_miRSCAPE_predicted_expression_summary.csv"),
  file.path(m2_dir, "cellxgene_Model2_global_FDR_summary.csv"),
  file.path(m2_dir, "cellxgene_Model2_primary_candidate_precursor_proxy_summary.csv")
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"))
}

dec_sum <- read.csv(file.path(dec_dir, "cellxgene_conservative_decoupleR_global_FDR_summary.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
dec_pri <- read.csv(file.path(dec_dir, "cellxgene_conservative_decoupleR_primary_candidates_all_filters.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

gtex_sum <- read.csv(file.path(gtex_dir, "cellxgene_GTEx_miRSCAPE_global_FDR_summary.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
gtex_pri <- read.csv(file.path(gtex_dir, "cellxgene_GTEx_primary_candidate_miRSCAPE_predicted_expression_summary.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

m2_sum <- read.csv(file.path(m2_dir, "cellxgene_Model2_global_FDR_summary.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
m2_pri <- read.csv(file.path(m2_dir, "cellxgene_Model2_primary_candidate_precursor_proxy_summary.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

# Method-level summary.
dec_method <- data.frame(
  method = paste0("decoupleR_", dec_sum$filter),
  contrast = dec_sum$contrast,
  n_tests = dec_sum$n_tests,
  n_global_FDR005 = dec_sum$n_global_FDR005,
  n_global_FDR010 = dec_sum$n_global_FDR010,
  top_signal = dec_sum$top_miRNA,
  top_logFC = dec_sum$top_logFC,
  top_P_value = dec_sum$top_P_value,
  top_global_FDR = dec_sum$top_global_FDR,
  top_direction = dec_sum$top_direction,
  interpretation_level = ifelse(dec_sum$n_global_FDR005 > 0, "significant_method_specific", "negative"),
  stringsAsFactors = FALSE
)

gtex_method <- data.frame(
  method = "miRSCAPE_Model1_GTEx_pancreas",
  contrast = gtex_sum$contrast,
  n_tests = gtex_sum$n_tests,
  n_global_FDR005 = gtex_sum$n_global_FDR005,
  n_global_FDR010 = gtex_sum$n_global_FDR010,
  top_signal = gtex_sum$top_miRNA,
  top_logFC = gtex_sum$top_logFC,
  top_P_value = gtex_sum$top_P_value,
  top_global_FDR = gtex_sum$top_global_FDR,
  top_direction = gtex_sum$top_direction,
  interpretation_level = ifelse(gtex_sum$n_global_FDR005 > 0, "significant_method_specific", "negative"),
  stringsAsFactors = FALSE
)

m2_method <- data.frame(
  method = "miRSCAPE_Model2_TCGA_PAAD",
  contrast = m2_sum$contrast,
  n_tests = m2_sum$n_tests,
  n_global_FDR005 = m2_sum$n_global_FDR005,
  n_global_FDR010 = m2_sum$n_global_FDR010,
  top_signal = m2_sum$top_miRNA_feature,
  top_logFC = m2_sum$top_logFC,
  top_P_value = m2_sum$top_P_value,
  top_global_FDR = m2_sum$top_global_FDR,
  top_direction = m2_sum$top_direction,
  interpretation_level = ifelse(m2_sum$n_global_FDR005 > 0, "significant_method_specific", "negative"),
  stringsAsFactors = FALSE
)

method_level <- rbind(dec_method, gtex_method, m2_method)
write.csv(method_level, file.path(out_dir, "cellxgene_miRNA_method_level_summary.csv"), row.names = FALSE)

# Primary candidate summary.
primary_mature <- c("hsa-miR-195-5p", "hsa-miR-16-5p", "hsa-miR-15a-5p", "hsa-miR-15b-5p", "hsa-miR-649", "hsa-miR-6838-5p")
primary_strict <- c("MIR195_5P", "MIR16_5P", "MIR15A_5P", "MIR15B_5P", "MIR649", "MIR6838_5P")
primary_precursor <- c("hsa-mir-195", "hsa-mir-16-1", "hsa-mir-16-2", "hsa-mir-15a", "hsa-mir-15b", "hsa-mir-649", "hsa-mir-6838")

primary_rows <- list()

# decoupleR primary.
if (nrow(dec_pri) > 0) {
  d <- dec_pri
  d$method <- paste0("decoupleR_", d$filter)
  d$feature <- d$miRNA
  d$FDR_used <- d$global_FDR_across_all_contrasts
  d$direction <- d$activity_direction
  d$is_primary_candidate <- d$feature %in% primary_strict
  d$is_significant <- d$FDR_used < 0.05
  primary_rows[["decoupleR"]] <- d[, c("method", "contrast", "feature", "logFC", "P.Value", "adj.P.Val", "FDR_used", "direction", "is_primary_candidate", "is_significant")]
}

# GTEx primary.
if (nrow(gtex_pri) > 0) {
  d <- gtex_pri
  d$method <- "miRSCAPE_Model1_GTEx_pancreas"
  d$feature <- d$miRNA
  d$FDR_used <- d$adj.P.Val
  d$direction <- d$predicted_expression_direction
  d$is_primary_candidate <- d$feature %in% primary_mature
  d$is_significant <- d$FDR_used < 0.05
  primary_rows[["GTEx"]] <- d[, c("method", "contrast", "feature", "logFC", "P.Value", "adj.P.Val", "FDR_used", "direction", "is_primary_candidate", "is_significant")]
}

# Model2 primary proxies.
if (nrow(m2_pri) > 0) {
  d <- m2_pri
  d$method <- "miRSCAPE_Model2_TCGA_PAAD"
  d$feature <- d$miRNA_feature
  d$FDR_used <- d$adj.P.Val
  d$direction <- d$predicted_expression_direction
  d$is_primary_candidate <- d$feature %in% primary_precursor
  d$is_significant <- d$FDR_used < 0.05
  primary_rows[["Model2"]] <- d[, c("method", "contrast", "feature", "logFC", "P.Value", "adj.P.Val", "FDR_used", "direction", "is_primary_candidate", "is_significant")]
}

primary_summary <- do.call(rbind, primary_rows)
rownames(primary_summary) <- NULL
write.csv(primary_summary, file.path(out_dir, "cellxgene_primary_candidate_summary.csv"), row.names = FALSE)

# Conservative interpretation.
n_dec_sig <- sum(dec_method$n_global_FDR005 > 0)
n_gtex_sig_rows <- sum(gtex_method$n_global_FDR005 > 0)
n_m2_sig_rows <- sum(m2_method$n_global_FDR005 > 0)
n_primary_sig <- sum(primary_summary$is_significant, na.rm = TRUE)

gtex_top_T2D_ND <- gtex_method$top_signal[gtex_method$contrast == "T2D_vs_ND"][1]
m2_top_T2D_ND <- m2_method$top_signal[m2_method$contrast == "T2D_vs_ND"][1]

interpretation <- c(
  "CELLxGENE beta-only miRNA robustness summary",
  "",
  paste0("Conservative decoupleR detected globally significant miRNA activity in ", n_dec_sig, " method-contrast rows."),
  paste0("GTEx-miRSCAPE detected globally significant predicted miRNA expression in ", n_gtex_sig_rows, " contrast rows."),
  paste0("TCGA-PAAD miRSCAPE Model2 detected globally significant predicted miRNA expression in ", n_m2_sig_rows, " contrast rows."),
  paste0("The top T2D vs ND signal in GTEx-miRSCAPE was ", gtex_top_T2D_ND, "."),
  paste0("The top T2D vs ND signal in TCGA-PAAD Model2 was ", m2_top_T2D_ND, "."),
  paste0("The number of significant primary-candidate/proxy rows was ", n_primary_sig, "."),
  "",
  "Conservative interpretation:",
  "The CELLxGENE beta-only analysis reproduced the negative decoupleR result and did not validate the original primary miRNA candidates.",
  "miRSCAPE Model1 and Model2 detected method-specific predicted miRNA signals, but their top signals differed.",
  "Therefore, the CELLxGENE beta-only analysis supports the same final conclusion: the miRNA layer remains exploratory and no robust miRNA regulator was identified."
)

writeLines(interpretation, file.path(out_dir, "cellxgene_miRNA_report_interpretation.txt"))

cat("\nMethod-level summary:\n")
print(method_level)

cat("\nPrimary candidate summary:\n")
print(primary_summary)

cat("\nInterpretation:\n")
cat(paste(interpretation, collapse = "\n"), "\n")

cat("\nDone.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
