out_dir <- "results/27_miRSCAPE_Model1_Model2_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model1_file <- "results/23_miRSCAPE_global_FDR/miRSCAPE_all_contrasts_global_FDR.csv"
model2_file <- "results/26_miRSCAPE_Model2_prediction/Model2_all_contrasts_global_FDR.csv"

model1 <- read.csv(model1_file, stringsAsFactors = FALSE)
model2 <- read.csv(model2_file, stringsAsFactors = FALSE)

to_precursor_key <- function(x) {
  y <- tolower(x)
  y <- gsub("^hsa-mir-", "hsa-mir-", y)
  y <- gsub("^hsa-mir", "hsa-mir", y)
  y <- gsub("-5p$", "", y)
  y <- gsub("-3p$", "", y)
  y
}

model1$key <- to_precursor_key(model1$miRNA)
model2$key <- to_precursor_key(model2$miRNA_feature)

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

model1_best_list <- list()
model2_best_list <- list()

for (contrast in contrasts) {
  x1 <- model1[model1$contrast == contrast, ]
  x1 <- x1[order(x1$key, x1$global_FDR_across_all_contrasts, x1$P.Value), ]
  x1 <- x1[!duplicated(x1$key), ]
  model1_best_list[[contrast]] <- x1
  
  x2 <- model2[model2$contrast == contrast, ]
  x2 <- x2[order(x2$key, x2$global_FDR_across_all_Model2_tests, x2$P.Value), ]
  x2 <- x2[!duplicated(x2$key), ]
  model2_best_list[[contrast]] <- x2
}

model1_best <- do.call(rbind, model1_best_list)
model2_best <- do.call(rbind, model2_best_list)

summary_list <- list()

for (contrast in contrasts) {
  m1 <- model1_best[model1_best$contrast == contrast, ]
  m2 <- model2_best[model2_best$contrast == contrast, ]
  
  m1_sig <- m1[m1$global_FDR_across_all_contrasts < 0.05, ]
  m2_sig <- m2[m2$global_FDR_across_all_Model2_tests < 0.05, ]
  
  common <- merge(
    m1[, c(
      "key", "miRNA", "logFC", "P.Value", "adj.P.Val",
      "global_FDR_across_all_contrasts",
      "predicted_expression_direction"
    )],
    m2[, c(
      "key", "miRNA_feature", "logFC", "P.Value", "adj.P.Val",
      "global_FDR_across_all_Model2_tests",
      "predicted_expression_direction"
    )],
    by = "key",
    all = FALSE,
    suffixes = c("_Model1_GTEx", "_Model2_TCGA")
  )
  
  names(common) <- gsub("P.Value_Model1_GTEx", "P_Model1_GTEx", names(common))
  names(common) <- gsub("P.Value_Model2_TCGA", "P_Model2_TCGA", names(common))
  names(common) <- gsub("adj.P.Val_Model1_GTEx", "within_FDR_Model1_GTEx", names(common))
  names(common) <- gsub("adj.P.Val_Model2_TCGA", "within_FDR_Model2_TCGA", names(common))
  
  common$direction_match <- sign(common$logFC_Model1_GTEx) == sign(common$logFC_Model2_TCGA)
  common$model1_global_sig <- common$global_FDR_across_all_contrasts < 0.05
  common$model2_global_sig <- common$global_FDR_across_all_Model2_tests < 0.05
  common$both_global_sig <- common$model1_global_sig & common$model2_global_sig
  
  write.csv(
    common,
    file.path(out_dir, paste0(contrast, "_Model1_Model2_common_key_comparison.csv")),
    row.names = FALSE
  )
  
  sig_overlap_keys <- intersect(m1_sig$key, m2_sig$key)
  sig_overlap <- common[common$key %in% sig_overlap_keys, ]
  
  write.csv(
    sig_overlap,
    file.path(out_dir, paste0(contrast, "_Model1_Model2_global_sig_overlap.csv")),
    row.names = FALSE
  )
  
  cor_logFC <- NA
  if (nrow(common) >= 3) {
    cor_logFC <- suppressWarnings(
      cor(common$logFC_Model1_GTEx, common$logFC_Model2_TCGA, method = "spearman")
    )
  }
  
  summary_list[[contrast]] <- data.frame(
    contrast = contrast,
    n_Model1_tests = nrow(m1),
    n_Model2_tests = nrow(m2),
    n_common_precursor_keys = nrow(common),
    n_Model1_global_FDR005 = nrow(m1_sig),
    n_Model2_global_FDR005 = nrow(m2_sig),
    n_overlap_global_FDR005 = length(sig_overlap_keys),
    n_direction_match_common_keys = sum(common$direction_match, na.rm = TRUE),
    direction_match_fraction = ifelse(nrow(common) > 0, mean(common$direction_match, na.rm = TRUE), NA),
    spearman_logFC_correlation = cor_logFC,
    Model1_top = m1[order(m1$global_FDR_across_all_contrasts, m1$P.Value), "miRNA"][1],
    Model2_top = m2[order(m2$global_FDR_across_all_Model2_tests, m2$P.Value), "miRNA_feature"][1],
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "Model1_Model2_miRSCAPE_comparison_summary.csv"),
  row.names = FALSE
)

print(summary_df)

primary_map <- data.frame(
  mature_candidate = c(
    "hsa-miR-195-5p",
    "hsa-miR-16-5p",
    "hsa-miR-15a-5p",
    "hsa-miR-15b-5p",
    "hsa-miR-649",
    "hsa-miR-6838-5p"
  ),
  model_key = c(
    "hsa-mir-195",
    "hsa-mir-16",
    "hsa-mir-15a",
    "hsa-mir-15b",
    "hsa-mir-649",
    "hsa-mir-6838"
  ),
  stringsAsFactors = FALSE
)

primary_list <- list()

for (contrast in contrasts) {
  m1 <- model1_best[model1_best$contrast == contrast, ]
  m2 <- model2_best[model2_best$contrast == contrast, ]
  
  for (i in seq_len(nrow(primary_map))) {
    key <- primary_map$model_key[i]
    
    m1_row <- m1[m1$key == key, ]
    
    if (key == "hsa-mir-16") {
      m2_row <- m2[grepl("^hsa-mir-16", m2$key), ]
    } else {
      m2_row <- m2[m2$key == key, ]
    }
    
    if (nrow(m1_row) == 0) {
      m1_row <- data.frame(
        miRNA = NA,
        logFC = NA,
        P.Value = NA,
        adj.P.Val = NA,
        global_FDR_across_all_contrasts = NA,
        predicted_expression_direction = NA
      )
    }
    
    if (nrow(m2_row) == 0) {
      m2_row <- data.frame(
        miRNA_feature = NA,
        logFC = NA,
        P.Value = NA,
        adj.P.Val = NA,
        global_FDR_across_all_Model2_tests = NA,
        predicted_expression_direction = NA
      )
    }
    
    m2_row <- m2_row[order(m2_row$global_FDR_across_all_Model2_tests, m2_row$P.Value), ]
    m2_row <- m2_row[1, ]
    
    primary_list[[paste(contrast, i, sep = "_")]] <- data.frame(
      contrast = contrast,
      mature_candidate = primary_map$mature_candidate[i],
      model_key = key,
      Model1_feature = m1_row$miRNA[1],
      Model1_logFC = m1_row$logFC[1],
      Model1_P = m1_row$P.Value[1],
      Model1_within_FDR = m1_row$adj.P.Val[1],
      Model1_global_FDR = m1_row$global_FDR_across_all_contrasts[1],
      Model1_direction = m1_row$predicted_expression_direction[1],
      Model2_feature = m2_row$miRNA_feature[1],
      Model2_logFC = m2_row$logFC[1],
      Model2_P = m2_row$P.Value[1],
      Model2_within_FDR = m2_row$adj.P.Val[1],
      Model2_global_FDR = m2_row$global_FDR_across_all_Model2_tests[1],
      Model2_direction = m2_row$predicted_expression_direction[1],
      stringsAsFactors = FALSE
    )
  }
}

primary_df <- do.call(rbind, primary_list)

write.csv(
  primary_df,
  file.path(out_dir, "primary_candidate_Model1_Model2_comparison.csv"),
  row.names = FALSE
)

print(primary_df)

cat("Done: Model1 vs Model2 miRSCAPE comparison completed.\n")
