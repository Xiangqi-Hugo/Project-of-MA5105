# Four-model miRNA robustness comparison
#
# Four models compared:
#   Model A: Conservative target-set GSEA
#   Model B: Conservative decoupleR activity inference
#   Model C: miRSCAPE Model1 trained on GTEx Pancreas
#   Model D: miRSCAPE Model2 trained on TCGA-PAAD
#
# Main outputs:
#   results/31_four_model_miRNA_comparison/four_model_method_level_summary.csv
#   results/31_four_model_miRNA_comparison/four_model_overlap_summary.csv
#   results/31_four_model_miRNA_comparison/primary_candidate_four_model_long_table.csv
#   results/31_four_model_miRNA_comparison/primary_candidate_four_model_wide_table.csv

out_dir <- "results/31_four_model_miRNA_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")
filters <- c("min15_max200", "min15_max300", "min15_max500")

gsea_dir <- "results/29_conservative_miRNA_GSEA_filtering"
dec_dir <- "results/30_decoupleR_conservative_miRNA_activity"
m1_dir <- "results/23_miRSCAPE_global_FDR"
m2_dir <- "results/26_miRSCAPE_Model2_prediction"
m12_dir <- "results/27_miRSCAPE_Model1_Model2_comparison"

primary_candidates <- data.frame(
  mature_candidate = c(
    "hsa-miR-195-5p",
    "hsa-miR-16-5p",
    "hsa-miR-15a-5p",
    "hsa-miR-15b-5p",
    "hsa-miR-649",
    "hsa-miR-6838-5p"
  ),
  gsea_decoupleR_ID = c(
    "MIR195_5P",
    "MIR16_5P",
    "MIR15A_5P",
    "MIR15B_5P",
    "MIR649",
    "MIR6838_5P"
  ),
  model2_proxy_1 = c(
    "hsa-mir-195",
    "hsa-mir-16-1",
    "hsa-mir-15a",
    "hsa-mir-15b",
    "hsa-mir-649",
    "hsa-mir-6838"
  ),
  model2_proxy_2 = c(
    NA,
    "hsa-mir-16-2",
    NA,
    NA,
    NA,
    NA
  ),
  key = c(
    "hsa-mir-195",
    "hsa-mir-16",
    "hsa-mir-15a",
    "hsa-mir-15b",
    "hsa-mir-649",
    "hsa-mir-6838"
  ),
  stringsAsFactors = FALSE
)

# --------------------------
# Helper functions
# --------------------------
safe_read <- function(file) {
  if (!file.exists(file)) {
    warning("Missing file: ", file)
    return(NULL)
  }
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

to_key_from_mature <- function(x) {
  y <- tolower(x)
  y <- gsub("-5p$", "", y)
  y <- gsub("-3p$", "", y)
  y
}

to_key_from_term <- function(x) {
  # MIR195_5P -> hsa-mir-195
  # MIR15A_5P -> hsa-mir-15a
  # MIR649 -> hsa-mir-649
  y <- toupper(x)
  y <- gsub("^MIR", "", y)
  y <- gsub("_5P$", "", y)
  y <- gsub("_3P$", "", y)
  y <- tolower(y)
  y <- paste0("hsa-mir-", y)
  y <- gsub("hsa-mir-let", "hsa-let", y)
  y
}

to_key_from_model2 <- function(x) {
  y <- tolower(x)
  y <- gsub("-5p$", "", y)
  y <- gsub("-3p$", "", y)
  # Combine hsa-mir-16-1 and hsa-mir-16-2 at family key level for overlap.
  y <- gsub("^hsa-mir-16-[12]$", "hsa-mir-16", y)
  y
}

get_first <- function(x) {
  if (length(x) == 0) return(NA)
  x[1]
}

direction_from_logfc <- function(logfc, first_group) {
  if (is.na(logfc)) return(NA_character_)
  if (logfc > 0) paste0("higher_in_", first_group) else paste0("lower_in_", first_group)
}

direction_from_nes <- function(nes, first_group) {
  if (is.na(nes)) return(NA_character_)
  if (nes > 0) paste0("targets_up_in_", first_group) else paste0("targets_down_in_", first_group)
}

# --------------------------
# Load Model1 and Model2 global tables once
# --------------------------
model1_all <- safe_read(file.path(m1_dir, "miRSCAPE_all_contrasts_global_FDR.csv"))
model2_all <- safe_read(file.path(m2_dir, "Model2_all_contrasts_global_FDR.csv"))

if (!is.null(model1_all)) {
  model1_all$key <- to_key_from_mature(model1_all$miRNA)
}

if (!is.null(model2_all)) {
  model2_all$key <- to_key_from_model2(model2_all$miRNA_feature)
}

# --------------------------
# Method-level summary and overlap
# --------------------------
method_summary_rows <- list()
overlap_rows <- list()

for (filter in filters) {
  for (contrast in contrasts) {
    first_group <- sub("_vs_.*$", "", contrast)

    # Conservative GSEA
    gsea_file <- file.path(gsea_dir, paste0(contrast, "_", filter, "_all_refiltered.csv"))
    gsea <- safe_read(gsea_file)

    if (!is.null(gsea) && nrow(gsea) > 0) {
      gsea$key <- to_key_from_term(gsea$ID)
      gsea_sig <- gsea[gsea$conservative_FDR < 0.05, ]
      gsea_top <- gsea[order(gsea$conservative_FDR, gsea$P_value_for_refilter), ][1, ]

      method_summary_rows[[paste(filter, contrast, "GSEA", sep = "_")]] <- data.frame(
        filter = filter,
        contrast = contrast,
        model = "A_conservative_GSEA",
        n_tests = nrow(gsea),
        n_global_or_conservative_FDR005 = nrow(gsea_sig),
        n_global_or_conservative_FDR010 = sum(gsea$conservative_FDR < 0.10, na.rm = TRUE),
        top_feature = gsea_top$ID,
        top_key = gsea_top$key,
        top_effect = gsea_top$NES,
        top_FDR = gsea_top$conservative_FDR,
        top_direction = direction_from_nes(gsea_top$NES, first_group),
        stringsAsFactors = FALSE
      )
    } else {
      gsea_sig <- data.frame(key = character())
    }

    # Conservative decoupleR
    dec_file <- file.path(dec_dir, filter, "decoupleR_all_contrasts_global_FDR.csv")
    dec <- safe_read(dec_file)

    if (!is.null(dec) && nrow(dec) > 0) {
      dec <- dec[dec$contrast == contrast, ]
      dec$key <- to_key_from_term(dec$miRNA)
      dec_sig <- dec[dec$global_FDR_across_all_contrasts < 0.05, ]
      dec_top <- dec[order(dec$global_FDR_across_all_contrasts, dec$P.Value), ][1, ]

      method_summary_rows[[paste(filter, contrast, "decoupleR", sep = "_")]] <- data.frame(
        filter = filter,
        contrast = contrast,
        model = "B_conservative_decoupleR",
        n_tests = nrow(dec),
        n_global_or_conservative_FDR005 = nrow(dec_sig),
        n_global_or_conservative_FDR010 = sum(dec$global_FDR_across_all_contrasts < 0.10, na.rm = TRUE),
        top_feature = dec_top$miRNA,
        top_key = dec_top$key,
        top_effect = dec_top$logFC,
        top_FDR = dec_top$global_FDR_across_all_contrasts,
        top_direction = dec_top$activity_direction,
        stringsAsFactors = FALSE
      )
    } else {
      dec_sig <- data.frame(key = character())
    }

    # miRSCAPE Model1
    if (!is.null(model1_all) && nrow(model1_all) > 0) {
      m1 <- model1_all[model1_all$contrast == contrast, ]
      # Collapse by key, keeping best mature arm.
      m1 <- m1[order(m1$key, m1$global_FDR_across_all_contrasts, m1$P.Value), ]
      m1 <- m1[!duplicated(m1$key), ]
      m1_sig <- m1[m1$global_FDR_across_all_contrasts < 0.05, ]
      m1_top <- m1[order(m1$global_FDR_across_all_contrasts, m1$P.Value), ][1, ]

      method_summary_rows[[paste(filter, contrast, "Model1", sep = "_")]] <- data.frame(
        filter = filter,
        contrast = contrast,
        model = "C_miRSCAPE_Model1_GTEx_pancreas",
        n_tests = nrow(m1),
        n_global_or_conservative_FDR005 = nrow(m1_sig),
        n_global_or_conservative_FDR010 = sum(m1$global_FDR_across_all_contrasts < 0.10, na.rm = TRUE),
        top_feature = m1_top$miRNA,
        top_key = m1_top$key,
        top_effect = m1_top$logFC,
        top_FDR = m1_top$global_FDR_across_all_contrasts,
        top_direction = m1_top$predicted_expression_direction,
        stringsAsFactors = FALSE
      )
    } else {
      m1_sig <- data.frame(key = character())
    }

    # miRSCAPE Model2
    if (!is.null(model2_all) && nrow(model2_all) > 0) {
      m2 <- model2_all[model2_all$contrast == contrast, ]
      # Collapse by key, keeping best precursor feature.
      m2 <- m2[order(m2$key, m2$global_FDR_across_all_Model2_tests, m2$P.Value), ]
      m2 <- m2[!duplicated(m2$key), ]
      m2_sig <- m2[m2$global_FDR_across_all_Model2_tests < 0.05, ]
      m2_top <- m2[order(m2$global_FDR_across_all_Model2_tests, m2$P.Value), ][1, ]

      method_summary_rows[[paste(filter, contrast, "Model2", sep = "_")]] <- data.frame(
        filter = filter,
        contrast = contrast,
        model = "D_miRSCAPE_Model2_TCGA_PAAD",
        n_tests = nrow(m2),
        n_global_or_conservative_FDR005 = nrow(m2_sig),
        n_global_or_conservative_FDR010 = sum(m2$global_FDR_across_all_Model2_tests < 0.10, na.rm = TRUE),
        top_feature = m2_top$miRNA_feature,
        top_key = m2_top$key,
        top_effect = m2_top$logFC,
        top_FDR = m2_top$global_FDR_across_all_Model2_tests,
        top_direction = m2_top$predicted_expression_direction,
        stringsAsFactors = FALSE
      )
    } else {
      m2_sig <- data.frame(key = character())
    }

    keys_gsea <- unique(gsea_sig$key)
    keys_dec <- unique(dec_sig$key)
    keys_m1 <- unique(m1_sig$key)
    keys_m2 <- unique(m2_sig$key)

    all_sig_keys <- unique(c(keys_gsea, keys_dec, keys_m1, keys_m2))

    if (length(all_sig_keys) > 0) {
      support_count <- data.frame(
        key = all_sig_keys,
        A_conservative_GSEA = all_sig_keys %in% keys_gsea,
        B_conservative_decoupleR = all_sig_keys %in% keys_dec,
        C_miRSCAPE_Model1_GTEx_pancreas = all_sig_keys %in% keys_m1,
        D_miRSCAPE_Model2_TCGA_PAAD = all_sig_keys %in% keys_m2,
        stringsAsFactors = FALSE
      )

      support_count$n_models_supported <- rowSums(support_count[, 2:5])

      write.csv(
        support_count[order(-support_count$n_models_supported, support_count$key), ],
        file.path(out_dir, paste0(contrast, "_", filter, "_four_model_significant_key_support.csv")),
        row.names = FALSE
      )

      n_at_least_2 <- sum(support_count$n_models_supported >= 2)
      n_at_least_3 <- sum(support_count$n_models_supported >= 3)
      n_all_4 <- sum(support_count$n_models_supported == 4)
    } else {
      n_at_least_2 <- 0
      n_at_least_3 <- 0
      n_all_4 <- 0
    }

    overlap_rows[[paste(filter, contrast, sep = "_")]] <- data.frame(
      filter = filter,
      contrast = contrast,
      n_sig_A_conservative_GSEA = length(keys_gsea),
      n_sig_B_conservative_decoupleR = length(keys_dec),
      n_sig_C_miRSCAPE_Model1 = length(keys_m1),
      n_sig_D_miRSCAPE_Model2 = length(keys_m2),
      overlap_A_B = length(intersect(keys_gsea, keys_dec)),
      overlap_A_C = length(intersect(keys_gsea, keys_m1)),
      overlap_A_D = length(intersect(keys_gsea, keys_m2)),
      overlap_B_C = length(intersect(keys_dec, keys_m1)),
      overlap_B_D = length(intersect(keys_dec, keys_m2)),
      overlap_C_D = length(intersect(keys_m1, keys_m2)),
      n_keys_supported_by_at_least_2_models = n_at_least_2,
      n_keys_supported_by_at_least_3_models = n_at_least_3,
      n_keys_supported_by_all_4_models = n_all_4,
      stringsAsFactors = FALSE
    )
  }
}

method_summary <- do.call(rbind, method_summary_rows)
overlap_summary <- do.call(rbind, overlap_rows)

write.csv(
  method_summary,
  file.path(out_dir, "four_model_method_level_summary.csv"),
  row.names = FALSE
)

write.csv(
  overlap_summary,
  file.path(out_dir, "four_model_overlap_summary.csv"),
  row.names = FALSE
)

# --------------------------
# Primary candidate long table
# --------------------------
primary_rows <- list()

for (filter in filters) {
  for (contrast in contrasts) {
    first_group <- sub("_vs_.*$", "", contrast)

    gsea <- safe_read(file.path(gsea_dir, paste0(contrast, "_", filter, "_all_refiltered.csv")))
    dec <- safe_read(file.path(dec_dir, filter, "decoupleR_all_contrasts_global_FDR.csv"))

    if (!is.null(dec)) dec <- dec[dec$contrast == contrast, ]
    if (!is.null(model1_all)) m1 <- model1_all[model1_all$contrast == contrast, ] else m1 <- NULL
    if (!is.null(model2_all)) m2 <- model2_all[model2_all$contrast == contrast, ] else m2 <- NULL

    for (i in seq_len(nrow(primary_candidates))) {
      mature <- primary_candidates$mature_candidate[i]
      term <- primary_candidates$gsea_decoupleR_ID[i]
      key <- primary_candidates$key[i]
      p1 <- primary_candidates$model2_proxy_1[i]
      p2 <- primary_candidates$model2_proxy_2[i]

      # A: Conservative GSEA
      if (!is.null(gsea)) {
        row <- gsea[gsea$ID == term, ]
      } else {
        row <- data.frame()
      }

      if (nrow(row) == 0) {
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "A_conservative_GSEA",
          feature_used = term,
          key = key,
          available = FALSE,
          effect_type = "NES",
          effect = NA,
          P_value = NA,
          FDR = NA,
          significant = FALSE,
          direction = NA,
          stringsAsFactors = FALSE
        )
      } else {
        row <- row[order(row$conservative_FDR, row$P_value_for_refilter), ][1, ]
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "A_conservative_GSEA",
          feature_used = term,
          key = key,
          available = TRUE,
          effect_type = "NES",
          effect = row$NES,
          P_value = row$P_value_for_refilter,
          FDR = row$conservative_FDR,
          significant = row$conservative_FDR < 0.05,
          direction = direction_from_nes(row$NES, first_group),
          stringsAsFactors = FALSE
        )
      }

      # B: Conservative decoupleR
      if (!is.null(dec)) {
        row <- dec[dec$miRNA == term, ]
      } else {
        row <- data.frame()
      }

      if (nrow(row) == 0) {
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "B_conservative_decoupleR",
          feature_used = term,
          key = key,
          available = FALSE,
          effect_type = "activity_logFC",
          effect = NA,
          P_value = NA,
          FDR = NA,
          significant = FALSE,
          direction = NA,
          stringsAsFactors = FALSE
        )
      } else {
        row <- row[order(row$global_FDR_across_all_contrasts, row$P.Value), ][1, ]
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "B_conservative_decoupleR",
          feature_used = term,
          key = key,
          available = TRUE,
          effect_type = "activity_logFC",
          effect = row$logFC,
          P_value = row$P.Value,
          FDR = row$global_FDR_across_all_contrasts,
          significant = row$global_FDR_across_all_contrasts < 0.05,
          direction = row$activity_direction,
          stringsAsFactors = FALSE
        )
      }

      # C: miRSCAPE Model1 exact mature
      if (!is.null(m1)) {
        row <- m1[m1$miRNA == mature, ]
      } else {
        row <- data.frame()
      }

      if (nrow(row) == 0) {
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "C_miRSCAPE_Model1_GTEx_pancreas",
          feature_used = mature,
          key = key,
          available = FALSE,
          effect_type = "predicted_expression_logFC",
          effect = NA,
          P_value = NA,
          FDR = NA,
          significant = FALSE,
          direction = NA,
          stringsAsFactors = FALSE
        )
      } else {
        row <- row[order(row$global_FDR_across_all_contrasts, row$P.Value), ][1, ]
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "C_miRSCAPE_Model1_GTEx_pancreas",
          feature_used = row$miRNA,
          key = key,
          available = TRUE,
          effect_type = "predicted_expression_logFC",
          effect = row$logFC,
          P_value = row$P.Value,
          FDR = row$global_FDR_across_all_contrasts,
          significant = row$global_FDR_across_all_contrasts < 0.05,
          direction = row$predicted_expression_direction,
          stringsAsFactors = FALSE
        )
      }

      # D: miRSCAPE Model2 precursor proxy
      if (!is.null(m2)) {
        proxies <- c(p1, p2)
        proxies <- proxies[!is.na(proxies)]
        row <- m2[m2$miRNA_feature %in% proxies, ]
      } else {
        row <- data.frame()
      }

      if (nrow(row) == 0) {
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "D_miRSCAPE_Model2_TCGA_PAAD",
          feature_used = paste(c(p1, p2)[!is.na(c(p1, p2))], collapse = ";"),
          key = key,
          available = FALSE,
          effect_type = "predicted_expression_logFC",
          effect = NA,
          P_value = NA,
          FDR = NA,
          significant = FALSE,
          direction = NA,
          stringsAsFactors = FALSE
        )
      } else {
        row <- row[order(row$global_FDR_across_all_Model2_tests, row$P.Value), ][1, ]
        primary_rows[[length(primary_rows) + 1]] <- data.frame(
          filter = filter,
          contrast = contrast,
          mature_candidate = mature,
          model = "D_miRSCAPE_Model2_TCGA_PAAD",
          feature_used = row$miRNA_feature,
          key = key,
          available = TRUE,
          effect_type = "predicted_expression_logFC",
          effect = row$logFC,
          P_value = row$P.Value,
          FDR = row$global_FDR_across_all_Model2_tests,
          significant = row$global_FDR_across_all_Model2_tests < 0.05,
          direction = row$predicted_expression_direction,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

primary_long <- do.call(rbind, primary_rows)

write.csv(
  primary_long,
  file.path(out_dir, "primary_candidate_four_model_long_table.csv"),
  row.names = FALSE
)

# Wide summary: one row per candidate/filter/contrast, support count across four models.
wide_rows <- list()

for (filter in filters) {
  for (contrast in contrasts) {
    for (cand in primary_candidates$mature_candidate) {
      x <- primary_long[
        primary_long$filter == filter &
          primary_long$contrast == contrast &
          primary_long$mature_candidate == cand,
      ]

      wide_rows[[length(wide_rows) + 1]] <- data.frame(
        filter = filter,
        contrast = contrast,
        mature_candidate = cand,
        n_models_available = sum(x$available, na.rm = TRUE),
        n_models_significant = sum(x$significant, na.rm = TRUE),
        significant_models = paste(x$model[x$significant], collapse = ";"),
        available_models = paste(x$model[x$available], collapse = ";"),
        directions_available = paste(
          paste(x$model[x$available], x$direction[x$available], sep = "="),
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
}

primary_wide <- do.call(rbind, wide_rows)

write.csv(
  primary_wide,
  file.path(out_dir, "primary_candidate_four_model_wide_table.csv"),
  row.names = FALSE
)

print(method_summary)
print(overlap_summary)
print(primary_wide)

cat("\nDone: four-model miRNA comparison completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
