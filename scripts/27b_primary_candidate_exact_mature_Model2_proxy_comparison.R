out_dir <- "results/27_miRSCAPE_Model1_Model2_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model1_file <- "results/23_miRSCAPE_global_FDR/miRSCAPE_all_contrasts_global_FDR.csv"
model2_file <- "results/26_miRSCAPE_Model2_prediction/Model2_all_contrasts_global_FDR.csv"

model1 <- read.csv(model1_file, stringsAsFactors = FALSE)
model2 <- read.csv(model2_file, stringsAsFactors = FALSE)

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

primary_map <- data.frame(
  mature_candidate = c(
    "hsa-miR-195-5p",
    "hsa-miR-16-5p",
    "hsa-miR-15a-5p",
    "hsa-miR-15b-5p",
    "hsa-miR-649",
    "hsa-miR-6838-5p"
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
  stringsAsFactors = FALSE
)

out_list <- list()

for (contrast in contrasts) {
  m1 <- model1[model1$contrast == contrast, ]
  m2 <- model2[model2$contrast == contrast, ]

  for (i in seq_len(nrow(primary_map))) {
    mature <- primary_map$mature_candidate[i]
    p1 <- primary_map$model2_proxy_1[i]
    p2 <- primary_map$model2_proxy_2[i]

    m1_row <- m1[m1$miRNA == mature, ]

    if (nrow(m1_row) == 0) {
      m1_row <- data.frame(
        miRNA = NA,
        logFC = NA,
        P.Value = NA,
        adj.P.Val = NA,
        global_FDR_across_all_contrasts = NA,
        predicted_expression_direction = NA
      )
    } else {
      m1_row <- m1_row[order(m1_row$global_FDR_across_all_contrasts, m1_row$P.Value), ]
      m1_row <- m1_row[1, ]
    }

    proxies <- c(p1, p2)
    proxies <- proxies[!is.na(proxies)]
    m2_row <- m2[m2$miRNA_feature %in% proxies, ]

    if (nrow(m2_row) == 0) {
      m2_row <- data.frame(
        miRNA_feature = NA,
        logFC = NA,
        P.Value = NA,
        adj.P.Val = NA,
        global_FDR_across_all_Model2_tests = NA,
        predicted_expression_direction = NA
      )
    } else {
      m2_row <- m2_row[order(m2_row$global_FDR_across_all_Model2_tests, m2_row$P.Value), ]
      m2_row <- m2_row[1, ]
    }

    out_list[[paste(contrast, i, sep = "_")]] <- data.frame(
      contrast = contrast,
      mature_candidate = mature,
      Model1_exact_mature_feature = m1_row$miRNA[1],
      Model1_logFC = m1_row$logFC[1],
      Model1_P = m1_row$P.Value[1],
      Model1_within_FDR = m1_row$adj.P.Val[1],
      Model1_global_FDR = m1_row$global_FDR_across_all_contrasts[1],
      Model1_direction = m1_row$predicted_expression_direction[1],
      Model2_precursor_proxy_feature = m2_row$miRNA_feature[1],
      Model2_logFC = m2_row$logFC[1],
      Model2_P = m2_row$P.Value[1],
      Model2_within_FDR = m2_row$adj.P.Val[1],
      Model2_global_FDR = m2_row$global_FDR_across_all_Model2_tests[1],
      Model2_direction = m2_row$predicted_expression_direction[1],
      stringsAsFactors = FALSE
    )
  }
}

primary_exact <- do.call(rbind, out_list)

primary_exact$Model1_global_significant <- !is.na(primary_exact$Model1_global_FDR) & primary_exact$Model1_global_FDR < 0.05
primary_exact$Model2_global_significant <- !is.na(primary_exact$Model2_global_FDR) & primary_exact$Model2_global_FDR < 0.05
primary_exact$direction_match <- !is.na(primary_exact$Model1_logFC) &
  !is.na(primary_exact$Model2_logFC) &
  sign(primary_exact$Model1_logFC) == sign(primary_exact$Model2_logFC)

write.csv(
  primary_exact,
  file.path(out_dir, "primary_candidate_Model1_exact_mature_Model2_proxy_comparison.csv"),
  row.names = FALSE
)

summary_exact <- do.call(
  rbind,
  lapply(contrasts, function(contrast) {
    x <- primary_exact[primary_exact$contrast == contrast, ]
    data.frame(
      contrast = contrast,
      n_primary_candidates = nrow(x),
      n_Model1_available_exact_mature = sum(!is.na(x$Model1_exact_mature_feature)),
      n_Model2_available_precursor_proxy = sum(!is.na(x$Model2_precursor_proxy_feature)),
      n_Model1_global_FDR005 = sum(x$Model1_global_significant, na.rm = TRUE),
      n_Model2_global_FDR005 = sum(x$Model2_global_significant, na.rm = TRUE),
      n_both_global_FDR005 = sum(x$Model1_global_significant & x$Model2_global_significant, na.rm = TRUE),
      n_direction_match_where_both_available = sum(x$direction_match, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  summary_exact,
  file.path(out_dir, "primary_candidate_exact_comparison_summary.csv"),
  row.names = FALSE
)

print(summary_exact)
print(primary_exact)

cat("Done: exact mature Model1 vs precursor proxy Model2 primary candidate comparison completed.\n")
