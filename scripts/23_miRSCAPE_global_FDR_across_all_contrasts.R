out_dir <- "results/23_miRSCAPE_global_FDR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mir_dir <- "results/22_miRSCAPE_prediction"

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

all_list <- list()

for (contrast in contrasts) {
  file <- file.path(mir_dir, paste0(contrast, "_miRSCAPE_predicted_miRNA_all.csv"))
  x <- read.csv(file, stringsAsFactors = FALSE)
  x$contrast <- contrast
  all_list[[contrast]] <- x
}

all_tests <- do.call(rbind, all_list)

all_tests$global_FDR_across_all_contrasts <- p.adjust(
  all_tests$P.Value,
  method = "BH"
)

all_tests <- all_tests[
  order(all_tests$global_FDR_across_all_contrasts, all_tests$P.Value),
]

write.csv(
  all_tests,
  file.path(out_dir, "miRSCAPE_all_contrasts_global_FDR.csv"),
  row.names = FALSE
)

summary_df <- do.call(
  rbind,
  lapply(contrasts, function(contrast) {
    x <- all_tests[all_tests$contrast == contrast, ]
    x <- x[order(x$global_FDR_across_all_contrasts, x$P.Value), ]
    
    data.frame(
      contrast = contrast,
      n_tests = nrow(x),
      n_global_FDR005 = sum(x$global_FDR_across_all_contrasts < 0.05),
      n_global_FDR010 = sum(x$global_FDR_across_all_contrasts < 0.10),
      top_miRNA = x$miRNA[1],
      top_logFC = x$logFC[1],
      top_P_value = x$P.Value[1],
      top_within_contrast_FDR = x$adj.P.Val[1],
      top_global_FDR = x$global_FDR_across_all_contrasts[1],
      top_direction = x$predicted_expression_direction[1],
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  summary_df,
  file.path(out_dir, "miRSCAPE_global_FDR_summary.csv"),
  row.names = FALSE
)

primary_candidates <- c(
  "hsa-miR-195-5p",
  "hsa-miR-16-5p",
  "hsa-miR-15a-5p",
  "hsa-miR-15b-5p",
  "hsa-miR-649",
  "hsa-miR-6838-5p"
)

primary <- all_tests[all_tests$miRNA %in% primary_candidates, ]

write.csv(
  primary,
  file.path(out_dir, "primary_candidates_miRSCAPE_global_FDR.csv"),
  row.names = FALSE
)

print(summary_df)
print(primary)

cat("Done: miRSCAPE global FDR completed.\n")
