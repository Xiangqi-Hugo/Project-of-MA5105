out_dir <- "results/17_decoupleR_global_FDR"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dec_dir <- "results/15_decoupleR_miRNA_activity"

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

all_list <- list()

for (contrast in contrasts) {
  file <- file.path(dec_dir, paste0(contrast, "_decoupleR_activity_all.csv"))
  x <- read.csv(file, stringsAsFactors = FALSE)
  x$contrast <- contrast
  all_list[[contrast]] <- x
}

all_tests <- do.call(rbind, all_list)

all_tests$global_FDR_across_all_contrasts <- p.adjust(
  all_tests$P.Value,
  method = "BH"
)

all_tests <- all_tests[order(all_tests$global_FDR_across_all_contrasts, all_tests$P.Value), ]

write.csv(
  all_tests,
  file.path(out_dir, "decoupleR_all_contrasts_global_FDR.csv"),
  row.names = FALSE
)

summary_df <- aggregate(
  global_FDR_across_all_contrasts ~ contrast,
  data = all_tests,
  FUN = function(x) sum(x < 0.05)
)

colnames(summary_df) <- c("contrast", "n_global_FDR005")

summary_df$n_global_FDR010 <- aggregate(
  global_FDR_across_all_contrasts ~ contrast,
  data = all_tests,
  FUN = function(x) sum(x < 0.10)
)$global_FDR_across_all_contrasts

summary_df$n_tests <- aggregate(
  miRNA ~ contrast,
  data = all_tests,
  FUN = length
)$miRNA

top_by_contrast <- do.call(
  rbind,
  lapply(contrasts, function(contrast) {
    x <- all_tests[all_tests$contrast == contrast, ]
    x <- x[order(x$global_FDR_across_all_contrasts, x$P.Value), ]
    data.frame(
      contrast = contrast,
      top_miRNA = x$miRNA[1],
      top_logFC = x$logFC[1],
      top_P_value = x$P.Value[1],
      top_within_contrast_FDR = x$adj.P.Val[1],
      top_global_FDR = x$global_FDR_across_all_contrasts[1],
      top_activity_direction = x$activity_direction[1],
      stringsAsFactors = FALSE
    )
  })
)

summary_df <- merge(
  summary_df,
  top_by_contrast,
  by = "contrast",
  all.x = TRUE
)

write.csv(
  summary_df,
  file.path(out_dir, "decoupleR_global_FDR_summary.csv"),
  row.names = FALSE
)

print(summary_df)

# Also export the six primary candidates with global FDR.
primary_candidates <- c(
  "MIR195_5P",
  "MIR16_5P",
  "MIR15B_5P",
  "MIR15A_5P",
  "MIR649",
  "MIR6838_5P"
)

primary <- all_tests[all_tests$miRNA %in% primary_candidates, ]

primary <- primary[order(primary$contrast, primary$P.Value), ]

write.csv(
  primary,
  file.path(out_dir, "primary_candidates_decoupleR_global_FDR.csv"),
  row.names = FALSE
)

print(primary)

message("Done: decoupleR global FDR across all contrasts completed.")
