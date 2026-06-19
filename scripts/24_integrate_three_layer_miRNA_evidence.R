out_dir <- "results/24_integrated_three_layer_miRNA_evidence"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strict_dir <- "results/14_strict_miRNA_target_GSEA"
dec_dir <- "results/17_decoupleR_global_FDR"
mirscape_dir <- "results/23_miRSCAPE_global_FDR"

to_hsa_name <- function(x) {
  y <- x
  
  y <- gsub("^MIRLET", "LET", y)
  y <- gsub("^MIR", "MIR", y)
  
  y <- tolower(y)
  y <- gsub("_", "-", y)
  
  y <- ifelse(
    grepl("^mir", y),
    paste0("hsa-", sub("^mir", "miR-", y)),
    y
  )
  
  y <- ifelse(
    grepl("^let", y),
    paste0("hsa-", y),
    y
  )
  
  y <- gsub("-5p$", "-5p", y)
  y <- gsub("-3p$", "-3p", y)
  
  y
}

strict_to_hsa <- function(x) {
  y <- x
  y <- gsub("^MIR", "", y)
  y <- tolower(y)
  y <- gsub("_", "-", y)
  y <- paste0("hsa-miR-", y)
  y <- gsub("hsa-miR-let", "hsa-let", y)
  y
}

contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

summary_list <- list()

for (contrast in contrasts) {
  message("Integrating contrast: ", contrast)
  
  strict_file <- file.path(strict_dir, paste0(contrast, "_strict_miRNA_GSEA_FDR005.csv"))
  dec_file <- file.path(dec_dir, "decoupleR_all_contrasts_global_FDR.csv")
  mir_file <- file.path(mirscape_dir, "miRSCAPE_all_contrasts_global_FDR.csv")
  
  strict <- read.csv(strict_file, stringsAsFactors = FALSE)
  dec <- read.csv(dec_file, stringsAsFactors = FALSE)
  mir <- read.csv(mir_file, stringsAsFactors = FALSE)
  
  dec <- dec[dec$contrast == contrast, ]
  mir <- mir[mir$contrast == contrast, ]
  
  if (nrow(strict) > 0) {
    strict_tbl <- data.frame(
      strict_ID = strict$ID,
      miRNA_hsa = strict_to_hsa(strict$ID),
      strict_NES = strict$NES,
      strict_FDR = strict$p.adjust,
      strict_target_direction = strict$target_direction,
      strict_inferred_activity = strict$inferred_miRNA_activity,
      stringsAsFactors = FALSE
    )
  } else {
    strict_tbl <- data.frame(
      strict_ID = character(),
      miRNA_hsa = character(),
      strict_NES = numeric(),
      strict_FDR = numeric(),
      strict_target_direction = character(),
      strict_inferred_activity = character(),
      stringsAsFactors = FALSE
    )
  }
  
  dec_tbl <- data.frame(
    miRNA_hsa = strict_to_hsa(dec$miRNA),
    decoupleR_logFC = dec$logFC,
    decoupleR_P = dec$P.Value,
    decoupleR_within_FDR = dec$adj.P.Val,
    decoupleR_global_FDR = dec$global_FDR_across_all_contrasts,
    decoupleR_direction = dec$activity_direction,
    stringsAsFactors = FALSE
  )
  
  mir_tbl <- data.frame(
    miRNA_hsa = mir$miRNA,
    miRSCAPE_logFC = mir$logFC,
    miRSCAPE_P = mir$P.Value,
    miRSCAPE_within_FDR = mir$adj.P.Val,
    miRSCAPE_global_FDR = mir$global_FDR_across_all_contrasts,
    miRSCAPE_direction = mir$predicted_expression_direction,
    stringsAsFactors = FALSE
  )
  
  # Union of strict significant hits and miRSCAPE globally significant hits.
  mir_sig <- mir_tbl[mir_tbl$miRSCAPE_global_FDR < 0.05, ]
  
  all_mirnas <- unique(c(strict_tbl$miRNA_hsa, mir_sig$miRNA_hsa))
  
  integrated <- data.frame(miRNA_hsa = all_mirnas, stringsAsFactors = FALSE)
  
  integrated <- merge(integrated, strict_tbl, by = "miRNA_hsa", all.x = TRUE)
  integrated <- merge(integrated, dec_tbl, by = "miRNA_hsa", all.x = TRUE)
  integrated <- merge(integrated, mir_tbl, by = "miRNA_hsa", all.x = TRUE)
  
  integrated$strict_GSEA_significant <- !is.na(integrated$strict_FDR) & integrated$strict_FDR < 0.05
  integrated$decoupleR_global_significant <- !is.na(integrated$decoupleR_global_FDR) & integrated$decoupleR_global_FDR < 0.05
  integrated$miRSCAPE_global_significant <- !is.na(integrated$miRSCAPE_global_FDR) & integrated$miRSCAPE_global_FDR < 0.05
  
  integrated$evidence_class <- ifelse(
    integrated$strict_GSEA_significant & integrated$decoupleR_global_significant & integrated$miRSCAPE_global_significant,
    "three_layer_global_support",
    ifelse(
      integrated$strict_GSEA_significant & integrated$miRSCAPE_global_significant,
      "strict_GSEA_plus_miRSCAPE_global",
      ifelse(
        integrated$strict_GSEA_significant & integrated$decoupleR_global_significant,
        "strict_GSEA_plus_decoupleR_global",
        ifelse(
          integrated$strict_GSEA_significant,
          "strict_GSEA_only",
          ifelse(
            integrated$miRSCAPE_global_significant,
            "miRSCAPE_global_only",
            "other"
          )
        )
      )
    )
  )
  
  integrated <- integrated[order(
    integrated$evidence_class,
    integrated$strict_FDR,
    integrated$miRSCAPE_global_FDR,
    integrated$decoupleR_global_FDR
  ), ]
  
  write.csv(
    integrated,
    file.path(out_dir, paste0(contrast, "_three_layer_miRNA_evidence.csv")),
    row.names = FALSE
  )
  
  primary <- c(
    "hsa-miR-195-5p",
    "hsa-miR-16-5p",
    "hsa-miR-15a-5p",
    "hsa-miR-15b-5p",
    "hsa-miR-649",
    "hsa-miR-6838-5p"
  )
  
  primary_tbl <- integrated[integrated$miRNA_hsa %in% primary, ]
  
  write.csv(
    primary_tbl,
    file.path(out_dir, paste0(contrast, "_primary_candidate_three_layer_evidence.csv")),
    row.names = FALSE
  )
  
  one <- data.frame(
    contrast = contrast,
    n_integrated_miRNAs = nrow(integrated),
    n_strict_GSEA_significant = sum(integrated$strict_GSEA_significant, na.rm = TRUE),
    n_decoupleR_global_significant = sum(integrated$decoupleR_global_significant, na.rm = TRUE),
    n_miRSCAPE_global_significant = sum(integrated$miRSCAPE_global_significant, na.rm = TRUE),
    n_three_layer_global_support = sum(integrated$evidence_class == "three_layer_global_support", na.rm = TRUE),
    n_strict_plus_miRSCAPE = sum(integrated$evidence_class == "strict_GSEA_plus_miRSCAPE_global", na.rm = TRUE),
    n_strict_only = sum(integrated$evidence_class == "strict_GSEA_only", na.rm = TRUE),
    n_miRSCAPE_only = sum(integrated$evidence_class == "miRSCAPE_global_only", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  summary_list[[contrast]] <- one
}

summary_df <- do.call(rbind, summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "three_layer_miRNA_evidence_summary.csv"),
  row.names = FALSE
)

print(summary_df)

cat("Done: three-layer miRNA evidence integration completed.\n")
