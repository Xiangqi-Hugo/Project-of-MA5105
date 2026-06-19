# Step 33: GO pathway threshold sensitivity analysis
#
# Purpose:
# Test whether the pathway-level interpretation is stable across beta-cell count thresholds.
#
# Inputs:
#   results/32_beta_threshold_sensitivity/beta_min50/*_edgeR_all.csv
#   results/32_beta_threshold_sensitivity/beta_min100/*_edgeR_all.csv
#   results/32_beta_threshold_sensitivity/beta_min200/*_edgeR_all.csv
#
# Outputs:
#   results/33_pathway_threshold_sensitivity/GO_GSEA_threshold_summary.csv
#   results/33_pathway_threshold_sensitivity/GO_GSEA_key_pathway_tracking.csv
#   results/33_pathway_threshold_sensitivity/GO_GSEA_top20_overlap_vs_100.csv
#   results/33_pathway_threshold_sensitivity/<threshold>/<contrast>_GO_BP_GSEA_all.csv
#   results/33_pathway_threshold_sensitivity/<threshold>/<contrast>_GO_BP_GSEA_FDR005.csv

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  pkgs <- c("clusterProfiler", "org.Hs.eg.db", "dplyr")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

in_dir <- "results/32_beta_threshold_sensitivity"
out_dir <- "results/33_pathway_threshold_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

thresholds <- c(50, 100, 200)
contrasts <- c("T2D_vs_ND", "PD_vs_ND", "T2D_vs_PD")

key_patterns <- c(
  "organic acid catabolic process",
  "small molecule catabolic process",
  "generation of precursor metabolites and energy",
  "cellular respiration",
  "aerobic respiration",
  "oxidation of organic compounds",
  "fatty acid metabolic process",
  "lipid transport",
  "hormone secretion",
  "peptide secretion",
  "hormone transport",
  "proton transmembrane transport",
  "proton motive force-driven mitochondrial ATP synthesis",
  "mitochondrial ATP synthesis",
  "ribosome biogenesis",
  "ribonucleoprotein complex biogenesis",
  "RNA splicing",
  "chromosome segregation"
)

run_one_gsea <- function(threshold, contrast_name) {
  threshold_label <- paste0("beta_min", threshold)
  threshold_dir <- file.path(out_dir, threshold_label)
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)

  de_file <- file.path(in_dir, threshold_label, paste0(contrast_name, "_edgeR_all.csv"))

  if (!file.exists(de_file)) {
    warning("Missing DE file: ", de_file)
    return(NULL)
  }

  de <- read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (!all(c("gene", "logFC", "PValue") %in% colnames(de))) {
    stop("DE file does not contain gene, logFC, and PValue columns: ", de_file)
  }

  de <- de[!is.na(de$gene) & !is.na(de$logFC) & !is.na(de$PValue), ]
  de <- de[de$gene != "", ]
  de$PValue[de$PValue == 0] <- min(de$PValue[de$PValue > 0], na.rm = TRUE)

  de$rank_metric <- sign(de$logFC) * (-log10(de$PValue))

  # Remove duplicated gene symbols by keeping the strongest absolute rank.
  de <- de[order(abs(de$rank_metric), decreasing = TRUE), ]
  de <- de[!duplicated(de$gene), ]

  # Map gene symbols to Entrez IDs.
  suppressMessages({
    map <- bitr(
      de$gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  })

  de2 <- merge(
    de,
    map,
    by.x = "gene",
    by.y = "SYMBOL",
    all = FALSE
  )

  # If multiple symbols map to same Entrez ID, keep strongest rank.
  de2 <- de2[order(abs(de2$rank_metric), decreasing = TRUE), ]
  de2 <- de2[!duplicated(de2$ENTREZID), ]

  gene_list <- de2$rank_metric
  names(gene_list) <- de2$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)

  if (length(gene_list) < 1000) {
    warning("Too few mapped genes for GSEA: ", threshold_label, " ", contrast_name)
    return(NULL)
  }

  message("Running gseGO: ", threshold_label, " ", contrast_name, " with ", length(gene_list), " genes")

  gsea <- suppressMessages(
    gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
  )

  res <- as.data.frame(gsea)

  if (nrow(res) == 0) {
    warning("No GSEA rows returned for: ", threshold_label, " ", contrast_name)
    return(NULL)
  }

  first_group <- sub("_vs_.*$", "", contrast_name)
  res$threshold <- threshold
  res$contrast <- contrast_name
  res$direction <- ifelse(
    res$NES > 0,
    paste0("higher_in_", first_group),
    paste0("lower_in_", first_group)
  )

  res <- res[order(res$p.adjust, res$pvalue), ]

  write.csv(
    res,
    file.path(threshold_dir, paste0(contrast_name, "_GO_BP_GSEA_all.csv")),
    row.names = FALSE
  )

  write.csv(
    res[res$p.adjust < 0.05, ],
    file.path(threshold_dir, paste0(contrast_name, "_GO_BP_GSEA_FDR005.csv")),
    row.names = FALSE
  )

  write.csv(
    res[res$p.adjust < 0.10, ],
    file.path(threshold_dir, paste0(contrast_name, "_GO_BP_GSEA_FDR010.csv")),
    row.names = FALSE
  )

  return(res)
}

all_results <- list()
summary_rows <- list()
key_rows <- list()

for (threshold in thresholds) {
  for (contrast_name in contrasts) {
    res <- run_one_gsea(threshold, contrast_name)

    if (is.null(res)) next

    key <- paste0("beta_min", threshold, "_", contrast_name)
    all_results[[key]] <- res

    sig005 <- res[res$p.adjust < 0.05, ]
    sig010 <- res[res$p.adjust < 0.10, ]

    top <- res[order(res$p.adjust, res$pvalue), ][1, ]

    summary_rows[[key]] <- data.frame(
      threshold = threshold,
      contrast = contrast_name,
      n_terms_tested = nrow(res),
      n_FDR005 = nrow(sig005),
      n_FDR010 = nrow(sig010),
      top_ID = top$ID,
      top_Description = top$Description,
      top_NES = top$NES,
      top_FDR = top$p.adjust,
      top_direction = top$direction,
      stringsAsFactors = FALSE
    )

    for (pattern in key_patterns) {
      hit <- res[grepl(pattern, res$Description, ignore.case = TRUE, fixed = TRUE), ]

      if (nrow(hit) == 0) {
        key_rows[[paste(key, pattern, sep = "__")]] <- data.frame(
          threshold = threshold,
          contrast = contrast_name,
          key_pattern = pattern,
          found = FALSE,
          best_ID = NA,
          best_Description = NA,
          NES = NA,
          FDR = NA,
          direction = NA,
          stringsAsFactors = FALSE
        )
      } else {
        hit <- hit[order(hit$p.adjust, hit$pvalue), ][1, ]
        key_rows[[paste(key, pattern, sep = "__")]] <- data.frame(
          threshold = threshold,
          contrast = contrast_name,
          key_pattern = pattern,
          found = TRUE,
          best_ID = hit$ID,
          best_Description = hit$Description,
          NES = hit$NES,
          FDR = hit$p.adjust,
          direction = hit$direction,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
key_df <- do.call(rbind, key_rows)

write.csv(
  summary_df,
  file.path(out_dir, "GO_GSEA_threshold_summary.csv"),
  row.names = FALSE
)

write.csv(
  key_df,
  file.path(out_dir, "GO_GSEA_key_pathway_tracking.csv"),
  row.names = FALSE
)

# Compare top 20 significant pathways against threshold 100.
overlap_rows <- list()

for (contrast_name in contrasts) {
  ref_key <- paste0("beta_min100_", contrast_name)
  ref <- all_results[[ref_key]]

  if (is.null(ref)) next

  ref_top <- head(ref[ref$p.adjust < 0.05, "Description"], 20)

  for (threshold in thresholds) {
    cur_key <- paste0("beta_min", threshold, "_", contrast_name)
    cur <- all_results[[cur_key]]

    if (is.null(cur)) next

    cur_top <- head(cur[cur$p.adjust < 0.05, "Description"], 20)

    inter <- intersect(ref_top, cur_top)
    uni <- union(ref_top, cur_top)

    overlap_rows[[paste(threshold, contrast_name, sep = "_")]] <- data.frame(
      comparison = paste0("beta_min", threshold, "_vs_beta_min100"),
      threshold = threshold,
      contrast = contrast_name,
      n_top20_threshold100 = length(ref_top),
      n_top20_current_threshold = length(cur_top),
      n_overlap = length(inter),
      jaccard = ifelse(length(uni) > 0, length(inter) / length(uni), NA),
      overlap_terms = paste(inter, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

overlap_df <- do.call(rbind, overlap_rows)

write.csv(
  overlap_df,
  file.path(out_dir, "GO_GSEA_top20_overlap_vs_100.csv"),
  row.names = FALSE
)

# Stable significant GO terms across all thresholds.
stable_rows <- list()

for (contrast_name in contrasts) {
  term_sets <- list()

  for (threshold in thresholds) {
    key <- paste0("beta_min", threshold, "_", contrast_name)
    res <- all_results[[key]]
    if (is.null(res)) {
      term_sets[[as.character(threshold)]] <- character()
    } else {
      term_sets[[as.character(threshold)]] <- res[res$p.adjust < 0.05, "Description"]
    }
  }

  stable_terms <- Reduce(intersect, term_sets)
  union_terms <- Reduce(union, term_sets)

  stable_rows[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    n_stable_GO_terms_all_thresholds = length(stable_terms),
    n_union_GO_terms_any_threshold = length(union_terms),
    stable_fraction_of_union = ifelse(length(union_terms) > 0, length(stable_terms) / length(union_terms), NA),
    stable_terms_top30 = paste(head(stable_terms, 30), collapse = ";"),
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(Description = stable_terms),
    file.path(out_dir, paste0(contrast_name, "_stable_GO_terms_all_thresholds.csv")),
    row.names = FALSE
  )
}

stable_df <- do.call(rbind, stable_rows)

write.csv(
  stable_df,
  file.path(out_dir, "stable_GO_terms_all_thresholds_summary.csv"),
  row.names = FALSE
)

print(summary_df)
print(overlap_df)
print(stable_df)

cat("\nDone: GO pathway threshold sensitivity analysis completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
