# Step 34: miRNA-target-pathway indirect effect analysis
#
# Purpose:
# Test whether inferred miRNA candidates could indirectly explain key pathway changes
# through their target genes.
#
# This is not a causal test.
# It is a direction-consistency and target-overlap analysis.
#
# Main question:
# Do miRNA target genes overlap with key disease pathways?
# Are those overlapping target genes higher or lower in the expected disease group?
#
# Inputs:
#   results/14_strict_miRNA_target_GSEA/strict_miRNA_TERM2GENE.csv
#   results/32_beta_threshold_sensitivity/beta_min100/*_edgeR_all.csv
#   results/33_pathway_threshold_sensitivity/GO_GSEA_key_pathway_tracking.csv
#   results/31_four_model_miRNA_comparison/primary_candidate_four_model_wide_table.csv  optional
#   results/31_four_model_miRNA_comparison/four_model_method_level_summary.csv          optional
#
# Outputs:
#   results/34_miRNA_indirect_pathway_effect/selected_miRNA_terms_used.csv
#   results/34_miRNA_indirect_pathway_effect/selected_key_GO_pathways_used.csv
#   results/34_miRNA_indirect_pathway_effect/miRNA_target_pathway_overlap_table.csv
#   results/34_miRNA_indirect_pathway_effect/miRNA_target_pathway_direction_consistency.csv
#   results/34_miRNA_indirect_pathway_effect/primary_candidate_indirect_effect_summary.csv
#   results/34_miRNA_indirect_pathway_effect/report_ready_key_indirect_findings.csv

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  pkgs <- c("org.Hs.eg.db", "AnnotationDbi", "dplyr", "stringr")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(stringr)
})

out_dir <- "results/34_miRNA_indirect_pathway_effect"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

term2gene_file <- "results/14_strict_miRNA_target_GSEA/strict_miRNA_TERM2GENE.csv"
key_pathway_file <- "results/33_pathway_threshold_sensitivity/GO_GSEA_key_pathway_tracking.csv"
primary_four_model_file <- "results/31_four_model_miRNA_comparison/primary_candidate_four_model_wide_table.csv"
four_model_method_file <- "results/31_four_model_miRNA_comparison/four_model_method_level_summary.csv"

if (!file.exists(term2gene_file)) {
  stop("Missing strict miRNA TERM2GENE file: ", term2gene_file)
}
if (!file.exists(key_pathway_file)) {
  stop("Missing key pathway tracking file: ", key_pathway_file)
}

find_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) return(hit[1])
  lower_names <- tolower(colnames(df))
  lower_candidates <- tolower(candidates)
  idx <- match(lower_candidates, lower_names)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(colnames(df)[idx[1]])
  if (required) stop("Cannot find required column. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

canonical_mirna <- function(x) {
  x <- as.character(x)
  x <- gsub("^hsa[-_]", "", x, ignore.case = TRUE)
  x <- gsub("^mirna_", "", x, ignore.case = TRUE)
  x <- gsub("^miRNA_", "", x, ignore.case = TRUE)
  x <- gsub("^MIRNA_", "", x, ignore.case = TRUE)
  x <- gsub("^miR-", "MIR", x, ignore.case = TRUE)
  x <- gsub("^mir-", "MIR", x, ignore.case = TRUE)
  x <- gsub("^MIR-", "MIR", x, ignore.case = TRUE)
  x <- gsub("-", "_", x)
  x <- gsub("\\.", "_", x)
  toupper(x)
}

pretty_mirna <- function(term) {
  term <- canonical_mirna(term)
  x <- tolower(term)
  x <- gsub("_", "-", x)
  x <- gsub("^mir", "hsa-miR-", x)
  x
}

first_group_from_contrast <- function(contrast) {
  sub("_vs_.*$", "", contrast)
}

opposite_direction <- function(direction, first_group) {
  if (is.na(direction) || direction == "") return(NA_character_)
  if (direction == paste0("higher_in_", first_group)) return(paste0("lower_in_", first_group))
  if (direction == paste0("lower_in_", first_group)) return(paste0("higher_in_", first_group))
  NA_character_
}

# ----------------------------
# 1. Load strict miRNA target network
# ----------------------------

term2gene <- read.csv(term2gene_file, stringsAsFactors = FALSE, check.names = FALSE)

term_col <- find_col(
  term2gene,
  c("term", "TERM", "gs_name", "gs", "miRNA", "mirna", "source", "ID")
)

gene_col <- find_col(
  term2gene,
  c("gene", "GENE", "target", "TARGET", "target_gene", "gene_symbol", "SYMBOL")
)

network <- data.frame(
  miRNA_term_raw = as.character(term2gene[[term_col]]),
  target_gene = as.character(term2gene[[gene_col]]),
  stringsAsFactors = FALSE
)

network <- network[!is.na(network$miRNA_term_raw) & !is.na(network$target_gene), ]
network <- network[network$miRNA_term_raw != "" & network$target_gene != "", ]
network$miRNA_term <- canonical_mirna(network$miRNA_term_raw)
network$target_gene <- toupper(network$target_gene)

# ----------------------------
# 2. Select miRNAs to test
# ----------------------------

primary_terms <- canonical_mirna(c(
  "MIR195_5P",
  "MIR16_5P",
  "MIR15A_5P",
  "MIR15B_5P",
  "MIR649",
  "MIR6838_5P",
  "MIR144_3P",
  "MIR302C_5P",
  "MIR6507_5P",
  "MIR381_3P",
  "MIR338_5P",
  "MIR300"
))

method_top_terms <- character()

if (file.exists(four_model_method_file)) {
  method_summary <- read.csv(four_model_method_file, stringsAsFactors = FALSE, check.names = FALSE)
  if ("top_feature" %in% colnames(method_summary)) {
    method_top_terms <- canonical_mirna(method_summary$top_feature)
  }
}

selected_terms <- unique(c(primary_terms, method_top_terms))
selected_terms <- selected_terms[selected_terms %in% unique(network$miRNA_term)]

selected_mirna_df <- data.frame(
  miRNA_term = selected_terms,
  pretty_name = pretty_mirna(selected_terms),
  n_targets_in_network = as.integer(table(network$miRNA_term)[selected_terms]),
  is_primary_candidate = selected_terms %in% primary_terms,
  stringsAsFactors = FALSE
)

selected_mirna_df <- selected_mirna_df[order(!selected_mirna_df$is_primary_candidate, selected_mirna_df$miRNA_term), ]

write.csv(
  selected_mirna_df,
  file.path(out_dir, "selected_miRNA_terms_used.csv"),
  row.names = FALSE
)

# ----------------------------
# 3. Load key stable GO pathways
# ----------------------------

key_pathways <- read.csv(key_pathway_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("threshold", "contrast", "key_pattern", "found", "best_ID", "best_Description", "NES", "FDR", "direction")
missing_cols <- setdiff(required_cols, colnames(key_pathways))
if (length(missing_cols) > 0) {
  stop("Missing columns in key pathway file: ", paste(missing_cols, collapse = ", "))
}

# Use the main threshold 100.
# Keep stable and report-relevant pathways only.
pathway_keep_patterns <- c(
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

selected_pathways <- key_pathways[
  key_pathways$threshold == 100 &
    key_pathways$found == TRUE &
    key_pathways$key_pattern %in% pathway_keep_patterns,
]

# Keep FDR < 0.10 to include near-significant beta-cell secretory terms.
selected_pathways <- selected_pathways[!is.na(selected_pathways$FDR) & selected_pathways$FDR < 0.10, ]

selected_pathways <- selected_pathways[!duplicated(
  paste(selected_pathways$contrast, selected_pathways$best_ID, sep = "__")
), ]

write.csv(
  selected_pathways,
  file.path(out_dir, "selected_key_GO_pathways_used.csv"),
  row.names = FALSE
)

# ----------------------------
# 4. Build GO pathway gene sets
# ----------------------------

go_ids <- unique(selected_pathways$best_ID)

go_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_ids,
    keytype = "GOALL",
    columns = c("SYMBOL", "GOALL", "ONTOLOGYALL")
  )
)

go_map <- go_map[!is.na(go_map$SYMBOL) & !is.na(go_map$GOALL), ]
go_map <- go_map[go_map$GOALL %in% go_ids, ]
go_map <- go_map[is.na(go_map$ONTOLOGYALL) | go_map$ONTOLOGYALL == "BP", ]
go_map$SYMBOL <- toupper(go_map$SYMBOL)

# ----------------------------
# 5. Optional model direction evidence
# ----------------------------

model_direction_df <- data.frame()

if (file.exists(primary_four_model_file)) {
  primary_wide <- read.csv(primary_four_model_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (all(c("filter", "contrast", "mature_candidate", "n_models_significant", "significant_models", "directions_available") %in% colnames(primary_wide))) {
    # Use min15_max500 because it gives the widest primary-candidate availability.
    primary_wide2 <- primary_wide[primary_wide$filter == "min15_max500", ]
    primary_wide2$miRNA_term <- canonical_mirna(primary_wide2$mature_candidate)

    model_direction_df <- primary_wide2[, c(
      "contrast",
      "mature_candidate",
      "miRNA_term",
      "n_models_available",
      "n_models_significant",
      "significant_models",
      "available_models",
      "directions_available"
    )]
  }
}

# ----------------------------
# 6. Run indirect effect table
# ----------------------------

de_dir <- "results/32_beta_threshold_sensitivity/beta_min100"
contrasts <- unique(selected_pathways$contrast)

rows <- list()

for (contrast_name in contrasts) {
  de_file <- file.path(de_dir, paste0(contrast_name, "_edgeR_all.csv"))

  if (!file.exists(de_file)) {
    warning("Missing DE file: ", de_file)
    next
  }

  de <- read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("gene", "logFC", "FDR", "PValue") %in% colnames(de))) {
    stop("DE file missing required columns: ", de_file)
  }

  de$gene <- toupper(de$gene)
  de <- de[!is.na(de$gene) & de$gene != "", ]
  de <- de[!duplicated(de$gene), ]

  universe <- unique(de$gene)
  de_logfc <- de$logFC
  names(de_logfc) <- de$gene

  first_group <- first_group_from_contrast(contrast_name)

  pathways_this <- selected_pathways[selected_pathways$contrast == contrast_name, ]

  for (i in seq_len(nrow(pathways_this))) {
    p <- pathways_this[i, ]
    go_id <- p$best_ID

    pathway_genes <- unique(go_map$SYMBOL[go_map$GOALL == go_id])
    pathway_genes <- intersect(pathway_genes, universe)

    if (length(pathway_genes) == 0) next

    for (mir in selected_terms) {
      target_genes <- unique(network$target_gene[network$miRNA_term == mir])
      target_genes <- intersect(target_genes, universe)

      if (length(target_genes) == 0) next

      overlap_genes <- intersect(target_genes, pathway_genes)
      n_overlap <- length(overlap_genes)

      a <- n_overlap
      b <- length(target_genes) - n_overlap
      c <- length(pathway_genes) - n_overlap
      d <- length(universe) - length(union(target_genes, pathway_genes))

      fisher_p <- NA_real_
      odds_ratio <- NA_real_

      if (all(c(a, b, c, d) >= 0)) {
        ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
        fisher_p <- ft$p.value
        odds_ratio <- unname(ft$estimate)
      }

      overlap_logfc <- de_logfc[overlap_genes]
      target_logfc <- de_logfc[target_genes]

      overlap_mean_logFC <- ifelse(n_overlap > 0, mean(overlap_logfc, na.rm = TRUE), NA)
      overlap_median_logFC <- ifelse(n_overlap > 0, median(overlap_logfc, na.rm = TRUE), NA)
      overlap_frac_up <- ifelse(n_overlap > 0, mean(overlap_logfc > 0, na.rm = TRUE), NA)
      overlap_frac_down <- ifelse(n_overlap > 0, mean(overlap_logfc < 0, na.rm = TRUE), NA)

      target_mean_logFC <- mean(target_logfc, na.rm = TRUE)
      target_median_logFC <- median(target_logfc, na.rm = TRUE)
      target_frac_up <- mean(target_logfc > 0, na.rm = TRUE)
      target_frac_down <- mean(target_logfc < 0, na.rm = TRUE)

      overlap_target_direction <- ifelse(
        is.na(overlap_mean_logFC),
        NA,
        ifelse(overlap_mean_logFC > 0, paste0("higher_in_", first_group), paste0("lower_in_", first_group))
      )

      whole_target_direction <- ifelse(
        target_mean_logFC > 0,
        paste0("higher_in_", first_group),
        paste0("lower_in_", first_group)
      )

      target_implied_miRNA_direction <- opposite_direction(overlap_target_direction, first_group)

      pathway_direction <- as.character(p$direction)
      pathway_target_direction_consistent <- !is.na(overlap_target_direction) &&
        overlap_target_direction == pathway_direction

      indirect_logic <- NA_character_

      if (n_overlap < 3) {
        indirect_logic <- "too_few_overlap_targets"
      } else if (pathway_target_direction_consistent && overlap_target_direction == paste0("lower_in_", first_group)) {
        indirect_logic <- paste0("target_genes_lower_with_pathway; compatible_with_higher_miRNA_in_", first_group)
      } else if (pathway_target_direction_consistent && overlap_target_direction == paste0("higher_in_", first_group)) {
        indirect_logic <- paste0("target_genes_higher_with_pathway; compatible_with_lower_miRNA_in_", first_group)
      } else if (!pathway_target_direction_consistent) {
        indirect_logic <- "target_direction_opposes_pathway_direction"
      }

      model_row <- model_direction_df[
        model_direction_df$contrast == contrast_name &
          model_direction_df$miRNA_term == mir,
      ]

      if (nrow(model_row) > 0) {
        model_directions <- paste(unique(model_row$directions_available), collapse = " | ")
        n_models_available <- model_row$n_models_available[1]
        n_models_significant <- model_row$n_models_significant[1]
        significant_models <- model_row$significant_models[1]
      } else {
        model_directions <- NA_character_
        n_models_available <- NA
        n_models_significant <- NA
        significant_models <- NA_character_
      }

      rows[[paste(contrast_name, go_id, mir, sep = "__")]] <- data.frame(
        contrast = contrast_name,
        first_group = first_group,
        miRNA_term = mir,
        miRNA_pretty_name = pretty_mirna(mir),
        is_primary_candidate = mir %in% primary_terms,
        pathway_key_pattern = p$key_pattern,
        pathway_GO_ID = go_id,
        pathway_Description = p$best_Description,
        pathway_NES = p$NES,
        pathway_FDR = p$FDR,
        pathway_direction = pathway_direction,
        universe_n = length(universe),
        miRNA_target_n_in_universe = length(target_genes),
        pathway_gene_n_in_universe = length(pathway_genes),
        overlap_target_pathway_n = n_overlap,
        overlap_genes = paste(overlap_genes, collapse = ";"),
        fisher_OR = odds_ratio,
        fisher_P = fisher_p,
        target_mean_logFC_all_targets = target_mean_logFC,
        target_median_logFC_all_targets = target_median_logFC,
        target_frac_up_all_targets = target_frac_up,
        target_frac_down_all_targets = target_frac_down,
        target_direction_all_targets = whole_target_direction,
        overlap_mean_logFC = overlap_mean_logFC,
        overlap_median_logFC = overlap_median_logFC,
        overlap_frac_up = overlap_frac_up,
        overlap_frac_down = overlap_frac_down,
        overlap_target_direction = overlap_target_direction,
        target_implied_miRNA_direction = target_implied_miRNA_direction,
        pathway_target_direction_consistent = pathway_target_direction_consistent,
        indirect_logic = indirect_logic,
        n_models_available_for_primary_candidate = n_models_available,
        n_models_significant_for_primary_candidate = n_models_significant,
        significant_models_for_primary_candidate = significant_models,
        model_directions_for_primary_candidate = model_directions,
        stringsAsFactors = FALSE
      )
    }
  }
}

overlap_table <- do.call(rbind, rows)

if (nrow(overlap_table) == 0) {
  stop("No miRNA-target-pathway overlap rows were generated.")
}

overlap_table$fisher_FDR <- p.adjust(overlap_table$fisher_P, method = "BH")

overlap_table <- overlap_table[order(
  overlap_table$contrast,
  overlap_table$fisher_FDR,
  -overlap_table$overlap_target_pathway_n
), ]

write.csv(
  overlap_table,
  file.path(out_dir, "miRNA_target_pathway_overlap_table.csv"),
  row.names = FALSE
)

direction_table <- overlap_table[
  overlap_table$overlap_target_pathway_n >= 3,
]

direction_table <- direction_table[order(
  direction_table$contrast,
  direction_table$is_primary_candidate,
  direction_table$fisher_FDR,
  -direction_table$overlap_target_pathway_n
), ]

write.csv(
  direction_table,
  file.path(out_dir, "miRNA_target_pathway_direction_consistency.csv"),
  row.names = FALSE
)

# ----------------------------
# 7. Summaries for primary candidates
# ----------------------------

primary_rows <- direction_table[direction_table$is_primary_candidate == TRUE, ]

if (nrow(primary_rows) > 0) {
  primary_summary <- primary_rows %>%
    group_by(contrast, miRNA_term, miRNA_pretty_name) %>%
    summarise(
      n_pathways_tested_with_overlap_ge3 = n(),
      n_pathways_target_enriched_FDR005 = sum(fisher_FDR < 0.05, na.rm = TRUE),
      n_pathways_target_enriched_FDR010 = sum(fisher_FDR < 0.10, na.rm = TRUE),
      n_pathways_target_direction_consistent = sum(pathway_target_direction_consistent, na.rm = TRUE),
      n_pathways_target_direction_opposes = sum(!pathway_target_direction_consistent, na.rm = TRUE),
      main_target_implied_miRNA_direction = names(sort(table(target_implied_miRNA_direction), decreasing = TRUE))[1],
      best_pathway = pathway_Description[which.min(fisher_FDR)],
      best_pathway_FDR = min(fisher_FDR, na.rm = TRUE),
      best_pathway_overlap_n = overlap_target_pathway_n[which.min(fisher_FDR)],
      model_directions = paste(unique(na.omit(model_directions_for_primary_candidate)), collapse = " | "),
      .groups = "drop"
    )
} else {
  primary_summary <- data.frame()
}

write.csv(
  primary_summary,
  file.path(out_dir, "primary_candidate_indirect_effect_summary.csv"),
  row.names = FALSE
)

# ----------------------------
# 8. Report-ready key findings
# ----------------------------

report_ready <- direction_table[
  direction_table$fisher_FDR < 0.10 &
    direction_table$pathway_target_direction_consistent == TRUE,
]

report_ready <- report_ready[order(
  report_ready$contrast,
  report_ready$fisher_FDR,
  -report_ready$overlap_target_pathway_n
), ]

report_ready <- report_ready[, c(
  "contrast",
  "miRNA_pretty_name",
  "miRNA_term",
  "is_primary_candidate",
  "pathway_Description",
  "pathway_direction",
  "overlap_target_pathway_n",
  "overlap_mean_logFC",
  "overlap_target_direction",
  "target_implied_miRNA_direction",
  "fisher_OR",
  "fisher_P",
  "fisher_FDR",
  "indirect_logic",
  "overlap_genes",
  "model_directions_for_primary_candidate",
  "n_models_significant_for_primary_candidate"
)]

write.csv(
  report_ready,
  file.path(out_dir, "report_ready_key_indirect_findings.csv"),
  row.names = FALSE
)

cat("\nSelected miRNA terms:\n")
print(selected_mirna_df)

cat("\nSelected pathways:\n")
print(selected_pathways[, c("threshold", "contrast", "key_pattern", "best_ID", "best_Description", "NES", "FDR", "direction")])

cat("\nPrimary candidate indirect summary:\n")
print(primary_summary)

cat("\nTop report-ready indirect findings:\n")
print(head(report_ready, 30))

cat("\nDone: miRNA-target-pathway indirect effect analysis completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
