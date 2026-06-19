# Step 34b: relaxed exploratory miRNA indirect effect tables
#
# Purpose:
# Create report-friendly tables from Step 34.
#
# This script separates:
#   1. high-confidence indirect effects
#   2. relaxed exploratory directional patterns
#   3. top miRNA-pathway overlaps by enrichment FDR
#   4. primary-candidate report summary
#
# Important:
# The relaxed table is exploratory.
# It must not be interpreted as causal evidence.
#
# Inputs:
#   results/34_miRNA_indirect_pathway_effect/miRNA_target_pathway_overlap_table.csv
#   results/34_miRNA_indirect_pathway_effect/primary_candidate_indirect_effect_summary.csv
#
# Outputs:
#   results/34b_miRNA_indirect_exploratory_tables/high_confidence_indirect_effects.csv
#   results/34b_miRNA_indirect_exploratory_tables/exploratory_directional_indirect_patterns.csv
#   results/34b_miRNA_indirect_exploratory_tables/top_miRNA_pathway_overlap_by_FDR.csv
#   results/34b_miRNA_indirect_exploratory_tables/primary_candidate_report_summary.csv
#   results/34b_miRNA_indirect_exploratory_tables/miRNA_indirect_effect_report_sentence.txt

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr", repos = "https://cloud.r-project.org")
  }
  library(dplyr)
})

in_dir <- "results/34_miRNA_indirect_pathway_effect"
out_dir <- "results/34b_miRNA_indirect_exploratory_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

overlap_file <- file.path(in_dir, "miRNA_target_pathway_overlap_table.csv")
primary_file <- file.path(in_dir, "primary_candidate_indirect_effect_summary.csv")

if (!file.exists(overlap_file)) {
  stop("Missing file: ", overlap_file)
}

overlap <- read.csv(overlap_file, stringsAsFactors = FALSE, check.names = FALSE)

required <- c(
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
  "pathway_target_direction_consistent",
  "fisher_OR",
  "fisher_P",
  "fisher_FDR",
  "indirect_logic",
  "overlap_genes"
)

missing <- setdiff(required, colnames(overlap))
if (length(missing) > 0) {
  stop("Missing columns in overlap table: ", paste(missing, collapse = ", "))
}

# Make logical robust.
overlap$pathway_target_direction_consistent <- as.logical(overlap$pathway_target_direction_consistent)
overlap$is_primary_candidate <- as.logical(overlap$is_primary_candidate)

# ----------------------------
# 1. High-confidence table
# ----------------------------
# Strict table. This is expected to be empty in current results.
high_conf <- overlap %>%
  filter(
    !is.na(fisher_FDR),
    fisher_FDR < 0.10,
    pathway_target_direction_consistent == TRUE,
    overlap_target_pathway_n >= 3
  ) %>%
  arrange(contrast, fisher_FDR, desc(overlap_target_pathway_n)) %>%
  select(
    contrast,
    miRNA_pretty_name,
    miRNA_term,
    is_primary_candidate,
    pathway_Description,
    pathway_direction,
    overlap_target_pathway_n,
    overlap_mean_logFC,
    overlap_target_direction,
    target_implied_miRNA_direction,
    fisher_OR,
    fisher_P,
    fisher_FDR,
    indirect_logic,
    overlap_genes
  )

write.csv(
  high_conf,
  file.path(out_dir, "high_confidence_indirect_effects.csv"),
  row.names = FALSE
)

# ----------------------------
# 2. Relaxed exploratory table
# ----------------------------
# Relaxed conditions:
#   - overlap >= 3 target genes
#   - target direction agrees with pathway direction
#   - Fisher OR > 1 means target genes are more common in the pathway than expected
#   - Fisher P < 0.25 is permissive and exploratory
#
# This is NOT a significance threshold.
# It is only used to show possible directional patterns.

exploratory <- overlap %>%
  filter(
    overlap_target_pathway_n >= 3,
    pathway_target_direction_consistent == TRUE,
    !is.na(fisher_OR),
    fisher_OR > 1,
    !is.na(fisher_P),
    fisher_P < 0.25
  ) %>%
  mutate(
    evidence_level = case_when(
      fisher_FDR < 0.10 ~ "FDR_supported",
      fisher_P < 0.05 ~ "nominal_overlap_only",
      TRUE ~ "weak_exploratory_only"
    ),
    interpretation = case_when(
      grepl("compatible_with_higher_miRNA", indirect_logic) ~
        "Targets and pathway are lower; this is compatible with higher miRNA-mediated repression.",
      grepl("compatible_with_lower_miRNA", indirect_logic) ~
        "Targets and pathway are higher; this is compatible with lower miRNA-mediated de-repression.",
      TRUE ~ "Directional pattern only."
    )
  ) %>%
  arrange(
    contrast,
    is_primary_candidate,
    evidence_level,
    fisher_P,
    desc(overlap_target_pathway_n)
  ) %>%
  select(
    contrast,
    miRNA_pretty_name,
    miRNA_term,
    is_primary_candidate,
    pathway_Description,
    pathway_direction,
    overlap_target_pathway_n,
    overlap_mean_logFC,
    overlap_target_direction,
    target_implied_miRNA_direction,
    fisher_OR,
    fisher_P,
    fisher_FDR,
    evidence_level,
    interpretation,
    overlap_genes
  )

write.csv(
  exploratory,
  file.path(out_dir, "exploratory_directional_indirect_patterns.csv"),
  row.names = FALSE
)

# ----------------------------
# 3. Top overlaps by FDR
# ----------------------------
# This table is useful when high-confidence results are empty.
# It shows the best overlaps, even if not significant.

top_by_fdr <- overlap %>%
  filter(
    overlap_target_pathway_n >= 3,
    !is.na(fisher_FDR)
  ) %>%
  arrange(contrast, fisher_FDR, fisher_P, desc(overlap_target_pathway_n)) %>%
  group_by(contrast) %>%
  slice_head(n = 30) %>%
  ungroup() %>%
  mutate(
    direction_status = ifelse(
      pathway_target_direction_consistent,
      "direction_consistent",
      "direction_opposed"
    ),
    report_use = ifelse(
      fisher_FDR < 0.10 & pathway_target_direction_consistent,
      "high_confidence",
      ifelse(
        fisher_P < 0.25 & pathway_target_direction_consistent,
        "exploratory_only",
        "not_report_as_effect"
      )
    )
  ) %>%
  select(
    contrast,
    miRNA_pretty_name,
    miRNA_term,
    is_primary_candidate,
    pathway_Description,
    pathway_direction,
    overlap_target_pathway_n,
    overlap_target_direction,
    target_implied_miRNA_direction,
    fisher_OR,
    fisher_P,
    fisher_FDR,
    direction_status,
    report_use,
    overlap_genes
  )

write.csv(
  top_by_fdr,
  file.path(out_dir, "top_miRNA_pathway_overlap_by_FDR.csv"),
  row.names = FALSE
)

# ----------------------------
# 4. Primary-candidate report summary
# ----------------------------

if (file.exists(primary_file)) {
  primary <- read.csv(primary_file, stringsAsFactors = FALSE, check.names = FALSE)

  primary_report <- primary %>%
    mutate(
      final_indirect_effect_call = case_when(
        n_pathways_target_enriched_FDR010 > 0 ~
          "possible_FDR_supported_indirect_effect",
        n_pathways_target_direction_consistent > 0 & n_pathways_target_enriched_FDR010 == 0 ~
          "directional_but_not_enriched",
        TRUE ~
          "no_indirect_support"
      ),
      report_interpretation = case_when(
        final_indirect_effect_call == "possible_FDR_supported_indirect_effect" ~
          "This candidate had at least one pathway with target enrichment at FDR < 0.10.",
        final_indirect_effect_call == "directional_but_not_enriched" ~
          "This candidate showed some target-direction consistency, but target-pathway enrichment was not significant.",
        TRUE ~
          "This candidate did not support a miRNA-mediated indirect pathway effect."
      )
    ) %>%
    arrange(contrast, final_indirect_effect_call, miRNA_term)

  write.csv(
    primary_report,
    file.path(out_dir, "primary_candidate_report_summary.csv"),
    row.names = FALSE
  )
} else {
  primary_report <- data.frame()
  write.csv(
    primary_report,
    file.path(out_dir, "primary_candidate_report_summary.csv"),
    row.names = FALSE
  )
}

# ----------------------------
# 5. Report sentence
# ----------------------------

n_high <- nrow(high_conf)
n_exploratory <- nrow(exploratory)
n_top <- nrow(top_by_fdr)

report_sentence <- c(
  "miRNA-target-pathway indirect effect analysis was used as an exploratory mechanism check.",
  paste0("High-confidence indirect effects were defined as target-pathway enrichment FDR < 0.10 with direction consistency. The number of high-confidence effects was ", n_high, "."),
  paste0("A relaxed exploratory table identified ", n_exploratory, " directional patterns using overlap >= 3 target genes, Fisher OR > 1, Fisher P < 0.25, and direction consistency."),
  "These relaxed patterns should not be interpreted as validated miRNA-mediated regulation.",
  "The final interpretation should remain conservative: no statistically robust miRNA-mediated indirect pathway effect was identified."
)

writeLines(
  report_sentence,
  con = file.path(out_dir, "miRNA_indirect_effect_report_sentence.txt")
)

cat("\nHigh-confidence effects:", n_high, "\n")
cat("Relaxed exploratory directional patterns:", n_exploratory, "\n")
cat("Top-overlap rows:", n_top, "\n")

cat("\nTop relaxed exploratory patterns:\n")
print(head(exploratory, 30))

cat("\nReport sentence:\n")
cat(paste(report_sentence, collapse = "\n"), "\n")

cat("\nDone: relaxed exploratory miRNA indirect tables completed.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
