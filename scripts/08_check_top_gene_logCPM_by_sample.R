logcpm_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"
sample_file <- "results/07_beta_DE_edgeR/samples_used_for_DE.csv"

out_dir <- "results/09_sanity_checks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

logcpm <- read.csv(
  logcpm_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

genes_to_check <- c(
  "ELFN1", "FOXA3", "OPLAH", "TBX2-AS1", "PHLDB2",
  "FAIM2", "OPRD1", "DKK3", "A1CF", "ASCL2",
  "INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8"
)

genes_to_check <- genes_to_check[genes_to_check %in% rownames(logcpm)]

expr_list <- list()

for (g in genes_to_check) {
  df <- data.frame(
    gene = g,
    sample_prefix = colnames(logcpm),
    logCPM = as.numeric(logcpm[g, ]),
    stringsAsFactors = FALSE
  )
  
  df <- merge(
    df,
    sample_info[, c("sample_prefix", "disease_group", "donor_or_pair_id", "n_beta_cells", "beta_fraction")],
    by = "sample_prefix",
    all.x = TRUE
  )
  
  expr_list[[g]] <- df
}

expr_df <- do.call(rbind, expr_list)

write.csv(
  expr_df,
  file.path(out_dir, "top_gene_and_beta_marker_logCPM_by_sample.csv"),
  row.names = FALSE
)

group_summary <- aggregate(
  logCPM ~ gene + disease_group,
  data = expr_df,
  FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x))
)

group_summary_clean <- data.frame(
  gene = group_summary$gene,
  disease_group = group_summary$disease_group,
  mean_logCPM = group_summary$logCPM[, "mean"],
  median_logCPM = group_summary$logCPM[, "median"],
  sd_logCPM = group_summary$logCPM[, "sd"]
)

write.csv(
  group_summary_clean,
  file.path(out_dir, "top_gene_and_beta_marker_group_summary.csv"),
  row.names = FALSE
)

print(group_summary_clean)
