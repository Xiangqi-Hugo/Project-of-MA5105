suppressPackageStartupMessages({
  library(decoupleR)
  library(limma)
  library(ggplot2)
})

expr_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"
sample_file <- "results/07_beta_DE_edgeR/samples_used_for_DE.csv"
term2gene_file <- "results/14_strict_miRNA_target_GSEA/strict_miRNA_TERM2GENE.csv"

out_dir <- "results/15_decoupleR_miRNA_activity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading expression matrix...")
expr <- read.csv(
  expr_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

message("Expression matrix dimension:")
print(dim(expr))

message("Reading sample metadata...")
sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

sample_info <- sample_info[match(colnames(expr), sample_info$sample_prefix), ]

if (any(is.na(sample_info$sample_prefix))) {
  stop("Sample metadata does not match expression matrix columns.")
}

sample_info$disease_group <- factor(
  sample_info$disease_group,
  levels = c("ND", "PD", "T2D")
)

message("Group counts:")
print(table(sample_info$disease_group))

message("Reading strict miRNA TERM2GENE...")
term2gene <- read.csv(term2gene_file, stringsAsFactors = FALSE)

colnames(term2gene) <- c("source", "target")

# Keep only genes present in expression matrix.
term2gene <- term2gene[term2gene$target %in% rownames(expr), ]

# Remove duplicated edges.
term2gene <- unique(term2gene)

# Require at least 10 detected targets per miRNA.
target_counts <- table(term2gene$source)
keep_sources <- names(target_counts[target_counts >= 10])
term2gene <- term2gene[term2gene$source %in% keep_sources, ]

# miRNAs usually repress target mRNAs.
# We set mor = -1 so lower target expression gives higher inferred miRNA activity.
net <- data.frame(
  source = term2gene$source,
  target = term2gene$target,
  mor = -1,
  stringsAsFactors = FALSE
)

message("Network summary:")
print(length(unique(net$source)))
print(nrow(net))

# Gene-wise scaling across samples.
# This makes each gene comparable before activity inference.
expr_z <- t(scale(t(expr)))
expr_z[is.na(expr_z)] <- 0

message("Running decoupleR weighted mean activity inference...")

act_long <- decoupleR::run_wmean(
  mat = expr_z,
  network = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 10,
  times = 100
)

write.csv(
  act_long,
  file.path(out_dir, "decoupleR_wmean_activity_long.csv"),
  row.names = FALSE
)

# Convert long table to miRNA × sample matrix.
activity_mat <- tapply(
  act_long$score,
  list(act_long$source, act_long$condition),
  mean
)

activity_mat <- as.matrix(activity_mat)

# Match sample order.
activity_mat <- activity_mat[, sample_info$sample_prefix, drop = FALSE]

write.csv(
  activity_mat,
  file.path(out_dir, "decoupleR_wmean_activity_matrix.csv"),
  row.names = TRUE
)

message("Activity matrix dimension:")
print(dim(activity_mat))

# Differential activity analysis using limma.
design <- model.matrix(~ 0 + disease_group, data = sample_info)
colnames(design) <- levels(sample_info$disease_group)

fit <- lmFit(activity_mat, design)

contrasts <- makeContrasts(
  T2D_vs_ND = T2D - ND,
  PD_vs_ND = PD - ND,
  T2D_vs_PD = T2D - PD,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

contrast_names <- colnames(contrasts)

summary_list <- list()

for (contrast in contrast_names) {
  message("Testing activity contrast: ", contrast)
  
  tab <- topTable(
    fit2,
    coef = contrast,
    number = Inf,
    sort.by = "P"
  )
  
  tab$miRNA <- rownames(tab)
  tab <- tab[, c("miRNA", setdiff(colnames(tab), "miRNA"))]
  
  first_group <- sub("_vs_.*$", "", contrast)
  
  tab$activity_direction <- ifelse(
    tab$logFC > 0,
    paste0("higher_in_", first_group),
    paste0("lower_in_", first_group)
  )
  
  write.csv(
    tab,
    file.path(out_dir, paste0(contrast, "_decoupleR_activity_all.csv")),
    row.names = FALSE
  )
  
  sig005 <- tab[tab$adj.P.Val < 0.05, ]
  sig010 <- tab[tab$adj.P.Val < 0.10, ]
  
  write.csv(
    sig005,
    file.path(out_dir, paste0(contrast, "_decoupleR_activity_FDR005.csv")),
    row.names = FALSE
  )
  
  write.csv(
    sig010,
    file.path(out_dir, paste0(contrast, "_decoupleR_activity_FDR010.csv")),
    row.names = FALSE
  )
  
  summary_one <- data.frame(
    contrast = contrast,
    n_miRNAs_tested = nrow(tab),
    n_FDR005 = nrow(sig005),
    n_higher_first_group_FDR005 = sum(sig005$logFC > 0),
    n_lower_first_group_FDR005 = sum(sig005$logFC < 0),
    n_FDR010 = nrow(sig010),
    top_miRNA = tab$miRNA[1],
    top_logFC = tab$logFC[1],
    top_FDR = tab$adj.P.Val[1],
    stringsAsFactors = FALSE
  )
  
  summary_list[[contrast]] <- summary_one
}

summary_df <- do.call(rbind, summary_list)

write.csv(
  summary_df,
  file.path(out_dir, "decoupleR_activity_summary.csv"),
  row.names = FALSE
)

print(summary_df)

# Primary candidate summary.
primary_candidates <- c(
  "MIR195_5P",
  "MIR16_5P",
  "MIR15B_5P",
  "MIR15A_5P",
  "MIR649",
  "MIR6838_5P"
)

candidate_list <- list()

for (contrast in contrast_names) {
  tab <- read.csv(
    file.path(out_dir, paste0(contrast, "_decoupleR_activity_all.csv")),
    stringsAsFactors = FALSE
  )
  
  subtab <- tab[tab$miRNA %in% primary_candidates, ]
  subtab$contrast <- contrast
  
  candidate_list[[contrast]] <- subtab
}

candidate_df <- do.call(rbind, candidate_list)

write.csv(
  candidate_df,
  file.path(out_dir, "primary_T2D_vs_ND_candidates_decoupleR_summary.csv"),
  row.names = FALSE
)

print(candidate_df)

# Long activity table for plots.
activity_long <- data.frame(
  miRNA = rep(rownames(activity_mat), times = ncol(activity_mat)),
  sample_prefix = rep(colnames(activity_mat), each = nrow(activity_mat)),
  activity = as.numeric(activity_mat),
  stringsAsFactors = FALSE
)

activity_long <- merge(
  activity_long,
  sample_info[, c("sample_prefix", "disease_group")],
  by = "sample_prefix",
  all.x = TRUE
)

write.csv(
  activity_long,
  file.path(out_dir, "decoupleR_activity_long_with_metadata.csv"),
  row.names = FALSE
)

# Plot primary candidates.
plot_df <- activity_long[activity_long$miRNA %in% primary_candidates, ]

plot_df$disease_group <- factor(
  plot_df$disease_group,
  levels = c("ND", "PD", "T2D")
)

p <- ggplot(
  plot_df,
  aes(x = disease_group, y = activity)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.8) +
  facet_wrap(~ miRNA, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 12) +
  labs(
    title = "decoupleR inferred miRNA activity for primary candidates",
    x = "Disease group",
    y = "decoupleR weighted mean activity"
  )

ggsave(
  filename = file.path(out_dir, "primary_candidate_miRNA_activity_boxplots.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

message("Done: decoupleR miRNA activity analysis completed.")
