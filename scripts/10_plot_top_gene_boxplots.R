logcpm_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"
sample_file <- "results/07_beta_DE_edgeR/samples_used_for_DE.csv"
out_dir <- "results/09_sanity_checks/top_gene_boxplots"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

logcpm <- read.csv(
  logcpm_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_info <- read.csv(sample_file, stringsAsFactors = FALSE)

genes <- c(
  "ELFN1", "FOXA3", "OPLAH", "TBX2-AS1", "PHLDB2",
  "FAIM2", "OPRD1", "DKK3", "A1CF", "ASCL2",
  "INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8"
)

genes <- genes[genes %in% rownames(logcpm)]

for (g in genes) {
  df <- data.frame(
    sample_prefix = colnames(logcpm),
    logCPM = as.numeric(logcpm[g, ]),
    stringsAsFactors = FALSE
  )
  
  df <- merge(
    df,
    sample_info[, c("sample_prefix", "disease_group")],
    by = "sample_prefix",
    all.x = TRUE
  )
  
  df$disease_group <- factor(df$disease_group, levels = c("ND", "PD", "T2D"))
  
  png(
    filename = file.path(out_dir, paste0(g, "_boxplot.png")),
    width = 900,
    height = 700,
    res = 150
  )
  
  boxplot(
    logCPM ~ disease_group,
    data = df,
    main = paste0(g, " logCPM by disease group"),
    xlab = "Disease group",
    ylab = "logCPM"
  )
  
  stripchart(
    logCPM ~ disease_group,
    data = df,
    vertical = TRUE,
    method = "jitter",
    pch = 16,
    add = TRUE
  )
  
  dev.off()
}

message("Done: top gene boxplots saved.")
