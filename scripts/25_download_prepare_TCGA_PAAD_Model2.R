suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    BiocManager::install("TCGAbiolinks", ask = FALSE, update = FALSE)
  }
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

out_dir <- "results/25_TCGA_PAAD_Model2_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

project <- "TCGA-PAAD"

message("Querying TCGA-PAAD mRNA data...")

query_mrna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

message("Downloading TCGA-PAAD mRNA data...")
GDCdownload(query_mrna)

message("Preparing TCGA-PAAD mRNA data...")
mrna_se <- GDCprepare(query_mrna)

saveRDS(
  mrna_se,
  file.path(out_dir, "TCGA_PAAD_mRNA_STAR_Counts_SE.rds")
)

message("mRNA SummarizedExperiment assays:")
print(assayNames(mrna_se))

# Prefer TPM if available.
assay_names <- assayNames(mrna_se)

if ("tpm_unstranded" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "tpm_unstranded")
} else if ("fpkm_unstrand" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "fpkm_unstrand")
} else if ("unstranded" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "unstranded")
} else {
  mrna_mat <- assay(mrna_se, 1)
}

mrna_meta <- as.data.frame(rowData(mrna_se))

# Get gene symbols.
symbol_col <- intersect(
  c("gene_name", "external_gene_name", "gene_symbol", "symbol"),
  colnames(mrna_meta)
)[1]

if (is.na(symbol_col)) {
  stop("Could not find gene symbol column in mRNA rowData.")
}

gene_symbol <- as.character(mrna_meta[[symbol_col]])

keep_gene <- !is.na(gene_symbol) &
  gene_symbol != "" &
  gene_symbol != "NA"

mrna_mat <- mrna_mat[keep_gene, , drop = FALSE]
gene_symbol <- gene_symbol[keep_gene]

rownames(mrna_mat) <- gene_symbol
storage.mode(mrna_mat) <- "numeric"

# Aggregate duplicated gene symbols by mean.
mrna_sum <- rowsum(mrna_mat, group = rownames(mrna_mat), reorder = FALSE)
gene_counts <- table(rownames(mrna_mat))
mrna_mat <- mrna_sum / as.numeric(gene_counts[rownames(mrna_sum)])
mrna_mat <- as.matrix(mrna_mat)

message("Final mRNA matrix dimension:")
print(dim(mrna_mat))

message("Querying TCGA-PAAD miRNA data...")

query_mirna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = "Primary Tumor"
)

message("Downloading TCGA-PAAD miRNA data...")
GDCdownload(query_mirna)

message("Preparing TCGA-PAAD miRNA data...")
mirna_se <- GDCprepare(query_mirna)

saveRDS(
  mirna_se,
  file.path(out_dir, "TCGA_PAAD_miRNA_SE.rds")
)

# TCGAbiolinks may return miRNA data as SummarizedExperiment or data.frame.
if (inherits(mirna_se, "SummarizedExperiment")) {
  message("miRNA SummarizedExperiment assays:")
  print(assayNames(mirna_se))
  
  mir_assay_names <- assayNames(mirna_se)
  
  if ("reads_per_million_miRNA_mapped" %in% mir_assay_names) {
    mirna_mat <- assay(mirna_se, "reads_per_million_miRNA_mapped")
  } else if ("read_count" %in% mir_assay_names) {
    mirna_mat <- assay(mirna_se, "read_count")
  } else {
    mirna_mat <- assay(mirna_se, 1)
  }
  
  mirna_meta <- as.data.frame(rowData(mirna_se))
  
  id_col <- intersect(
    c("miRNA_ID", "miRNA", "mirna_id", "name"),
    colnames(mirna_meta)
  )[1]
  
  if (!is.na(id_col)) {
    rownames(mirna_mat) <- as.character(mirna_meta[[id_col]])
  }
  
} else {
  message("miRNA object class:")
  print(class(mirna_se))
  
  mirna_df <- as.data.frame(mirna_se)
  
  print(head(mirna_df))
  print(colnames(mirna_df))
  
  id_col <- intersect(
    c("miRNA_ID", "miRNA", "mirna_id", "miRNA_IDs"),
    colnames(mirna_df)
  )[1]
  
  if (is.na(id_col)) {
    id_col <- colnames(mirna_df)[1]
  }
  
  sample_cols <- grep("^TCGA-", colnames(mirna_df), value = TRUE)
  
  if (length(sample_cols) == 0) {
    stop("Could not detect TCGA miRNA sample columns.")
  }
  
  mirna_mat <- as.matrix(mirna_df[, sample_cols, drop = FALSE])
  storage.mode(mirna_mat) <- "numeric"
  rownames(mirna_mat) <- as.character(mirna_df[[id_col]])
}

storage.mode(mirna_mat) <- "numeric"

# Remove invalid miRNA names.
keep_mirna <- !is.na(rownames(mirna_mat)) &
  rownames(mirna_mat) != "" &
  grepl("hsa-", rownames(mirna_mat), ignore.case = TRUE)

mirna_mat <- mirna_mat[keep_mirna, , drop = FALSE]

# Aggregate duplicated miRNAs by mean.
mirna_sum <- rowsum(mirna_mat, group = rownames(mirna_mat), reorder = FALSE)
mirna_counts <- table(rownames(mirna_mat))
mirna_mat <- mirna_sum / as.numeric(mirna_counts[rownames(mirna_sum)])
mirna_mat <- as.matrix(mirna_mat)

message("Final miRNA matrix dimension:")
print(dim(mirna_mat))

# Pair by TCGA sample core: first 16 characters = participant + sample vial.
# Example: TCGA-XX-YYYY-01A
get_tcga_sample_core <- function(x) substr(x, 1, 16)

mrna_core <- get_tcga_sample_core(colnames(mrna_mat))
mirna_core <- get_tcga_sample_core(colnames(mirna_mat))

mrna_df <- data.frame(
  mrna_col = colnames(mrna_mat),
  core = mrna_core,
  stringsAsFactors = FALSE
)

mirna_df <- data.frame(
  mirna_col = colnames(mirna_mat),
  core = mirna_core,
  stringsAsFactors = FALSE
)

# Use one mRNA and one miRNA column per core.
mrna_df <- mrna_df[!duplicated(mrna_df$core), ]
mirna_df <- mirna_df[!duplicated(mirna_df$core), ]

pair_df <- merge(mrna_df, mirna_df, by = "core", all = FALSE)

message("Paired TCGA-PAAD mRNA-miRNA samples:")
print(nrow(pair_df))

write.csv(
  pair_df,
  file.path(out_dir, "TCGA_PAAD_Model2_paired_sample_table.csv"),
  row.names = FALSE
)

mrna_paired <- mrna_mat[, pair_df$mrna_col, drop = FALSE]
mirna_paired <- mirna_mat[, pair_df$mirna_col, drop = FALSE]

colnames(mrna_paired) <- pair_df$core
colnames(mirna_paired) <- pair_df$core

if (!identical(colnames(mrna_paired), colnames(mirna_paired))) {
  stop("Paired mRNA and miRNA columns are not aligned.")
}

# Match genes to target β-cell pseudo-bulk matrix.
target_file <- "results/21_GTEx_miRSCAPE_training_input/GSE221156_beta_pseudobulk_logCPM_common_genes_for_miRSCAPE.rds"

if (!file.exists(target_file)) {
  target_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"
  target <- read.csv(target_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  target <- as.matrix(target)
} else {
  target <- readRDS(target_file)
}

storage.mode(target) <- "numeric"

common_genes <- intersect(rownames(mrna_paired), rownames(target))

message("Common genes between TCGA-PAAD mRNA and β-cell target:")
print(length(common_genes))

mrna_common <- mrna_paired[common_genes, , drop = FALSE]
target_common <- target[common_genes, , drop = FALSE]

saveRDS(
  mrna_common,
  file.path(out_dir, "TCGA_PAAD_Model2_bulk_mRNA_common_genes.rds")
)

saveRDS(
  mirna_paired,
  file.path(out_dir, "TCGA_PAAD_Model2_bulk_miRNA.rds")
)

saveRDS(
  target_common,
  file.path(out_dir, "GSE221156_beta_pseudobulk_common_genes_for_Model2_miRSCAPE.rds")
)

write.csv(
  mrna_common,
  file.path(out_dir, "TCGA_PAAD_Model2_bulk_mRNA_common_genes.csv"),
  row.names = TRUE
)

write.csv(
  mirna_paired,
  file.path(out_dir, "TCGA_PAAD_Model2_bulk_miRNA.csv"),
  row.names = TRUE
)

summary_df <- data.frame(
  item = c(
    "TCGA_PAAD_paired_samples",
    "TCGA_mRNA_genes_common_with_target",
    "TCGA_miRNAs",
    "GSE221156_beta_samples"
  ),
  value = c(
    ncol(mrna_common),
    nrow(mrna_common),
    nrow(mirna_paired),
    ncol(target_common)
  )
)

write.csv(
  summary_df,
  file.path(out_dir, "TCGA_PAAD_Model2_input_summary.csv"),
  row.names = FALSE
)

print(summary_df)

cat("\nPrimary candidate availability:\n")
primary <- c(
  "hsa-miR-195-5p",
  "hsa-miR-16-5p",
  "hsa-miR-15a-5p",
  "hsa-miR-15b-5p",
  "hsa-miR-649",
  "hsa-miR-6838-5p"
)

print(data.frame(
  miRNA = primary,
  available_in_Model2_TCGA_PAAD = primary %in% rownames(mirna_paired)
))

cat("\nDone: TCGA-PAAD Model2 input prepared.\n")
