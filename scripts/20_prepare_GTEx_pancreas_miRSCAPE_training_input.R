suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  library(data.table)
})

out_dir <- "results/21_GTEx_miRSCAPE_training_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mrna_file <- "data/reference/GTEx/GTEx_mRNA_expression.gct.gz"
mirna_file <- "data/reference/GTEx/GTEx_miRNA_TPM_matrix.txt.gz"
pair_file <- "results/19_GTEx_miRSCAPE_reference_audit/GTEx_mRNA_miRNA_TPM_paired_pancreas_samples.csv"
target_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"

mean_by_rownames <- function(mat) {
  ids <- rownames(mat)
  sum_mat <- rowsum(mat, group = ids, reorder = FALSE)
  counts <- table(ids)
  sum_mat <- sum_mat / as.numeric(counts[rownames(sum_mat)])
  return(as.matrix(sum_mat))
}

message("Reading pancreas pair table...")
pairs <- read.csv(pair_file, stringsAsFactors = FALSE)

if (nrow(pairs) < 30) {
  stop("Too few pancreas paired samples. Pancreas-only model is not recommended.")
}

message("Pancreas paired samples: ", nrow(pairs))

# ---------- Read GTEx mRNA TPM ----------
message("Reading GTEx mRNA header...")
con <- gzfile(mrna_file, "rt")
readLines(con, n = 2)
mrna_header <- readLines(con, n = 1)
close(con)

mrna_fields <- strsplit(mrna_header, "\t")[[1]]

needed_mrna_cols <- c("Name", "Description", pairs$mrna_sample_id)
missing_mrna <- setdiff(needed_mrna_cols, mrna_fields)

if (length(missing_mrna) > 0) {
  stop("Missing mRNA columns: ", paste(head(missing_mrna, 10), collapse = ", "))
}

message("Reading selected GTEx mRNA columns...")
mrna_dt <- fread(
  mrna_file,
  skip = 2,
  select = needed_mrna_cols,
  data.table = FALSE
)

message("mRNA table dimension:")
print(dim(mrna_dt))

gene_symbol <- as.character(mrna_dt$Description)

keep_gene <- !is.na(gene_symbol) &
  gene_symbol != "" &
  gene_symbol != "NA" &
  gene_symbol != "Description"

mrna_dt <- mrna_dt[keep_gene, ]
gene_symbol <- gene_symbol[keep_gene]

mrna_mat <- as.matrix(mrna_dt[, pairs$mrna_sample_id, drop = FALSE])
storage.mode(mrna_mat) <- "numeric"
rownames(mrna_mat) <- gene_symbol
colnames(mrna_mat) <- pairs$core_id

message("Aggregating duplicated mRNA gene symbols...")
mrna_mat <- mean_by_rownames(mrna_mat)

message("Final GTEx pancreas mRNA matrix dimension:")
print(dim(mrna_mat))

# ---------- Read GTEx miRNA TPM ----------
message("Reading GTEx miRNA header...")
con <- gzfile(mirna_file, "rt")
mirna_header <- readLines(con, n = 1)
close(con)

mirna_fields <- strsplit(mirna_header, "\t")[[1]]

message("First 10 miRNA header fields:")
print(head(mirna_fields, 10))

mirna_sample_idx <- match(pairs$mirna_sample_id, mirna_fields)

if (any(is.na(mirna_sample_idx))) {
  missing_ids <- pairs$mirna_sample_id[is.na(mirna_sample_idx)]
  stop("Missing miRNA sample columns: ", paste(head(missing_ids, 10), collapse = ", "))
}

# The first column contains miRNA IDs.
# Its header may be blank, so we select it by position.
mirna_id_idx <- 1
needed_mirna_idx <- c(mirna_id_idx, mirna_sample_idx)

message("Reading selected GTEx miRNA columns by position...")
mirna_dt <- fread(
  mirna_file,
  select = needed_mirna_idx,
  data.table = FALSE
)

colnames(mirna_dt) <- c("miRNA_id", pairs$mirna_sample_id)

message("miRNA table dimension:")
print(dim(mirna_dt))

mirna_id <- as.character(mirna_dt$miRNA_id)

keep_mirna <- !is.na(mirna_id) &
  mirna_id != "" &
  mirna_id != "miRNA_id"

mirna_dt <- mirna_dt[keep_mirna, ]
mirna_id <- mirna_id[keep_mirna]

mirna_mat <- as.matrix(mirna_dt[, pairs$mirna_sample_id, drop = FALSE])
storage.mode(mirna_mat) <- "numeric"
rownames(mirna_mat) <- mirna_id
colnames(mirna_mat) <- pairs$core_id

message("Aggregating duplicated miRNA IDs...")
mirna_mat <- mean_by_rownames(mirna_mat)

message("Final GTEx pancreas miRNA matrix dimension:")
print(dim(mirna_mat))

# ---------- Match columns ----------
common_core <- intersect(colnames(mrna_mat), colnames(mirna_mat))
common_core <- intersect(common_core, pairs$core_id)

mrna_mat <- mrna_mat[, common_core, drop = FALSE]
mirna_mat <- mirna_mat[, common_core, drop = FALSE]

if (!identical(colnames(mrna_mat), colnames(mirna_mat))) {
  stop("mRNA and miRNA columns are not aligned.")
}

message("Matched GTEx pancreas training samples:")
print(length(common_core))

# ---------- Read target beta-cell pseudo-bulk mRNA ----------
message("Reading GSE221156 beta-cell pseudo-bulk target matrix...")
target <- read.csv(
  target_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

target <- as.matrix(target)
storage.mode(target) <- "numeric"

message("Target beta-cell pseudo-bulk matrix dimension before gene matching:")
print(dim(target))

# ---------- Gene overlap ----------
common_genes <- intersect(rownames(mrna_mat), rownames(target))

message("Common genes between GTEx mRNA and beta pseudo-bulk target:")
print(length(common_genes))

mrna_mat_common <- mrna_mat[common_genes, , drop = FALSE]
target_common <- target[common_genes, , drop = FALSE]

# ---------- Save outputs ----------
saveRDS(
  mrna_mat,
  file.path(out_dir, "GTEx_pancreas_bulk_mRNA_TPM_all_genes.rds")
)

saveRDS(
  mirna_mat,
  file.path(out_dir, "GTEx_pancreas_bulk_miRNA_TPM.rds")
)

saveRDS(
  mrna_mat_common,
  file.path(out_dir, "GTEx_pancreas_bulk_mRNA_TPM_common_genes.rds")
)

saveRDS(
  target_common,
  file.path(out_dir, "GSE221156_beta_pseudobulk_logCPM_common_genes_for_miRSCAPE.rds")
)

write.csv(
  mrna_mat_common,
  file.path(out_dir, "GTEx_pancreas_bulk_mRNA_TPM_common_genes.csv"),
  row.names = TRUE
)

write.csv(
  mirna_mat,
  file.path(out_dir, "GTEx_pancreas_bulk_miRNA_TPM.csv"),
  row.names = TRUE
)

write.csv(
  target_common,
  file.path(out_dir, "GSE221156_beta_pseudobulk_logCPM_common_genes_for_miRSCAPE.csv"),
  row.names = TRUE
)

write.csv(
  data.frame(
    core_id = common_core,
    mrna_sample_id = pairs$mrna_sample_id[match(common_core, pairs$core_id)],
    mirna_sample_id = pairs$mirna_sample_id[match(common_core, pairs$core_id)],
    tissue = "Pancreas"
  ),
  file.path(out_dir, "GTEx_pancreas_training_sample_pairs_used.csv"),
  row.names = FALSE
)

summary_df <- data.frame(
  item = c(
    "GTEx_pancreas_paired_samples",
    "GTEx_mRNA_genes_all",
    "GTEx_miRNAs",
    "GSE221156_beta_samples",
    "common_mRNA_genes_used"
  ),
  value = c(
    length(common_core),
    nrow(mrna_mat),
    nrow(mirna_mat),
    ncol(target_common),
    length(common_genes)
  )
)

write.csv(
  summary_df,
  file.path(out_dir, "GTEx_pancreas_miRSCAPE_training_input_summary.csv"),
  row.names = FALSE
)

print(summary_df)

cat("\nExample miRNAs:\n")
print(head(rownames(mirna_mat), 20))

cat("\nExample common genes:\n")
print(head(common_genes, 20))

cat("\nDone: GTEx pancreas miRSCAPE training input prepared.\n")
