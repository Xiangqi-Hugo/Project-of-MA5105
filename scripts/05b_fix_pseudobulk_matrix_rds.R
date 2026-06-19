suppressPackageStartupMessages({
  library(Matrix)
})

csv_file <- "results/06_beta_pseudobulk/beta_pseudobulk_counts_matrix.csv"
fixed_rds <- "results/06_beta_pseudobulk/beta_pseudobulk_counts_matrix_FIXED.rds"

message("Reading pseudo-bulk CSV...")

pb_df <- read.csv(
  csv_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

message("Data frame dimension:")
print(dim(pb_df))

pb_mat <- as.matrix(pb_df)
storage.mode(pb_mat) <- "numeric"

message("Matrix dimension:")
print(dim(pb_mat))

pb_sparse <- Matrix(
  pb_mat,
  sparse = TRUE
)

message("Sparse matrix dimension:")
print(dim(pb_sparse))
print(class(pb_sparse))
print(isS4(pb_sparse))

saveRDS(pb_sparse, fixed_rds)

message("Saved fixed RDS:")
message(fixed_rds)
