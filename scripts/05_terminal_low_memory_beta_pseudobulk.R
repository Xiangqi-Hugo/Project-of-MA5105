suppressPackageStartupMessages({
  library(Matrix)
})

metadata_file <- "results/03_metadata/sample_metadata_54.csv"
raw_dir <- "data/raw"
out_dir <- "results/06_beta_pseudobulk"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_cols <- c(
  "sample_prefix",
  "barcodes",
  "features",
  "matrix",
  "disease_group",
  "donor_or_pair_id"
)

missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop(paste("Missing columns in metadata:", paste(missing_cols, collapse = ", ")))
}

marker_list <- list(
  beta = c("INS", "IAPP", "PCSK1", "PCSK2", "MAFA", "PDX1", "NKX6-1", "SLC30A8"),
  alpha = c("GCG", "ARX", "IRX2"),
  delta = c("SST", "HHEX"),
  pp = c("PPY"),
  epsilon = c("GHRL"),
  acinar = c("PRSS1", "PRSS2", "CPA1", "CTRB2", "CELA3A"),
  ductal = c("KRT19", "SOX9", "CFTR", "KRT8"),
  endothelial = c("PECAM1", "VWF", "KDR"),
  immune = c("PTPRC", "LST1", "CD74")
)

read_10x_sparse <- function(mtx_file, barcode_file, feature_file) {
  message("Reading matrix: ", mtx_file)
  
  counts <- Matrix::readMM(gzfile(mtx_file))
  counts <- as(counts, "dgCMatrix")
  
  barcodes <- read.delim(
    gzfile(barcode_file),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  features <- read.delim(
    gzfile(feature_file),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(features) >= 2) {
    gene_names <- features[[2]]
  } else {
    gene_names <- features[[1]]
  }
  
  gene_names <- make.unique(gene_names)
  
  rownames(counts) <- gene_names
  colnames(counts) <- barcodes[[1]]
  
  counts
}

score_celltype <- function(counts, marker_list) {
  lib_size <- Matrix::colSums(counts)
  lib_size[lib_size == 0] <- 1
  
  score_mat <- matrix(
    0,
    nrow = ncol(counts),
    ncol = length(marker_list)
  )
  
  rownames(score_mat) <- colnames(counts)
  colnames(score_mat) <- names(marker_list)
  
  for (ct in names(marker_list)) {
    markers <- marker_list[[ct]]
    markers <- markers[markers %in% rownames(counts)]
    
    if (length(markers) == 0) {
      next
    }
    
    marker_counts <- counts[markers, , drop = FALSE]
    norm_marker_counts <- t(t(marker_counts) / lib_size) * 10000
    score_mat[, ct] <- Matrix::colMeans(log1p(norm_marker_counts))
  }
  
  score_mat
}

pseudobulk_list <- list()
summary_list <- list()

for (i in seq_len(nrow(metadata))) {
  sample_id <- metadata$sample_prefix[i]
  
  message("====================================")
  message("Processing sample ", i, " / ", nrow(metadata), ": ", sample_id)
  message("Disease group: ", metadata$disease_group[i])
  message("====================================")
  
  mtx_file <- file.path(raw_dir, metadata$matrix[i])
  barcode_file <- file.path(raw_dir, metadata$barcodes[i])
  feature_file <- file.path(raw_dir, metadata$features[i])
  
  if (!file.exists(mtx_file)) stop(paste("Missing matrix file:", mtx_file))
  if (!file.exists(barcode_file)) stop(paste("Missing barcode file:", barcode_file))
  if (!file.exists(feature_file)) stop(paste("Missing feature file:", feature_file))
  
  counts <- read_10x_sparse(
    mtx_file = mtx_file,
    barcode_file = barcode_file,
    feature_file = feature_file
  )
  
  nFeature_RNA <- Matrix::colSums(counts > 0)
  nCount_RNA <- Matrix::colSums(counts)
  
  mt_genes <- grepl("^MT-", rownames(counts))
  
  if (sum(mt_genes) > 0) {
    percent_mt <- Matrix::colSums(counts[mt_genes, , drop = FALSE]) / nCount_RNA * 100
  } else {
    percent_mt <- rep(0, ncol(counts))
  }
  
  percent_mt[is.na(percent_mt)] <- 0
  
  keep_cells <- nFeature_RNA >= 200 &
    nFeature_RNA <= 7500 &
    nCount_RNA >= 500 &
    percent_mt <= 25
  
  counts_qc <- counts[, keep_cells, drop = FALSE]
  
  if (ncol(counts_qc) == 0) {
    warning("No cells passed QC in sample: ", sample_id)
    
    pb <- rep(0, nrow(counts))
    names(pb) <- rownames(counts)
    
    summary_one <- data.frame(
      sample_prefix = sample_id,
      disease_group = metadata$disease_group[i],
      donor_or_pair_id = metadata$donor_or_pair_id[i],
      use_for_conservative_DE = metadata$disease_group[i] != "MIXED",
      n_cells_raw = ncol(counts),
      n_cells_after_QC = 0,
      n_beta_cells = 0,
      beta_fraction = 0,
      total_beta_UMI = 0,
      stringsAsFactors = FALSE
    )
    
    pseudobulk_list[[sample_id]] <- pb
    summary_list[[sample_id]] <- summary_one
    
    rm(counts, counts_qc)
    gc()
    next
  }
  
  score_mat <- score_celltype(counts_qc, marker_list)
  best_type <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
  
  INS_expr <- if ("INS" %in% rownames(counts_qc)) {
    as.numeric(counts_qc["INS", ])
  } else {
    rep(0, ncol(counts_qc))
  }
  
  IAPP_expr <- if ("IAPP" %in% rownames(counts_qc)) {
    as.numeric(counts_qc["IAPP", ])
  } else {
    rep(0, ncol(counts_qc))
  }
  
  beta_cells <- best_type == "beta" & (INS_expr > 0 | IAPP_expr > 0)
  
  if (sum(beta_cells) > 0) {
    beta_counts <- counts_qc[, beta_cells, drop = FALSE]
    pb <- Matrix::rowSums(beta_counts)
  } else {
    pb <- rep(0, nrow(counts_qc))
    names(pb) <- rownames(counts_qc)
  }
  
  names(pb) <- rownames(counts_qc)
  
  pseudobulk_list[[sample_id]] <- pb
  
  summary_one <- data.frame(
    sample_prefix = sample_id,
    disease_group = metadata$disease_group[i],
    donor_or_pair_id = metadata$donor_or_pair_id[i],
    use_for_conservative_DE = metadata$disease_group[i] != "MIXED",
    n_cells_raw = ncol(counts),
    n_cells_after_QC = ncol(counts_qc),
    n_beta_cells = sum(beta_cells),
    beta_fraction = round(sum(beta_cells) / ncol(counts_qc), 4),
    total_beta_UMI = sum(pb),
    stringsAsFactors = FALSE
  )
  
  summary_list[[sample_id]] <- summary_one
  
  rm(counts, counts_qc, score_mat, beta_cells, INS_expr, IAPP_expr, pb)
  if (exists("beta_counts")) rm(beta_counts)
  gc()
}

message("Combining sample-level pseudo-bulk vectors...")

all_genes <- unique(unlist(lapply(pseudobulk_list, names)))

pb_matrix <- matrix(
  0,
  nrow = length(all_genes),
  ncol = length(pseudobulk_list)
)

rownames(pb_matrix) <- all_genes
colnames(pb_matrix) <- names(pseudobulk_list)

for (sample_id in names(pseudobulk_list)) {
  x <- pseudobulk_list[[sample_id]]
  pb_matrix[names(x), sample_id] <- x
}

pb_matrix <- as(pb_matrix, "dgCMatrix")

sample_summary <- do.call(rbind, summary_list)

write.csv(
  sample_summary,
  file.path(out_dir, "beta_pseudobulk_sample_summary.csv"),
  row.names = FALSE
)

saveRDS(
  pb_matrix,
  file.path(out_dir, "beta_pseudobulk_counts_matrix.rds")
)

write.csv(
  as.matrix(pb_matrix),
  file.path(out_dir, "beta_pseudobulk_counts_matrix.csv"),
  row.names = TRUE
)

overall_summary <- data.frame(
  n_samples = ncol(pb_matrix),
  n_genes = nrow(pb_matrix),
  total_beta_cells = sum(sample_summary$n_beta_cells),
  median_beta_cells_per_sample = median(sample_summary$n_beta_cells),
  min_beta_cells_per_sample = min(sample_summary$n_beta_cells),
  max_beta_cells_per_sample = max(sample_summary$n_beta_cells),
  conservative_DE_samples = sum(sample_summary$use_for_conservative_DE),
  stringsAsFactors = FALSE
)

write.csv(
  overall_summary,
  file.path(out_dir, "beta_pseudobulk_overall_summary.csv"),
  row.names = FALSE
)

print(overall_summary)

message("Done.")
message("Saved:")
message(file.path(out_dir, "beta_pseudobulk_counts_matrix.rds"))
message(file.path(out_dir, "beta_pseudobulk_sample_summary.csv"))
message(file.path(out_dir, "beta_pseudobulk_overall_summary.csv"))
