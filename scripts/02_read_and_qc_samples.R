suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

manifest_file <- "results/01_manifest/sample_manifest.csv"
raw_dir <- "data/raw"
out_obj_dir <- "objects/02_seurat_samples"
out_qc_dir <- "results/02_qc"

dir.create(out_obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_qc_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)

if (!all(manifest$complete_10x_triplet)) {
  stop("Some samples are incomplete. Please check sample_manifest.csv first.")
}

qc_list <- list()

for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample_prefix[i]
  
  message("====================================")
  message("Reading sample ", i, " / ", nrow(manifest))
  message(sample_id)
  message("====================================")
  
  mtx_file <- file.path(raw_dir, manifest$matrix[i])
  barcode_file <- file.path(raw_dir, manifest$barcodes[i])
  feature_file <- file.path(raw_dir, manifest$features[i])
  
  counts <- ReadMtx(
    mtx = mtx_file,
    cells = barcode_file,
    features = feature_file,
    feature.column = 2,
    unique.features = TRUE
  )
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  seu$sample_id <- sample_id
  
  # mitochondrial percentage
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # ribosomal percentage, optional but useful
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  # raw QC numbers before filtering
  n_cells_before <- ncol(seu)
  median_nFeature_before <- median(seu$nFeature_RNA)
  median_nCount_before <- median(seu$nCount_RNA)
  median_percent_mt_before <- median(seu$percent.mt)
  
  # Conservative first-pass QC.
  # These thresholds are intentionally broad.
  # We only remove clear low-quality cells here.
  seu_qc <- subset(
    seu,
    subset =
      nFeature_RNA >= 200 &
      nFeature_RNA <= 7500 &
      nCount_RNA >= 500 &
      percent.mt <= 25
  )
  
  n_cells_after <- ncol(seu_qc)
  median_nFeature_after <- median(seu_qc$nFeature_RNA)
  median_nCount_after <- median(seu_qc$nCount_RNA)
  median_percent_mt_after <- median(seu_qc$percent.mt)
  
  qc_one <- data.frame(
    sample_id = sample_id,
    n_cells_before_QC = n_cells_before,
    n_cells_after_QC = n_cells_after,
    cells_removed = n_cells_before - n_cells_after,
    percent_cells_removed = round(100 * (n_cells_before - n_cells_after) / n_cells_before, 2),
    median_nFeature_before = median_nFeature_before,
    median_nFeature_after = median_nFeature_after,
    median_nCount_before = median_nCount_before,
    median_nCount_after = median_nCount_after,
    median_percent_mt_before = median_percent_mt_before,
    median_percent_mt_after = median_percent_mt_after,
    stringsAsFactors = FALSE
  )
  
  qc_list[[sample_id]] <- qc_one
  
  saveRDS(
    seu_qc,
    file = file.path(out_obj_dir, paste0(sample_id, ".seurat_qc.rds"))
  )
  
  rm(counts, seu, seu_qc)
  gc()
}

qc_summary <- do.call(rbind, qc_list)

write.csv(
  qc_summary,
  file.path(out_qc_dir, "sample_qc_summary.csv"),
  row.names = FALSE
)

overall_summary <- data.frame(
  n_samples = nrow(qc_summary),
  total_cells_before_QC = sum(qc_summary$n_cells_before_QC),
  total_cells_after_QC = sum(qc_summary$n_cells_after_QC),
  total_cells_removed = sum(qc_summary$cells_removed),
  percent_cells_removed = round(
    100 * sum(qc_summary$cells_removed) / sum(qc_summary$n_cells_before_QC),
    2
  ),
  median_cells_per_sample_after_QC = median(qc_summary$n_cells_after_QC),
  min_cells_per_sample_after_QC = min(qc_summary$n_cells_after_QC),
  max_cells_per_sample_after_QC = max(qc_summary$n_cells_after_QC)
)

write.csv(
  overall_summary,
  file.path(out_qc_dir, "overall_qc_summary.csv"),
  row.names = FALSE
)

print(overall_summary)
message("Done: sample-level Seurat objects and QC summaries were saved.")
