suppressPackageStartupMessages({
  library(Seurat)
})

metadata_file <- "results/03_metadata/sample_metadata_54.csv"
obj_dir <- "objects/02_seurat_samples"
out_obj_dir <- "objects/03_merged"
out_res_dir <- "results/04_merge"

dir.create(out_obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_res_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_cols <- c(
  "sample_prefix",
  "GSM",
  "islet_title",
  "is_mixed_islet",
  "islet_pair",
  "disease_group",
  "donor_or_pair_id",
  "n_cells_after_QC"
)

missing_cols <- setdiff(required_cols, colnames(metadata))

if (length(missing_cols) > 0) {
  stop(paste("Missing metadata columns:", paste(missing_cols, collapse = ", ")))
}

metadata$use_for_conservative_DE <- metadata$disease_group != "MIXED"

message("Metadata group summary:")
print(table(metadata$disease_group, useNA = "ifany"))

merged_obj <- NULL
merge_log <- list()

for (i in seq_len(nrow(metadata))) {
  sample_id <- metadata$sample_prefix[i]
  rds_file <- file.path(obj_dir, paste0(sample_id, ".seurat_qc.rds"))
  
  message("====================================")
  message("Loading sample ", i, " / ", nrow(metadata))
  message(sample_id)
  message("====================================")
  
  if (!file.exists(rds_file)) {
    stop(paste("Missing Seurat object:", rds_file))
  }
  
  seu <- readRDS(rds_file)
  
  n_cells <- ncol(seu)
  
  # Add sample-level metadata to every cell.
  seu$sample_prefix <- sample_id
  seu$GSM <- metadata$GSM[i]
  seu$islet_title <- metadata$islet_title[i]
  seu$is_mixed_islet <- metadata$is_mixed_islet[i]
  seu$islet_pair <- metadata$islet_pair[i]
  seu$disease_group <- metadata$disease_group[i]
  seu$donor_or_pair_id <- metadata$donor_or_pair_id[i]
  seu$use_for_conservative_DE <- metadata$use_for_conservative_DE[i]
  
  # Keep original cell barcode before merge.
  seu$original_barcode <- colnames(seu)
  
  # Make cell names globally unique.
  seu <- RenameCells(seu, add.cell.id = sample_id)
  
  merge_log[[sample_id]] <- data.frame(
    sample_prefix = sample_id,
    disease_group = metadata$disease_group[i],
    donor_or_pair_id = metadata$donor_or_pair_id[i],
    use_for_conservative_DE = metadata$use_for_conservative_DE[i],
    n_cells = n_cells,
    stringsAsFactors = FALSE
  )
  
  if (is.null(merged_obj)) {
    merged_obj <- seu
  } else {
    merged_obj <- merge(
      x = merged_obj,
      y = seu,
      project = "GSE221156_54samples"
    )
  }
  
  rm(seu)
  gc()
}

merge_summary <- do.call(rbind, merge_log)

write.csv(
  merge_summary,
  file.path(out_res_dir, "merge_sample_summary.csv"),
  row.names = FALSE
)

cell_group_summary <- as.data.frame(table(merged_obj$disease_group))
colnames(cell_group_summary) <- c("disease_group", "n_cells")

write.csv(
  cell_group_summary,
  file.path(out_res_dir, "merged_cell_group_summary.csv"),
  row.names = FALSE
)

sample_group_summary <- as.data.frame(table(merge_summary$disease_group))
colnames(sample_group_summary) <- c("disease_group", "n_samples")

write.csv(
  sample_group_summary,
  file.path(out_res_dir, "merged_sample_group_summary.csv"),
  row.names = FALSE
)

overall_summary <- data.frame(
  n_samples = length(unique(merged_obj$sample_prefix)),
  n_cells = ncol(merged_obj),
  n_genes = nrow(merged_obj),
  n_ND_samples = sum(merge_summary$disease_group == "ND"),
  n_PD_samples = sum(merge_summary$disease_group == "PD"),
  n_T2D_samples = sum(merge_summary$disease_group == "T2D"),
  n_MIXED_samples = sum(merge_summary$disease_group == "MIXED"),
  n_conservative_DE_samples = sum(merge_summary$use_for_conservative_DE),
  stringsAsFactors = FALSE
)

write.csv(
  overall_summary,
  file.path(out_res_dir, "merged_overall_summary.csv"),
  row.names = FALSE
)

saveRDS(
  merged_obj,
  file.path(out_obj_dir, "GSE221156_54samples_QC_merged_raw.rds")
)

print(overall_summary)
message("Done: merged 54-sample Seurat object saved.")
