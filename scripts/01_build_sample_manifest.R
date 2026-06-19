raw_dir <- "data/raw"
out_dir <- "results/01_manifest"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(raw_dir)

files <- files[grepl("barcodes|features|genes|matrix", files)]

file_type <- ifelse(grepl("barcodes", files), "barcodes",
             ifelse(grepl("features|genes", files), "features",
             ifelse(grepl("matrix", files), "matrix", "other")))

sample_prefix <- files
sample_prefix <- sub("_barcodes.*$", "", sample_prefix)
sample_prefix <- sub("_features.*$", "", sample_prefix)
sample_prefix <- sub("_genes.*$", "", sample_prefix)
sample_prefix <- sub("_matrix.*$", "", sample_prefix)

file_table <- data.frame(
  sample_prefix = sample_prefix,
  file_type = file_type,
  file = files,
  stringsAsFactors = FALSE
)

samples <- sort(unique(file_table$sample_prefix))

manifest <- data.frame(
  sample_prefix = samples,
  barcodes = NA_character_,
  features = NA_character_,
  matrix = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(samples)) {
  s <- samples[i]
  subtab <- file_table[file_table$sample_prefix == s, ]
  
  if (any(subtab$file_type == "barcodes")) {
    manifest$barcodes[i] <- subtab$file[subtab$file_type == "barcodes"][1]
  }
  
  if (any(subtab$file_type == "features")) {
    manifest$features[i] <- subtab$file[subtab$file_type == "features"][1]
  }
  
  if (any(subtab$file_type == "matrix")) {
    manifest$matrix[i] <- subtab$file[subtab$file_type == "matrix"][1]
  }
}

manifest$has_barcodes <- !is.na(manifest$barcodes)
manifest$has_features <- !is.na(manifest$features)
manifest$has_matrix <- !is.na(manifest$matrix)
manifest$complete_10x_triplet <- manifest$has_barcodes &
                                 manifest$has_features &
                                 manifest$has_matrix

write.csv(
  manifest,
  file.path(out_dir, "sample_manifest.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  n_samples_detected = nrow(manifest),
  n_complete_samples = sum(manifest$complete_10x_triplet),
  n_incomplete_samples = sum(!manifest$complete_10x_triplet)
)

write.csv(
  summary_table,
  file.path(out_dir, "sample_manifest_summary.csv"),
  row.names = FALSE
)

print(summary_table)

if (any(!manifest$complete_10x_triplet)) {
  warning("Some samples are incomplete. Please check results/01_manifest/sample_manifest.csv")
}
