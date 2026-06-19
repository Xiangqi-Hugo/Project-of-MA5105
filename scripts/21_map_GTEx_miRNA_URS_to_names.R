suppressPackageStartupMessages({
  library(data.table)
})

in_dir <- "results/21_GTEx_miRSCAPE_training_input"
out_dir <- "results/21_GTEx_miRSCAPE_training_input"
map_file <- "data/reference/RNAcentral/rnacentral_mirbase_mappings.tsv"

mirna_rds <- file.path(in_dir, "GTEx_pancreas_bulk_miRNA_TPM.rds")

message("Reading GTEx pancreas miRNA matrix...")
mirna_mat <- readRDS(mirna_rds)

message("Current miRNA matrix dimension:")
print(dim(mirna_mat))

message("Current rowname examples:")
print(head(rownames(mirna_mat), 20))

message("Reading RNAcentral-miRBase mapping...")
map <- fread(
  map_file,
  header = FALSE,
  sep = "\t",
  data.table = FALSE
)

message("Mapping file dimension:")
print(dim(map))

message("First rows of mapping file:")
print(head(map))

# Detect URS column.
is_urs_col <- sapply(map, function(x) any(grepl("^URS[0-9A-F]+", x)))

# Detect hsa-miR / hsa-let column.
is_hsa_mir_col <- sapply(map, function(x) {
  any(grepl("^hsa-(miR|let)", x, ignore.case = TRUE))
})

urs_col <- which(is_urs_col)[1]
name_col <- which(is_hsa_mir_col)[1]

if (is.na(urs_col)) {
  stop("Cannot detect RNAcentral URS column.")
}

if (is.na(name_col)) {
  stop("Cannot detect hsa-miR name column.")
}

message("Detected URS column: ", urs_col)
message("Detected miRNA name column: ", name_col)

map_small <- data.frame(
  URS = sub("_.*$", "", as.character(map[[urs_col]])),
  miRNA_name = as.character(map[[name_col]]),
  stringsAsFactors = FALSE
)

map_small <- map_small[
  grepl("^URS", map_small$URS) &
    grepl("^hsa-(miR|let)", map_small$miRNA_name, ignore.case = TRUE),
]

map_small <- unique(map_small)

write.csv(
  map_small,
  file.path(out_dir, "RNAcentral_URS_to_hsa_miRNA_mapping_detected.csv"),
  row.names = FALSE
)

current_ids <- rownames(mirna_mat)

mapped_names <- map_small$miRNA_name[match(current_ids, map_small$URS)]

mapping_result <- data.frame(
  original_URS = current_ids,
  mapped_miRNA_name = mapped_names,
  stringsAsFactors = FALSE
)

write.csv(
  mapping_result,
  file.path(out_dir, "GTEx_pancreas_miRNA_URS_mapping_result.csv"),
  row.names = FALSE
)

n_mapped <- sum(!is.na(mapped_names))

message("Mapped miRNA rows: ", n_mapped, " / ", length(current_ids))

keep <- !is.na(mapped_names)

mirna_mapped <- mirna_mat[keep, , drop = FALSE]
rownames(mirna_mapped) <- mapped_names[keep]

# Aggregate duplicated miRNA names by mean.
sum_mat <- rowsum(mirna_mapped, group = rownames(mirna_mapped), reorder = FALSE)
counts <- table(rownames(mirna_mapped))
mirna_mapped <- sum_mat / as.numeric(counts[rownames(sum_mat)])
mirna_mapped <- as.matrix(mirna_mapped)

saveRDS(
  mirna_mapped,
  file.path(out_dir, "GTEx_pancreas_bulk_miRNA_TPM_miRBase_names.rds")
)

write.csv(
  mirna_mapped,
  file.path(out_dir, "GTEx_pancreas_bulk_miRNA_TPM_miRBase_names.csv"),
  row.names = TRUE
)

summary_df <- data.frame(
  item = c(
    "miRNA_rows_before_mapping",
    "miRNA_rows_mapped_to_hsa_miRBase",
    "miRNA_rows_after_name_aggregation"
  ),
  value = c(
    nrow(mirna_mat),
    n_mapped,
    nrow(mirna_mapped)
  )
)

write.csv(
  summary_df,
  file.path(out_dir, "GTEx_pancreas_miRNA_mapping_summary.csv"),
  row.names = FALSE
)

print(summary_df)

cat("\nMapped miRNA examples:\n")
print(head(rownames(mirna_mapped), 30))

cat("\nPrimary candidate check:\n")
candidate_patterns <- c(
  "hsa-miR-195-5p",
  "hsa-miR-16-5p",
  "hsa-miR-15a-5p",
  "hsa-miR-15b-5p"
)
print(candidate_patterns %in% rownames(mirna_mapped))

cat("\nDone: GTEx miRNA URS-to-name mapping completed.\n")
