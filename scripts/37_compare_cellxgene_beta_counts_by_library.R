# Step 37: Compare our marker-based beta-cell counts with CELLxGENE beta-only metadata
#
# This script is designed for the CELLxGENE beta-cell-only H5AD file.
# In that file, cell_type is "type B pancreatic cell" for all cells.
#
# It compares:
#   our selected beta cells per sample
#   vs
#   CELLxGENE beta cells per LibraryID
#
# Main inputs:
#   results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv
#   data/reference/CELLxGENE/cell_metadata.csv
#
# Outputs:
#   results/37_cellxgene_beta_count_comparison/cellxgene_vs_our_beta_counts_by_library.csv
#   results/37_cellxgene_beta_count_comparison/cellxgene_vs_our_beta_count_summary.csv
#   results/37_cellxgene_beta_count_comparison/largest_beta_count_differences.csv

out_dir <- "results/37_cellxgene_beta_count_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

our_file <- "results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv"
cxg_file <- "data/reference/CELLxGENE/cell_metadata.csv"

if (!file.exists(our_file)) stop("Missing file: ", our_file)
if (!file.exists(cxg_file)) stop("Missing file: ", cxg_file)

find_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) return(hit[1])
  lower_names <- tolower(colnames(df))
  lower_candidates <- tolower(candidates)
  idx <- match(lower_candidates, lower_names)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(colnames(df)[idx[1]])
  if (required) stop("Cannot find required column. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

extract_library_id <- function(x) {
  x <- as.character(x)
  sub("^.*_(MS[0-9]+)$", "\\1", x)
}

our <- read.csv(our_file, stringsAsFactors = FALSE, check.names = FALSE)
cxg <- read.csv(cxg_file, stringsAsFactors = FALSE, check.names = FALSE)

our_sample_col <- find_col(our, c("sample_prefix", "sample", "Sample", "sample_id", "Sample_ID", "library", "LibraryID"))
our_beta_col <- find_col(our, c("n_beta_cells", "beta_cells", "Beta cells", "n_beta", "beta_cell_count", "selected_beta_cells"))
our_group_col <- find_col(our, c("disease_group", "group", "Group", "disease", "Disease", "condition", "DiseaseState"), required = FALSE)

cxg_library_col <- find_col(cxg, c("LibraryID", "library", "library_id", "sample", "sample_id"))
cxg_celltype_col <- find_col(cxg, c("cell_type", "celltype", "CellType", "author_cell_type", "annotation"))

our2 <- data.frame(
  our_sample = as.character(our[[our_sample_col]]),
  LibraryID = extract_library_id(our[[our_sample_col]]),
  our_beta_cells = as.numeric(our[[our_beta_col]]),
  stringsAsFactors = FALSE
)

if (!is.na(our_group_col)) {
  our2$group <- as.character(our[[our_group_col]])
}

cxg2 <- data.frame(
  LibraryID = as.character(cxg[[cxg_library_col]]),
  cell_type = as.character(cxg[[cxg_celltype_col]]),
  stringsAsFactors = FALSE
)

cxg2$is_beta <- grepl("type B pancreatic cell|beta|β", cxg2$cell_type, ignore.case = TRUE)

cxg_counts <- aggregate(
  cxg2$is_beta,
  by = list(LibraryID = cxg2$LibraryID),
  FUN = function(z) sum(z, na.rm = TRUE)
)
colnames(cxg_counts)[2] <- "cellxgene_beta_cells"

cmp <- merge(our2, cxg_counts, by = "LibraryID", all = TRUE)

cmp$difference_cellxgene_minus_ours <- cmp$cellxgene_beta_cells - cmp$our_beta_cells
cmp$absolute_difference <- abs(cmp$difference_cellxgene_minus_ours)
cmp$ratio_cellxgene_to_ours <- cmp$cellxgene_beta_cells / cmp$our_beta_cells

cmp <- cmp[order(-cmp$absolute_difference), ]

write.csv(cmp, file.path(out_dir, "cellxgene_vs_our_beta_counts_by_library.csv"), row.names = FALSE)

matched <- !is.na(cmp$our_beta_cells) & !is.na(cmp$cellxgene_beta_cells)
total_our <- sum(cmp$our_beta_cells[matched], na.rm = TRUE)
total_cxg <- sum(cmp$cellxgene_beta_cells[matched], na.rm = TRUE)
diff_total <- sum(cmp$difference_cellxgene_minus_ours[matched], na.rm = TRUE)

summary_df <- data.frame(
  metric = c(
    "n_libraries_in_our_summary",
    "n_libraries_in_cellxgene_metadata",
    "n_matched_libraries",
    "total_our_beta_cells_matched",
    "total_cellxgene_beta_cells_matched",
    "total_difference_cellxgene_minus_ours",
    "percent_difference_vs_cellxgene",
    "median_absolute_difference_per_library",
    "max_absolute_difference_per_library",
    "correlation_our_vs_cellxgene"
  ),
  value = c(
    length(unique(our2$LibraryID)),
    length(unique(cxg_counts$LibraryID)),
    sum(matched),
    total_our,
    total_cxg,
    diff_total,
    100 * diff_total / total_cxg,
    median(cmp$absolute_difference[matched], na.rm = TRUE),
    max(cmp$absolute_difference[matched], na.rm = TRUE),
    suppressWarnings(cor(cmp$our_beta_cells[matched], cmp$cellxgene_beta_cells[matched], use = "complete.obs"))
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, file.path(out_dir, "cellxgene_vs_our_beta_count_summary.csv"), row.names = FALSE)

largest <- head(cmp, 20)
write.csv(largest, file.path(out_dir, "largest_beta_count_differences.csv"), row.names = FALSE)

cat("\nCELLxGENE vs our beta-count summary:\n")
print(summary_df)

cat("\nLargest differences:\n")
print(largest)

cat("\nDone.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
