out_dir <- "results/19_GTEx_miRSCAPE_reference_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mrna_file <- "data/reference/GTEx/GTEx_mRNA_expression.gct.gz"
mirna_file <- "data/reference/GTEx/GTEx_miRNA_TPM_matrix.txt.gz"
sample_file <- "data/reference/GTEx/GTEx_sample_attributes.tsv"

get_gct_sample_ids <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con))
  readLines(con, n = 2)
  header <- readLines(con, n = 1)
  fields <- strsplit(header, "\t")[[1]]
  sample_ids <- fields[-c(1, 2)]
  return(sample_ids)
}

get_mirna_sample_ids <- function(file) {
  con <- gzfile(file, "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1)
  fields <- strsplit(header, "\t")[[1]]
  sample_ids <- fields[grepl("^GTEX-", fields)]
  return(sample_ids)
}

get_core_id <- function(x) {
  sub("-SM-.*$", "", x)
}

message("Reading mRNA header...")
mrna_samples <- get_gct_sample_ids(mrna_file)

message("Reading miRNA header...")
mirna_samples <- get_mirna_sample_ids(mirna_file)

message("Example mRNA sample IDs:")
print(head(mrna_samples, 10))

message("Example miRNA sample IDs:")
print(head(mirna_samples, 10))

exact_paired_samples <- intersect(mrna_samples, mirna_samples)

mrna_df <- data.frame(
  mrna_sample_id = mrna_samples,
  core_id = get_core_id(mrna_samples),
  stringsAsFactors = FALSE
)

mirna_df <- data.frame(
  mirna_sample_id = mirna_samples,
  core_id = get_core_id(mirna_samples),
  stringsAsFactors = FALSE
)

core_paired_ids <- intersect(mrna_df$core_id, mirna_df$core_id)

message("mRNA full sample IDs: ", length(mrna_samples))
message("miRNA full sample IDs: ", length(mirna_samples))
message("exact full-ID paired samples: ", length(exact_paired_samples))
message("core-ID paired biological samples: ", length(core_paired_ids))

write.csv(
  data.frame(sample_id = exact_paired_samples),
  file.path(out_dir, "GTEx_exact_fullID_paired_sample_ids.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(core_id = core_paired_ids),
  file.path(out_dir, "GTEx_coreID_paired_sample_ids.csv"),
  row.names = FALSE
)

# Check duplicate core IDs.
mrna_core_counts <- as.data.frame(table(mrna_df$core_id))
colnames(mrna_core_counts) <- c("core_id", "n_mRNA_ids")

mirna_core_counts <- as.data.frame(table(mirna_df$core_id))
colnames(mirna_core_counts) <- c("core_id", "n_miRNA_ids")

core_count_table <- merge(
  mrna_core_counts,
  mirna_core_counts,
  by = "core_id",
  all = TRUE
)

core_count_table$n_mRNA_ids[is.na(core_count_table$n_mRNA_ids)] <- 0
core_count_table$n_miRNA_ids[is.na(core_count_table$n_miRNA_ids)] <- 0

write.csv(
  core_count_table,
  file.path(out_dir, "GTEx_coreID_duplicate_check.csv"),
  row.names = FALSE
)

# Keep one mRNA and one miRNA sample ID per core ID.
# This is for audit and later training matrix construction.
mrna_unique <- mrna_df[!duplicated(mrna_df$core_id), ]
mirna_unique <- mirna_df[!duplicated(mirna_df$core_id), ]

pair_table <- merge(
  mrna_unique,
  mirna_unique,
  by = "core_id",
  all = FALSE
)

message("Reading sample attributes...")
sample_info <- read.delim(
  sample_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"SAMPID" %in% colnames(sample_info)) {
  stop("SAMPID column was not found in sample attributes.")
}

if (!"SMTSD" %in% colnames(sample_info)) {
  stop("SMTSD column was not found in sample attributes.")
}

sample_info$core_id <- get_core_id(sample_info$SAMPID)

# Add tissue annotation using mRNA full sample ID first.
idx_full <- match(pair_table$mrna_sample_id, sample_info$SAMPID)

pair_table$SMTS <- NA
pair_table$SMTSD <- NA
pair_table$ANALYTE_TYPE <- NA

pair_table$SMTS[!is.na(idx_full)] <- sample_info$SMTS[idx_full[!is.na(idx_full)]]
pair_table$SMTSD[!is.na(idx_full)] <- sample_info$SMTSD[idx_full[!is.na(idx_full)]]
pair_table$ANALYTE_TYPE[!is.na(idx_full)] <- sample_info$ANALYTE_TYPE[idx_full[!is.na(idx_full)]]

# If some tissue labels are still missing, fill by core ID.
missing_tissue <- is.na(pair_table$SMTSD)

if (any(missing_tissue)) {
  core_first <- sample_info[!duplicated(sample_info$core_id), ]
  idx_core <- match(pair_table$core_id[missing_tissue], core_first$core_id)
  
  pair_table$SMTS[missing_tissue] <- core_first$SMTS[idx_core]
  pair_table$SMTSD[missing_tissue] <- core_first$SMTSD[idx_core]
  pair_table$ANALYTE_TYPE[missing_tissue] <- core_first$ANALYTE_TYPE[idx_core]
}

write.csv(
  pair_table,
  file.path(out_dir, "GTEx_mRNA_miRNA_TPM_coreID_pair_table.csv"),
  row.names = FALSE
)

if (nrow(pair_table) == 0) {
  stop("No core-ID paired samples found. Need to inspect header format.")
}

tissue_counts <- as.data.frame(table(pair_table$SMTSD))
colnames(tissue_counts) <- c("tissue", "n_paired_samples")
tissue_counts <- tissue_counts[order(tissue_counts$n_paired_samples, decreasing = TRUE), ]

write.csv(
  tissue_counts,
  file.path(out_dir, "GTEx_mRNA_miRNA_TPM_paired_sample_counts_by_tissue.csv"),
  row.names = FALSE
)

pancreas_pairs <- pair_table[
  grepl("pancreas", pair_table$SMTSD, ignore.case = TRUE),
]

write.csv(
  pancreas_pairs,
  file.path(out_dir, "GTEx_mRNA_miRNA_TPM_paired_pancreas_samples.csv"),
  row.names = FALSE
)

summary_df <- data.frame(
  metric = c(
    "mRNA_full_sample_IDs",
    "miRNA_full_sample_IDs",
    "exact_full_ID_pairs",
    "core_ID_pairs",
    "pancreas_core_ID_pairs"
  ),
  value = c(
    length(mrna_samples),
    length(mirna_samples),
    length(exact_paired_samples),
    nrow(pair_table),
    nrow(pancreas_pairs)
  )
)

write.csv(
  summary_df,
  file.path(out_dir, "GTEx_mRNA_miRNA_pairing_summary.csv"),
  row.names = FALSE
)

cat("\nSummary:\n")
print(summary_df)

cat("\nTop paired tissues:\n")
print(head(tissue_counts, 30))

cat("\nPancreas paired samples:", nrow(pancreas_pairs), "\n")

cat("\nDone: GTEx mRNA-miRNA TPM core-ID pairing audit completed.\n")
