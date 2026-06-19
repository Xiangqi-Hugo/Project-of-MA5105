suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    BiocManager::install("TCGAbiolinks", ask = FALSE, update = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
})

out_dir <- "results/25_TCGA_PAAD_Model2_input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

project <- "TCGA-PAAD"

gdc_mrna_dir <- "data/reference/TCGA_PAAD_Model2/GDC_mRNA"
gdc_mirna_dir <- "data/reference/TCGA_PAAD_Model2/GDC_miRNA_clean"
dir.create(gdc_mrna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gdc_mirna_dir, recursive = TRUE, showWarnings = FALSE)

mrna_se_file <- file.path(out_dir, "TCGA_PAAD_mRNA_STAR_Counts_SE.rds")

mean_by_rownames <- function(mat) {
  sum_mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  counts <- table(rownames(mat))
  sum_mat <- sum_mat / as.numeric(counts[rownames(sum_mat)])
  as.matrix(sum_mat)
}

get_tcga_sample_core <- function(x) {
  substr(x, 1, 16)
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

if (file.exists(mrna_se_file)) {
  message("Using existing mRNA SummarizedExperiment:")
  message(mrna_se_file)
  mrna_se <- readRDS(mrna_se_file)
} else {
  message("Querying TCGA-PAAD mRNA data...")

  query_mrna <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )

  message("Downloading TCGA-PAAD mRNA data...")
  GDCdownload(
    query_mrna,
    method = "api",
    files.per.chunk = 10,
    directory = gdc_mrna_dir
  )

  message("Preparing TCGA-PAAD mRNA data...")
  mrna_se <- GDCprepare(query_mrna, directory = gdc_mrna_dir)
  saveRDS(mrna_se, mrna_se_file)
}

message("mRNA assay names:")
print(assayNames(mrna_se))

assay_names <- assayNames(mrna_se)

if ("tpm_unstranded" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "tpm_unstranded")
} else if ("fpkm_unstranded" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "fpkm_unstranded")
} else if ("fpkm_unstrand" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "fpkm_unstrand")
} else if ("unstranded" %in% assay_names) {
  mrna_mat <- assay(mrna_se, "unstranded")
} else {
  mrna_mat <- assay(mrna_se, 1)
}

mrna_meta <- as.data.frame(rowData(mrna_se))

symbol_col <- find_col(
  mrna_meta,
  c("gene_name", "external_gene_name", "gene_symbol", "symbol")
)

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
mrna_mat <- mean_by_rownames(mrna_mat)

message("Final TCGA-PAAD mRNA matrix dimension:")
print(dim(mrna_mat))

message("Querying TCGA-PAAD miRNA data...")

query_mirna <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = "Primary Tumor"
)

mirna_results <- getResults(query_mirna)

write.csv(
  mirna_results,
  file.path(out_dir, "TCGA_PAAD_miRNA_GDCquery_results.csv"),
  row.names = FALSE
)

message("miRNA files in query:")
print(nrow(mirna_results))

existing_files <- list.files(gdc_mirna_dir, recursive = TRUE, full.names = TRUE)
existing_quant_files <- existing_files[
  grepl("mirna|mirnas|quantification", basename(existing_files), ignore.case = TRUE) &
    grepl("\\.txt$|\\.tsv$", basename(existing_files), ignore.case = TRUE)
]

if (length(existing_quant_files) == 0) {
  message("No local miRNA quantification txt files found. Downloading with files.per.chunk = 1...")

  GDCdownload(
    query_mirna,
    method = "api",
    files.per.chunk = 1,
    directory = gdc_mirna_dir
  )
} else {
  message("Existing miRNA quantification files found. Reusing local files.")
  message(length(existing_quant_files))
}

files <- list.files(gdc_mirna_dir, recursive = TRUE, full.names = TRUE)

files <- files[
  grepl("mirna|mirnas|quantification", basename(files), ignore.case = TRUE) &
    grepl("\\.txt$|\\.tsv$", basename(files), ignore.case = TRUE)
]

if (length(files) == 0) {
  stop("No miRNA quantification txt files were found after download.")
}

message("Direct miRNA txt files found:")
print(length(files))
print(head(files))

file_col <- find_col(mirna_results, c("file_name", "filename", "File.Name"))
case_col <- find_col(mirna_results, c("cases", "case_submitter_id", "submitter_id", "sample_submitter_id"))

if (is.na(file_col)) {
  stop("Could not find file_name column in GDC query results.")
}

if (is.na(case_col)) {
  stop("Could not find cases/sample ID column in GDC query results.")
}

query_map <- data.frame(
  file_name = as.character(mirna_results[[file_col]]),
  sample_id = as.character(mirna_results[[case_col]]),
  stringsAsFactors = FALSE
)

query_map$sample_id <- sub("[,;].*$", "", query_map$sample_id)

file_map <- data.frame(
  path = files,
  file_name = basename(files),
  stringsAsFactors = FALSE
)

idx <- match(file_map$file_name, query_map$file_name)
file_map$sample_id <- query_map$sample_id[idx]

missing_sample <- is.na(file_map$sample_id)

if (any(missing_sample)) {
  message("Some files did not match exactly by file_name. Trying partial matching.")

  for (i in which(missing_sample)) {
    hit <- which(query_map$file_name == file_map$file_name[i])

    if (length(hit) == 0) {
      hit <- grep(file_map$file_name[i], query_map$file_name, fixed = TRUE)
    }

    if (length(hit) == 1) {
      file_map$sample_id[i] <- query_map$sample_id[hit]
    }
  }
}

file_map <- file_map[!is.na(file_map$sample_id), ]

if (nrow(file_map) == 0) {
  stop("No miRNA files could be mapped to TCGA sample IDs.")
}

write.csv(
  file_map,
  file.path(out_dir, "TCGA_PAAD_miRNA_file_to_sample_map.csv"),
  row.names = FALSE
)

one <- fread(file_map$path[1], data.table = FALSE)

message("Example miRNA file columns:")
print(colnames(one))
print(head(one))

id_col <- find_col(
  one,
  c("miRNA_ID", "miRNA", "mirna_id", "miRNA_IDs")
)

if (is.na(id_col)) {
  id_col <- colnames(one)[1]
}

value_col <- find_col(
  one,
  c("reads_per_million_miRNA_mapped", "RPM", "rpm", "read_count", "ReadCount")
)

if (is.na(value_col)) {
  numeric_cols <- colnames(one)[sapply(one, is.numeric)]
  if (length(numeric_cols) == 0) {
    stop("No numeric miRNA expression column was detected.")
  }
  value_col <- numeric_cols[1]
}

message("Using miRNA ID column: ", id_col)
message("Using miRNA value column: ", value_col)

value_list <- list()

for (i in seq_len(nrow(file_map))) {
  d <- fread(file_map$path[i], data.table = FALSE)

  ids <- as.character(d[[id_col]])
  vals <- suppressWarnings(as.numeric(d[[value_col]]))

  keep <- !is.na(ids) & ids != "" & !is.na(vals)
  ids <- ids[keep]
  vals <- vals[keep]

  names(vals) <- ids
  value_list[[file_map$sample_id[i]]] <- vals
}

all_ids <- sort(unique(unlist(lapply(value_list, names))))

mirna_mat <- matrix(
  0,
  nrow = length(all_ids),
  ncol = length(value_list),
  dimnames = list(all_ids, names(value_list))
)

for (sample_id in names(value_list)) {
  vals <- value_list[[sample_id]]
  mirna_mat[names(vals), sample_id] <- vals
}

storage.mode(mirna_mat) <- "numeric"

keep_mirna <- !is.na(rownames(mirna_mat)) &
  rownames(mirna_mat) != "" &
  grepl("^hsa-", rownames(mirna_mat), ignore.case = TRUE)

mirna_mat <- mirna_mat[keep_mirna, , drop = FALSE]
mirna_mat <- mean_by_rownames(mirna_mat)

message("Final TCGA-PAAD miRNA matrix dimension:")
print(dim(mirna_mat))

mrna_df <- data.frame(
  mrna_col = colnames(mrna_mat),
  core = get_tcga_sample_core(colnames(mrna_mat)),
  stringsAsFactors = FALSE
)

mirna_df <- data.frame(
  mirna_col = colnames(mirna_mat),
  core = get_tcga_sample_core(colnames(mirna_mat)),
  stringsAsFactors = FALSE
)

mrna_df <- mrna_df[!duplicated(mrna_df$core), ]
mirna_df <- mirna_df[!duplicated(mirna_df$core), ]

pair_df <- merge(mrna_df, mirna_df, by = "core", all = FALSE)

message("Paired TCGA-PAAD mRNA-miRNA samples:")
print(nrow(pair_df))

if (nrow(pair_df) < 30) {
  stop("Too few paired TCGA-PAAD mRNA-miRNA samples. Check file mapping.")
}

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

target_file <- "results/07_beta_DE_edgeR/beta_pseudobulk_logCPM_TMM.csv"

target <- read.csv(
  target_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

target <- as.matrix(target)
storage.mode(target) <- "numeric"

common_genes <- intersect(rownames(mrna_paired), rownames(target))

message("Common genes between TCGA-PAAD mRNA and beta-cell target:")
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

primary <- c(
  "hsa-miR-195-5p",
  "hsa-miR-16-5p",
  "hsa-miR-15a-5p",
  "hsa-miR-15b-5p",
  "hsa-miR-649",
  "hsa-miR-6838-5p"
)

availability <- data.frame(
  miRNA = primary,
  available_in_Model2_TCGA_PAAD = primary %in% rownames(mirna_paired),
  stringsAsFactors = FALSE
)

write.csv(
  availability,
  file.path(out_dir, "TCGA_PAAD_Model2_primary_candidate_availability.csv"),
  row.names = FALSE
)

print(availability)

cat("\nDone: direct-parse TCGA-PAAD Model2 input prepared.\n")
