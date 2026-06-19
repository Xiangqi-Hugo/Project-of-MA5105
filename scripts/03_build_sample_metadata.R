manifest_file <- "results/01_manifest/sample_manifest.csv"
qc_file <- "results/02_qc/sample_qc_summary.csv"

out_dir <- "results/03_metadata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
qc <- read.csv(qc_file, stringsAsFactors = FALSE)

# Extract GSM accession from local sample_prefix.
manifest$GSM <- sub("^(GSM[0-9]+).*", "\\1", manifest$sample_prefix)

# GEO sample title mapping from GSE221156 GEO page.
geo_map <- data.frame(
  GSM = c(
    "GSM6846488","GSM6846489","GSM6846490","GSM6846491","GSM6846492","GSM6846493",
    "GSM6846494","GSM6846495","GSM6846496","GSM6846497","GSM6846498","GSM6846499",
    "GSM6846500","GSM6846501","GSM6846502","GSM6846503","GSM6846504","GSM6846505",
    "GSM6846506","GSM6846507","GSM6846508","GSM6846509","GSM6846510","GSM6846511",
    "GSM6846512","GSM6846513","GSM6846514","GSM6846515","GSM6846516","GSM6846517",
    "GSM6846518","GSM6846519","GSM6846520","GSM6846521","GSM6846522","GSM6846523",
    "GSM6846524","GSM6846525","GSM6846526","GSM6846527","GSM6846528","GSM6846529",
    "GSM6846530","GSM6846531","GSM6846532","GSM6846533","GSM6846534","GSM6846535",
    "GSM6846536","GSM6846537","GSM6846538","GSM6846539","GSM6846540","GSM6846541"
  ),
  islet_title = c(
    "Islet28","Islet29","Islet30","Islet31","Islet32","Islet33",
    "Islet34","Islet37","Islet38","Islet39","Islet40","Islet41",
    "Islet42","Islet44","Islet45","Islet47","Islet48","Islet47/48_MS19006",
    "Islet47/48_MS19007","Islet50","Islet52","Islet53","Islet54","Islet55",
    "Islet56","Islet57","Islet58","Islet57/58_MS19016","Islet57/58_MS19017","Islet59",
    "Islet60","Islet59/60_MS19020","Islet59/60_MS19021","Islet61","Islet62","Islet63",
    "Islet68","Islet70/71_MS19034","Islet70/71_MS19035","Islet73","Islet80","Islet84/85_MS19042",
    "Islet84/85_MS19043","Islet89","Islet91","Islet101","Islet104","Islet107",
    "Islet118/119_MS20086","Islet118/119_MS20087","Islet123","Islet126","Islet127","Islet128"
  ),
  stringsAsFactors = FALSE
)

ND_islets <- c(
  "Islet29","Islet34","Islet37","Islet38","Islet42","Islet45","Islet47",
  "Islet54","Islet57","Islet59","Islet61","Islet68","Islet73","Islet80",
  "Islet104","Islet107","Islet123"
)

PD_islets <- c(
  "Islet40","Islet50","Islet53","Islet55","Islet56","Islet58","Islet62",
  "Islet63","Islet70","Islet85","Islet89","Islet101","Islet118","Islet127"
)

T2D_islets <- c(
  "Islet28","Islet30","Islet31","Islet32","Islet33","Islet39","Islet41",
  "Islet44","Islet48","Islet52","Islet60","Islet71","Islet84","Islet91",
  "Islet126","Islet119","Islet128"
)

# Extract clean Islet IDs.
# For simple samples: Islet28 -> Islet28
# For mixed samples: Islet47/48_MS19006 -> Islet47/48
geo_map$is_mixed_islet <- grepl("/", geo_map$islet_title)

geo_map$islet_pair <- sub("_MS.*$", "", geo_map$islet_title)

assign_single_disease <- function(x) {
  if (x %in% ND_islets) return("ND")
  if (x %in% PD_islets) return("PD")
  if (x %in% T2D_islets) return("T2D")
  return(NA_character_)
}

assign_pair_disease <- function(pair) {
  # pair example: Islet47/48
  if (!grepl("/", pair)) {
    return(assign_single_disease(pair))
  }
  
  nums <- unlist(strsplit(sub("^Islet", "", pair), "/"))
  islets <- paste0("Islet", nums)
  groups <- sapply(islets, assign_single_disease)
  groups <- unique(groups)
  
  if (length(groups) == 1 && !is.na(groups[1])) {
    return(groups[1])
  }
  
  return("MIXED")
}

geo_map$disease_group <- sapply(geo_map$islet_pair, assign_pair_disease)
geo_map$donor_or_pair_id <- geo_map$islet_pair

metadata <- merge(
  manifest,
  geo_map,
  by = "GSM",
  all.x = TRUE
)

metadata <- merge(
  metadata,
  qc,
  by.x = "sample_prefix",
  by.y = "sample_id",
  all.x = TRUE
)

metadata <- metadata[order(metadata$GSM), ]

write.csv(
  metadata,
  file.path(out_dir, "sample_metadata_54.csv"),
  row.names = FALSE
)

group_summary <- as.data.frame(table(metadata$disease_group, useNA = "ifany"))
colnames(group_summary) <- c("disease_group", "n_samples")

write.csv(
  group_summary,
  file.path(out_dir, "sample_metadata_group_summary.csv"),
  row.names = FALSE
)

mixed_samples <- metadata[metadata$disease_group == "MIXED", ]

write.csv(
  mixed_samples,
  file.path(out_dir, "mixed_or_ambiguous_samples.csv"),
  row.names = FALSE
)

print(group_summary)

cat("\nMixed or ambiguous samples:\n")
print(mixed_samples[, c("sample_prefix", "islet_title", "donor_or_pair_id", "disease_group", "n_cells_after_QC")])
