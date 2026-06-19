out_dir <- "results/18_miRSCAPE_setup"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(xgboost)
  library(Seurat)
})

source("tools/miRSCAPE/code/miRSCAPE.R")

setup <- data.frame(
  item = c(
    "xgboost_loaded",
    "Seurat_loaded",
    "miRSCAPE_function_exists",
    "bulkTransform_function_exists",
    "modifySeuratObject_function_exists"
  ),
  value = c(
    TRUE,
    TRUE,
    exists("miRSCAPE"),
    exists("bulkTransform"),
    exists("modifySeuratObject")
  )
)

write.csv(
  setup,
  file.path(out_dir, "miRSCAPE_setup_check.csv"),
  row.names = FALSE
)

sink(file.path(out_dir, "miRSCAPE_sessionInfo.txt"))
print(sessionInfo())
sink()

print(setup)

cat("Done: miRSCAPE setup check completed.\n")
