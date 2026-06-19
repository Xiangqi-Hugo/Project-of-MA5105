# T2D beta-cell pseudo-bulk analysis

This repository contains the analysis code for a human islet single-cell RNA-seq project on beta-cell dysfunction in type 2 diabetes.

## Study design

The primary analysis used GSE221156 human islet scRNA-seq data. Beta cells were selected by marker-based annotation and aggregated into sample-level pseudo-bulk profiles. This design avoided cell-level pseudoreplication.

A second CELLxGENE label-based sensitivity analysis used a beta-cell-only h5ad object. CELLxGENE beta cells were aggregated by LibraryID and compared with the marker-based beta-cell analysis.

## Main workflow

1. Build sample manifest.
2. Run per-sample quality control.
3. Define disease groups and MIXED samples.
4. Select beta cells by marker scoring.
5. Generate sample-level beta-cell pseudo-bulk counts.
6. Run edgeR differential expression analysis.
7. Run GO ORA and GO GSEA.
8. Run beta-cell threshold sensitivity analysis.
9. Compare marker-based beta-cell selection with CELLxGENE beta-cell labels.
10. Run miRNA inference using conservative GSEA, decoupleR, and miRSCAPE.

## Main conclusion

The beta-cell transcriptomic disease signal was robust. T2D beta cells showed strong gene-level and pathway-level changes. This signal was stable across beta-cell thresholds and was reproduced using CELLxGENE beta-cell labels.

The miRNA layer remained exploratory. No primary miRNA regulator was consistently validated across conservative GSEA, decoupleR, miRSCAPE Model 1, miRSCAPE Model 2, and indirect target-pathway analysis.

## Repository structure

scripts/ contains analysis scripts.

results_summary/ contains small summary tables.

tools/ contains external tool files if included.

## Data availability

Large files are not included. Raw 10x matrices, h5ad files, RDS objects, full pseudo-bulk matrices, and predicted miRNA expression matrices should be downloaded or generated locally.

Primary data source: GEO GSE221156.

## Notes

This repository documents the analysis workflow and provides reproducible code. It does not include large raw or intermediate data files.
