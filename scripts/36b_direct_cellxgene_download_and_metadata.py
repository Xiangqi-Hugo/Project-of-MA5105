#!/usr/bin/env python3
"""
Step 36b: Direct CELLxGENE H5AD download and metadata export.

Use this when you already have a direct CELLxGENE H5AD link, for example:
https://datasets.cellxgene.cziscience.com/faacd98e-af54-459d-8a7a-74427766bc81.h5ad

Outputs:
  data/reference/CELLxGENE/cellxgene_direct_download.h5ad
  data/reference/CELLxGENE/cell_metadata.csv
  data/reference/CELLxGENE/cell_metadata_columns.txt
  data/reference/CELLxGENE/cellxgene_celltype_counts.csv
  data/reference/CELLxGENE/cellxgene_sample_celltype_summary.csv
"""

from __future__ import annotations

import sys
import subprocess
from pathlib import Path
from urllib.request import urlretrieve


URL = "https://datasets.cellxgene.cziscience.com/faacd98e-af54-459d-8a7a-74427766bc81.h5ad"
OUT_DIR = Path("data/reference/CELLxGENE")
H5AD_PATH = OUT_DIR / "cellxgene_direct_download.h5ad"


def install_if_missing(package: str, import_name: str | None = None) -> None:
    import_name = import_name or package
    try:
        __import__(import_name)
    except ImportError:
        print(f"[INFO] Installing missing package: {package}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])


def find_col(columns, candidates):
    lower_map = {str(c).lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    for col in columns:
        low = str(col).lower()
        for cand in candidates:
            if cand.lower() in low:
                return col
    return None


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    install_if_missing("anndata", "anndata")
    install_if_missing("pandas", "pandas")

    import anndata as ad
    import pandas as pd

    if not H5AD_PATH.exists():
        print(f"[INFO] Downloading H5AD from:\n{URL}")
        print(f"[INFO] Output:\n{H5AD_PATH}")
        print("[INFO] This file may be large. Please keep the network stable.")
        urlretrieve(URL, H5AD_PATH)
    else:
        print(f"[INFO] H5AD already exists. Skipping download:\n{H5AD_PATH}")

    print("[INFO] Reading H5AD obs metadata in backed mode...")
    adata = ad.read_h5ad(H5AD_PATH, backed="r")
    obs = adata.obs.reset_index()

    print(f"[INFO] Number of cells: {obs.shape[0]}")
    print(f"[INFO] Number of metadata columns: {obs.shape[1]}")

    columns_path = OUT_DIR / "cell_metadata_columns.txt"
    with open(columns_path, "w", encoding="utf-8") as f:
        for col in obs.columns:
            f.write(str(col) + "\n")
    print(f"[INFO] Wrote metadata columns: {columns_path}")

    meta_path = OUT_DIR / "cell_metadata.csv"
    obs.to_csv(meta_path, index=False)
    print(f"[INFO] Wrote full cell metadata: {meta_path}")

    celltype_col = find_col(
        obs.columns,
        [
            "cell_type",
            "celltype",
            "CellType",
            "cell_type_ontology_term_label",
            "author_cell_type",
            "annotation",
            "cell_type_original",
            "cell_type_label"
        ],
    )

    sample_col = find_col(
        obs.columns,
        [
            "sample",
            "sample_id",
            "Sample",
            "Sample_ID",
            "library",
            "LibraryID",
            "library_id",
            "donor_id",
            "donor",
            "assay"
        ],
    )

    print("\n[INFO] Candidate cell-type column:", celltype_col)
    print("[INFO] Candidate sample column:", sample_col)

    if celltype_col is not None:
        counts = (
            obs[celltype_col]
            .astype(str)
            .value_counts(dropna=False)
            .reset_index()
        )
        counts.columns = [celltype_col, "n_cells"]
        counts.to_csv(OUT_DIR / "cellxgene_celltype_counts.csv", index=False)
        print("[INFO] Wrote cellxgene_celltype_counts.csv")
        print("\n[INFO] Top cell types:")
        print(counts.head(30).to_string(index=False))
    else:
        print("[WARN] No obvious cell-type column found. Inspect cell_metadata_columns.txt.")

    if sample_col is not None and celltype_col is not None:
        summary = (
            obs.groupby([sample_col, celltype_col], dropna=False)
            .size()
            .reset_index(name="n_cells")
            .sort_values([sample_col, "n_cells"], ascending=[True, False])
        )
        summary.to_csv(OUT_DIR / "cellxgene_sample_celltype_summary.csv", index=False)
        print("[INFO] Wrote cellxgene_sample_celltype_summary.csv")
    else:
        print("[WARN] Sample-celltype summary was not generated because one required column was missing.")

    print("\n[DONE] Direct CELLxGENE H5AD metadata export completed.")
    print("\n[NEXT] Run:")
    print("Rscript scripts/35_beta_identity_audit_and_cellxgene_template.R")
    print("\nThen send:")
    print("cat results/35_beta_identity_audit/optional_cellxgene_beta_count_comparison.csv")


if __name__ == "__main__":
    main()
