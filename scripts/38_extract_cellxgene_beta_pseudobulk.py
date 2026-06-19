#!/usr/bin/env python3
"""
Step 38: Build CELLxGENE beta-cell pseudo-bulk matrices.

Input:
  data/reference/CELLxGENE/cellxgene_direct_download.h5ad
  results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv

Main outputs:
  results/38_cellxgene_beta_pseudobulk/
    cellxgene_beta_pseudobulk_counts_matrix_all_libraries.csv
    cellxgene_beta_pseudobulk_sample_summary_all_libraries.csv
    cellxgene_beta_pseudobulk_counts_matrix_matched41.csv
    cellxgene_beta_pseudobulk_sample_summary_matched41.csv
    cellxgene_beta_pseudobulk_overall_summary.csv

Design:
  - Use CELLxGENE beta-only H5AD.
  - Aggregate beta cells by LibraryID.
  - Use LibraryID as the sample-level pseudo-bulk unit.
  - Create a matched41 matrix for the same conservative sample set used in the main report:
      clear ND/PD/T2D samples
      excluding MIXED samples
      excluding MS17003 low-beta outlier
"""

from __future__ import annotations

from pathlib import Path
import re
import sys
import subprocess


H5AD_PATH = Path("data/reference/CELLxGENE/cellxgene_direct_download.h5ad")
OUR_SUMMARY = Path("results/06_beta_pseudobulk/beta_pseudobulk_sample_summary.csv")
OUT_DIR = Path("results/38_cellxgene_beta_pseudobulk")


def install_if_missing(package: str, import_name: str | None = None) -> None:
    import_name = import_name or package
    try:
        __import__(import_name)
    except ImportError:
        print(f"[INFO] Installing missing package: {package}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])


def map_disease(x: str) -> str:
    x0 = str(x).lower()
    if "pre" in x0:
        return "PD"
    if "type 2" in x0 or "t2d" in x0 or "diabetes mellitus" in x0 or "diabetic" in x0:
        return "T2D"
    if "normal" in x0 or "non-diabetic" in x0 or "healthy" in x0 or x0 == "nd":
        return "ND"
    return "UNKNOWN"


def extract_ms_id(x: str) -> str:
    x = str(x)
    m = re.search(r"(MS[0-9]+)", x)
    return m.group(1) if m else x


def main() -> None:
    if not H5AD_PATH.exists():
        raise FileNotFoundError(f"Missing H5AD file: {H5AD_PATH}")
    if not OUR_SUMMARY.exists():
        raise FileNotFoundError(f"Missing original sample summary: {OUR_SUMMARY}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    install_if_missing("anndata", "anndata")
    install_if_missing("pandas", "pandas")
    install_if_missing("scipy", "scipy")
    install_if_missing("numpy", "numpy")

    import numpy as np
    import pandas as pd
    import anndata as ad
    import scipy.sparse as sp

    print("[INFO] Reading H5AD in backed mode...")
    adata = ad.read_h5ad(H5AD_PATH, backed="r")
    obs = adata.obs.copy()
    var = adata.var.copy()

    if "LibraryID" not in obs.columns:
        raise RuntimeError("CELLxGENE metadata does not contain LibraryID.")
    if "disease" not in obs.columns:
        raise RuntimeError("CELLxGENE metadata does not contain disease.")

    # Determine gene names.
    if "feature_name" in var.columns:
        gene_names = var["feature_name"].astype(str).values
    elif "gene_symbols" in var.columns:
        gene_names = var["gene_symbols"].astype(str).values
    else:
        gene_names = var.index.astype(str).values

    # De-duplicate gene names by summing later if needed.
    library_ids = obs["LibraryID"].astype(str).values
    unique_libs = sorted(pd.unique(library_ids))

    print(f"[INFO] Cells in CELLxGENE beta H5AD: {adata.n_obs:,}")
    print(f"[INFO] Genes/features: {adata.n_vars:,}")
    print(f"[INFO] Unique LibraryID units: {len(unique_libs)}")

    pseudobulk_rows = []
    summary_rows = []

    # Matrix will be genes x libraries.
    # Use per-library slicing to avoid loading the full matrix into dense memory.
    for i, lib in enumerate(unique_libs, start=1):
        idx = np.where(library_ids == lib)[0]
        print(f"[INFO] Aggregating {i}/{len(unique_libs)}: {lib} ({len(idx)} cells)")

        X = adata.X[idx, :]
        if sp.issparse(X):
            summed = np.asarray(X.sum(axis=0)).ravel()
        else:
            summed = np.asarray(X).sum(axis=0).ravel()

        # Round only if values are very close to integers.
        near_integer_fraction = np.mean(np.isclose(summed, np.round(summed), atol=1e-6))
        if near_integer_fraction > 0.999:
            summed = np.round(summed).astype(np.int64)

        pseudobulk_rows.append(summed)

        sub = obs.iloc[idx]
        disease_values = sub["disease"].astype(str).value_counts(dropna=False)
        disease_main = disease_values.index[0]
        group_main = map_disease(disease_main)

        islet_values = (
            ";".join(sorted(pd.unique(sub["Islet"].astype(str))))
            if "Islet" in sub.columns else ""
        )

        summary_rows.append({
            "LibraryID": lib,
            "MS_ID": extract_ms_id(lib),
            "n_cellxgene_beta_cells": len(idx),
            "disease_raw_main": disease_main,
            "disease_group": group_main,
            "n_disease_labels_in_unit": len(disease_values),
            "disease_label_counts": ";".join([f"{k}:{v}" for k, v in disease_values.items()]),
            "Islet": islet_values,
            "is_split_islet_library": bool("_Islet" in lib),
        })

    # Build matrix genes x libraries.
    mat = np.vstack(pseudobulk_rows).T
    pb = pd.DataFrame(mat, index=gene_names, columns=unique_libs)

    # If duplicated gene names exist, sum duplicates.
    if pb.index.duplicated().any():
        print("[WARN] Duplicated gene names detected. Summing duplicated genes.")
        pb = pb.groupby(pb.index).sum()

    summary = pd.DataFrame(summary_rows)

    # Save all libraries.
    all_counts_path = OUT_DIR / "cellxgene_beta_pseudobulk_counts_matrix_all_libraries.csv"
    all_summary_path = OUT_DIR / "cellxgene_beta_pseudobulk_sample_summary_all_libraries.csv"
    pb.to_csv(all_counts_path)
    summary.to_csv(all_summary_path, index=False)

    # Build matched41 set from original conservative sample summary.
    our = pd.read_csv(OUR_SUMMARY)
    sample_col = None
    for c in ["sample_prefix", "sample", "Sample", "sample_id", "LibraryID"]:
        if c in our.columns:
            sample_col = c
            break
    if sample_col is None:
        raise RuntimeError("Could not find sample column in original sample summary.")

    group_col = None
    for c in ["disease_group", "group", "Group", "disease"]:
        if c in our.columns:
            group_col = c
            break
    if group_col is None:
        raise RuntimeError("Could not find disease group column in original sample summary.")

    beta_col = None
    for c in ["n_beta_cells", "beta_cells", "selected_beta_cells", "n_beta"]:
        if c in our.columns:
            beta_col = c
            break
    if beta_col is None:
        raise RuntimeError("Could not find beta-cell count column in original sample summary.")

    our["MS_ID"] = our[sample_col].astype(str).apply(extract_ms_id)
    our["group_original"] = our[group_col].astype(str)
    our["n_beta_cells_original"] = our[beta_col].astype(float)

    target = our[
        (our["group_original"].isin(["ND", "PD", "T2D"])) &
        (our["MS_ID"] != "MS17003") &
        (our["n_beta_cells_original"] >= 100)
    ].copy()

    target_libs = [x for x in target["MS_ID"].tolist() if x in pb.columns]
    missing = [x for x in target["MS_ID"].tolist() if x not in pb.columns]

    print(f"[INFO] Matched conservative libraries found in CELLxGENE matrix: {len(target_libs)}")
    if missing:
        print("[WARN] Conservative libraries missing from CELLxGENE matrix:", ",".join(missing))

    matched_counts = pb.loc[:, target_libs]
    matched_summary = summary[summary["LibraryID"].isin(target_libs)].copy()
    matched_summary = matched_summary.merge(
        target[["MS_ID", "group_original", "n_beta_cells_original"]],
        on="MS_ID",
        how="left"
    )
    # Use original disease group for strict comparability.
    matched_summary["disease_group_for_DE"] = matched_summary["group_original"]

    matched_counts_path = OUT_DIR / "cellxgene_beta_pseudobulk_counts_matrix_matched41.csv"
    matched_summary_path = OUT_DIR / "cellxgene_beta_pseudobulk_sample_summary_matched41.csv"
    matched_counts.to_csv(matched_counts_path)
    matched_summary.to_csv(matched_summary_path, index=False)

    overall = pd.DataFrame([
        {"metric": "n_all_cellxgene_beta_cells", "value": int(summary["n_cellxgene_beta_cells"].sum())},
        {"metric": "n_all_library_units", "value": int(summary.shape[0])},
        {"metric": "n_matched41_library_units", "value": int(matched_summary.shape[0])},
        {"metric": "n_matched41_cellxgene_beta_cells", "value": int(matched_summary["n_cellxgene_beta_cells"].sum())},
        {"metric": "n_matched41_original_marker_beta_cells", "value": int(matched_summary["n_beta_cells_original"].sum())},
    ])
    overall.to_csv(OUT_DIR / "cellxgene_beta_pseudobulk_overall_summary.csv", index=False)

    print("\n[DONE] CELLxGENE beta pseudo-bulk extraction completed.")
    print("[INFO] Outputs written to:", OUT_DIR)
    print("\nOverall summary:")
    print(overall.to_string(index=False))


if __name__ == "__main__":
    main()
