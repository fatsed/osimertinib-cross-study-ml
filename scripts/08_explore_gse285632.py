"""Exploratory PCA of GSE285632 expression values (NOT a classification model).

Requirements: pandas, numpy, matplotlib, scikit-learn.

Run after step 07:
    python scripts/08_explore_gse285632.py

Input is the CSV exported by step 07. We preserve the workbook's original
abbreviated sample labels; the source-to-GSM mapping and measurement units
are not yet independently verified. Do not infer biological significance or
cross-study performance from this exploratory plot.
"""

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/processed/GSE285632_expression_as_provided.csv"
RESULTS = ROOT / "results/GSE285632/overview"
FIGURES = ROOT / "figures/GSE285632/overview"
LABELS = [f"{group}-{rep}" for group in ("P", "P+C", "D", "D+C", "R", "R+C") for rep in (1, 2, 3)]


def main():
    if not INPUT.is_file():
        raise FileNotFoundError("Run script 07 first. Missing: " + str(INPUT))
    df = pd.read_csv(INPUT)
    required = ["source_excel_row", "gene_symbol"] + LABELS
    if list(df.columns) != required:
        raise ValueError(f"Unexpected matrix columns: {list(df.columns)}")
    if df["gene_symbol"].isna().any():
        raise ValueError("Missing gene names in gene_symbol")

    data = df[LABELS].to_numpy(dtype=float)
    if not np.isfinite(data).all() or (data < 0).any():
        raise ValueError("The expression matrix contains missing, nonfinite, or negative values")

    # For this exploration only: remove *all* rows with ambiguous repeated gene
    # symbols rather than silently treating them as identical features.
    repeated = df["gene_symbol"].duplicated(keep=False).to_numpy()
    unique_data = data[~repeated]
    variance = unique_data.var(axis=1)
    eligible = np.flatnonzero(variance > 1e-10)
    if eligible.size < 2:
        raise ValueError("Not enough varying gene features for PCA")
    top = eligible[np.argsort(variance[eligible])[-min(500, eligible.size):]]

    # Samples x genes; the PCA is fit ONLY to make an exploratory plot.
    samples_by_genes = unique_data[top].T
    scaled = StandardScaler().fit_transform(samples_by_genes)
    pca = PCA(n_components=2, svd_solver="full")
    coords = pca.fit_transform(scaled)

    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    record = pd.DataFrame({
        "workbook_sample": LABELS,
        "group_as_written": [s.rsplit("-", 1)[0] for s in LABELS],
        "PC1": coords[:, 0], "PC2": coords[:, 1],
    })
    record.to_csv(RESULTS / "PCA_sample_coordinates.csv", index=False)

    fig, ax = plt.subplots(figsize=(9, 6.6))
    groups = ("P", "P+C", "D", "D+C", "R", "R+C")
    markers = ("o", "s", "^", "D", "v", "P")
    for group, marker in zip(groups, markers):
        inds = [i for i, name in enumerate(LABELS) if name.rsplit("-", 1)[0] == group]
        ax.scatter(coords[inds, 0], coords[inds, 1], marker=marker, s=78, label=group, alpha=0.86)

    pc1, pc2 = pca.explained_variance_ratio_ * 100
    ax.set(xlabel=f"PC1 ({pc1:.1f}% explained variance)",
           ylabel=f"PC2 ({pc2:.1f}% explained variance)",
           title="GSE285632 | Exploratory PCA of 18 workbook samples")
    ax.legend(title="Workbook labels (3 samples each)", frameon=False, ncol=3, loc="best")
    ax.margins(0.18)
    fig.text(0.5, 0.018,
             "Exploratory only; sample-level coordinates are saved separately. +C and GSM mapping unverified.",
             ha="center", fontsize=8, color="0.3")
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    image_file = FIGURES / "GSE285632_PCA_exploratory.png"
    fig.savefig(image_file, dpi=190)
    plt.close(fig)

    report = {
        "dataset": "GSE285632", "analysis": "exploratory PCA, no predictive modeling",
        "samples": len(LABELS), "total_gene_rows": len(df),
        "rows_excluded_due_to_repeated_gene_names": int(repeated.sum()),
        "genes_used_for_PCA": len(top),
        "pca_variance_explained_PC1": float(pca.explained_variance_ratio_[0]),
        "pca_variance_explained_PC2": float(pca.explained_variance_ratio_[1]),
        "preprocessing": "Top 500 variable nonduplicated gene symbols; gene-wise standardization; PCA",
        "limitations": [
            "Workbook column names are not yet confirmed against individual GEO GSM identifiers.",
            "The meaning of +C and the expression units are not yet verified.",
            "The plot uses all samples and cannot estimate prediction performance.",
            "All 18 samples belong to the same study; no cross-study validation was performed.",
            "The PCA is descriptive and gene selection here is not an ML feature-selection result.",
        ],
    }
    (RESULTS / "PCA_exploratory_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("GSE285632 EXPLORATORY PCA COMPLETED (NO MODEL TRAINED)")
    print(f"Samples: {len(LABELS)} | input rows: {len(df)} | PCA genes: {len(top)}")
    print(f"Explained variance: PC1={pc1:.1f}%, PC2={pc2:.1f}%")
    print("Figure:", image_file.relative_to(ROOT))
    print("Report:", (RESULTS / "PCA_exploratory_report.json").relative_to(ROOT))
    print("WARNING: +C and sample-to-GSM mapping remain unverified.")


if __name__ == "__main__":
    main()
