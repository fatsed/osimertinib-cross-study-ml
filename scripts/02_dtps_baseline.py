"""Sample audit and a small DTP-versus-control baseline for GSE193258.

This is a within-study experiment, with one cell line held out at a time.
It does not classify acquired resistance or establish cross-study generalization.
"""

import csv
import gzip
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (ConfusionMatrixDisplay, accuracy_score,
                             balanced_accuracy_score, roc_auc_score)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = ROOT / "data/raw/GSE193258_RNAseq_log2TPM_abundance.tsv.gz"
RESULTS = ROOT / "results" / "GSE193258" / "dtp_vs_control"
FIGURES = ROOT / "figures" / "GSE193258" / "dtp_vs_control"
DATA_CHECK = ROOT / "results" / "GSE193258" / "data_check"
CELL_LINES = ("H1975", "HCC2935", "HCC827", "PC9")


def sample_info(name):
    """Use the sample naming scheme from the official GSE193258 record."""
    match = re.fullmatch(r"(H1975|HCC2935|HCC827|PC9)_(.+)_([123])", name)
    if match is None:
        raise ValueError(f"Unrecognized sample column: {name}")
    cell_line, treatment, replicate = match.groups()
    if treatment == "DMSO":
        state = "control"
    elif treatment == "osi_DTP":
        state = "DTP"
    elif treatment == "osi_acute":
        state = "acute"
    elif treatment.startswith("short_wash"):
        state = "short_washout"
    elif treatment.startswith("long_wash"):
        state = "long_washout"
    else:
        raise ValueError(f"Treatment needs manual verification: {name}")
    return {"sample": name, "cell_line": cell_line,
            "condition": state, "replicate": int(replicate)}


def write_csv(path, columns, records):
    with path.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(records)


def main():
    if not DATA_FILE.is_file():
        raise FileNotFoundError(f"Missing matrix: {DATA_FILE}\nRun 01_inspect_gse193258.py first.")
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA_CHECK.mkdir(parents=True, exist_ok=True)

    with gzip.open(DATA_FILE, "rt", encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file, delimiter="\t")
        header = next(reader)
        if header[0] != "gene":
            raise ValueError(f"Unexpected first column: {header[0]}")
        metadata = [sample_info(name) for name in header[1:]]
        if len(metadata) != 60 or len({item["sample"] for item in metadata}) != 60:
            raise ValueError("Expected 60 distinct GSE193258 sample columns.")

        group_counts = Counter((m["cell_line"], m["condition"]) for m in metadata)
        expected_conditions = ("control", "DTP", "acute", "short_washout", "long_washout")
        for cell_line in CELL_LINES:
            for condition in expected_conditions:
                if group_counts[cell_line, condition] != 3:
                    raise ValueError(f"Check metadata: {cell_line} / {condition} has "
                                     f"{group_counts[cell_line, condition]} samples, not 3.")
        write_csv(DATA_CHECK / "GSE193258_sample_metadata.csv",
                  ["sample", "cell_line", "condition", "replicate"], metadata)

        selected = [i for i, m in enumerate(metadata)
                    if m["condition"] in ("control", "DTP")]
        selected_info = [metadata[i] for i in selected]
        genes = []
        values = []
        for row in reader:
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"Malformed expression row for gene {row[0]}")
            genes.append(row[0])
            values.append([float(row[i + 1]) for i in selected])

    X = np.asarray(values, dtype=float).T   # samples x genes
    y = np.array([int(m["condition"] == "DTP") for m in selected_info])
    groups = np.array([m["cell_line"] for m in selected_info])
    if X.shape[0] != 24 or not np.isfinite(X).all():
        raise ValueError(f"Unexpected matrix shape or invalid entries: {X.shape}")

    print(f"Matrix: {X.shape[0]} selected samples x {X.shape[1]} genes")
    print("Groups checked: 4 cell lines x 5 conditions x 3 replicates")
    print("Baseline samples: 12 control, 12 DTP")

    fold_results = []
    predictions = []
    for cell_line in CELL_LINES:
        train = groups != cell_line
        test = groups == cell_line
        # All transformations are fitted on the training cell lines only.
        model = make_pipeline(
            SimpleImputer(strategy="median"),
            VarianceThreshold(threshold=1e-10),
            SelectKBest(score_func=f_classif, k=50),
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42),
        )
        model.fit(X[train], y[train])
        pred = model.predict(X[test])
        prob = model.predict_proba(X[test])[:, 1]
        score = {
            "held_out_cell_line": cell_line,
            "train_n": int(train.sum()),
            "test_n": int(test.sum()),
            "accuracy": round(float(accuracy_score(y[test], pred)), 4),
            "balanced_accuracy": round(float(balanced_accuracy_score(y[test], pred)), 4),
            "roc_auc": round(float(roc_auc_score(y[test], prob)), 4),
        }
        fold_results.append(score)
        for meta, true, p, probability in zip(
            [m for m, chosen in zip(selected_info, test) if chosen],
            y[test], pred, prob,
        ):
            predictions.append({"sample": meta["sample"],
                                "held_out_cell_line": cell_line,
                                "true_label": int(true), "predicted_label": int(p),
                                "probability_DTP": round(float(probability), 6)})
        print(f"{cell_line}: balanced accuracy={score['balanced_accuracy']:.3f}, "
              f"ROC-AUC={score['roc_auc']:.3f}")

    write_csv(RESULTS / "GSE193258_baseline_folds.csv",
              list(fold_results[0]), fold_results)
    write_csv(RESULTS / "GSE193258_baseline_predictions.csv",
              list(predictions[0]), predictions)

    actual = np.array([p["true_label"] for p in predictions])
    predicted = np.array([p["predicted_label"] for p in predictions])
    display = ConfusionMatrixDisplay.from_predictions(
        actual, predicted, labels=[0, 1], display_labels=["Control", "DTP"],
        values_format="d", colorbar=False,
    )
    display.ax_.set_title("Leave-one-cell-line-out baseline (24 samples)")
    display.figure_.tight_layout()
    display.figure_.savefig(FIGURES / "GSE193258_baseline_confusion_matrix.png", dpi=180)
    plt.close(display.figure_)

    print("Mean balanced accuracy:",
          round(float(np.mean([r["balanced_accuracy"] for r in fold_results])), 3))
    print("Outputs: results/GSE193258/data_check/GSE193258_sample_metadata.csv,")
    print("         results/GSE193258/dtp_vs_control/GSE193258_baseline_folds.csv,")
    print("         results/GSE193258/dtp_vs_control/GSE193258_baseline_predictions.csv,")
    print("         figures/GSE193258/dtp_vs_control/GSE193258_baseline_confusion_matrix.png")
    print("Note: within-study cell-line holdout, NOT cross-study evaluation.")


if __name__ == "__main__":
    main()
