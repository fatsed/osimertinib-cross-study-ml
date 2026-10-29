"""Compare osimertinib DTP samples with acute-treatment samples in GSE193258.

Leave one cell line out for testing. All preprocessing and gene selection
are fitted on training samples. This is a within-study analysis.
"""

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (ConfusionMatrixDisplay, balanced_accuracy_score,
                             roc_auc_score)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
EXPRESSION_FILE = ROOT / "data" / "raw" / "GSE193258_RNAseq_log2TPM_abundance.tsv.gz"
RESULTS_DIR = ROOT / "results" / "GSE193258" / "dtp_vs_acute"
FIGURES_DIR = ROOT / "figures" / "GSE193258" / "dtp_vs_acute"
CELL_LINES = ("H1975", "HCC2935", "HCC827", "PC9")


def sample_metadata(name):
    match = re.fullmatch(r"(H1975|HCC2935|HCC827|PC9)_osi_(DTP|acute)_([123])", name)
    if not match:
        return None
    cell_line, condition, replicate = match.groups()
    return {"sample": name, "cell_line": cell_line,
            "condition": condition, "replicate": int(replicate)}


def load_data():
    if not EXPRESSION_FILE.is_file():
        raise FileNotFoundError(f"Missing {EXPRESSION_FILE}. Run script 01 first.")

    with gzip.open(EXPRESSION_FILE, "rt", encoding="utf-8-sig", newline="") as file:
        rows = csv.reader(file, delimiter="\t")
        header = next(rows)
        if header[0] != "gene":
            raise ValueError("Unexpected expression-matrix format")

        selected = [(i, sample_metadata(name)) for i, name in enumerate(header[1:])]
        selected = [(i, item) for i, item in selected if item is not None]
        samples = [item for _, item in selected]
        indices = [i + 1 for i, _ in selected]
        if len(samples) != 24 or len({s["sample"] for s in samples}) != 24:
            raise ValueError("Expected 24 distinct acute/DTP samples")

        for cell_line in CELL_LINES:
            for condition in ("acute", "DTP"):
                if sum(s["cell_line"] == cell_line and s["condition"] == condition
                       for s in samples) != 3:
                    raise ValueError(f"Check samples: {cell_line} {condition}")

        values = []
        genes = []
        for row in rows:
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(f"Wrong column count for {row[0]}")
            genes.append(row[0])
            values.append([float(row[i]) for i in indices])

    X = np.asarray(values, dtype=float).T
    if not np.isfinite(X).all():
        raise ValueError("Missing or non-finite expression values")
    y = np.asarray([int(s["condition"] == "DTP") for s in samples])
    groups = np.asarray([s["cell_line"] for s in samples])
    return X, y, groups, samples, genes


def build_model():
    return make_pipeline(
        VarianceThreshold(threshold=1e-10),
        SelectKBest(score_func=f_classif, k=50),
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42),
    )


def grouped_evaluation(X, y, groups, samples=None):
    folds = []
    predictions = []
    for held_out in CELL_LINES:
        train = groups != held_out
        test = groups == held_out
        model = build_model()
        model.fit(X[train], y[train])
        pred = model.predict(X[test])
        prob = model.predict_proba(X[test])[:, 1]
        folds.append({
            "held_out_cell_line": held_out,
            "train_n": int(train.sum()), "test_n": int(test.sum()),
            "balanced_accuracy": round(float(balanced_accuracy_score(y[test], pred)), 4),
            "roc_auc": round(float(roc_auc_score(y[test], prob)), 4),
        })
        if samples is not None:
            test_samples = [s for s, is_test in zip(samples, test) if is_test]
            for s, real, predicted, probability in zip(test_samples, y[test], pred, prob):
                predictions.append({
                    "sample": s["sample"], "cell_line": held_out,
                    "true_label": int(real), "predicted_label": int(predicted),
                    "probability_DTP": round(float(probability), 6),
                })
    return folds, predictions


def shuffled_labels(y, groups, rng):
    # Keep three acute and three DTP labels per cell line.
    result = y.copy()
    for cell_line in CELL_LINES:
        positions = np.flatnonzero(groups == cell_line)
        result[positions] = rng.permutation(y[positions])
    return result


def save_csv(path, records):
    with path.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=30,
                        help="Number of within-cell-line label-shuffle checks (default: 30)")
    args = parser.parse_args()
    if args.permutations < 1:
        parser.error("--permutations must be at least 1")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    X, y, groups, samples, genes = load_data()
    folds, predictions = grouped_evaluation(X, y, groups, samples)
    observed = float(np.mean([f["balanced_accuracy"] for f in folds]))
    print(f"Samples: {len(samples)}, genes: {len(genes)}")
    print("Comparison: DTP (21 days) vs acute (24 hours), both osimertinib-treated")
    for f in folds:
        print(f"{f['held_out_cell_line']}: balanced accuracy={f['balanced_accuracy']:.3f}, "
              f"ROC-AUC={f['roc_auc']:.3f}")
    print(f"Mean balanced accuracy: {observed:.3f}")

    rng = np.random.default_rng(42)
    null_scores = []
    for i in range(args.permutations):
        y_shuffled = shuffled_labels(y, groups, rng)
        shuffled_folds, _ = grouped_evaluation(X, y_shuffled, groups)
        null_scores.append(float(np.mean([f["balanced_accuracy"] for f in shuffled_folds])))
    p_value = (1 + sum(score >= observed for score in null_scores)) / (args.permutations + 1)
    print(f"Label-shuffle check: {args.permutations} repeats, "
          f"mean null score={np.mean(null_scores):.3f}, empirical p={p_value:.3f}")

    prefix = "GSE193258_DTP_vs_acute"
    save_csv(RESULTS_DIR / f"{prefix}_folds.csv", folds)
    save_csv(RESULTS_DIR / f"{prefix}_predictions.csv", predictions)
    report = {
        "dataset": "GSE193258",
        "comparison": "DTP vs acute osimertinib treatment",
        "total_samples": len(samples), "genes": len(genes),
        "train_test_design": "leave-one-cell-line-out, same study",
        "mean_balanced_accuracy": observed,
        "number_of_label_permutations": args.permutations,
        "mean_permuted_balanced_accuracy": float(np.mean(null_scores)),
        "permutation_p_value": float(p_value),
        "interpretation_limits": [
            "This does not assess acquired resistance.",
            "The groups differ in treatment duration (21 days vs 24 hours).",
            "All samples originate from one study; this is not cross-study validation.",
            "With six test samples per fold, performance estimates are imprecise.",
        ],
    }
    (RESULTS_DIR / f"{prefix}_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    true = [p["true_label"] for p in predictions]
    predicted = [p["predicted_label"] for p in predictions]
    display = ConfusionMatrixDisplay.from_predictions(
        true, predicted, labels=[0, 1], display_labels=["Acute", "DTP"],
        values_format="d", colorbar=False,
    )
    display.ax_.set_title("DTP vs acute: cell-line-held-out predictions")
    display.figure_.tight_layout()
    display.figure_.savefig(FIGURES_DIR / f"{prefix}_confusion_matrix.png", dpi=180)
    plt.close(display.figure_)
    print(f"Saved report, folds and predictions to {RESULTS_DIR.relative_to(ROOT)}/{prefix}_*")
    print(f"Saved confusion matrix to {FIGURES_DIR.relative_to(ROOT)}/{prefix}_confusion_matrix.png")
    print("NOTE: This experiment does not test DTP versus acquired resistance.")


if __name__ == "__main__":
    main()
