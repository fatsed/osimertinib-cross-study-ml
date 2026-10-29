"""Explore PC9 response states and run a small DTP-versus-resistant pilot.

Dataset: GSE255958. This is a within-study check, not external validation.
"""

import json
import warnings
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import ConfusionMatrixDisplay, accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "processed" / "GSE255958_PC9_TPM.csv"
RESULTS = ROOT / "results" / "GSE255958" / "dtp_vs_resistant"
FIGURES = ROOT / "figures" / "GSE255958"

# Sample labels taken from the GSE255958 GEO record, not inferred from expression values.
SAMPLES = {
    **{f"GSM808337{i}": "Control" for i in (4, 5, 6)},
    **{f"GSM808337{i}": "Acute" for i in (7, 8, 9)},
    **{f"GSM808338{i}": "DTP" for i in (0, 1, 2)},
    **{f"GSM808338{i}": "Acquired_resistant" for i in (3, 4, 5)},
}


def load_data():
    if not INPUT.is_file():
        raise FileNotFoundError(f"Run script 04 first. Missing: {INPUT}")
    df = pd.read_csv(INPUT)
    if "gene_id" not in df.columns:
        raise ValueError("The expression matrix needs a gene_id column")
    missing = [sample for sample in SAMPLES if sample not in df.columns]
    if missing:
        raise ValueError(f"Missing sample columns: {missing}")
    if df["gene_id"].isna().any() or df["gene_id"].duplicated().any():
        raise ValueError("Missing or duplicate gene identifiers: inspect before modeling")
    expression = df[list(SAMPLES)].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    if not np.isfinite(expression).all() or (expression < 0).any():
        raise ValueError("TPM should have finite, nonnegative values")
    # Genes are rows in the input, but samples must be rows for scikit-learn.
    X = np.log2(1 + expression.T)
    return X, df["gene_id"].astype(str).to_numpy()


def plot_pca(X, names):
    folder = FIGURES / "overview"
    folder.mkdir(parents=True, exist_ok=True)
    # This unsupervised figure uses all 12 samples, not the classification folds.
    variance = X.var(axis=0)
    top = np.argsort(variance)[-min(500, X.shape[1]):]
    model = PCA(n_components=2)
    coords = model.fit_transform(X[:, top])
    fig, ax = plt.subplots(figsize=(9, 6))
    # Offset the three replicate labels around each cluster so they do not overlap.
    label_offsets = [(-30, 18), (30, 18), (0, -27)]
    for state in dict.fromkeys(SAMPLES.values()):
        indices = [i for i, name in enumerate(names) if SAMPLES[name] == state]
        ax.scatter(coords[indices, 0], coords[indices, 1], label=state.replace("_", " "), s=65, zorder=3)
        for replicate, i in enumerate(indices):
            dx, dy = label_offsets[replicate]
            ax.annotate(
                names[i][-2:],
                xy=(coords[i, 0], coords[i, 1]),
                xytext=(dx, dy),
                textcoords="offset points",
                ha="right" if dx < 0 else ("left" if dx > 0 else "center"),
                va="bottom" if dy > 0 else "top",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.9),
                arrowprops=dict(arrowstyle="-", color="0.55", linewidth=0.7, shrinkA=2, shrinkB=4),
            )
    ax.margins(x=0.22, y=0.22)
    ax.set_xlabel(f"PC1 ({model.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({model.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("GSE255958: PC9 treatment states (exploratory PCA)")
    ax.legend(frameon=False)
    fig.tight_layout()
    path = folder / "PC9_treatment_states_PCA.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def predict_out_of_fold(X, y):
    predicted = np.zeros(len(y), dtype=int)
    probabilities = np.zeros(len(y), dtype=float)
    for train, test in LeaveOneOut().split(X):
        pipeline = make_pipeline(
            VarianceThreshold(threshold=1e-10),
            SelectKBest(score_func=f_classif, k=10),
            StandardScaler(),
            LogisticRegression(C=0.1, class_weight="balanced", max_iter=1500),
        )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            warnings.filterwarnings("ignore", category=UserWarning, module="sklearn.feature_selection")
            pipeline.fit(X[train], y[train])
        predicted[test] = pipeline.predict(X[test])
        probabilities[test] = pipeline.predict_proba(X[test])[:, 1]
    return predicted, probabilities


def main():
    X_all, genes = load_data()
    sample_names = list(SAMPLES)
    print(f"Loaded {len(sample_names)} PC9 samples and {len(genes)} gene IDs")
    print("Groups:", pd.Series(list(SAMPLES.values())).value_counts().to_dict())
    pca_path = plot_pca(X_all, sample_names)
    print("PCA:", pca_path.relative_to(ROOT))

    # PC9 DTP samples and PC9-R resistant samples, three biological replicates each.
    selected = [i for i, name in enumerate(sample_names)
                if SAMPLES[name] in ("DTP", "Acquired_resistant")]
    names = [sample_names[i] for i in selected]
    X = X_all[selected]
    y = np.array([int(SAMPLES[name] == "Acquired_resistant") for name in names])
    predictions, probabilities = predict_out_of_fold(X, y)
    accuracy = accuracy_score(y, predictions)
    balanced = balanced_accuracy_score(y, predictions)
    auc = roc_auc_score(y, probabilities)

    # Enumerate all 20 possible ways to label 3 of the 6 samples resistant.
    # At this sample size, permutation resolution is limited and labels may not be exchangeable.
    null_scores = []
    for positives in combinations(range(len(y)), 3):
        shuffled = np.zeros(len(y), dtype=int)
        shuffled[list(positives)] = 1
        perm_predictions, _ = predict_out_of_fold(X, shuffled)
        null_scores.append(balanced_accuracy_score(shuffled, perm_predictions))
    exact_p = sum(score >= balanced - 1e-10 for score in null_scores) / len(null_scores)

    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "sample_id": names,
        "state": [SAMPLES[name] for name in names],
        "true_label_resistant": y,
        "predicted_label_resistant": predictions,
        "probability_resistant": np.round(probabilities, 6),
    }).to_csv(RESULTS / "PC9_DTP_vs_resistant_predictions.csv", index=False)

    fig, ax = plt.subplots(figsize=(5, 4.5))
    ConfusionMatrixDisplay.from_predictions(
        y, predictions, labels=[0, 1], display_labels=["DTP", "Resistant"],
        ax=ax, colorbar=False, values_format="d")
    ax.set_title("PC9 pilot: leave-one-sample-out")
    fig.tight_layout()
    cm_path = RESULTS / "PC9_DTP_vs_resistant_confusion_matrix.png"
    fig.savefig(cm_path, dpi=180)
    plt.close(fig)

    report = {
        "dataset": "GSE255958",
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255958",
        "expression_unit": "RSEM TPM transformed with log2(TPM + 1)",
        "total_PC9_samples": len(sample_names),
        "genes": len(genes),
        "comparison": "PC9 DTP 9 days vs PC9-R acquired resistant treated 48 hours",
        "comparison_samples": len(y),
        "validation": "leave-one-sample-out; all preprocessing and gene selection within each training fold",
        "model": "10 selected genes + scaled regularized logistic regression (pilot only)",
        "accuracy": float(accuracy),
        "balanced_accuracy": float(balanced),
        "roc_auc_out_of_fold": float(auc),
        "exact_label_permutations": len(null_scores),
        "exploratory_permutation_p": float(exact_p),
        "limitations": [
            "Only three replicates per treatment state, all from one study.",
            "PC9 persisters and PC9-R resistant cells have different selection histories.",
            "A model can learn cell-line derivative or treatment-context differences.",
            "Leave-one-sample-out is not external validation; results are highly unstable.",
            "The permutation comparison assumes sample labels are exchangeable, which may not hold.",
            "This cannot establish cross-study generalization or a validated biomarker.",
        ],
    }
    (RESULTS / "PC9_DTP_vs_resistant_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")

    print("\nPILOT ONLY — 3 DTP vs 3 acquired-resistant samples")
    print(f"Leave-one-out accuracy: {accuracy:.3f}")
    print(f"Balanced accuracy: {balanced:.3f}")
    print(f"Out-of-fold ROC-AUC: {auc:.3f}")
    print(f"Exploratory label permutation: p={exact_p:.3f} ({len(null_scores)} assignments)")
    print("Saved results to:", RESULTS.relative_to(ROOT))
    print("Important: this is NOT cross-study validation.")


if __name__ == "__main__":
    main()
