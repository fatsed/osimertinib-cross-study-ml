# Research implementation plan

**Working title:** Confounding-Aware Cross-Study Evaluation of Osimertinib Response States

## Research question
Can transcriptomic models distinguish osimertinib drug-tolerant persister (DTP) states from acquired-resistant states in independent studies without relying mostly on study or cell-line identity?

## Status
Planning and data audit only. No classification performance has been established.

## Implementation stages
1. **Source audit:** register GEO accession, source publication, sample IDs, cell line, exposure, response label, experimental replicate and file format. Confirm labels from original metadata.
2. **Data integrity:** download into `data/raw/` (excluded from Git), record provenance and checksums; check gene identifiers, units, missingness, sample counts and duplicate samples.
3. **Feasibility assessment:** cross-tabulate study × cell line × response state. If state and origin are perfectly confounded, do not claim independent cross-study classification. Consider a revised estimand or dataset collection.
4. **First baseline:** only if samples and labels support it, build a carefully separated baseline (e.g., regularized logistic regression); fit feature selection and scaling on training folds only.
5. **Robust evaluation:** compare random splits with study-held-out or grouped splits only where scientifically feasible; report uncertainty and baselines. Avoid leakage and pseudo-replication.
6. **Confounding checks:** study-origin predictions, label imbalance, cell-line stratification, negative controls and sensitivity analyses.
7. **Scientific reporting:** generate genuine Table 2 data preview, flowchart, run logs, evaluation tables and reproducible figures; then draft Methodology and Results.

## Candidate GEO datasets
| GEO accession | Preliminary relevance | Notes |
| --- | --- | --- |
| [GSE193258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193258) | Osimertinib tolerance, withdrawal and controls; bulk RNA-seq | 60 samples. Do **not** assign acquired-resistant labels to washout states. Initial inspection dataset. |
| [GSE255958](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255958) | Multiple treatment states, including DTP/resistant in selected cell-line settings | Requires sample-level audit before use. |
| [GSE285632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285632) | HCC827 parental, DTP and resistant, with experimental perturbation groups | 24 samples; EV vs OC must not be mixed blindly. |

## Data interpretation
The `log2TPM` values in GSE193258 are normalized expression measurements, **not** raw read counts. Do not feed log2TPM into DESeq2 as raw counts. Use the appropriately documented matrix for each downstream step.

## First deliverable
Run `scripts/01_inspect_gse193258.py`, inspect the output, and confirm what rows, columns and treatment labels mean **before** training a model.
