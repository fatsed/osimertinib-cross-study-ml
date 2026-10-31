# GSE285632 — Uploaded supplementary workbook audit

**File:** `GSE285632_All.DEG.Expression.xlsx` (user-supplied).  
**GEO record:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285632  
**Audit date:** 2026-10-30

## Verified file contents

| Worksheet | Feature rows (excluding title and headers) | Sample-value columns | Additional contents |
|---|---:|---|---|
| `mRNA` | 57,764 unique nonempty gene IDs | `P1`, `P2`, `P3`, `D1`, `D2`, `D3` | No resistant-labeled columns or differential-statistic columns in this worksheet |
| `circRNA` | 2,997 unique nonempty circRNA IDs | `D1`, `D2`, `D3`, `O1`, `O2`, `O3` | `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `diffState` |

All 346,584 mRNA numeric entries (57,764 genes x 6 samples) were present. No duplicate or empty gene IDs were found in the mRNA worksheet. The mRNA values are exactly consistent (within floating-point rounding) with a base-10 transformation of nonnegative integer values: `value = log10(integer + 1)`. **This identifies the numerical transform, not the provenance or normalization of the original integers; do not assert they are raw counts without verification.** GEO sample pages describe some processed data as FPKM, which makes a units/provenance cross-check necessary before scientific analysis.

The `circRNA` worksheet includes differential-statistic results, and its `diffState` labels should not be interpreted as per-sample treatment-state labels. Of 2,997 circRNA entries, 146 are marked Up, 45 Down and 2,806 Not significant; only **4** rows have adjusted p-value < 0.05. This disagreement means the threshold defining `diffState` must be established from the original pipeline before using this classification.

## Mapping and feasibility caveats

- The worksheet headings alone do not provide an explicit one-to-one `P1`/`D1` to GEO sample accession mapping. `P` likely denotes Parental and `D` likely denotes DTP, but confirm this rather than assume.
- The circRNA `D` and `O` codes must also be mapped against the source metadata before use. The identical letter `D` in two different worksheets does not prove identical samples.
- The parent GEO study lists **24** samples across eight groups: DMSO, Osi, Parental-EV, Parental-OC, DTP-EV, DTP-OC, Resistant-EV, Resistant-OC, with three replicates per group. **Only 6 per-sample columns per worksheet occur in this uploaded workbook**. This does not constitute a complete expression matrix for the 24 samples.
- **No directly labeled acquired-resistant sample expression columns appear in the uploaded mRNA worksheet.** Do not use it yet as a standalone DTP-vs-acquired-resistant modeling dataset.
- Selection/filtering and any batch correction must be fitted on training data only. The circRNA differential-test statistics must not be used as precomputed predictive features in a cross-validation pipeline.

## Next steps

1. Inspect the second GEO supplementary file `GSE285632_differentially-expressed.xlsx` for sample-level expression columns for resistant-EV and DTP-EV groups.
2. Obtain and verify exact column-to-GSM IDs and expression-value units from the source.
3. If the supplementary data still lack a full matrix, evaluate whether source sequencing reads from SRA can be processed into a consistent count matrix (a substantially larger task), or use a different appropriately labeled public dataset.
4. Preserve raw original files and audit decisions; do not merge different scales/experimental arms blindly.

**Conclusion:** This uploaded Excel workbook is useful for an initial P-versus-D style investigation and candidate-gene checks, but it does **not yet** support the intended cross-study DTP-versus-acquired-resistance classifier. No independent validation or new model was run in this audit.
