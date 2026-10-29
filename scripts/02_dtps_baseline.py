"""Download and audit the public GSE193258 expression matrix (standard library only).

This script does NOT train a model or label acquired-resistant cells.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
RESULTS_DIR = ROOT / "results" / "GSE193258" / "data_check"
FILENAME = "GSE193258_RNAseq_log2TPM_abundance.tsv.gz"
SOURCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193258/suppl/"
    + FILENAME
)
SOURCE_PAGE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193258"


def download_if_needed(path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        print(f"Using existing file: {path}")
        return
    print(f"Downloading public processed expression file from:\n  {SOURCE_URL}")
    temp_path = path.with_name(path.name + ".part")
    request = urllib.request.Request(
        SOURCE_URL, headers={"User-Agent": "ResearchDataAudit/1.0 (academic)"}
    )
    try:
        with urllib.request.urlopen(request, timeout=90) as response, temp_path.open("wb") as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
        temp_path.replace(path)
    except (OSError, urllib.error.URLError) as exc:
        temp_path.unlink(missing_ok=True)
        print("\nDownload failed:", exc)
        print("Open the official GEO page and download the log2TPM file manually:")
        print("  " + SOURCE_PAGE)
        print(f"Save the original file as:\n  {path}")
        sys.exit(1)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_matrix(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 3:
            raise ValueError("Unexpected header: expected gene ID and multiple sample columns")
        # The first column usually contains gene identifiers. Verify it manually before analysis.
        sample_names = header[1:]
        first_five_gene_rows = []
        total_gene_rows = 0
        first_col_examples = []
        for row in reader:
            if not row or len(row) < 2 or row[0].startswith("#"):
                continue
            total_gene_rows += 1
            if len(first_five_gene_rows) < 5:
                first_five_gene_rows.append(row)
                first_col_examples.append(row[0])
    if len(first_five_gene_rows) < 5:
        raise ValueError("Not enough gene-expression rows to create a preview")

    # Five samples by five features: a format illustration only.
    preview_path = RESULTS_DIR / "GSE193258_table2_preview.csv"
    gene_ids = [r[0] for r in first_five_gene_rows]
    selected_sample_names = sample_names[:5]
    with preview_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample_column"] + gene_ids)
        for sample_idx, name in enumerate(selected_sample_names):
            values = []
            for gene_row in first_five_gene_rows:
                values.append(gene_row[sample_idx + 1] if len(gene_row) > sample_idx + 1 else "")
            writer.writerow([name] + values)

    return {
        "geo_accession": "GSE193258",
        "source_page": SOURCE_PAGE,
        "source_file": SOURCE_URL,
        "local_file": str(path.relative_to(ROOT)),
        "sha256": sha256(path),
        "expression_type": "submitter-provided log2TPM abundance; not raw counts",
        "number_of_expression_rows": total_gene_rows,
        "number_of_matrix_columns_excluding_first_id_column": len(sample_names),
        "first_column_header": header[0],
        "first_five_first_column_values": first_col_examples,
        "first_five_sample_column_headers": selected_sample_names,
        "preview_file": str(preview_path.relative_to(ROOT)),
        "limitations": [
            "Sample columns and row identifiers must be validated against GEO metadata.",
            "This matrix alone is not evidence of acquired-resistant labels.",
            "The preview is not a classification result.",
            "Different studies may use incompatible gene IDs and expression units.",
        ],
    }


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / FILENAME
    download_if_needed(path)
    try:
        report = inspect_matrix(path)
    except (OSError, UnicodeError, ValueError, EOFError) as exc:
        print("Unable to inspect downloaded file:", exc)
        print("Verify the file format against the official GEO record:", SOURCE_PAGE)
        sys.exit(1)
    report_path = RESULTS_DIR / "GSE193258_dataset_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nSUCCESS: inspection completed (not model training)")
    print("Gene-expression rows:", report["number_of_expression_rows"])
    print("Columns after gene-ID column:", report["number_of_matrix_columns_excluding_first_id_column"])
    print("Example sample columns:", report["first_five_sample_column_headers"])
    print("Preview:", report["preview_file"])
    print("Report:", report_path.relative_to(ROOT))


if __name__ == "__main__":
    main()
