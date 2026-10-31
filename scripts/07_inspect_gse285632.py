"""Audit the sample-level table in GSE285632_differentially-expressed.xlsx.

This script does not download data, harmonize expression scales, train a model,
or assume a GEO sample ID for any of the abbreviated workbook columns.

Put the original workbook in data/raw and run:
    python scripts/07_inspect_gse285632.py

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285632
"""
import argparse
import csv
import hashlib
import json
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
XMLNS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
EXPECTED = [f"{group}-{rep}" for group in ("P", "P+C", "D", "D+C", "R", "R+C") for rep in (1, 2, 3)]


def read_shared_strings(archive):
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(element.itertext()) for element in root.findall(f"{XMLNS}si")]


def iter_rows(archive, shared):
    with archive.open("xl/worksheets/sheet1.xml") as stream:
        for _, elem in ET.iterparse(stream, events=("end",)):
            if elem.tag != f"{XMLNS}row":
                continue
            values = {}
            for cell in elem.findall(f"{XMLNS}c"):
                item = cell.find(f"{XMLNS}v")
                if item is None:
                    continue
                value = item.text
                if cell.attrib.get("t") == "s":
                    value = shared[int(value)]
                col = re.match(r"[A-Z]+", cell.attrib["r"]).group()
                values[col] = value
            yield int(elem.attrib["r"]), values
            elem.clear()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "data/raw/GSE285632_differentially-expressed.xlsx")
    parser.add_argument("--output", type=Path, default=ROOT / "results/GSE285632/data_check")
    parser.add_argument("--matrix", type=Path, default=ROOT / "data/processed/GSE285632_expression_as_provided.csv")
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Input workbook not found: {args.input}")
    args.output.mkdir(parents=True, exist_ok=True)
    args.matrix.parent.mkdir(parents=True, exist_ok=True)

    gene_counts = Counter()
    seen_profiles = defaultdict(set)
    sample_zero = Counter()
    sample_min, sample_max = {}, {}
    n_genes = n_values = n_noninteger_log = 0
    preview = []
    sha256 = hashlib.sha256(args.input.read_bytes()).hexdigest()
    col_letters = [chr(c) for c in range(ord("B"), ord("S") + 1)]

    with zipfile.ZipFile(args.input) as archive, args.matrix.open("w", newline="", encoding="utf-8") as matrix_file:
        shared = read_shared_strings(archive)
        writer = csv.writer(matrix_file)
        writer.writerow(["source_excel_row", "gene_symbol"] + EXPECTED)
        for row_number, row in iter_rows(archive, shared):
            if row_number == 1:
                continue
            if row_number == 2:
                observed = [row.get(col) for col in col_letters]
                if row.get("A") != "Symbol" or observed != EXPECTED:
                    raise ValueError(f"Unexpected header. Got: {row.get('A')}, {observed}")
                sample_min = {label: math.inf for label in EXPECTED}
                sample_max = {label: -math.inf for label in EXPECTED}
                continue
            gene = row.get("A")
            if not gene:
                raise ValueError(f"Blank gene symbol in row {row_number}")
            raw = [row.get(col) for col in col_letters]
            if any(v is None for v in raw):
                raise ValueError(f"Missing sample value in row {row_number}")
            vals = [float(v) for v in raw]
            if not all(math.isfinite(v) and v >= 0 for v in vals):
                raise ValueError(f"Invalid expression value in row {row_number}")
            n_genes += 1
            gene_counts[gene] += 1
            seen_profiles[gene].add(tuple(vals))
            n_values += len(vals)
            for label, value in zip(EXPECTED, vals):
                if value == 0:
                    sample_zero[label] += 1
                sample_min[label] = min(sample_min[label], value)
                sample_max[label] = max(sample_max[label], value)
                raw_estimate = 10**value - 1
                if abs(raw_estimate - round(raw_estimate)) > 1e-5:
                    n_noninteger_log += 1
            writer.writerow([row_number, gene] + vals)
            if gene in {"A1BG", "A1CF", "A2M", "A2ML1", "A3GALT2"} and len(preview) < 5:
                preview.append([gene, vals[6], vals[7], vals[8], vals[12], vals[13], vals[14]])

    repeated = sorted([(gene, count, len(seen_profiles[gene])) for gene, count in gene_counts.items() if count > 1])
    with (args.output / "duplicate_gene_symbols.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["gene_symbol", "total_rows", "different_numeric_profiles"])
        writer.writerows(repeated)

    with (args.output / "expression_preview_D_vs_R.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["gene_symbol", "D-1", "D-2", "D-3", "R-1", "R-2", "R-3"])
        writer.writerows(preview)

    with (args.output / "sample_columns_unverified.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(["excel_sample", "provisional_state", "condition_marker", "geo_sample_id", "mapping_verified"])
        for label in EXPECTED:
            group = label.rsplit("-", 1)[0]
            state = {"P": "Parental", "D": "DTP", "R": "Resistant"}[group.replace("+C", "")]
            writer.writerow([label, state, "+C" if "+C" in group else "without +C", "", "NO"])

    report = {
        "dataset": "GSE285632", "input_filename": args.input.name,
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285632",
        "source_sha256": sha256,
        "sheet_name": "Sheet1", "expression_rows": n_genes,
        "distinct_gene_symbols": len(gene_counts), "repeated_gene_symbol_names": len(repeated),
        "extra_rows_due_to_repeated_symbols": sum(x-1 for x in gene_counts.values() if x > 1),
        "sample_columns": len(EXPECTED), "sample_names": EXPECTED,
        "finite_nonnegative_numeric_entries": n_values,
        "entries_not_compatible_with_log10_integer_plus_one": n_noninteger_log,
        "zero_fraction_by_sample": {s: round(sample_zero[s]/n_genes, 5) for s in EXPECTED},
        "numeric_range_by_sample": {s: [sample_min[s], sample_max[s]] for s in EXPECTED},
        "cautions": [
            "P, D and R are provisional interpretations; verify the exact mapping to GEO GSM samples.",
            "The +C condition is NOT assigned to EV/OC until provenance is confirmed.",
            "The values are mathematically compatible with log10(integer+1); original measurement units are unverified.",
            "Duplicate gene symbols may represent different underlying features; do not automatically average or discard them.",
            "All 18 samples belong to ONE GEO study; they are NOT an independent cross-study test.",
            "This script does not perform classification or differential expression analysis.",
        ],
    }
    (args.output / "workbook_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("GSE285632 WORKBOOK AUDIT COMPLETE (NO MODEL TRAINED)")
    print(f"Gene rows: {n_genes} | Unique gene symbols: {len(gene_counts)} | Samples: {len(EXPECTED)}")
    print(f"Nonmissing numeric values: {n_values}")
    print(f"Repeated symbols: {len(repeated)} (extra rows: {report['extra_rows_due_to_repeated_symbols']})")
    print(f"Unexpected log10(integer+1) values: {n_noninteger_log}")
    print("Audit:", args.output)
    print("As-provided expression matrix:", args.matrix)
    print("CAUTION: Sample-to-GSM mapping and expression units are not yet verified.")

if __name__ == "__main__":
    main()
