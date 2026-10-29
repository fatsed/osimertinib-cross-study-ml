"""Inspect PC9 RNA-seq samples in GEO series GSE255958.

Reads the submitter's RSEM TPM files; this script does not train a model.
"""

import csv
import gzip
import hashlib
import io
import json
import tarfile
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "results" / "GSE255958" / "data_check"
ARCHIVE = RAW_DIR / "GSE255958_RAW.tar"
SOURCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/"
    "GSE255958/suppl/GSE255958_RAW.tar"
)

# Sample titles and accessions are taken from the official GEO series record.
GROUPS = [
    ("Control", "PC9", "DMSO, 48 h", range(8083374, 8083377)),
    ("Acute", "PC9", "Osimertinib, 48 h", range(8083377, 8083380)),
    ("DTP", "PC9", "Osimertinib, 9 days", range(8083380, 8083383)),
    ("Acquired_resistant", "PC9-R", "Osimertinib, 48 h", range(8083383, 8083386)),
]


def download_archive():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if ARCHIVE.exists() and ARCHIVE.stat().st_size:
        print("Using existing archive:", ARCHIVE)
        return
    temp = ARCHIVE.with_suffix(".tar.part")
    print("Downloading the official GEO archive (about 50 MB)...")
    print(SOURCE_URL)
    request = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temp.open("wb") as f:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
        temp.replace(ARCHIVE)
    except Exception as exc:
        temp.unlink(missing_ok=True)
        raise RuntimeError(
            f"Download failed: {exc}\n"
            "Download GSE255958_RAW.tar from the GEO page and place it in data/raw/.\n"
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255958"
        ) from exc


def sample_definitions():
    rows = []
    for state, line, treatment, ids in GROUPS:
        for rep, accession in enumerate(ids, start=1):
            rows.append({
                "sample_id": f"GSM{accession}",
                "cell_line": line,
                "response_state": state,
                "treatment": treatment,
                "replicate": rep,
            })
    return rows


def read_tpm(archive, member):
    fileobj = archive.extractfile(member)
    if fileobj is None:
        raise ValueError(f"Cannot open {member.name}")
    with fileobj, gzip.GzipFile(fileobj=fileobj) as unpacked:
        with io.TextIOWrapper(unpacked, encoding="utf-8-sig") as text:
            table = csv.DictReader(text, delimiter="\t")
            if not table.fieldnames or not {"gene_id", "TPM"}.issubset(table.fieldnames):
                raise ValueError(f"Missing gene_id or TPM columns in {member.name}: {table.fieldnames}")
            values = {}
            for row in table:
                gene = row["gene_id"].strip()
                if not gene:
                    continue
                if gene in values:
                    raise ValueError(f"Repeated gene_id {gene} in {member.name}")
                values[gene] = float(row["TPM"])
            if not values:
                raise ValueError(f"No genes read from {member.name}")
            return values


def write_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    download_archive()
    samples = sample_definitions()
    by_id = {s["sample_id"]: s for s in samples}
    expression = {}

    try:
        with tarfile.open(ARCHIVE, "r:") as archive:
            members = {}
            for member in archive.getmembers():
                if not member.isfile() or not member.name.endswith(".genes.results.gz"):
                    continue
                for sample_id in by_id:
                    if Path(member.name).name.startswith(sample_id + "_"):
                        if sample_id in members:
                            raise ValueError(f"Duplicate file for {sample_id}")
                        members[sample_id] = member
                        break
            missing = sorted(set(by_id) - set(members))
            if missing:
                raise ValueError("Missing sample files in archive: " + ", ".join(missing))
            for sample in samples:
                sid = sample["sample_id"]
                member = members[sid]
                sample["source_filename"] = Path(member.name).name
                expression[sid] = read_tpm(archive, member)
                print(sid, sample["response_state"], len(expression[sid]), "genes")
    except (tarfile.TarError, EOFError, OSError) as exc:
        raise RuntimeError(
            f"The archive may be incomplete or damaged: {exc}\n"
            f"Check or re-download: {ARCHIVE}"
        ) from exc

    shared_genes = set.intersection(*(set(v) for v in expression.values()))
    first_id = samples[0]["sample_id"]
    ordered_genes = [gene for gene in expression[first_id] if gene in shared_genes]
    if len(ordered_genes) < 1000:
        raise ValueError("Too few shared genes across selected samples")

    ids = [s["sample_id"] for s in samples]
    metadata_columns = ["sample_id", "cell_line", "response_state", "treatment", "replicate", "source_filename"]
    write_csv(OUTPUT_DIR / "PC9_sample_metadata.csv", metadata_columns,
              ([s[col] for col in metadata_columns] for s in samples))

    # Intermediate processed data, excluded from Git. Gene identifiers are left unchanged.
    matrix_path = PROCESSED_DIR / "GSE255958_PC9_TPM.csv"
    write_csv(matrix_path, ["gene_id"] + ids,
              ([gene] + [expression[sid][gene] for sid in ids] for gene in ordered_genes))

    # Five actual samples and five actual genes, suitable for a small Table 2 preview.
    preview_genes = ordered_genes[:5]
    write_csv(OUTPUT_DIR / "PC9_expression_preview.csv", ["sample_id", "response_state"] + preview_genes,
              ([s["sample_id"], s["response_state"]] +
               [expression[s["sample_id"]][gene] for gene in preview_genes] for s in samples[:5]))

    digest = hashlib.sha256()
    with ARCHIVE.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    report = {
        "dataset": "GSE255958",
        "source_url": SOURCE_URL,
        "archive_sha256": digest.hexdigest(),
        "total_samples_in_geo_series": 82,
        "selected_pc9_samples": len(samples),
        "group_counts": dict(Counter(s["response_state"] for s in samples)),
        "common_gene_ids": len(ordered_genes),
        "data_type": "submitter-provided RSEM TPM (not raw counts)",
        "processed_matrix": "data/processed/GSE255958_PC9_TPM.csv",
        "limitations": [
            "Only PC9 samples from one GEO study are included.",
            "PC9-R is a separately derived resistant population; treatment histories differ.",
            "Three replicates per response state; this is not independent cross-study validation.",
            "Gene IDs and TPM scale must be harmonized before combining with GSE193258 log2TPM.",
        ],
    }
    (OUTPUT_DIR / "PC9_inspection_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\nSUCCESS: data inspection completed (no model trained)")
    print("Samples:", len(samples), "| Shared gene IDs:", len(ordered_genes))
    print("Metadata, preview, and report in:", OUTPUT_DIR.relative_to(ROOT))
    print("Matrix in:", matrix_path.relative_to(ROOT))


if __name__ == "__main__":
    main()
