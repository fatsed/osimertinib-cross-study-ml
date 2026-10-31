"""Summarize the available osimertinib datasets before cross-study modeling.

Run from the repository root:
    python scripts/06_dataset_audit.py

The GSE285632 records below are transcribed from the official GEO sample list.
They are metadata only; expression data have not been downloaded or checked.
"""

import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "cross_study" / "data_audit"

GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
FIELDS = [
    "study", "sample_id", "cell_line", "response_state", "treatment",
    "perturbation", "replicate", "expression_available_locally", "source_url",
]


def read_rows(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path.relative_to(ROOT)}. Run the earlier inspection scripts first.")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def gse193258():
    path = ROOT / "results/GSE193258/data_check/GSE193258_sample_metadata.csv"
    rows = []
    for item in read_rows(path):
        state = {
            "control": "Control", "acute": "Acute", "DTP": "DTP",
            "short_washout": "Short_washout", "long_washout": "Long_washout",
        }.get(item["condition"])
        if state is None:
            raise ValueError(f"Unknown GSE193258 condition: {item['condition']}")
        rows.append({
            "study": "GSE193258", "sample_id": item["sample"],
            "cell_line": item["cell_line"], "response_state": state,
            "treatment": "DMSO" if state == "Control" else "Osimertinib / washout as applicable",
            "perturbation": "None", "replicate": item["replicate"],
            "expression_available_locally": "yes",
            "source_url": GEO + "GSE193258",
        })
    if len(rows) != 60:
        raise ValueError(f"Expected 60 GSE193258 records; found {len(rows)}")
    return rows


def gse255958():
    path = ROOT / "results/GSE255958/data_check/PC9_sample_metadata.csv"
    rows = []
    for item in read_rows(path):
        state = item["response_state"]
        if state not in {"Control", "Acute", "DTP", "Acquired_resistant"}:
            raise ValueError(f"Unknown GSE255958 state: {state}")
        rows.append({
            "study": "GSE255958", "sample_id": item["sample_id"],
            "cell_line": item["cell_line"], "response_state": state,
            "treatment": item["treatment"], "perturbation": "None",
            "replicate": item["replicate"], "expression_available_locally": "yes",
            "source_url": GEO + "GSE255958",
        })
    if len(rows) != 12:
        raise ValueError(f"Expected 12 selected GSE255958 PC9 samples; found {len(rows)}")
    return rows


def gse285632():
    # Each group has three GEO samples. EV and OC must not be pooled blindly.
    groups = [
        ("DMSO", "Control", "None", "DMSO"),
        ("Osi", "Osi_treated", "None", "Osimertinib; precise state to verify"),
        ("Parental-EV", "Parental", "EV", "Experimental context to verify"),
        ("Parental-OC", "Parental", "OC", "cATF6 overexpression"),
        ("DTP-EV", "DTP", "EV", "Osimertinib DTP"),
        ("DTP-OC", "DTP", "OC", "Osimertinib DTP + cATF6 overexpression"),
        ("Resistant-EV", "Acquired_resistant", "EV", "Resistant state"),
        ("Resistant-OC", "Acquired_resistant", "OC", "Resistant + cATF6 overexpression"),
    ]
    rows = []
    for group_idx, (group, state, perturbation, treatment) in enumerate(groups):
        for replicate in (1, 2, 3):
            gsm = 8706066 + group_idx * 3 + replicate - 1
            rows.append({
                "study": "GSE285632", "sample_id": f"GSM{gsm}",
                "cell_line": "HCC827", "response_state": state,
                "treatment": treatment, "perturbation": perturbation,
                "replicate": str(replicate),
                "expression_available_locally": "no - metadata only",
                "source_url": GEO + "GSE285632",
            })
    assert len(rows) == 24
    return rows


def write_table(path, fields, records):
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def main():
    rows = gse193258() + gse255958() + gse285632()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    write_table(OUTPUT / "sample_registry.csv", FIELDS, rows)

    grouped = Counter((r["study"], r["cell_line"], r["response_state"], r["perturbation"]) for r in rows)
    counts = [
        dict(study=study, cell_line=cell, response_state=state,
             perturbation=perturbation, samples=n)
        for (study, cell, state, perturbation), n in sorted(grouped.items())
    ]
    write_table(OUTPUT / "group_counts.csv",
                ["study", "cell_line", "response_state", "perturbation", "samples"], counts)

    by_study = Counter(r["study"] for r in rows)
    candidate = Counter(r["study"] for r in rows if
                        r["response_state"] in {"DTP", "Acquired_resistant"} and
                        (r["study"] != "GSE285632" or r["perturbation"] == "EV"))
    report = {
        "metadata_records_total": len(rows),
        "metadata_records_by_study": dict(by_study),
        "DTP_or_resistant_candidates_by_study": dict(candidate),
        "important_notes": [
            "GSE193258 has DTP samples but no confirmed acquired-resistant group; washout is not resistance.",
            "Only 12 PC9-related samples were selected from the 82 samples in GSE255958.",
            "GSE255958 DTP uses PC9; its resistant samples use the PC9-R derivative.",
            "GSE285632 contributes only GEO sample metadata for now; expression has not been audited.",
            "For GSE285632, EV-only DTP/resistant is the proposed first comparison; OC is an intervention group.",
            "GSE255958 vs GSE285632 changes cell line, study protocol and possibly measurement scale together.",
            "These records are not 96 independent samples for one valid classifier.",
            "No cross-study classification performance is established by this audit.",
        ],
        "decision": "Do not claim cross-study generalization until expression availability, gene IDs, labels and confounding are checked.",
    }
    (OUTPUT / "feasibility_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print("DATASET AUDIT COMPLETED")
    for study, n in sorted(by_study.items()):
        print(f"  {study}: {n} metadata records")
    print("Output:", OUTPUT.relative_to(ROOT))
    print("Reminder: GSE285632 expression data have not been loaded or validated.")


if __name__ == "__main__":
    main()
