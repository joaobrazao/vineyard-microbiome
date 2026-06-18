#!/usr/bin/env python3
"""
Build the curated vineyard amplicon set from the raw harvest + manual review.

Takes the raw ENA harvest (1_harvest_ena.py) and the manual curation verdicts
(curation/inspection.tsv), keeps only studies marked KEEP, and writes:
  1. metadata/curated_metadata.tsv  — run-level table, KEEP studies only
  2. metadata/project_summary.tsv   — one row per study (status, n_runs, etc.)
  3. metadata/fastq_manifest.tsv    — one row per downloadable FASTQ file

inspection.tsv must have columns: study_accession, verdict, reason, study_title
where verdict is one of KEEP / DISCARD / REVIEW. Only KEEP is carried forward;
DISCARD and REVIEW studies are excluded from the curated set.

Requires: pandas
    pip install pandas
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def build_summary(df: pd.DataFrame, verdict: dict[str, str]) -> pd.DataFrame:
    summary = (
        df.groupby("study_accession")
          .agg(
            n_runs=("run_accession", "count"),
            n_samples=("sample_accession", lambda s: s[s != ""].nunique()),
            study_title=("study_title", "first"),
            first_public_year=("first_public_year", "first"),
            platforms=("instrument_platform",
                       lambda s: ", ".join(sorted({x for x in s if x}))),
            targets=("guessed_target",
                     lambda s: ", ".join(sorted({x for x in s if x}))),
            countries=("country",
                       lambda s: ", ".join(sorted({x for x in s if x}))),
          )
          .reset_index()
    )
    summary.insert(1, "curation_status",
                   summary["study_accession"].map(verdict).fillna("UNKNOWN"))
    return summary.sort_values(["first_public_year", "n_runs"],
                              ascending=[False, False])


def build_manifest(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        urls = [u for u in r["fastq_ftp"].split(";") if u]
        md5s = r["fastq_md5"].split(";") if r["fastq_md5"] else []
        sizes = r["fastq_bytes"].split(";") if r["fastq_bytes"] else []
        for i, url in enumerate(urls):
            rows.append({
                "run_accession": r["run_accession"],
                "study_accession": r["study_accession"],
                "sample_accession": r["sample_accession"],
                "url": "https://" + url if not url.startswith("http") else url,
                "md5": md5s[i] if i < len(md5s) else "",
                "bytes": sizes[i] if i < len(sizes) else "",
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", type=Path, default=Path("metadata/raw_harvest.tsv"),
                    help="Raw harvest TSV from 1_harvest_ena.py")
    ap.add_argument("--inspection", type=Path,
                    default=Path("curation/inspection.tsv"),
                    help="Manual verdicts (study_accession, verdict, ...)")
    ap.add_argument("--outdir", type=Path, default=Path("metadata"))
    ap.add_argument("--keep", nargs="+", default=["KEEP"],
                    help="Verdicts to carry into the curated set (default: KEEP)")
    args = ap.parse_args()

    if not args.raw.exists():
        sys.exit(f"[error] {args.raw} not found. Run 1_harvest_ena.py first.")
    if not args.inspection.exists():
        sys.exit(f"[error] {args.inspection} not found.")

    raw = pd.read_csv(args.raw, sep="\t", dtype=str).fillna("")
    ins = pd.read_csv(args.inspection, sep="\t", dtype=str).fillna("")
    verdict = dict(zip(ins["study_accession"], ins["verdict"]))

    keep_studies = {s for s, v in verdict.items() if v in set(args.keep)}
    cur = raw[raw["study_accession"].isin(keep_studies)].copy()

    args.outdir.mkdir(parents=True, exist_ok=True)

    cur_path = args.outdir / "curated_metadata.tsv"
    cur.to_csv(cur_path, sep="\t", index=False)
    print(f"[ok] {cur_path}: {len(cur)} runs / "
          f"{cur['study_accession'].nunique()} studies", file=sys.stderr)

    summ = build_summary(cur, verdict)
    summ_path = args.outdir / "project_summary.tsv"
    summ.to_csv(summ_path, sep="\t", index=False)
    print(f"[ok] {summ_path}: {len(summ)} studies", file=sys.stderr)

    man = build_manifest(cur)
    man_path = args.outdir / "fastq_manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False)
    gb = pd.to_numeric(man["bytes"], errors="coerce").sum() / 1e9
    print(f"[ok] {man_path}: {len(man)} files, ~{gb:.1f} GB", file=sys.stderr)

    # Provenance: report what was dropped
    dropped = pd.Series(verdict).value_counts()
    print("\n[info] verdict breakdown across all inspected studies:",
          file=sys.stderr)
    print(dropped.to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
