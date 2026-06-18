#!/usr/bin/env python3
"""
Vineyard/wine amplicon data harvester.

Searches ENA (which mirrors NCBI SRA) for AMPLICON libraries matching
vineyard/wine/grape-related keywords, filters by year, and produces:
  1. A merged, cleaned metadata TSV
  2. A FASTQ download manifest (URLs + md5)
  3. A short summary of projects found

Requires: pandas, requests
    pip install pandas requests
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from io import StringIO
from pathlib import Path

import pandas as pd
import requests

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"

# Fields we want for every run. Add/remove as needed.
FIELDS = [
    "run_accession", "study_accession", "sample_accession", "experiment_accession",
    "scientific_name", "tax_id", "library_strategy", "library_source",
    "library_selection", "library_layout", "instrument_platform",
    "instrument_model", "read_count", "base_count",
    "first_public", "last_updated", "country", "collection_date",
    "location", "environment_biome", "environment_feature",
    "environment_material", "host", "host_scientific_name",
    "isolation_source", "sample_title", "study_title", "study_alias",
    "sample_alias", "library_name", "experiment_title",
    "fastq_ftp", "fastq_md5", "fastq_bytes", "submitted_ftp",
]

# Keyword strategy:
#  - PRECISE keywords: matched as whole-word tokens (no wildcards).
#    Safe for short/ambiguous terms.
#  - WILDCARD keywords: matched as substrings (*term*).
#    Use ONLY for unambiguous long compounds where you want morphological
#    variants ("vineyard" -> "vineyards", "vineyard-soil").
DEFAULT_PRECISE_KEYWORDS = [
    "terroir",
]
DEFAULT_WILDCARD_KEYWORDS = [
    "vineyard",      # -> vineyards, vineyard-soil, vineyard_rhizosphere
]

# Blocklist: if ANY of these substrings appears in study/sample/host/source
# fields, the run is dropped. Tuned for vineyard work — extend as needed.
DEFAULT_BLOCKLIST = [
    # mammals that contain "ine"/"vine" or use "wine"-adjacent terms
    "swine", "porcine", "bovine", "ovine", "equine", "feline", "canine",
    "murine", "piglet", "pig ", "pigs", "cattle", "calf", "calves",
    "sheep", "goat", "buffalo", "horse",
    # human/clinical
    "human gut", "patient", "clinical", "hospital", "infant", "fecal",
    "stool", "feces", "vaginal", "oral microbiome", "skin microbiome",
    # other off-topic but keyword-adjacent
    "grapeseed oil", "grapefruit",
    "ravine", "vineland",  # placenames that contain "vine"
    # virology / non-microbiome studies sometimes called "AMPLICON"
    "h5n1", "influenza", "sars-cov", "covid", "hiv",
]


def build_query(precise: list[str], wildcard: list[str],
                year_from: int, year_to: int,
                library_strategies: list[str] | None = None) -> str:
    """Build an ENA query string with two-tier keyword matching."""
    # Fields where on-topic terms tend to appear.
    text_fields = [
        "study_title", "experiment_title", "sample_title",
        "isolation_source", "host", "scientific_name",
        "environment_material", "environment_feature", "environment_biome",
        "sample_alias", "library_name",  # extra hooks for sparse metadata
    ]

    kw_clauses = []

    # Precise: whole-word match via quoted phrase
    for kw in precise:
        kw_escaped = kw.replace('"', '').strip()
        if not kw_escaped:
            continue
        per_field = [f'{f}="{kw_escaped}"' for f in text_fields]
        kw_clauses.append("(" + " OR ".join(per_field) + ")")

    # Wildcard: substring match for morphological variants
    for kw in wildcard:
        kw_escaped = kw.replace('"', '').strip()
        if not kw_escaped:
            continue
        per_field = [f'{f}="*{kw_escaped}*"' for f in text_fields]
        kw_clauses.append("(" + " OR ".join(per_field) + ")")

    if not kw_clauses:
        raise ValueError("No keywords supplied")
    kw_block = "(" + " OR ".join(kw_clauses) + ")"

    # library_strategy: by default include AMPLICON + OTHER (catches misannotated)
    if library_strategies is None:
        library_strategies = ["AMPLICON", "OTHER"]
    strat_clause = "(" + " OR ".join(
        f'library_strategy="{s}"' for s in library_strategies) + ")"

    parts = [
        strat_clause,
        f'first_public>={year_from}-01-01',
        f'first_public<={year_to}-12-31',
        kw_block,
    ]
    return " AND ".join(parts)


def apply_blocklist(df: pd.DataFrame, blocklist: list[str]) -> pd.DataFrame:
    """Drop rows whose text fields contain any blocklisted substring."""
    if df.empty or not blocklist:
        return df
    cols = ["study_title", "experiment_title", "sample_title",
            "isolation_source", "host", "scientific_name",
            "environment_material", "environment_feature",
            "environment_biome", "sample_alias", "library_name"]
    cols = [c for c in cols if c in df.columns]
    blob = df[cols].fillna("").agg(" ".join, axis=1).str.lower()

    pattern = "|".join(re.escape(b.lower()) for b in blocklist)
    mask_bad = blob.str.contains(pattern, regex=True, na=False)
    n_drop = mask_bad.sum()
    if n_drop:
        print(f"[info] blocklist removed {n_drop} runs", file=sys.stderr)
        # Show a small sample of what got dropped so you can sanity-check
        dropped_studies = (df.loc[mask_bad, ["study_accession", "study_title"]]
                             .drop_duplicates().head(10))
        for _, row in dropped_studies.iterrows():
            print(f"  dropped: {row['study_accession']}  {row['study_title'][:80]}",
                  file=sys.stderr)
    return df.loc[~mask_bad].copy()


def query_ena(query: str, fields: list[str], limit: int = 0) -> pd.DataFrame:
    """Run an ENA portal query and return a DataFrame."""
    params = {
        "result": "read_run",
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "limit": limit,  # 0 = no limit
    }
    print(f"[info] querying ENA…", file=sys.stderr)
    r = requests.get(ENA_PORTAL, params=params, timeout=120)
    r.raise_for_status()
    if not r.text.strip():
        return pd.DataFrame(columns=fields)
    df = pd.read_csv(StringIO(r.text), sep="\t", dtype=str).fillna("")
    return df


def clean_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Light harmonization. Expand as you discover quirks in your own data."""
    if df.empty:
        return df

    # Normalize collection_date to year where possible.
    def to_year(s: str) -> str:
        s = s.strip()
        if not s:
            return ""
        for fmt in ("%Y-%m-%d", "%Y-%m", "%Y/%m/%d", "%d/%m/%Y", "%Y"):
            try:
                return str(pd.to_datetime(s, format=fmt).year)
            except (ValueError, TypeError):
                continue
        try:
            return str(pd.to_datetime(s, errors="coerce").year)
        except Exception:
            return ""

    df = df.copy()
    df["collection_year"] = df["collection_date"].apply(to_year)
    df["first_public_year"] = df["first_public"].str[:4]

    # Heuristic: tag likely amplicon target region from titles.
    def guess_region(row) -> str:
        blob = " ".join([
            row.get("study_title", ""), row.get("experiment_title", ""),
            row.get("sample_title", ""), row.get("library_name", ""),
        ]).lower()
        # crude but useful
        markers = [
            ("its2", "ITS2"), ("its1", "ITS1"), ("its ", "ITS"),
            ("v3-v4", "16S V3-V4"), ("v3v4", "16S V3-V4"),
            ("v4-v5", "16S V4-V5"), ("v1-v3", "16S V1-V3"),
            ("v4", "16S V4"), ("16s", "16S (unspecified)"),
            ("18s", "18S"), ("d1-d2", "26S/28S D1-D2"),
            ("26s", "26S"), ("28s", "28S"),
        ]
        for needle, label in markers:
            if needle in blob:
                return label
        return ""

    df["guessed_target"] = df.apply(guess_region, axis=1)

    # Sort: newest first, then by study
    df = df.sort_values(["first_public_year", "study_accession"],
                        ascending=[False, True])
    return df


def write_outputs(df: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    meta_path = outdir / "raw_harvest.tsv"
    df.to_csv(meta_path, sep="\t", index=False)
    print(f"[ok] wrote {meta_path} ({len(df)} runs)", file=sys.stderr)

    # Download manifest: one row per FASTQ file (ENA gives ';'-separated lists for paired).
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
    manifest = pd.DataFrame(rows)
    manifest_path = outdir / "raw_fastq_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print(f"[ok] wrote {manifest_path} ({len(manifest)} files)", file=sys.stderr)

    # Project-level summary
    if not df.empty:
        summary = (
            df.groupby("study_accession")
              .agg(n_runs=("run_accession", "count"),
                   study_title=("study_title", "first"),
                   first_public_year=("first_public_year", "first"),
                   platforms=("instrument_platform",
                              lambda s: ", ".join(sorted(set(s)))),
                   targets=("guessed_target",
                            lambda s: ", ".join(sorted({x for x in s if x}))))
              .sort_values("first_public_year", ascending=False)
        )
        summary_path = outdir / "raw_project_summary.tsv"
        summary.to_csv(summary_path, sep="\t")
        print(f"[ok] wrote {summary_path} ({len(summary)} projects)",
              file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year-from", type=int, default=2015,
                    help="Earliest first_public year (default: 2015)")
    ap.add_argument("--year-to", type=int, default=2026,
                    help="Latest first_public year (default: 2026)")
    ap.add_argument("--precise", nargs="+", default=DEFAULT_PRECISE_KEYWORDS,
                    help="Whole-word keywords (safe for short/ambiguous terms)")
    ap.add_argument("--wildcard", nargs="+", default=DEFAULT_WILDCARD_KEYWORDS,
                    help="Substring keywords (catches morphological variants)")
    ap.add_argument("--blocklist", nargs="+", default=DEFAULT_BLOCKLIST,
                    help="Substrings that disqualify a run (post-filter)")
    ap.add_argument("--blocklist-file", type=Path, default=None,
                    help="Path to file with one blocklist term per line "
                         "(overrides --blocklist, no quoting needed)")
    ap.add_argument("--no-blocklist", action="store_true",
                    help="Disable blocklist filtering")
    ap.add_argument("--library-strategies", nargs="+",
                    default=["AMPLICON", "OTHER"],
                    help="Allowed library_strategy values (default: AMPLICON OTHER). "
                         "Add WGS if you want to include shotgun-bundled amplicons.")
    ap.add_argument("--outdir", type=Path, default=Path("metadata"))
    ap.add_argument("--limit", type=int, default=0,
                    help="Max rows from ENA (0 = unlimited)")
    args = ap.parse_args()

    print(f"[info] precise keywords:  {args.precise}", file=sys.stderr)
    print(f"[info] wildcard keywords: {args.wildcard}", file=sys.stderr)
    print(f"[info] year range: {args.year_from}–{args.year_to}", file=sys.stderr)

    query = build_query(args.precise, args.wildcard,
                        args.year_from, args.year_to,
                        library_strategies=args.library_strategies)
    print(f"[info] query:\n  {query}", file=sys.stderr)

    df = query_ena(query, FIELDS, limit=args.limit)
    print(f"[info] ENA returned {len(df)} runs before blocklist",
          file=sys.stderr)

    if not args.no_blocklist:
        blocklist = args.blocklist
        if args.blocklist_file:
            blocklist = [line.strip() for line in args.blocklist_file.read_text().splitlines()
                         if line.strip() and not line.strip().startswith("#")]
            print(f"[info] loaded {len(blocklist)} blocklist terms from "
                  f"{args.blocklist_file}", file=sys.stderr)
        df = apply_blocklist(df, blocklist)

    df = clean_metadata(df)
    write_outputs(df, args.outdir)

    if not df.empty:
        n_projects = df["study_accession"].nunique()
        print(f"\n[done] {len(df)} runs across {n_projects} projects",
              file=sys.stderr)


if __name__ == "__main__":
    main()