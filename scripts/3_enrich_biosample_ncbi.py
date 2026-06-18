#!/usr/bin/env python3
"""
Enrich vineyard amplicon metadata with full NCBI BioSample attributes.

Reads the metadata TSV from vineyard_amplicon_harvest.py and queries NCBI
BioSample for each unique sample_accession, pulling ALL submitter-provided
attributes (cultivar, soil_pH, vineyard, treatment, etc.) that ENA's flat
schema doesn't expose.

Output: a wide TSV where every BioSample attribute becomes its own column.

Requires: pandas, requests
    pip install pandas requests

NCBI rate limit: 3 requests/sec without API key, 10/sec with one.
Get a free key: https://www.ncbi.nlm.nih.gov/account/  -> set NCBI_API_KEY env var.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "")
SLEEP = 0.12 if API_KEY else 0.35  # respect rate limits


def esearch_biosample(accession: str) -> str | None:
    """Convert SAMN/SAMEA accession -> internal BioSample UID."""
    params = {"db": "biosample", "term": accession, "retmode": "json"}
    if API_KEY:
        params["api_key"] = API_KEY
    r = requests.get(f"{EUTILS}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else None


def efetch_biosample_xml(uids: list[str]) -> str:
    """Batch fetch BioSample XML for up to ~200 UIDs."""
    params = {"db": "biosample", "id": ",".join(uids), "rettype": "xml"}
    if API_KEY:
        params["api_key"] = API_KEY
    r = requests.get(f"{EUTILS}/efetch.fcgi", params=params, timeout=120)
    r.raise_for_status()
    return r.text


def parse_biosample_xml(xml_text: str) -> list[dict]:
    """Extract accession + all <Attribute> elements into dicts."""
    out = []
    root = ET.fromstring(xml_text)
    for bs in root.findall(".//BioSample"):
        rec = {"biosample_accession": bs.get("accession", "")}
        # Primary ID (sometimes differs from query accession)
        for idnode in bs.findall(".//Id"):
            if idnode.get("is_primary") == "1":
                rec["biosample_accession"] = idnode.text or rec["biosample_accession"]
        # All attributes — prefer harmonized_name as column key
        for attr in bs.findall(".//Attribute"):
            key = (attr.get("harmonized_name")
                   or attr.get("attribute_name")
                   or "unknown_attr")
            key = re.sub(r"[^a-zA-Z0-9_]+", "_", key).strip("_").lower()
            val = (attr.text or "").strip()
            if val and val.lower() not in {"not applicable", "not collected",
                                            "missing", "na", "n/a", ""}:
                # Handle duplicate keys by appending
                if key in rec:
                    rec[key] = f"{rec[key]} | {val}"
                else:
                    rec[key] = val
        out.append(rec)
    return out


def enrich(samples: list[str], batch_size: int = 100) -> pd.DataFrame:
    """For a list of BioSample accessions, return a wide attribute table."""
    print(f"[info] resolving {len(samples)} BioSample accessions -> UIDs",
          file=sys.stderr)
    acc_to_uid = {}
    for i, acc in enumerate(samples, 1):
        try:
            uid = esearch_biosample(acc)
            if uid:
                acc_to_uid[acc] = uid
        except requests.RequestException as e:
            print(f"[warn] {acc}: {e}", file=sys.stderr)
        time.sleep(SLEEP)
        if i % 25 == 0:
            print(f"  resolved {i}/{len(samples)}", file=sys.stderr)

    print(f"[info] fetching {len(acc_to_uid)} BioSample XML records",
          file=sys.stderr)
    all_records = []
    uids = list(acc_to_uid.values())
    for i in range(0, len(uids), batch_size):
        chunk = uids[i:i + batch_size]
        try:
            xml = efetch_biosample_xml(chunk)
            all_records.extend(parse_biosample_xml(xml))
        except (requests.RequestException, ET.ParseError) as e:
            print(f"[warn] batch {i}: {e}", file=sys.stderr)
        time.sleep(SLEEP)
        print(f"  fetched {min(i + batch_size, len(uids))}/{len(uids)}",
              file=sys.stderr)

    return pd.DataFrame(all_records).fillna("")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", type=Path,
                    default=Path("metadata/curated_metadata.tsv"),
                    help="Input TSV (curated set from 2_build_curated_set.py)")
    ap.add_argument("--outdir", type=Path,
                    default=Path("metadata"))
    ap.add_argument("--sample-col", default="sample_accession")
    ap.add_argument("--study-col", default="study_accession")
    ap.add_argument("--project", nargs="+", default=None,
                    help="Only enrich samples from these study/project "
                         "accessions, e.g. --project PRJNA1311317. "
                         "Default: all projects in the metadata.")
    args = ap.parse_args()

    if not args.metadata.exists():
        sys.exit(f"[error] {args.metadata} not found. Run the harvester first.")

    df = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")

    suffix = ""
    if args.project:
        if args.study_col not in df.columns:
            sys.exit(f"[error] column '{args.study_col}' not in {args.metadata}")
        wanted = set(args.project)
        df = df[df[args.study_col].isin(wanted)]
        if df.empty:
            sys.exit(f"[error] no rows match project(s): {', '.join(wanted)}")
        found = sorted(set(df[args.study_col]))
        missing = wanted - set(found)
        if missing:
            print(f"[warn] no rows for: {', '.join(sorted(missing))}",
                  file=sys.stderr)
        print(f"[info] filtered to {len(found)} project(s): {', '.join(found)}",
              file=sys.stderr)
        suffix = "." + "_".join(found)

    samples = sorted(set(df[args.sample_col]) - {""})
    print(f"[info] {len(samples)} unique BioSamples to enrich", file=sys.stderr)
    if not API_KEY:
        print("[hint] set NCBI_API_KEY env var for 3x faster fetching",
              file=sys.stderr)

    enriched = enrich(samples)

    # Merge back onto the run-level table
    merged = df.merge(enriched,
                      left_on=args.sample_col,
                      right_on="biosample_accession",
                      how="left")

    out_attrs = args.outdir / f"biosample_attributes{suffix}.tsv"
    enriched.to_csv(out_attrs, sep="\t", index=False)
    print(f"[ok] wrote {out_attrs} ({len(enriched)} samples, "
          f"{enriched.shape[1]} attribute columns)", file=sys.stderr)

    out_merged = args.outdir / f"curated_metadata_enriched{suffix}.tsv"
    merged.to_csv(out_merged, sep="\t", index=False)
    print(f"[ok] wrote {out_merged}", file=sys.stderr)

    # Quick view of the most populated attribute columns
    fill = (enriched.replace("", pd.NA).notna().sum()
                    .sort_values(ascending=False))
    print("\nMost populated BioSample attributes:", file=sys.stderr)
    print(fill.head(25).to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
