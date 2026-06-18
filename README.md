# Vineyard Microbiome — Public Amplicon Dataset Harvest

A reproducible pipeline that discovers, curates, and packages **vineyard / wine
amplicon (metabarcoding) sequencing datasets** from public archives (ENA, which
mirrors NCBI SRA), ready for download and downstream microbiome analysis.

## Current snapshot

| Stage | Studies | Runs | Notes |
|-------|--------:|-----:|-------|
| Raw harvest (ENA) | 166 | 25,316 | keyword + library-strategy search |
| **Curated (KEEP only)** | **130** | **22,777** | manually reviewed, on-topic |
| FASTQ manifest | — | 41,111 files | ≈ **322 GB** total |

- Samples in curated set: **21,488** unique BioSamples
- Platforms: Illumina (99%), plus PacBio, Nanopore, 454, Element
- Marker regions (guessed from titles): ITS / ITS2, 16S (V4 and others), 18S
- First-public years span **2019–2026**

## Directory layout

```
files/
├── README.md                       ← you are here
├── scripts/                        ← the pipeline, run in order
│   ├── 1_harvest_ena.py            ← search ENA, write raw harvest
│   ├── 2_build_curated_set.py      ← apply manual verdicts → curated set
│   └── 3_enrich_biosample_ncbi.py  ← (optional) pull NCBI BioSample attributes
├── metadata/                       ← the data products
│   ├── raw_harvest.tsv             ← all 166 studies, pre-curation (run-level)
│   ├── curated_metadata.tsv        ← 130 KEEP studies (run-level) — the main table
│   ├── project_summary.tsv         ← one row per study (status, n_runs, year, …)
│   └── fastq_manifest.tsv          ← one row per FASTQ file (url + md5 + bytes)
├── curation/                       ← the human review trail
│   ├── inspection.tsv              ← verdict for every study: KEEP / DISCARD / REVIEW
│   ├── discard_list.tsv            ← 25 rejects, with reasons
│   └── review_list.tsv             ← 11 ambiguous studies set aside (excluded)
└── fastq/                          ← (empty) download target for FASTQ files
```

## Pipeline

### 1. Harvest from ENA — `scripts/1_harvest_ena.py`
Queries the ENA portal for AMPLICON (and OTHER) libraries whose text fields
match vineyard keywords (`terroir`, `vineyard*`), filtered by first-public year.
Applies a blocklist (swine/clinical/grapefruit/virome/placenames/etc.), lightly
harmonizes metadata, and guesses the amplicon target region from titles.

```bash
cd files
python3 scripts/1_harvest_ena.py            # → metadata/raw_harvest.tsv (+ raw_*)
# tweak the search:
python3 scripts/1_harvest_ena.py --year-from 2015 --year-to 2026 \
        --wildcard vineyard --precise terroir
```
Writes `raw_harvest.tsv`, `raw_fastq_manifest.tsv`, `raw_project_summary.tsv`.

### 2. Build the curated set — `scripts/2_build_curated_set.py`
Each study in the raw harvest was manually inspected and given a verdict in
`curation/inspection.tsv`:
- **KEEP** (130) — confirmed vineyard/wine microbiome amplicon studies → kept
- **DISCARD** (25) — off-topic (insect gut, apples/cider, virome, teabags, …)
- **REVIEW** (11) — ambiguous (wine-region soils with no explicit confirmation)
  → **excluded** from the curated set

This script reads the raw harvest + verdicts, keeps only KEEP studies, and
regenerates the curated products. It is idempotent.

```bash
python3 scripts/2_build_curated_set.py
# → metadata/curated_metadata.tsv, project_summary.tsv, fastq_manifest.tsv
```

### 3. (Optional) Enrich with NCBI BioSample — `scripts/3_enrich_biosample_ncbi.py`
ENA's flat schema hides submitter attributes (cultivar, soil pH, treatment, …).
This step pulls the full BioSample record for each sample from NCBI.

```bash
export NCBI_API_KEY=...        # strongly recommended: 3× faster, fewer 429s
python3 scripts/3_enrich_biosample_ncbi.py
# → metadata/biosample_attributes.tsv, metadata/curated_metadata_enriched.tsv
```
> ⚠️ Heavy: ~21,500 unique BioSamples, one `esearch` request each. Without an API
> key this takes **2+ hours**. Consider batching `esearch` if re-running often.

Use `--project` to enrich only specific studies — it filters the metadata
*before* querying NCBI, so a single study takes seconds instead of hours:

```bash
python3 scripts/3_enrich_biosample_ncbi.py --project PRJNA1311317
python3 scripts/3_enrich_biosample_ncbi.py --project PRJNA1311317 PRJEB31962
```
Filtered runs are written with a project suffix so they never overwrite the
full enrichment, e.g. `biosample_attributes.PRJNA1311317.tsv` and
`curated_metadata_enriched.PRJNA1311317.tsv`.

### 4. Download FASTQs (not yet run)
`metadata/fastq_manifest.tsv` lists 41,111 files (~322 GB) with md5 checksums.
Example download into `fastq/` with verification:

```bash
tail -n +2 metadata/fastq_manifest.tsv | while IFS=$'\t' read -r run study sample url md5 bytes; do
  out="fastq/$(basename "$url")"
  [ -f "$out" ] || curl -sSL "$url" -o "$out"
  echo "$md5  $out" | md5sum -c -
done
```

### 5. Download fastq and create samplesheet (e.g., only from PRJEB40350)
```bash
mkdir -p fastq/PRJEB40350

# 1) download
awk -F'\t' '$2=="PRJEB40350"{print $4}' metadata/fastq_manifest.tsv | \
while read -r url; do
  out="fastq/PRJEB40350/$(basename "$url")"
  [ -f "$out" ] || curl -fsSL "$url" -o "$out"
done

# 2) md5 check
awk -F'\t' '$2=="PRJEB40350"{print $5"  fastq/PRJEB40350/"$5; sub(/.*\//,"",$4); print $5"  fastq/PRJEB40350/"$4}' metadata/fastq_manifest.tsv >/dev/null
awk -F'\t' '$2=="PRJEB40350"{print $5, "fastq/PRJEB40350/" substr($4, match($4,/[^/]+$/))}' metadata/fastq_manifest.tsv | \
  while read -r md5 f; do echo "$md5  $f"; done | md5sum -c -


DIR=vineyard-microbiome/fastq/PRJEB40350

awk -F'\t' -v dir="$DIR" '
  FNR==NR { if($2=="PRJEB40350") plat[$1]=$11; next }
  $2=="PRJEB40350"{
    bn=$4; sub(/.*\//,"",bn); run=$1
    if (bn ~ /_1\.fastq\.gz$/) r1[run]=dir"/"bn
    else if (bn ~ /_2\.fastq\.gz$/) r2[run]=dir"/"bn
    else se[run]=dir"/"bn
    if (!(run in seen)) { order[++n]=run; seen[run]=1 }
  }
  END{
    print "sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta"
    for(i=1;i<=n;i++){ run=order[i]
      f1=(run in r1)?r1[run]:se[run]; f2=(run in r2)?r2[run]:""
      print i","run","plat[run]","f1","f2"," }
  }
' metadata/curated_metadata.tsv metadata/fastq_manifest.tsv > samplesheet_PRJEB40350.csv

```

## Requirements
```bash
pip install pandas requests
```

## Provenance notes
- `metadata/raw_harvest.tsv` is the unfiltered ENA pull; `curated_metadata.tsv`
  is derived from it by `2_build_curated_set.py` — never edit the curated files
  by hand, edit `curation/inspection.tsv` and re-run step 2.
- The previous flat output (`vineyard_amplicon_out/`) and an empty scratch dir
  (`vineyard_amplicon_out_0/`) have been retired. If `vineyard_amplicon_out_0/`
  still appears, it was locked by another process at cleanup time and can be
  deleted manually.
```
