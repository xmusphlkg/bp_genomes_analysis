# Bordetella pertussis: NCBI metadata step1 (generated 2026-03-04 20:33)

This folder contains step-1 outputs for *Bordetella pertussis* genome assembly metadata (NCBI Datasets CLI).

## Pipeline steps
1) Fetch genome assembly report (JSONL) from NCBI Datasets.
2) Export min/extended TSV (including BioSample collection date & geographic location).
3) Clean metadata: derive country/year/month_key/week_key; create count tables.

## Key statistics
- Total assemblies in report: 3370
- Missing geographic location (raw): 84 (2.49%)
- Missing collection date (raw): 84 (2.49%)

### Date resolution breakdown
- year: 2528
- full_date: 670
- missing: 126
- year_month: 28
- unparsed: 18

- Year range (year available): 1920–2025

## Modeling note (important)
- Many samples only have **year-level** dates. For global spread + national weekly/monthly cases,
  use a **two-track strategy**:
  1) Global analysis at **annual** resolution: use `*_country_year_counts.csv`.
  2) High-resolution (monthly/weekly) analysis on a subset with full dates: use `*_samples_month_ready.csv`.

## Top countries by genomes (year available)
- USA: 1488
- China: 886
- Japan: 98
- New Zealand: 96
- Australia: 80
- South Africa: 63
- United Kingdom: 62
- Netherlands: 49
- France: 40
- Sweden: 39
- Taiwan: 34
- Finland: 30
- Italy: 29
- Canada: 28
- Czech Republic: 28
- Brazil: 26
- Argentina: 22
- India: 22
- Kenya: 21
- Denmark: 18

## Top countries by genomes (month-ready subset)
- China: 572
- Japan: 59
- Spain: 18
- Czech Republic: 10
- USA: 10
- United Kingdom: 10
- Mexico: 6
- Australia: 6
- Brazil: 2
- India: 2
- Italy: 2
- France: 1

## Output files
- `bp_genome_report.jsonl`: assembly report (JSONL)
- `bp_min_metadata.tsv`: small TSV for quick inspection
- `bp_extended_metadata.tsv`: extended TSV with BioSample date/geo
- `bp_metadata_clean.csv`: per-assembly cleaned metadata with derived time keys
- `bp_date_resolution_summary.csv`: date resolution counts
- `bp_country_year_counts.csv`: country×year genome counts
- `bp_country_month_counts.csv`: country×month counts (where month exists)
- `bp_country_week_counts.csv`: country×ISO-week counts (where full date exists)
- `bp_samples_month_ready.csv`: subset usable for monthly analysis
- `bp_country_month_counts_hires.csv`: month counts from month-ready subset
- `assembly_accessions_month_ready.txt`: accessions for downloading genomes (month-ready)
- `assembly_accessions_all.txt`: accessions for downloading genomes (all)

## Download genomes (optional)
To download genomes referenced by the accessions lists, run:

- `bash bin/run_step1.sh --download-genomes`

### Speed tips
- Use parallel *datasets* downloads by splitting the accession list (recommended):
  - `PARALLEL=4 bash bin/run_step1.sh --download-genomes`
- aria2 mode requires a URL list file (one URL per line). If `ARIA2_URLS` is not set, the pipeline will automatically fall back to `datasets` download.
  - Example: `USE_ARIA2=1 ARIA2_URLS=urls.txt ARIA2_JOBS=8 bash bin/run_step1.sh --download-genomes`

## Distributed raw-read expansion (3 servers)

The repository now includes a distributed implementation to accelerate the
SRA/ENA raw-read expansion stage.

### 1) Build a run-level download plan

```bash
python3 bin/raw_reads/10_build_download_plan.py
```

Default inputs are:
- `outputs/bp_cohort_D_validation.tsv` (highest priority)
- `outputs/bp_cohort_C_country_year.tsv` (second priority)

Default output is:
- `../step4_prn_validation/inputs/bp_raw_reads_download_plan.tsv`

The planner now enriches each run with ENA metadata when possible:
- `ena_library_layout`
- `ena_instrument_platform`
- `ena_fastq_ftp`
- `estimated_total_bytes`
- `download_strategy`

This lets downstream steps skip single-end / non-Illumina runs early and prefer direct gzipped FASTQ downloads over `prefetch + fasterq-dump`.

### 2) Optional split plan for multiple servers

```bash
python3 bin/raw_reads/11_split_download_plan.py --servers server1,server2,server3
```

Outputs are written to:
- `../step4_prn_validation/inputs/shards/bp_raw_reads_download_plan.<server>.tsv`
- `../step4_prn_validation/inputs/shards/bp_raw_reads_runs.<server>.txt`
- `../step4_prn_validation/inputs/shards/run_<server>.sh`
- `../step4_prn_validation/inputs/shards/bp_raw_reads_shard_summary.tsv`

Shards are now balanced greedily by estimated FASTQ bytes when metadata is available, instead of plain accession hashing.
These shard files are runtime acceleration artifacts and are intentionally not tracked in the release repository. A single-host release can call the same worker with a local plan, local `--threads`, and local `--jobs` settings.

### 3) Run one shard on each server

Use the generated launcher on each machine, or call the worker script directly:

```bash
bash bin/raw_reads/12_run_shard_to_assembly.sh \
  --plan-tsv ../step4_prn_validation/inputs/shards/bp_raw_reads_download_plan.server1.tsv \
  --run-list ../step4_prn_validation/inputs/shards/bp_raw_reads_runs.server1.txt \
  --workdir ../step4_prn_validation/work/server1 \
  --outdir ../step4_prn_validation/outputs/assemblies/server1 \
  --threads 12 \
  --jobs 2
```

Notes:
- This stage is CPU/IO bound; GPU is typically not used by `prefetch`,
  `fasterq-dump`, or `shovill`.
- The worker prefers direct ENA FASTQ downloads whenever paired-end Illumina files are available.
- Single-end and non-Illumina runs are skipped early and recorded in `run_status.tsv`.
- `fasterq-dump` fallback work is isolated under each shard workdir so temporary files no longer spill into the repository root.
- The shard workdir is now the canonical state entrypoint for distributed launches:
  `launcher.pid`, `launcher.log`, and `launcher.command.sh` are written there by the step4 launcher.

### Troubleshooting: worker appears stuck

If a shard launcher seems to hang or exits immediately, run the preflight check:

```bash
bash bin/raw_reads/13_preflight_env.sh
```

Current required tools:
- `prefetch` and `fasterq-dump` (from `sra-tools`)
- `shovill`
- `pigz` (optional but recommended)

Quick setup on each server:

```bash
CONDA_NO_PLUGINS=true CONDA_SOLVER=classic conda env create -f environment.yml
conda activate pertussis-prn-global
```

## One-command distributed launch with env

If you prefer central orchestration from one control machine:

1) Create a local-only `env` file and fill server values (prefer SSH key). This file is ignored by git and must not be committed:

```bash
cat > env <<'EOF'
HOST1=<server1-hostname-or-ip>
USER1=<ssh-user>
HOST2=<server2-hostname-or-ip>
USER2=<ssh-user>
SSH_KEY=~/.ssh/id_ed25519
CONDA_ENV=pertussis-prn-global
JOBS=20
THREADS=6
DATA_ROOT=<remote-or-shared-step4-data-root>
EOF
```

2) Launch all three servers from the control machine:

```bash
bash modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh --env-file env
```

3) Collect status summaries at any time:

```bash
bash modules/step4_prn_validation/bin/step4_00_collect_distributed_status.sh --env-file env
```

This reports both launcher state and `run_status.tsv` progress for each shard, using the same resolved shard paths as launch/sync.

4) Sync remote assemblies and logs back to the control machine after shards finish:

```bash
bash modules/step4_prn_validation/bin/step4_00_sync_distributed_outputs.sh --env-file env
```

Security note: avoid storing plaintext passwords when possible. If password mode is
used, keep `env` local-only and rotate credentials after long runs.

## External Raw-Read-Only Gap Fill

The distributed validation/assembly queue above is built from the current retained assembly-backed manifest.
To expand beyond the current assembly universe and move toward the larger public raw-read universe:

### 1) Fetch the full ENA taxon-level read-run catalog

```bash
python3 bin/raw_reads/15_fetch_taxon_read_run_catalog.py
```

Default output:
- `outputs/bp_ena_taxon_read_run_catalog.tsv`

This pulls the public `read_run` catalog for `tax_tree(520)` from ENA.

### 2) Subtract the current assembly-backed manifest and build an external-only gap-fill plan

```bash
python3 bin/raw_reads/16_build_external_gapfill.py
```

Default outputs:
- `outputs/bp_external_raw_reads_only_samples.tsv`
- `outputs/bp_external_raw_reads_only_plan.tsv`

Design notes:
- rows already covered by the retained assembly-backed manifest are removed by `run_accession`, `biosample_accession`, and linked sample accession overlap;
- the sample-level output summarizes all external-only samples;
- the run-level plan emits only paired-end Illumina runs with direct FASTQ links by default, so it can feed the existing shard worker without reintroducing incompatible long-read or single-end runs.

### 3) Split the external-only plan without touching the current validation shards

Use the existing splitter, but point it at a separate output directory:

```bash
python3 bin/raw_reads/11_split_download_plan.py \
  --plan outputs/bp_external_raw_reads_only_plan.tsv \
  --servers server1,server2,server3 \
  --outdir ../step4_prn_validation/inputs/external_gapfill_shards
```

This preserves the current `../step4_prn_validation/inputs/shards/` directory that may already be in active use.

## Post-Assembly QC And Merge

After raw-read runs finish assembling, collect the completed `contigs.fa` outputs:

```bash
python3 bin/raw_reads/17_collect_assembled_genomes.py
```

Default output:
- `outputs/bp_raw_read_assembly_manifest.tsv`

Run post-assembly QC. The script will use `quast.py` / `quast` and `checkm` if they are installed, and will mark rows as `pending_checkm` if CheckM metrics are unavailable.

```bash
python3 bin/raw_reads/18_qc_assembled_genomes.py
```

Default outputs:
- `outputs/bp_raw_read_assembly_qc.tsv`
- `outputs/bp_raw_read_assembly_qc_pass.tsv`

Finally merge QC-passed raw-read assemblies with the retained public assembly manifest:

```bash
python3 bin/raw_reads/19_merge_qc_passed_genomes.py
```

Default outputs:
- `outputs/bp_combined_public_plus_raw_read_manifest.tsv`
- `outputs/bp_raw_read_merge_exclusions.tsv`
