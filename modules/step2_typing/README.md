# Bordetella pertussis: step2 (QC + typing scaffolding)

This folder contains step-2 scripts that build on outputs from `modules/step1_ingest/`.

## Inputs (from step1)

- `../step1_ingest/outputs/bp_metadata_clean.csv`
- `../step1_ingest/outputs/bp_genomes_month_ready/ncbi_dataset/data/` (downloaded genomes)

> Note: some assemblies may be missing locally if the download was interrupted. Step2 will record missing genomes and proceed.

## What step2 produces

- QC-filtered sample table + accession list
- Genome FASTA path index (resolves GCA/GCF naming)
- Optional: MLST calls (if `mlst` is installed)
- Optional: marker scans via BLAST (if `blastn` is installed and you provide query FASTAs)

Outputs are written under `outputs/`.

## Quick start

From this folder:

```bash
bash bin/run_step2.sh
```

If you downloaded genomes to a different location (e.g. `../../pertussis_data/bp_genomes_qc/ncbi_dataset/data`), point step2 at it:

```bash
DATA_ROOT=../../pertussis_data/bp_genomes_qc/ncbi_dataset/data bash bin/run_step2.sh
```

## Optional tool installation (recommended)

These are not required for QC/path indexing, but are needed for typing/marker scans.

- `mlst` (Torsten Seemann)
- `ncbi-blast+` (`blastn`, `makeblastdb`)

If you are using conda:

```bash
conda install -c bioconda -c conda-forge mlst ncbi-blast
```

## Marker scan references

Step2 supports two levels of marker work:

1) `bin/step2_04_marker_scan_blast.py` (hit-level scan)
- BLAST query FASTAs against each assembly and report hit stats.

2) **Step6** (recommended for antigen/AMR analyses)
- Extract marker allele sequences (as `md5` hashes) for easy frequency tracking.
- Call 23S rRNA A2047G (macrolide resistance marker) by mapping reference position 2047.

### Provide marker FASTAs

Put query nucleotide FASTA files under:

- `references/markers/` (e.g. `prn_maker.fasta`, `ptxP_promoter.fasta`, `fim2.fasta`, `fim3.fasta`)

For 23S rRNA calling, also provide:

- `references/23S_rRNA.fasta`

See `references/README.md` for details.

### Run Step6

```bash
bash bin/run_step6.sh
```

Tune parallelism:

```bash
JOBS=40 bash bin/run_step6.sh
```
