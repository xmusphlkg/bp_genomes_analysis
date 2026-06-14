# Bordetella pertussis prn Disruption Analysis

This repository contains the analysis code, small curated inputs, frozen figure data,
and manuscript-facing reproducibility package for the study:

**Recurrent structural routes to pertactin gene disruption in Bordetella pertussis**

The public release is designed to support two levels of reuse:

- Reader-level inspection of the frozen manuscript figures, source-data workbooks,
  supplementary tables, and validation ledgers.
- Code-level reproduction of manuscript-facing extracts and figures from the
  tracked workflow inputs and public genome accessions.

Large recomputable workflow products, raw reads, local runtime overrides, submission
letters, author metadata, editorial checklist notes, and publisher full-text PDFs are
not part of the public release.

## Repository Layout

- `config/`: workflow settings, module contracts, and environment specifications.
- `workflow/`: cross-stage orchestration, Snakemake rules, shared helpers, and tests.
- `modules/`: stage-owned analysis code, schemas, references, and small curated inputs.
- `state/`: tracked manifests, checkpoints, ledgers, and small canonical registries.
- `docs/`: technical documentation, policies, and release-boundary notes.
- `manuscript/`: frozen manuscript-facing figures, figure data, supplementary tables,
  source-data workbooks, and scripts for rebuilding those release artifacts.
- `apps/`: optional interactive delivery scaffolds.

## Environment

The workflow uses Conda/Mamba environments defined under `config/env/`.

```bash
$HOME/miniforge3/bin/mamba env create -f config/env/environment_tool.yml
$HOME/miniforge3/bin/mamba env create -f config/env/environment_python.yml
$HOME/miniforge3/bin/mamba env create -f config/env/environment_r.yml
```

These files default to the official `conda-forge` and `bioconda` channels. If you
use a mirror, switch channels only after confirming that the mirror metadata is in
sync with those upstream channels.

## Reproducibility Modes

Lightweight manuscript reproduction uses the frozen package under `manuscript/`:

- `manuscript/figure_data/`
- `manuscript/supplementary/`
- `manuscript/submission_data/source_data/`
- `manuscript/figures/scripts/`
- `manuscript/figures/bin/`

Full recomputation uses the tracked workflow and public accessions to regenerate
large intermediate outputs outside Git. Generated workflow outputs are intentionally
not committed to this repository.

## Main Entry Points

```bash
bash workflow/bin/m0_foundation.sh
bash workflow/bin/m1_m2_qc.sh
bash workflow/bin/m3_snippy.sh --dry-run
bash workflow/bin/m4_phylogeny.sh --dry-run
bash workflow/bin/m5_asr.sh --dry-run
bash workflow/bin/run_full_workflow.sh --dry-run

bash modules/step4_prn_validation/bin/run_read_validation.sh --help
bash modules/step6_epi_transmission/bin/run_re_estimation.sh --help
bash modules/step6_epi_transmission/bin/run_transmission_models.sh --help

bash manuscript/bin/refresh_submission.sh
Rscript manuscript/figures/bin/render_main.R
Rscript manuscript/figures/bin/render_extended_data.R
```

## Public Release Boundary

The public package includes code, curated manifests, manuscript figure data,
rendered figures, supplementary tables, and source-data workbooks. It excludes:

- `manuscript/text/`
- cover letters, reporting-summary notes, editorial checklist notes, and author
  metadata templates
- `archive/`
- `outputs/`, `logs/`, `.snakemake/`, caches, virtual environments, and local
  runtime overrides
- raw reads, BAM/VCF/intermediate alignment files, and other large generated files

See `docs/release_package.md` and `manuscript/PUBLIC_RELEASE.md` for the detailed
release boundary.
