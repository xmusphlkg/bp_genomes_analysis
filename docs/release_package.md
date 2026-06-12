# Release Package Boundary

This document defines the repository surface intended for paper code and data sharing.

## Include In The Repository Release

- `README.md`, `LICENSE`, `.gitattributes`, `.gitignore`, and the top-level `Snakefile` symlink.
- `config/`, except untracked local overrides such as `config/runtime/runtime_envs.local.env`.
- `workflow/`, including `workflow/bin/`, `workflow/lib/`, `workflow/rules/`, and `workflow/tests/`.
- `modules/`, including active stage-owned code, references, schemas, tests, and small curated inputs.
- `state/`, including tracked manifests, checkpoints, ledgers, and schemas used as small canonical registries.
- `manuscript/PUBLIC_RELEASE.md`, `manuscript/figure_data/`, `manuscript/supplementary/`, and `manuscript/submission_data/` as the frozen submission-facing data package.
- `manuscript/figures/scripts/` and `manuscript/figures/bin/` as figure rendering code.
- `manuscript/scripts/freeze/`, `manuscript/scripts/sidecars/`, `manuscript/scripts/diagnostics/`, `manuscript/scripts/source_data/`, and `manuscript/bin/` only when the public package should support rebuilding frozen manuscript-facing extracts from workflow outputs and reviewer-facing diagnostic ledgers.
- `manuscript/figures/outputs/main/` and `manuscript/figures/outputs/extended_data/` when final figure renderings are being shared.
- `archive/docs/references/refer/` may be retained as a small human-readable literature cache when it contains redistributable notes only; do not redistribute publisher full-text PDFs publicly unless their licences explicitly permit it.
- `apps/` if the interactive Shiny scaffold is part of the release target.

## Exclude From The Repository Release

- `outputs/`, `logs/`, `.snakemake/`, `.pytest_cache/`, `__pycache__/`, `fastp.html`, and `fastp.json`.
- Local runtime files: root `env`, `.env`, virtualenv directories, and `config/runtime/runtime_envs.local.env`.
- Generated Step4 distributed shard folders under `modules/step4_prn_validation/inputs/shards/` and `modules/step4_prn_validation/inputs/external_gapfill_shards*/`.
- Module work/output folders such as `modules/*/work/` and `modules/*/outputs/`.
- Audit-only or historical folders named `archive/` or `_archive/`, except the optional literature cache noted above.
- Internal manuscript submission files: `manuscript/text/`, `manuscript/scripts/review/`, `manuscript/scripts/validate_commsbio_submission.py`, and `manuscript/focal_country_dynamics_audit.md`.
- Repository-maintenance tools and reports such as `bfg.jar` and `..bfg-report/`.
- Raw reads, BAM/VCF/intermediate alignment files, and large recomputable workflow products; keep these in external public archives, Git LFS, or the NAS-backed `pertussis_data/` data release.

## Reproducibility Boundary

- Reader-level reproduction is supported from the frozen manuscript package: `manuscript/figure_data/`, `manuscript/supplementary/`, `manuscript/submission_data/source_data/`, and the figure rendering scripts.
- Reviewer-risk and robustness ledgers are rebuilt from `manuscript/scripts/diagnostics/`; include that directory whenever claiming that the public code release can regenerate the manuscript-facing diagnostic tables.
- End-to-end recomputation is supported by the active workflow, but it is intentionally heavier: it rebuilds `outputs/workflow/` from tracked manifests, public accessions, configured environments, and external bioinformatics tools.
- Generated `outputs/workflow/` files are not required in Git for reproducibility; they are rebuild products. If exact frozen intermediates are needed for audit, publish them in an external data archive rather than the code repository.
- The historical three-server raw-read execution was an acceleration strategy. The release code should remain runnable on one host with local `threads`/`jobs`; distributed shard files are generated runtime plans, not source artifacts.

## Rebuild Entry Points

- The active M0-M5 execution path is `workflow/bin/run_full_workflow.sh`.
- Step4 read validation is owned by `modules/step4_prn_validation/bin/run_read_validation.sh`.
- Step6 transmission utilities are owned by `modules/step6_epi_transmission/bin/run_re_estimation.sh` and `modules/step6_epi_transmission/bin/run_transmission_models.sh`.
- The manuscript package refresh entry point is `manuscript/bin/refresh_submission.sh`.
- Main and supplementary figure assets are rendered with `manuscript/figures/bin/render_main.R` and `manuscript/figures/bin/render_extended_data.R`; the latter retains the historical renderer name while the Communications Biology submission labels the outputs as Supplementary Figures.

## Data Availability Notes

- Public genome/read accessions and cohort-level metadata are retained in `state/manifest/` and `manuscript/figure_data/project_genome_metadata_manifest.tsv`.
- Figure-specific source data workbooks live in `manuscript/submission_data/source_data/`.
- Plot-ready TSV extracts live in `manuscript/figure_data/` and are the preferred lightweight data-sharing layer for reviewers and readers.
- Large recomputable workflow outputs should be regenerated from the documented workflow or synchronized from the external project data root, not committed to the release repository.
