# Public Health Script Staging Area

These scripts define the ingestion and normalization workflow for public-health variables.

## Current Status

The core public-health refresh helpers are now active enough to support a documented manuscript-facing refresh order.

`ph_09_assess_vaccine_variable_coverage.py` is active and is called by `workflow/lib/run_foundation_checks.py` for vaccine-variable coverage reporting.

## Promotion Rule

Before promoting this module to active status, add:

- a documented execution order,
- required input freeze rules,
- output contracts, and
- a single runnable entry point.

## Active Coverage Helper

- python modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py --out outputs/workflow/checkpoints/vaccine_variable_coverage_report.json

## Minimal Refresh Order

For the manuscript-facing ecology layer, refresh these steps in order rather than in parallel:

1. `python modules/public_health/bin/ph_01_build_source_inventory.py`
2. `python modules/public_health/bin/ph_05_clean_vaccine_programs.py`
3. `python modules/public_health/bin/ph_11_clean_reporting_era_indicators.py`
4. `python modules/public_health/bin/ph_08_build_country_year_master.py`
5. `python modules/public_health/bin/ph_10_clean_highres_cases.py`
6. `python modules/public_health/bin/ph_12_audit_reporting_era_coverage.py`
7. `python modules/public_health/bin/ph_13_build_reporting_era_resolution_worklist.py`
8. `python modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py --out outputs/workflow/checkpoints/vaccine_variable_coverage_report.json`
9. `python workflow/lib/build_ap_exposure_index.py ...`
10. `python workflow/lib/panel_model.py ...`
11. `python manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py`
12. `Rscript manuscript/figures/bin/render_main.R`

`workflow/lib/panel_model.py` depends on the exposure-index output, and the manuscript/figure staging depends on the refreshed panel-model outputs.
On multi-core hosts, `workflow/lib/panel_model.py` can use process workers via `--max-workers` or `PANEL_MAX_WORKERS`.

`ph_01_build_source_inventory.py` now materializes three provenance layers:

- `ph_source_inventory.tsv`: planned source classes from the documentation inventory
- `ph_source_registry.tsv`: canonical per-document public-health source registry
- `ph_source_citation_map.tsv`: row-level citation map linking curated inputs to registry entries

`ph_11_clean_reporting_era_indicators.py` produces the reporting-era sidecar now used by `ph_08_build_country_year_master.py` and high-resolution focal-country surveillance outputs, and validates every cited URL against the canonical source registry.

## High-Resolution Surveillance Helper

For targeted monthly/weekly validation in the small overlap subset:

- `python modules/public_health/bin/ph_10_clean_highres_cases.py`

Outputs:

- `modules/public_health/outputs/ph_highres_cases.tsv`
- `modules/public_health/outputs/ph_highres_overlap_summary.tsv`

## Reporting-Era Helper

- `python modules/public_health/bin/ph_11_clean_reporting_era_indicators.py`
- `python modules/public_health/bin/ph_12_audit_reporting_era_coverage.py`
- `python modules/public_health/bin/ph_13_build_reporting_era_resolution_worklist.py`

Outputs:

- `modules/public_health/outputs/ph_reporting_era_indicators.tsv`
- `modules/public_health/outputs/ph_reporting_era_coverage_audit.tsv`
- `modules/public_health/outputs/ph_reporting_era_resolution_worklist.tsv`

## Optional Prn Curation Overlay

The coverage check can ingest an optional manual curation TSV to fill missing prn_in_vaccine values.

- Default path:
  - modules/public_health/inputs/curation/vaccine_formulation_curation.tsv
- Run example:
  - python modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py --prn-curation modules/public_health/inputs/curation/vaccine_formulation_curation.tsv --out outputs/workflow/checkpoints/vaccine_variable_coverage_report.json

Accepted input formats:

- Rich formulation curation:
  - country_iso3
  - prn_in_vaccine_curated
  - year_start
  - year_end
- Legacy overlay:
  - country_iso3
  - prn_in_vaccine

Optional columns:

- year_start
- year_end
- notes

`modules/public_health/inputs/curation/vaccine_formulation_curation.tsv` is now the single maintained source-of-truth for both the coverage check and formulation-aware manuscript analyses. The simpler legacy overlay remains parser-compatible if passed manually, but it is no longer tracked as a maintained repository input.
