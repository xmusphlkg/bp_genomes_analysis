# Vaccine Formulation Curation Inputs

This folder stores the auditable manual curation used by the public-health readiness checks and the manuscript-facing formulation-aware exposure builder.

## Files

- `vaccine_formulation_curation.tsv`: primary source-of-truth for formulation-aware curation used by `workflow/lib/build_ap_exposure_index.py`
- `vaccine_product_metadata.tsv`: role-specific product/programme metadata used to estimate PRN-positive vs PRN-free exposure shares for routine primary, routine booster, and maternal contexts
- `PRN vaccine status.md`: human-readable product/programme evidence log used to support TSV updates

## `vaccine_formulation_curation.tsv` Required Columns

- `country_iso3`
- `country_name`
- `year_start`
- `year_end`
- `ap_timing_anchor_year`
- `primary_series_formulation`
- `booster_formulation`
- `prn_in_vaccine_curated`
- `prn_in_vaccine_source_class`
- `formulation_confidence`
- `source_name`
- `source_url`
- `source_release_date`
- `notes`

## `vaccine_product_metadata.tsv` Required Columns

- `country_iso3`
- `country_name`
- `year_start`
- `year_end`
- `exposure_role`
- `region_scope`
- `product_name`
- `manufacturer`
- `product_platform`
- `ap_prn_positive_fraction`
- `population_share`
- `share_basis`
- `evidence_confidence`
- `source_name`
- `source_url`
- `source_release_date`
- `notes`

### `vaccine_product_metadata.tsv` Conventions

- `exposure_role` uses `routine_primary`, `routine_booster`, or `maternal`.
- `product_platform` uses `wp`, `ap_prn_positive`, `ap_prn_negative`, or `ap_mixed`.
- `ap_prn_positive_fraction` is the PRN-positive share within that row's acellular exposure. Use `1` for uniformly PRN-positive products, `0` for PRN-free acellular products, and fractional values such as `0.6` or `14/19`-derived decimals for documented mixed exposure summaries.
- `population_share` is the share of the national role-specific exposed population represented by that row. It can be `1` for a national summary row or split across multiple products when equal-split or region-share assumptions are needed.
- `share_basis` should state whether the share comes from a national single product, subnational region count, historical brand count, equal-split assumption, or another auditable rule.

## Usage

- `workflow/lib/build_ap_exposure_index.py` reads `vaccine_formulation_curation.tsv` by default when constructing `aPExposure V2`.
- `workflow/lib/build_ap_exposure_index.py` also reads `vaccine_product_metadata.tsv` by default when constructing role-specific product summaries and `aPExposure V3`.
- `workflow/lib/run_foundation_checks.py` and `modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py` now also read `vaccine_formulation_curation.tsv` by default and can derive the simpler Prn overlay directly from it.
- `manuscript/scripts/review/ms_15_build_selected_country_review_report.py` now also derives the selected-country program-history manifest from this same canonical TSV, so the review report no longer maintains a separate hand-edited formulation history table.
- In practice, `vaccine_formulation_curation.tsv` and `vaccine_product_metadata.tsv` now form the maintained manual curation layer.
- The coverage check still accepts the older two-column overlay format if one is supplied ad hoc, but that format is no longer tracked in this repository.

## Suggested Sources

- National immunization programme schedules or ministry dashboards
- National or regional regulatory product pages and package inserts
- WHO or PAHO programme documentation that identifies vaccine class or introduction milestones
- Product-level composition summaries from official reports

Keep every row auditable with explicit provenance in `source_name`, `source_url`, and `notes`.
