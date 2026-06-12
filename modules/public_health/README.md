# Public Health Integration

This module ingests, normalizes, and versions country-year public-health variables used to contextualize the *Bordetella pertussis* genomic analyses.

The goal is not to replace pathogen genomics with ecological correlations. The goal is to build a transparent, auditable country-year context layer that can be joined to genomic summaries in later steps.

## Scope

This module will manage:

- Official pertussis case counts and incidence
- Diagnosis/reporting-era milestones for surveillance comparability
- Vaccine coverage indicators
- Vaccine program descriptors
- Antimicrobial-use indicators
- Country-year genomic summary joins

## Core Principles

- Prefer official sources over secondary summaries
- Freeze all source data by release date
- Preserve raw values separately from normalized analytic values
- Keep country normalization explicit and reversible
- Record manual curation in plain-text sidecar files
- Separate provenance into planned inventory, canonical registry, and row-level citation map

## Main Inputs

- WHO Immunization Data portal exports
- CDC pertussis surveillance tables and provisional reports
- ECDC pertussis country-year data or report-derived values
- WHO/UNICEF WUENIC indicators
- WHO GLASS antimicrobial-use exports
- ECDC ESAC-Net exports for Europe
- Internal genomic country-year summaries from `bp_step6`
- Curated diagnosis/reporting-era milestone tables for focal surveillance systems

## Manual Download Checkpoints

These are the places where user-side manual download or export may be needed before an ingestion task can run cleanly.

- `PH-03` / WHO reported cases:
  expected input is a frozen official WHO Immunization Data portal export.
- `PH-04` / WUENIC coverage:
  expected input is a frozen official WUENIC export from the WHO Immunization Data portal.
- `PH-06` / WHO GLASS AMU:
  expected input is a frozen official GLASS AMU export or equivalent official release table.
- `PH-07` / ECDC ESAC-Net:
  expected input is a frozen official ESAC-Net export for Europe.
- `PH-05` / vaccine-program metadata:
  this may require manual extraction or curation from official schedule documents even if no single export table exists.
- `PH-11` / diagnosis and reporting era indicators:
  this requires manual curation from official surveillance guidance, case-definition pages, or ministry/public-health agency summaries.

Reminder for future task runs:

- If the frozen official source snapshot is already on disk, proceed with the script task.
- If it is not on disk yet, stop and ask the user to download or export it first.
- Record the access date, export date, release date, and any manual curation notes alongside the cleaned output.

## Planned Workflow

1. Fetch or manually export official source tables.
2. Store immutable source snapshots under a dated acquisition convention.
3. Normalize country names and map to ISO3.
4. Standardize year fields and incidence units.
5. Create separate clean tables for:
   - country-year surveillance values
   - diagnosis/reporting-era indicators
   - country-program metadata
   - antimicrobial-use variables
6. Merge clean tables into a master country-year context table.
7. Join the master context table to genomic country-year summaries in `bp_step6`.

## Expected Directory Layout

- `bin/`: ingestion and normalization entry points
- `outputs/`: clean analysis-ready tables and schema templates

## Planned Deliverables

- `outputs/ph_source_inventory.tsv`
- `outputs/ph_source_registry.tsv`
- `outputs/ph_source_citation_map.tsv`
- `outputs/ph_country_year_master.tsv`
- `outputs/ph_reporting_era_indicators.tsv`
- `outputs/ph_reporting_era_coverage_audit.tsv`
- `outputs/ph_reporting_era_resolution_worklist.tsv`
- `outputs/ph_country_program_metadata.tsv`
- `outputs/ph_country_year_master.schema.tsv`
- `outputs/ph_reporting_era_indicators.schema.tsv`
- `outputs/ph_country_program_metadata.schema.tsv`

## Recommended First Entrypoints

These entry points are the active implementations, and the following order keeps refreshes predictable:

- `bin/ph_01_build_source_inventory.py`
- `bin/ph_02_normalize_country_names.py`
- `bin/ph_03_clean_who_cases.py`
- `bin/ph_04_clean_wuenic.py`
- `bin/ph_05_clean_vaccine_programs.py`
- `bin/ph_06_clean_glass_amu.py`
- `bin/ph_07_clean_esacnet_amu.py`
- `bin/ph_08_build_country_year_master.py`
- `bin/ph_10_clean_highres_cases.py`
- `bin/ph_11_clean_reporting_era_indicators.py`
- `bin/ph_12_audit_reporting_era_coverage.py`
- `bin/ph_13_build_reporting_era_resolution_worklist.py`

## Joining Rules

The canonical join key for analytic outputs is:

- `country_iso3`
- `year`

Program metadata may require interval joins using:

- `country_iso3`
- `year_start`
- `year_end`

## Important Modeling Note

Public-health variables in this repository are intended for ecological association analyses. They must not be used to claim patient-level or isolate-level causation.
