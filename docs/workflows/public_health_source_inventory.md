# Public Health Source Inventory

This file records the official data sources planned for the country-year context layer.
It now also carries the human-facing public-health variable dictionary and the project-wide source freeze ledger, so related documentation does not need to be split across separate TSV files under `docs/`.

The machine-readable provenance outputs generated from this design now live in:

- `modules/public_health/outputs/ph_source_inventory.tsv`
- `modules/public_health/outputs/ph_source_registry.tsv`
- `modules/public_health/outputs/ph_source_citation_map.tsv`
- `modules/public_health/outputs/ph_reporting_era_resolution_worklist.tsv`

## Required Sources

| Domain | Source | URL | Expected Unit | Planned Freeze Metadata |
| --- | --- | --- | --- | --- |
| Reported cases | WHO Immunization Data portal | <https://immunizationdata.who.int/dashboard> | Country-year reported counts | export date, portal release note |
| Vaccine coverage | WHO/UNICEF WUENIC | <https://immunizationdata.who.int/dashboard> | Country-year percent coverage | export date, WUENIC release year |
| US surveillance detail | CDC pertussis surveillance | <https://www.cdc.gov/pertussis/php/surveillance/index.html> | US annual counts and updates | access date |
| US 2024 update | CDC provisional report | <https://www.cdc.gov/pertussis/media/pdfs/2025/01/pertuss-surv-report-2024_PROVISIONAL-508.pdf> | national and state counts | publication date and access date |
| EU surveillance detail | ECDC pertussis annual reports | <https://www.ecdc.europa.eu/en/publications-data/pertussis-annual-epidemiological-report-2023> | country-year cases or rates | publication date and access date |
| Antimicrobial use | WHO GLASS AMU | <https://www.who.int/news/item/25-09-2025-updated-who-dashboard-offers-new-insights-on-antimicrobial-resistance-and-use> | country-year DDD-based metrics | dashboard export date |
| European antimicrobial use | ECDC ESAC-Net | <https://www.ecdc.europa.eu/en/about-us/partnerships-and-networks/disease-and-laboratory-networks/esac-net> | country-year DDD-based metrics | export date |
| Bordetella typing framework | BIGSdb-Pasteur Bordetella | <https://bigsdb.pasteur.fr/bordetella/> | allele and lineage metadata | access date |

## Manual Curation Targets

The following variables may require manual curation from official documents or national schedules:

- `vaccine_program_type`
- `acellular_vs_whole_cell`
- `prn_in_vaccine`
- `booster_schedule`
- `program_change_year`

## Freeze Rules

- Record the access or export date for every table.
- Keep raw source snapshots immutable after freezing.
- If a source is updated later, create a new dated snapshot rather than overwriting the old one.
- Any manually curated value must cite the source document in a sidecar file.

## Priority Order for Conflicts

1. Official downloadable tables
2. Official dashboard exports
3. Official PDF reports
4. Official narrative pages
5. Literature only when official data are unavailable

## Row-Level Link Storage

The project-level inventory above is not the only place where source links are stored. The row-level curation files below now carry most of the country-specific links used in the current public-health layer:

- `modules/public_health/inputs/raw/report_cases/pertussis_diagnosis_reporting_era_indicators.csv`
  Columns: `primary_source_url`, `secondary_source_url`
  Scope: country or regional diagnosis/reporting-era milestones used by `PH-11`
- `modules/public_health/inputs/raw/vaccine_program_docs/vaccine_program.csv`
  Columns: `VaccinePregnantSource`, `VaccinePregnantIntroSource`
  Scope: maternal pertussis recommendation timing and introduction-year curation
- `modules/public_health/inputs/raw/vaccine_program_docs/source_meta.tsv`
  Column: `source_url`
  Scope: file-level provenance for vaccine-program raw inputs
- `modules/public_health/inputs/raw/report_cases/source_meta.tsv`
  Column: `source_url`
  Scope: file-level provenance for reporting-era and focal high-resolution case workbooks
- `modules/public_health/inputs/curation/PRN vaccine status.md`
  Final table column: `Official label / program source (primary)`
  Scope: product-level PRN status, label links, and program links used in formulation curation

The normalized reporting-era sidecar `modules/public_health/outputs/ph_reporting_era_indicators.tsv` also preserves the row-level `primary_source_url` and `secondary_source_url` fields after cleaning.

`modules/public_health/outputs/ph_reporting_era_resolution_worklist.tsv` now joins the cleaned country coverage status back to raw-country and proxy-source context so blocked countries can be resolved in batch rather than by one-off manual re-discovery.

## Provenance Layers

The public-health provenance design is now split deliberately into three layers:

1. Inventory
   `docs/public_health_source_inventory.md` and `modules/public_health/outputs/ph_source_inventory.tsv`
   Purpose: planned source classes and freeze expectations
2. Registry
   `modules/public_health/inputs/curation/public_health_source_registry.tsv` and `modules/public_health/outputs/ph_source_registry.tsv`
   Purpose: canonical per-document source registry used to validate row-level citations
3. Citation map
   `modules/public_health/outputs/ph_source_citation_map.tsv`
   Purpose: auditable linkage from each curated row or raw-sidecar record to the canonical registry

## Public Health Variable Dictionary

### Shared Join And Audit Fields

| Variable | Dataset | Level | Type | Required | Encoding | Source Priority | Description | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `country_iso3` | `all` | country-year | string | yes | ISO3 | derived | Canonical country key used for joins | Normalize before merge |
| `country_name` | `all` | country-year | string | yes | free text | official | Human-readable country name | Retain original source name separately if needed |
| `year` | `all` | country-year | integer | yes | Gregorian year | official | Observation year | No partial-year encoding in final tables |
| `source_release_date` | `all` | record | date | yes | `YYYY-MM-DD` | official | Release or publication date of source | Versioning field |
| `data_freeze_date` | `all` | record | date | yes | `YYYY-MM-DD` | derived | Date when the repository froze the source snapshot | Do not overwrite between freezes |
| `notes` | `all` | record | string | no | free text | derived | Free-text curation or caveat field | Keep short and factual |

### Country-Year Public Health Layer

| Variable | Dataset | Level | Type | Required | Encoding | Source Priority | Description | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `reported_cases` | `ph_country_year_master` | country-year | integer | yes | case count | WHO/CDC/ECDC | Official reported pertussis case count | Prefer official yearly totals |
| `incidence_per_100k` | `ph_country_year_master` | country-year | float | yes | cases per 100000 | WHO/CDC/ECDC | Standardized incidence value | Convert from per million when needed |
| `dtp3_coverage` | `ph_country_year_master` | country-year | float | yes | percent | WUENIC | DTP3 coverage indicator | Store as numeric percent not fraction |
| `booster_coverage` | `ph_country_year_master` | country-year | float | no | percent | WUENIC or official national source | Booster coverage when available | May be sparse |
| `post_covid_period` | `ph_country_year_master` | country-year | integer | yes | `0/1` | derived | Indicator for post-May-2023 analytic period | Define exactly in codebook |
| `macrolide_use_ddd_per_1000_per_day` | `ph_country_year_master` | country-year | float | no | DDD per 1000 inhabitants per day | GLASS or ESAC-Net | Country-level macrolide use metric | Keep source-specific origin |
| `total_antibiotic_use_ddd_per_1000_per_day` | `ph_country_year_master` | country-year | float | no | DDD per 1000 inhabitants per day | GLASS or ESAC-Net | Country-level total systemic antibiotic use metric | Secondary covariate |
| `genomes_count` | `ph_country_year_master` | country-year | integer | yes | genome count | internal genomic summary | Number of genomes retained in that country-year | Derived from genomic pipeline |
| `genomes_per_case` | `ph_country_year_master` | country-year | float | yes | ratio | derived | Genomic sampling density relative to reported burden | Key bias metric |
| `surveillance_source` | `ph_country_year_master` | country-year | string | yes | free text | official | Name of surveillance source used for cases | Audit column |

### Country-Program Metadata Layer

| Variable | Dataset | Level | Type | Required | Encoding | Source Priority | Description | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `vaccine_program_type` | `ph_country_program_metadata` | country-period | string | yes | categorical | official program documents | High-level program label | Examples: `wP-dominant`, `aP-dominant`, `mixed` |
| `acellular_vs_whole_cell` | `ph_country_program_metadata` | country-period | string | yes | categorical | official program documents | Whether routine program is primarily aP or wP | Allow `mixed` and `unknown` |
| `prn_in_vaccine` | `ph_country_program_metadata` | country-period | string | yes | `yes/no/partial/unknown` | official program documents | Whether PRN is included in routine vaccine products | Use `unknown` if unresolved |
| `booster_schedule` | `ph_country_program_metadata` | country-period | string | no | free text | official program documents | Short description of booster timing | Keep normalized summary plus note field |
| `program_change_year` | `ph_country_program_metadata` | country-period | integer | no | year | official program documents | Year of relevant program change | Can be interval start |

### Genomics-Bridge Analysis Layer

| Variable | Dataset | Level | Type | Required | Encoding | Source Priority | Description | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `n_genomes_prn_interpretable` | `bp_country_year_analysis_input` | country-year | integer | yes | genome count | internal genomic summary | Genomes with interpretable `prn` status | Used as denominator |
| `n_prn_disrupted` | `bp_country_year_analysis_input` | country-year | integer | yes | genome count | internal genomic summary | Number of disrupted `prn` genomes | Response count |
| `frac_prn_disrupted` | `bp_country_year_analysis_input` | country-year | float | yes | fraction | internal genomic summary | Proportion disrupted among interpretable genomes | Primary response |
| `n_read_supported_prn_disrupted` | `bp_country_year_analysis_input` | country-year | integer | no | genome count | internal genomic summary | Number of read-supported disrupted events | Validation-aware summary |
| `n_mr_marked` | `bp_country_year_analysis_input` | country-year | integer | no | genome count | internal genomic summary | Number of genomes carrying the 23S resistance marker | Secondary descriptive variable |
| `frac_23s_A2047G` | `bp_country_year_analysis_input` | country-year | float | no | fraction | internal genomic summary | Fraction carrying the 23S `A2047G` marker | Secondary descriptive variable |

## Project Source Freeze Ledger

This ledger keeps the human-facing freeze rules for both public-health sources and project-wide external repositories.

| Source | URL | Freeze Date | Scope | Access Method | Update Policy | Owner Module | Domain | Release Note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `NCBI Datasets` | <https://www.ncbi.nlm.nih.gov/datasets/> | `2026-03-21` | `project_source_ledger_snapshot` | datasets metadata refresh or manifest export | record accession snapshot date at ingest and create a new dated snapshot for each metadata refresh | `bp_step1` | assemblies | required in the technical spec for genome acquisition and metadata backbone |
| `GenBank` | <https://www.ncbi.nlm.nih.gov/genbank/> | `2026-03-21` | `project_source_ledger_snapshot` | linked assembly metadata via NCBI accessions | record accession snapshot date at ingest and create a new dated snapshot for each metadata refresh | `bp_step1` | assemblies | required in the technical spec as part of the public assembly backbone |
| `RefSeq` | <https://www.ncbi.nlm.nih.gov/refseq/> | `2026-03-21` | `project_source_ledger_snapshot` | linked assembly metadata via NCBI accessions | record accession snapshot date at ingest and create a new dated snapshot for each metadata refresh | `bp_step1` | assemblies | required in the technical spec as part of the public assembly backbone |
| `NCBI SRA` | <https://www.ncbi.nlm.nih.gov/sra> | `2026-03-21` | `project_source_ledger_snapshot` | run metadata lookup and accession recovery | record accession snapshot date at ingest and create a new dated snapshot for each metadata refresh | `bp_step1` | raw_reads | required where available for read-backed validation of `prn` events |
| `ENA` | <https://www.ebi.ac.uk/ena/browser/home> | `2026-03-21` | `project_source_ledger_snapshot` | run metadata lookup and accession recovery | record accession snapshot date at ingest and create a new dated snapshot for each metadata refresh | `bp_step1` | raw_reads | required where available for read-backed validation of `prn` events |
| `BIGSdb-Pasteur Bordetella` | <https://bigsdb.pasteur.fr/bordetella/> | `2026-03-21` | `project_source_ledger_snapshot` | manual lookup or exported allele and lineage tables | record access date or export date at ingest and create a new dated snapshot when allele resources change | `bp_step1` | typing_framework | required allele naming and lineage harmonization resource |
| `Public supplementary tables and accession lists` | `TBD_MULTIPLE_PUBLICATION_URLS` | `2026-03-21` | `project_source_ledger_snapshot` | manual extraction from publication supplements and accession lists | record the exact publication URL or DOI in a sidecar and create a new dated snapshot for each added supplement | `bp_step1` | literature_supplements | recommended for recovering newly reported strains and metadata |
| `WHO Immunization Data portal` | <https://immunizationdata.who.int/dashboard> | `2026-03-21` | `project_source_ledger_snapshot` | dashboard export or manual download | freeze by export date and portal release note and create a new dated snapshot if the portal content changes | `public_health` | reported_cases | technical spec notes the portal currently exposes reported VPD case trends through 2024 |
| `WHO/UNICEF WUENIC` | <https://immunizationdata.who.int/dashboard> | `2026-03-21` | `project_source_ledger_snapshot` | dashboard export or manual download | freeze by export date and WUENIC release year and create a new dated snapshot for each new WUENIC release | `public_health` | vaccine_coverage | technical spec notes the 2024 WUENIC release was published in 2025 |
| `CDC Pertussis Surveillance` | <https://www.cdc.gov/pertussis/php/surveillance/index.html> | `2026-03-21` | `project_source_ledger_snapshot` | manual extraction from official surveillance page | freeze by access date and create a new dated snapshot when CDC updates annual or provisional counts | `public_health` | us_surveillance_detail | required in the technical spec for US detail |
| `CDC provisional pertussis report` | <https://www.cdc.gov/pertussis/media/pdfs/2025/01/pertuss-surv-report-2024_PROVISIONAL-508.pdf> | `2026-03-21` | `project_source_ledger_snapshot` | manual extraction from official PDF | freeze by publication date and access date and create a new dated snapshot when CDC issues a superseding provisional or final report | `public_health` | us_surveillance_detail | technical spec identifies this PDF as the official one-page summary for US 2024 counts |
| `ECDC pertussis annual reports or data portals` | <https://www.ecdc.europa.eu/en/publications-data/pertussis-annual-epidemiological-report-2023> | `2026-03-21` | `project_source_ledger_snapshot` | manual extraction from annual reports or portal tables | freeze by publication date and access date and create a new dated snapshot when ECDC publishes a new report or portal update | `public_health` | eu_surveillance_detail | technical spec notes the referenced annual report was published in 2025 using data for 2024 retrieved on 2025-03-12 |
| `WHO portal and national schedules` | `TBD_MULTIPLE_OFFICIAL_URLS` | `2026-03-21` | `project_source_ledger_snapshot` | manual curation from WHO portal pages and official national immunization schedules | record exact cited URLs in a sidecar and create a new dated snapshot when program guidance changes | `public_health` | vaccine_program_characteristics | required for vaccine program type, acellular versus whole-cell usage, PRN inclusion, booster schedule, and program change year |
| `WHO GLASS AMU dashboard` | <https://www.who.int/news/item/25-09-2025-updated-who-dashboard-offers-new-insights-on-antimicrobial-resistance-and-use> | `2026-03-21` | `project_source_ledger_snapshot` | dashboard export or manual download | freeze by dashboard export date and create a new dated snapshot for each GLASS AMU update | `public_health` | antimicrobial_use | technical spec notes the September 25 2025 update reports validated AMU data through 2023 |
| `ECDC ESAC-Net` | <https://www.ecdc.europa.eu/en/about-us/partnerships-and-networks/disease-and-laboratory-networks/esac-net> | `2026-03-21` | `project_source_ledger_snapshot` | portal export or manual download | freeze by export date and create a new dated snapshot when ESAC-Net releases updated consumption tables | `public_health` | antimicrobial_use | recommended higher-resolution antimicrobial consumption source for Europe |
