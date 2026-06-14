# Audit Ledgers

This directory keeps historical wide TSV ledgers and sensitivity-audit inputs that
support the Communications Biology manuscript package.

These files are not the official numbered Supplementary Tables. The official
reader-facing tables are the 12 TSV files in `manuscript/supplementary/`, named
`Supplementary_Table_1_...tsv` through `Supplementary_Table_12_...tsv`.

If the journal requests additional machine-readable Supplementary Data files,
package the relevant audit ledgers by analysis role and upload them as
Supplementary Data, not as additional Supplementary Tables.

`epydemix_snapshot_manifest.tsv` records the pinned epydemix-data v1.1.0 files
used for focal-country contact priors. The raw CSV snapshot lives under
`modules/public_health/inputs/raw/epydemix-data/v1.1.0/`, which is ignored by
the repository and must be included in the release data archive or restored via
`MS05_EPYDEMIX_SNAPSHOT_DIR` before rebuilding `ms_05_build_focal_country_dynamics.py`.
