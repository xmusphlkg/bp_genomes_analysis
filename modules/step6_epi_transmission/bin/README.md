# Step6 Script Inventory

## Wrapper-Backed Entry Points

- `step6_06_estimate_reproduction_numbers.py`
- `step6_06_run_re_estimation.sh`
- `step6_07_fit_transmission_models.py`
- `step6_07_run_transmission_models.sh`
- `step6_10_build_focal_country_dynamics.py`

## Module Scripts Retained In Place

- `step6_01_build_country_year_genomic_summaries.py`
- `step6_02_join_public_health.py`
- `step6_03_fit_primary_models.py`
- `step6_04_run_sensitivity_models.py`
- `step6_05_run_amu_exploratory_sensitivity.py`
- `step6_08_cross_validation.py`

These files remain part of step6, but they are not currently documented as standalone top-level run commands.

## Layout Rule

- Wrapper shell scripts stay here when they are part of a documented execution path.
- Test helpers and synthetic fixtures belong in `../dev/`.
- Historical milestone notes should stay out of the release branch unless they are promoted into active documentation.
