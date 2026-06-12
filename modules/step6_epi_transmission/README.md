# Transmission Dynamics Workflow

This guide covers the supported transmission utilities in step6.

## Scope

- Estimate country-year effective reproduction number trajectories.
- Fit transmission models against the estimated trajectories.

## Default Input

The wrapper-backed workflow now starts from the workflow-native country-year panel:

- `workflow/epi/panel_model_country_year_dataset.tsv`

The R_e estimator requires `country`, `year`, and one of `cases`, `reported_cases`, or `case_count`.

## Commands

```bash
bash modules/step6_epi_transmission/bin/run_re_estimation.sh
bash modules/step6_epi_transmission/bin/run_transmission_models.sh
```

You can also pass explicit input and output paths to either wrapper.

## Default Outputs

- `bp_country_year_re_trajectories.tsv`
- `bp_re_summary_statistics.tsv`
- `bp_re_run_metadata.json`
- `bp_transmission_model_results.json`
- `bp_transmission_model_summary.tsv`
- `bp_transmission_model_diagnostics.json`
- `bp_transmission_models_metadata.json`

These files are written under the module data root used by the wrappers.

## Development-Only Helpers

- `tests/dev/generate_synthetic_data.py`
- `tests/dev/test_transmission_dynamics.py`

These files support smoke testing and development only; they are not part of the standard execution path.

## Current Limits

- `run_re_estimation.sh` and `run_transmission_models.sh` are the supported wrapper-backed transmission entry points.
- `step6_08_cross_validation.py` is retained as a utility script, not a standard workflow runner.
- Historical figure code contains a defensive synthetic-data fallback, but manuscript-facing Figure 4 runs should be treated as invalid if that branch is reached. The intended manuscript workflow uses the real step6 and public-health outputs listed above.
