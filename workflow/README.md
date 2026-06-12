# Root Scripts

This directory is reserved for repository-level orchestration and cross-step helpers.

## Active Runners

- `bin/m0_foundation.sh`: builds the unified manifest, runs the readiness checks, and records a version snapshot.
- `bin/m1_m2_qc.sh`: convenience runner for reads planning, assembly QC, and missingness modeling.
- `bin/m3_snippy.sh`: staged Snippy contig-mode bootstrap runner.
- `bin/m4_phylogeny.sh`: recombination filtering and IQ-TREE2 wrapper.
- `bin/m5_asr.sh`: rooted-tree dual-track ASR wrapper.
- `bin/m5_asr_sensitivity.sh`: ASR robustness reruns for support and composition scenarios.
- `bin/run_full_workflow.sh`: canonical M0-M5 repository-level orchestration wrapper.

## Active Cross-Step Helpers

- `lib/build_analysis_manifest.py`: builds `outputs/workflow/manifest/manifest.tsv` from stage outputs.
- `lib/run_foundation_checks.py`: runs reads availability, vaccine-variable coverage, and validation-feasibility checks, writing consolidated reports to `outputs/workflow/checkpoints/`.
- `lib/aggregate_assembly_qc.py`: merges per-sample QC reports for the root Snakemake contract.
- `lib/missingness_model.py`: workflow-native missingness diagnostics.
- `lib/build_ap_exposure_index.py`: ecology exposure-index builder used by `rules/epi_models.smk`.
- `lib/ipw_prevalence.py`: IPW prevalence layer used by `rules/epi_models.smk`.
- `lib/panel_model.py`: workflow-native ecology model runner.
- `lib/build_programme_country_period_panel.py`: programme-surveillance country-period panel builder retained as the ecology sidecar bridge.
- `lib/run_programme_surveillance_models.py`: programme-surveillance grouped-binomial model runner for the sidecar suite.
- `lib/run_programme_two_stage_uncertainty.py`: bootstrap uncertainty propagation for the programme-surveillance sidecar.
- `lib/its_feasibility.py`: interrupted time-series feasibility helper.
- `lib/mask_recombination.py`, `lib/compare_trees.py`, `lib/filter_alignment_by_missingness.py`, `lib/extract_iqtree_composition_report.py`: M4 phylogeny helpers.
- `lib/root_tree_on_tip.py`, `lib/asr_parsimony.py`, `lib/asr_pastml.py`, `lib/origin_events.py`, `lib/prune_tree_by_tips.py`: M5 ASR helpers.
- `manuscript/scripts/freeze/generate_supplementary_table_1.py`: manuscript-facing Supplementary Table 1 exporter.

## Parallel Runtime Notes

- `lib/panel_model.py` supports `--max-workers` or `PANEL_MAX_WORKERS` for leave-one-country-out fits.
- `lib/run_programme_surveillance_models.py` supports `--max-workers` or `PROGRAMME_MODELS_MAX_WORKERS` for sidecar leave-one-country-out fits.
- `lib/run_programme_two_stage_uncertainty.py` supports `--max-workers` or `PROGRAMME_UNCERTAINTY_MAX_WORKERS` for bootstrap replicates.
- `manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py` supports `--max-workers` or `MS05_MAX_WORKERS` for focal-country mechanistic multi-start fitting.

## Maintenance Helpers

- `bin/find-big-git-objects.sh`: repository maintenance helper for locating large git objects.

## Layout Rule

If a script is tied to a workflow stage, it should live under that stage's own `bin/` directory. If a helper is no longer part of routine project execution, remove it from the release branch rather than leaving it alongside active entry points.

## Environment Notes

- Runtime env paths are configured in `config/runtime/runtime_envs.env` plus optional untracked `config/runtime/runtime_envs.local.env`.
- Use `bash bin/bootstrap_runtime_envs.sh --check` to validate the configured prefixes.
- Use `bash bin/run_with_project_env.sh --script <path>` when you want the launcher to infer the correct env from the script header.
- Current Conda on this host is configured with a broken `libmamba` solver path.
- Use `CONDA_NO_PLUGINS=true CONDA_SOLVER=classic` for environment creation or `conda run` invocations tied to the M3/M4 workflow.
- M3 Snippy wrappers can fall back to Docker image `quay.io/biocontainers/snippy:4.6.0--hdfd78af_6` if Docker is available and the Conda env is absent.
- M4 can fall back to Docker image `quay.io/biocontainers/gubbins:3.4.3--py310hfc0ef84_1` for `run_gubbins.py` while using system `iqtree2` or the `pertussis-prn-global-bio` Conda env.
- On this host `pertussis_data` resolves to a shared external data mount; the M3 Docker wrappers auto-mount the resolved target so Snippy can see the reference and assembly FASTA files.
