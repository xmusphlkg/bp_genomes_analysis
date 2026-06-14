# Model Audit Remediation Todos

## Purpose

This tracker converts the model-code audit findings into executable remediation work.
The default policy is:

- fail closed when a manuscript-facing method is not statistically supportable
- prefer removing invalid certainty over preserving legacy outputs
- require an explicit development-only flag for smoke-test or synthetic branches

## Critical

- [x] Step6 `R_e` / transmission path
  - Files: `modules/step6_epi_transmission/bin/run_re_estimation.sh`, `step6_06_estimate_reproduction_numbers_v2.py`, `step6_07_fit_transmission_models.py`, `run_transmission_models.sh`, `modules/step6_epi_transmission/README.md`
  - Fix: refuse annual country-year case panels for renewal-based `R_e` unless an explicit development override is set; propagate metadata so downstream transmission models reject unsupported `R_e` inputs.
  - Verify: tests cover fail-closed behavior on annual inputs and still allow synthetic/dev fixtures with an explicit override.

- [x] Workflow ecology / programme IPW inference
  - Files: `workflow/lib/panel_model.py`, `run_programme_surveillance_models.py`, `run_programme_two_stage_uncertainty.py`
  - Fix: stop treating IPW pseudo-denominators as binomial trial counts; use a safer fractional-response / two-stage uncertainty workflow and manuscript-facing notes that reflect the new inferential scale.
  - Verify: regression tests cover serial/parallel parity and no longer rely on pseudo-binomial weights.

- [x] Full mechanistic selected-country model
  - Files: `manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py`
  - Fix: jointly estimate the main temporal basis coefficients instead of freezing them at least-squares initials; propagate full-parameter uncertainty into PPC and recovery.
  - Verify: tests cover optimizer task shapes, PPC perturbation scope, and recovery metrics on small fixtures.

- [x] ASR stochastic mapping defaults
  - Files: `manuscript/scripts/sidecars/ms_21_run_asr_stochastic_mapping.R`
  - Fix: default to empirically estimated transition rates, require an explicit smoke-test override, and validate resampling scenario registries before use.
  - Verify: smoke/full mode behavior is explicit in metadata and cannot silently overwrite manuscript-facing outputs.

## High

- [x] Step6 primary/sensitivity ecology models
  - Files: `modules/step6_epi_transmission/bin/step6_03_fit_primary_models.py`, `step6_04_run_sensitivity_models.py`
  - Fix: move from naive grouped-binomial independence to clustered/robust country-level inference and tighten output notes accordingly.
  - Verify: coefficient tables still render, and diagnostics report covariance type plus cluster counts.

- [x] Step6 mixed/clustered robustness script
  - Files: `modules/step6_epi_transmission/bin/step6_09_mixed_effects_models.py`
  - Fix: remove invalid grouped-binomial-plus-`var_weights` combinations, correct method labels, and keep only statistically coherent robustness paths.
  - Verify: script outputs no longer advertise HC3 when it is not used and do not normalize grouped-binomial likelihood weights arbitrarily.

- [x] Selected-country evidence summary pseudo-replication
  - Files: `manuscript/scripts/diagnostics/ms_12_build_submission_evidence_summary.py`
  - Fix: keep primary ecology robustness at the grouped country-year scale, relabel it as frequentist GLM diagnostics rather than posterior or mixed-effect output, and fail/diagnostic-label unstable rows.
  - Verify: no pseudo-row expansion remains in primary ecology robustness outputs; generated tables use `estimate`/`standard_error` and diagnostic status fields.

- [x] Selected-country DR / AIPW sidecar
  - Files: `manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py`
  - Fix: use out-of-fold nuisance predictions where feasible, relabel exploratory estimators clearly, and fail-close augmented estimates that leave the probability scale.
  - Verify: tests cover bounded AIPW helper behavior and generated missingness summaries carry estimator-status fields.

## Medium

- [x] Sensitivity-script provenance
  - Files: `modules/step6_epi_transmission/bin/step6_04_run_sensitivity_models.py`
  - Fix: separate immutable analysis-input construction from model-fitting output writes.
  - Verify: sensitivity runs cannot silently rebuild manuscript inputs without explicit invocation.

- [x] Missingness-model optimism
  - Files: `workflow/lib/missingness_model.py`, `manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py`
  - Fix: distinguish training diagnostics from out-of-fold calibration and ensure downstream summaries do not present in-sample metrics as generalization evidence.
  - Verify: tests still pass and diagnostics tables label metric provenance.

- [x] Two-stage uncertainty summary construction
  - Files: `workflow/lib/run_programme_two_stage_uncertainty.py`
  - Fix: replace ad hoc min/max CI fusion with a bootstrap-native interval summary and guard empty replicate tables.
  - Verify: empty or partially failed bootstrap runs write coherent diagnostic outputs instead of crashing.

- [x] aP exposure index transparency
  - Files: `workflow/lib/build_ap_exposure_index.py`
  - Fix: tighten availability flags, clarify arbitrary scoring components, and surface parameterization coverage in outputs.
  - Verify: V3 availability matches all required metadata components.

- [x] ASR representativeness resampling metadata
  - Files: `workflow/lib/run_m5_asr_resampling.py`
  - Fix: stop depending on archived inventory metadata and make stratified resampling explicitly stress-test-only.
  - Verify: active manifests drive block assignments.

- [x] Step6 CV utility
  - Files: `modules/step6_epi_transmission/bin/step6_08_cross_validation.py`
  - Fix: remove silent synthetic fallback from manuscript-facing behavior, fix TSV/CSV expectations, and align prediction/actual vectors before scoring.
  - Verify: fixture tests cover missing-data failure and alignment behavior.

- [x] Exploratory AMU sensitivity guardrails
  - Files: `modules/step6_epi_transmission/bin/step6_05_run_amu_exploratory_sensitivity.py`
  - Fix: do not fit standard GLMs below the exploratory sample threshold; emit explicit not-fit status instead.
  - Verify: low-overlap fixtures never produce misleading standard-fit diagnostics.

## Low / Reporting

- [x] Multiple-testing labels
  - Files: `workflow/lib/panel_model.py`, `run_programme_surveillance_models.py`, `modules/step6_epi_transmission/bin/step6_03_fit_primary_models.py`
  - Fix: avoid implying manuscript-wide FDR control when adjustment is only local to a term family.

- [x] Remote fallback reproducibility
  - Files: `manuscript/scripts/sidecars/ms_05_build_focal_country_dynamics.py`
  - Fix: prefer pinned local snapshots for contact/population priors and label any network fallback as non-canonical.

- [x] Study-dependence audit wording
  - Files: `manuscript/scripts/diagnostics/ms_18_build_study_dependence_audit.py`
  - Fix: keep heuristic calls clearly labeled as heuristic and not formal inferential tests.

## Completion Notes

- Update this file as each block lands.
- A block is only complete once code, tests, and output/metadata language all agree.
- Landed so far:
  - `R_e` / transmission now fail closed on unsupported annual disaggregation unless an explicit development override is set.
  - IPW programme models no longer treat pseudo-denominators as binomial trial counts.
  - `ms_05` jointly estimates basis and auxiliary parameters and propagates full-parameter uncertainty.
  - `ms_21` now defaults to empirical ARD stochastic mapping, requires explicit smoke-mode opt-in, and writes smoke outputs to dev-only paths.
  - `ms_16` now fail-closes out-of-range AIPW prevalence estimates and exposes raw augmented estimates only as diagnostics.
  - `ms_12` ecology robustness output now uses frequentist estimate fields and diagnostic/failure status labels.
  - Step6 ecology / sensitivity outputs now use country-cluster-robust covariance with HC1 fallback and explicit covariance metadata in notes.
  - Invalid weighted grouped-binomial robustness variants were removed from Step6 mixed-effects outputs.
  - Step6 CV no longer silently fabricates synthetic data or mismatches prediction/actual vectors.
  - AMU exploratory standard GLM fitting is now skipped below the minimum exploratory overlap threshold.
  - Two-stage uncertainty summaries now report bootstrap-native propagated intervals and tolerate empty replicate tables.
  - Missingness diagnostics now prefer out-of-fold probabilities when available and label training-scope metrics explicitly.
  - aP exposure V2/V3 availability now requires complete required components and reports parameter-grid/provenance fields.
  - ASR resampling outputs now carry stress-test-only inference scope and active-manifest block-assignment provenance through figure data and audit ledgers.
  - Model result tables now carry structured `q_value_scope` fields for local BH adjustment scope, including workflow panel, programme, Step6 primary, and Step6 sensitivity outputs.
  - `ms_05` contact priors now prefer a pinned local `epydemix-data v1.1.0` snapshot, record `source_access_mode` / `source_canonicality`, and require explicit opt-in for non-canonical network fallback.
  - Study-dependence outputs now carry `diagnostic_inference_scope` fields, and heuristic bootstrap/permutation wording no longer presents sensitivity screens as formal hypothesis tests.
