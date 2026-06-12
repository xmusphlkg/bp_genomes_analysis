# Selected-Country Rebuild Guide

This file is the active human-facing guide for the selected-country rebuild package. It replaces the older multi-file planning scaffold and should be enough to understand the current selected-country design without opening legacy design memos.
This directory should normally contain only this `README.md`.

## Current Role In The Project

The manuscript now uses a structural-evolution main frame. Selected-country curation still matters, but it now supports the programme-amplification and boundary-condition layers rather than leading the paper.

Current high-level country roles:

- Stage 1 primary epidemiologic set: `USA`, `NZL`, `AUS`, `GBR`, `JPN`
- Stage 2 mechanistically triangulated subset: `USA`, `NZL`, `AUS`, `JPN`
- context-only countries under the default screen: `CHN`, `FRA`, `BRA`, `CZE`

## What This Package Needs To Explain

- how country-program epochs are defined from the shared public-health curation
- why the primary epidemiologic set is asymmetric rather than a balanced country panel
- how Stage 1 eligibility is decoupled from Stage 2 mechanistic triangulation
- how Figure 4 uses country-program histories as heterogeneous amplification environments
- how Extended Data preserves Stage 1 / Stage 2 screening, missingness and block-dependence audits
- what remains bounded because of missingness, sparse epochs, or limited local rooted support

## Current Comparative Backbone

Figure-level logic:

- Figure 1 defines the public-genome atlas and recoverable-locus boundary
- Figure 2 makes the constrained *prn* structural solution space the core claim
- Figure 3 shows repeated origins across the rooted genome tree
- Figure 4 uses country-program histories as amplification environments
- Figure 5 bounds comparator, missingness, study-block and validation alternatives
- Figure 5 closes with the validation, constrained-evolution and identifiability synthesis

Results-level logic:

- recoverable-locus prevalence is the estimand, not population truth
- country-program epochs are amplification environments and audit strata, not the first narrative entry point
- `USA` and `NZL` provide the clearest positive PRN-exposed contrasts
- `AUS` is informative but directionally discordant and bounded
- `JPN` remains triangulated but bounded PRN-free context
- `GBR` is the clearest Stage-1-primary but non-triangulated comparator
- repeated emergence and constrained structural reuse remain the mechanism bridge
- local rooted reruns are retained as bounded audit context rather than as standalone decisive anchors
- multiverse readiness screening remains in the audit atlas, separating stage-1 epidemiologic eligibility from stage-2 mechanistic triangulation

## Curation And Methods Boundary

Country-program history is now generated from the canonical shared curation in:

- `modules/public_health/inputs/curation/vaccine_formulation_curation.tsv`

Source hierarchy:

- official national programme-history or immunisation guidance
- official product information for PRN content when needed
- country-specific peer-reviewed programme-history papers when they provide the clearest breakpoint

## Regeneration

Refresh the selected-country package with:

```bash
python manuscript/scripts/review/ms_15_build_selected_country_review_report.py
python manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py
```

These steps regenerate the figure-data tables used by the current selected-country comparison framework together with the multiverse screening, missingness, and ASR-audit sidecars cited in the revised manuscript.

## Practical Rule

If someone only reads one file in this directory, it should be this `README.md`.
