# bp_step5

Step 5 places `prn` states into a global phylogenetic framework and estimates repeated emergence, clade structure, and temporal spread.

## Purpose

This module transforms the study from a descriptive frequency analysis into an evolutionary analysis.

## Inputs

- `../modules/step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv`
- Harmonized marker and metadata tables from `bp_step2`
- Sample manifests defining the balanced and full phylogeny cohorts

## Main Questions

- How many times has `prn` disruption likely emerged independently?
- Which mechanism classes are phylogenetically clustered?
- Are disrupted clades geographically restricted or widely disseminated?
- Do disrupted clades expand differently from intact sister groups?

## Planned Analyses

- Balanced global phylogeny
- Full-data sensitivity phylogeny
- Ancestral-state reconstruction
- Independent-origin event counting
- Clade-level time summaries

## Planned Outputs

- `outputs/bp_phylogeny_manifest_balanced.tsv`
- `outputs/bp_phylogeny_manifest_full.tsv`
- `outputs/bp_global_phylogeny.nwk`
- `outputs/bp_prn_ancestral_states.tsv`
- `outputs/bp_prn_independent_origins.tsv`
- `outputs/bp_prn_clade_summary.tsv`

Schema templates are already provided in `outputs/`.

## Recommended Script Plan

- `bin/step5_01_build_phylogeny_manifest.py`
- `bin/step5_02_build_global_phylogeny.sh`
- `bin/step5_03_reconstruct_prn_states.py`
- `bin/step5_04_count_independent_origins.py`
- `bin/step5_05_summarize_clades.py`
- `bin/step5_06_build_missingness_model.py`

## Deliverable Standard

Main findings from this module should be reproducible under both:

- a balanced country-year-aware cohort, and
- a larger full-data sensitivity cohort.

`PHY-01` now materializes both manifest layers with explicit audit fields describing:

- manifest type and tree role,
- country-year rank and cap,
- deterministic selection rule,
- per-row inclusion reason for the balanced main tree versus the full sensitivity tree.

`PHY-06` now provides a cross-step missingness model used to support ASR sensitivity work. Its current outputs are written under `outputs/workflow/missingness_model/`.
