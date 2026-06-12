# `prn` Call Confidence Rules

## Purpose

`PRN-04` defines deterministic confidence labels for the assembly-based `prn` mechanism calls produced by `PRN-02`. The labels are written into `modules/step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv` as `prn_call_confidence` and are intentionally separate from later raw-read validation fields.

## Confidence Labels

- `assembly_high`: strong assembly-only evidence for the assigned `prn` state or mechanism.
- `assembly_moderate`: coherent structural evidence is present, but one supporting channel is weaker or non-orthogonal.
- `assembly_low`: disruption is plausible, but the assembly context or mechanism support is weak enough that downstream summaries should keep it separate from stronger calls.
- `insufficient_evidence`: the current pipeline does not contain enough locus evidence to support a mechanism claim.

## Inputs Used

- `prn_call_initial`
- `prn_mechanism_call`
- `prn_query_cov_pct`
- `prn_best_single_cov_pct`
- `prn_hsp_n`
- `bp_category`
- `modules/step4_prn_validation/outputs/bp_prn_is_hits.tsv` best-hit fields:
  `hit_support_tier`, `passes_prn02_support_rule`, `supports_assigned_mechanism`

## Rule Order

Rules are applied in order, one row at a time.

| `prn_call_confidence` | `prn04_rule` | Applies when |
| --- | --- | --- |
| `insufficient_evidence` | `no_current_step3_prn_input` | `prn_mechanism_call=insufficient_data` and the row has no current `step3` `prn` input |
| `insufficient_evidence` | `assembly_sequence_unavailable` | `prn_mechanism_call=insufficient_data` because `prn_call_initial=missing_fasta` |
| `insufficient_evidence` | `partial_prn_alignment_no_structural_upgrade` | `prn_mechanism_call=insufficient_data` because `prn_call_initial=partial` |
| `insufficient_evidence` | `other_insufficient_signal` | fallback for any other `insufficient_data` row |
| `assembly_high` | `intact_single_hsp_ge95cov` | `prn_mechanism_call=intact`, `prn_call_initial=intact`, `prn_hsp_n=1`, and both coverage fields are `>=95` |
| `assembly_moderate` | `intact_but_suboptimal_alignment` | fallback for any other `intact` row |
| `assembly_high` | `is481_supported_strong` | `coding_disrupted_is481` with `bp_category=insertion_like`, `prn_hsp_n>=2`, `prn_query_cov_pct>=95`, and the `PRN-03` best hit both passes the `PRN-02` support rule and is a `strong` match for the assigned mechanism |
| `assembly_moderate` | `is481_supported_moderate` | same as above, but the best hit is only `moderate` |
| `assembly_low` | `is481_label_without_supported_best_hit` | fallback for any other `coding_disrupted_is481` row |
| `assembly_moderate` | `other_is_supported_strong` | `coding_disrupted_other_is` with insertion-like structure and a `strong` mechanism-consistent non-`IS481` hit |
| `assembly_low` | `other_is_supported_moderate` | `coding_disrupted_other_is` with insertion-like structure and a `moderate` mechanism-consistent non-`IS481` hit |
| `assembly_low` | `other_is_incomplete_support_profile` | fallback for any other `coding_disrupted_other_is` row |
| `assembly_moderate` | `within_contig_split_alignment` | `coding_disrupted_inversion_or_rearrangement` with `bp_category=within_contig`, `prn_hsp_n>=2`, and `prn_query_cov_pct>=95` |
| `assembly_low` | `rearrangement_but_incomplete_structural_support` | fallback for any other `coding_disrupted_inversion_or_rearrangement` row |
| `assembly_low` | `insertion_like_without_supported_is` | `coding_disrupted_other` with insertion-like structure, `prn_hsp_n>=2`, and `prn_query_cov_pct>=95` |
| `assembly_low` | `other_disruption_limited_support` | fallback for any other `coding_disrupted_other` row |
| `assembly_low` | `fragmented_or_near_contig_end` | `uncertain_fragmented_assembly` rows |
| `insufficient_evidence` | `fallback_unclassified` | final safeguard if a future row falls outside the expected state space |

## Interpretation Notes

- `assembly_high` does not mean read-backed confirmation. It means the current assembly-side evidence is internally consistent and strong enough for descriptive summaries before `VAL-02`.
- `assembly_moderate` is the deliberate ceiling for within-contig rearrangement calls in this phase. Those rows have structural split evidence, but they do not yet have orthogonal read support.
- `assembly_low` keeps ambiguous assembly patterns visible instead of collapsing them into stronger classes.
- `mapped_from_legacy_mirror_or_alt_accession` is not itself a confidence penalty when the canonical sample and accession root still match.
- `read_validation_status` and `read_validation_support` remain orthogonal fields. They should be updated by validation tasks, not overwritten by `PRN-04`.

## Current Example Rows

- `SAMN02436203` (`GCA_000479575.2`) is `intact` with `99.60%` union coverage and `1` HSP, so it scores `assembly_high` via `intact_single_hsp_ge95cov`.
- `SAMN02436205` (`GCA_000479535.2`) is `coding_disrupted_is481` with an insertion-like breakpoint and a `strong` `IS481` best hit, so it scores `assembly_high` via `is481_supported_strong`.
- `SAMN11539676` (`GCA_027271095.1`) is `coding_disrupted_is481`, but its best `IS481` hit is only `moderate` (`86.28%` identity, `100.00%` query coverage), so it scores `assembly_moderate`.
- `SAMN07137591` (`GCA_002240475.1`) is `coding_disrupted_inversion_or_rearrangement` with `bp_category=within_contig` and `95.02%` union coverage, so it scores `assembly_moderate`.
- `SAMN02436242` (`GCA_000479775.2`) is insertion-like but has no supported IS hit, so it remains `assembly_low`.
- `SAMN02436206` (`GCA_000479615.2`) carries a fragmentation signal, so it remains `assembly_low`.
- `SAMN12523755` (`GCA_027271015.1`) has `prn_call_initial=missing_fasta`, so it is `insufficient_evidence`.
- `SAMN02470861` (`GCA_000193515.2`) has no current `step3` `prn` input and therefore remains `insufficient_evidence`.
