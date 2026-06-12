# bp_step4

Step 4 upgrades the current `prn` disruption workflow from coarse assembly-based calls to mechanism-aware structural classification with optional raw-read validation.

## Purpose

This module is the main technical differentiator of the redesigned study. Its output should answer not only whether `prn` appears disrupted, but also:

- what kind of disruption is present,
- how strong the evidence is,
- whether the event is supported by raw reads,
- whether the event is likely biological or artifactual.

## Inputs

- `../modules/step2_typing/outputs/bp_qc_merged_mlst_markers.tsv`
- `../modules/step3_prn_scan/outputs/bp_prn_disruption_calls.tsv`
- `../modules/step3_prn_scan/outputs/bp_prn_breakpoint_evidence.tsv`
- Assembly FASTA files indexed in previous steps
- Raw reads for the validation subset when available

## Target Output Classes

- `intact`
- `promoter_disrupted`
- `coding_disrupted_is481`
- `coding_disrupted_other_is`
- `coding_disrupted_deletion`
- `coding_disrupted_inversion_or_rearrangement`
- `coding_disrupted_other`
- `uncertain_fragmented_assembly`
- `insufficient_data`

## Evidence Channels

- BLAST alignment structure
- Query coverage
- Contig context
- Subject-gap characteristics
- Insertion-sequence matches
- Flanking-sequence coherence
- Raw-read support

## Planned Outputs

- `outputs/bp_prn_mechanism_calls.tsv`
- `outputs/bp_prn_is_hits.tsv`
- `outputs/bp_prn_event_catalog.tsv`
- `outputs/bp_prn_mechanism_summary.tsv`
- `outputs/bp_prn_confidence_summary.tsv`
- `outputs/bp_prn_breakpoint_summary.tsv`
- `outputs/bp_prn_unresolved_summary.tsv`
- `outputs/bp_prn_country_year_summary.tsv`
- `outputs/bp_prn_validation_subset.tsv`
- `outputs/bp_prn_read_validation.tsv`
- `outputs/bp_prn_read_validation_is_calls.tsv`
- `outputs/bp_prn_read_validation_tsd.tsv`
- `outputs/bp_prn_validation_summary.tsv`
- `outputs/bp_prn_is_hotspot_results.tsv`
- `outputs/bp_prn_is_hotspot_density.pdf`

Schema templates are already provided in `outputs/`.

## Recommended Script Plan

- `bin/step4_01_build_is_reference.py`
- `bin/step4_02_scan_prn_mechanisms.py`
- `bin/step4_02b_summarize_is_hits.py`
- `bin/step4_02c_score_prn_calls.py`
- `bin/step4_02d_build_prn_summary_tables.py`
- `bin/step4_03a_build_validation_subset.py`
- `bin/step4_03b_assess_validation_feasibility.py`
- `bin/step4_03c_prepare_is_reference.py`
- `bin/step4_03d_build_read_validation_batch.py`
- `bin/step4_03e_run_is_read_validation.sh`
- `bin/step4_03f_hotspot_test.py`
- `bin/step4_03_validate_prn_with_reads.py`
- `bin/step4_04_summarize_prn_validation.py`

## IS Reference Layer

`PRN-01` materializes a small curated insertion-sequence seed set for downstream `prn` screening.

Current reference artifacts:

- `references/is_elements/bp_is_reference.fasta`
- `references/is_elements/bp_is_reference_metadata.tsv`
- `references/is_elements/bp_is_reference.shell_safe.fasta`
- `references/is_elements/bp_is_reference_shell_safe_map.tsv`

Build them with:

```bash
python bin/step4_01_build_is_reference.py
```

Design notes:

- The reference layer is versioned and provenance-backed.
- Each row records accession-level source metadata and sequence hashes.
- The output is intentionally a scan-ready seed library, not the full mechanism engine.

## Deliverable Standard

Every final `prn` mechanism call should carry:

- a final mechanism label,
- a confidence score,
- a machine-readable evidence flag set,
- a raw-read validation status if reads are available.

`PRN-05` summary tables now provide:

- mechanism-level totals for figure-ready composition summaries,
- confidence-aware breakdowns,
- breakpoint and insertion-size summaries,
- unresolved-call rollups,
- country-year genomic aggregates for later integration steps.

`VAL-01` now provides:

- an explicit raw-read validation subset manifest with selection strata,
- per-sample inclusion reasons,
- accession provenance that distinguishes SRA-backed and ENA-backed read links.

## Read Validation Batch Artifacts

step4_03d_build_read_validation_batch.py now writes two manifests:

- work/read_validation/<label>/bp_prn_read_validation_batch.tsv
  - selected rows that can run read validation immediately
- work/read_validation/<label>/bp_prn_read_validation_missing_inputs.tsv
  - blocked or deferred rows with explicit missing input reasons
