# Output Contract Checklist

## Purpose

This checklist defines the minimum repository-wide completion standard for any implementation task that creates or changes scripts, tables, manifests, summaries, or model outputs.

A task should not be marked complete until every applicable item below is satisfied.

## When To Use This Checklist

Use this checklist for work in:

- `modules/public_health/`
- `modules/step1_ingest/`
- `modules/step4_prn_validation/`
- `modules/step5_phylogeny_asr/`
- `modules/step6_epi_transmission/`
- `manuscript/`
- `visualization/` when it depends on frozen tables

If an item is not applicable, the task report should say so explicitly instead of leaving the status ambiguous.

## Completion Gate

| Gate ID | Requirement | Pass condition |
| --- | --- | --- |
| `C01` | Write-scope compliance | All edits stay inside the task brief write scope. |
| `C02` | Deliverables present | Every declared script, table, schema, note, or README deliverable exists at the documented path. |
| `C03` | Machine-readable outputs | Any required data artifact is emitted as TSV, JSON, CSV, or another declared machine-readable format. |
| `C04` | Schema conformance | Output headers match the declared schema template exactly, or the task explicitly creates the new schema contract. |
| `C05` | Stable keys and missingness | Required identifiers, join keys, and missing values are explicit rather than silently dropped or imputed. |
| `C06` | Provenance columns | Every derived table includes source columns, traceable keys, or both. |
| `C07` | Freeze and release metadata | Outputs record the source release date, access date, export date, or data freeze date expected by the module. |
| `C08` | Audit trail for manual work | Any manual curation, override, or review-only decision is logged in a sidecar file or audit table. |
| `C09` | Script usability | New or changed scripts expose concise usage notes or CLI help and avoid hard-coded local-only paths. |
| `C10` | README or module notes | A module README is updated or created when behavior, inputs, outputs, or assumptions change materially. |
| `C11` | Basic validation checks | The task runs at least lightweight checks on headers, row counts, required fields, or expected statuses. |
| `C12` | Caveats recorded | Known limitations, unresolved cases, and uncertainty are written down rather than hidden. |
| `C13` | Reproducibility metadata | Major output directories or task notes record software version, environment note, or reference version when relevant. |
| `C14` | Task report completeness | The task report includes files changed, outputs generated, checks run, and remaining risks or follow-ups. |

## Practical Interpretation

### 1. Schema Conformance

- If a schema file already exists, the output header must match it exactly.
- If a task introduces a new output contract, the task should also create or document the schema.
- Schema changes must be deliberate and documented, not accidental drift.

Examples of current schema-driven outputs include:

- [ph_country_year_master.schema.tsv](modules/public_health/outputs/ph_country_year_master.schema.tsv)
- [ph_country_program_metadata.schema.tsv](modules/public_health/outputs/ph_country_program_metadata.schema.tsv)
- [bp_prn_mechanism_calls.schema.tsv](modules/step4_prn_validation/outputs/bp_prn_mechanism_calls.schema.tsv)
- [bp_prn_read_validation.schema.tsv](modules/step4_prn_validation/outputs/bp_prn_read_validation.schema.tsv)
- [bp_prn_independent_origins.schema.tsv](modules/step5_phylogeny_asr/outputs/bp_prn_independent_origins.schema.tsv)
- [bp_prn_clade_summary.schema.tsv](modules/step5_phylogeny_asr/outputs/bp_prn_clade_summary.schema.tsv)
- [bp_country_year_analysis_input.schema.tsv](modules/step6_epi_transmission/outputs/bp_country_year_analysis_input.schema.tsv)
- [bp_country_year_association_models.schema.tsv](modules/step6_epi_transmission/outputs/bp_country_year_association_models.schema.tsv)

### 2. Provenance And Freeze Metadata

At minimum, outputs should preserve the metadata needed to answer:

- what inputs were used,
- what source release or access window was assumed,
- what freeze date the task depends on,
- how a downstream module can trace a row back to upstream records.

Common fields include:

- `source_name`
- `source_url`
- `source_release_date`
- `access_date`
- `export_date`
- `data_freeze_date`
- accession fields or canonical sample identifiers

### 3. README And Usage Notes

Update the relevant module README when a task materially changes:

- expected inputs,
- output files,
- command-line usage,
- key assumptions,
- caveats or interpretation notes.

If the module already has a README, extend it. If not, add a short usage note in the most appropriate task-facing document.

### 4. Basic Validation Checks

Minimum acceptable checks usually include one or more of:

- header or schema check,
- required-column presence check,
- no empty values in mandatory fields,
- row-count sanity check,
- category/status coverage check,
- successful script help output,
- joinability or key-uniqueness check where relevant.

The goal is not exhaustive testing for every task. The goal is to avoid shipping malformed or non-auditable outputs.

### 5. Manual Review And Uncertainty

If a task cannot resolve all rows deterministically:

- keep unresolved cases explicit,
- label them with a status such as `manual_review_required` or `unknown`,
- document the rule boundary that prevented automatic resolution.

Never hide uncertainty by silently dropping rows or forcing unsupported classifications.

## Suggested Report Template

Every task close-out should report:

1. Files changed
2. Outputs generated
3. Checks run
4. Remaining risks or follow-ups

## Definition Of Done

A task is ready to mark `done` only when:

- all declared deliverables exist,
- applicable checklist gates pass,
- the task report records what was checked,
- any residual risk is visible to the next module rather than buried in code.
