# Documentation Map

This file is the main `docs/` guide. Keep it focused on active repository documentation rather than retired submission targets.

## Read First

- `architecture/project_status.md`: current repository state and active sources of truth.
- `../README.md`: repository-level navigation and execution boundary.
- `architecture/technical_spec_zh.md`: consolidated Chinese technical specification for the active code and script boundary.
- `architecture/technical_spec_en.md`: consolidated English technical specification for the active code and script boundary.
- `release_package.md`: code/data release boundary for the paper-sharing repository snapshot.
- `manuscript_rebuild/selected_country/README.md`: selected-country rebuild guide retained as analysis documentation.
- `../manuscript/README.md`: Communications Biology manuscript-facing deliverables.

## Active Subdirectories

- `architecture/`: project status and technical specifications.
- `manuscript_rebuild/selected_country/`: selected-country analysis notes and curation context.
- `policies/`: reusable project policies and source inventories.

## Submission Boundary

The active submission target is Communications Biology. Submission-facing text files live under `manuscript/text/commsbio_*.md`; source data and supplementary tables remain under `manuscript/submission_data/` and `manuscript/supplementary/`.

Retired prior-journal submission notes and manuscript drafts have been removed from the active tree.

## Do Not Move Without Updating Consumers

These paths are referenced by scripts or build steps and should stay stable unless the consuming code is updated too:

- `policies/public_health_source_inventory.md`
- `policies/public_health_variable_dictionary.tsv`
- `../state/ledgers/source_freeze_ledger.tsv`
- `../manuscript/submission_data/cohort/master_cohort_decision_log.tsv`
- `../manuscript/submission_data/cohort/master_cohort_flow_summary.tsv`
- `../manuscript/submission_data/validation/validation_priority_ledger.tsv`
- `../manuscript/submission_data/source_data/`

## Maintenance Rules

- Keep only current status, design, policy, contract and reference docs.
- Remove retired venue-switch notes after their content has been reflected in the active manuscript package.
- Do not duplicate module usage guides here if a stage README is already the active source of truth.
