# Duplicate Resolution Policy

## Purpose

This document converts the narrative duplicate hierarchy in the technical specification into a deterministic policy for future `bp_step1` manifest and duplicate-resolution scripts.

The policy is intended to drive:

- duplicate-group construction,
- representative-record selection within each group,
- audit output fields for `GC-03`,
- conservative handling of likely replicates that should not be auto-collapsed.

## Scope

This policy applies to public *Bordetella pertussis* genome records used to build the retained-sample manifest.

It does not do the following:

- replace genome QC thresholds,
- decide analytic cohort inclusion,
- infer biological identity from weak metadata alone,
- collapse records that only look similar by country, year, or generic isolate text.

`GC-03` should implement this policy after `GC-02` has made `raw_reads_available` explicit for every candidate record.

## Source Hierarchy From The Technical Spec

The governing order from the technical specification is:

1. Prefer a sample with raw reads over an assembly-only sample.
2. Prefer a more complete assembly over a lower-quality assembly.
3. Prefer richer metadata over sparse metadata.
4. Prefer primary records over repackaged records.
5. Retain replicate genomes only for dedicated sensitivity analyses.

The remainder of this document makes those rules executable.

## Normalized Fields

Future scripts should normalize the following policy fields before duplicate resolution.

| Policy field | Current `bp_step1` source | Notes |
| --- | --- | --- |
| `assembly_accession` | `Assembly Accession` | Keep full accession with prefix and version. |
| `assembly_accession_root` | derived from `Assembly Accession` | Strip `GCA_` or `GCF_` prefix and version suffix so mirror pairs can join. |
| `biosample_accession` | `Assembly BioSample Accession` | Strong duplicate identifier when non-empty. |
| `bioproject_accession` | `Assembly BioProject Accession` | Used for conservative review grouping only. |
| `source_database` | `Source Database` | Current values are mainly `SOURCE_DATABASE_GENBANK` and `SOURCE_DATABASE_REFSEQ`. |
| `raw_reads_available` | future `GC-02` field | Must be explicit `true` or `false` before automatic resolution. |
| `sra_run_accession` | future `GC-02` field | Strong duplicate identifier when non-empty. |
| `ena_run_accession` | future `GC-02` field | Strong duplicate identifier when non-empty. |
| `assembly_level` | `Assembly Level` | Used in QC ranking. |
| `n_contigs` | `Assembly Stats Number of Contigs` | Lower is better. |
| `contig_n50` | `Assembly Stats Contig N50` | Higher is better. |
| `assembly_release_date` | `Assembly Release Date` | Late tie-breaker only. |
| `country` | `country` | Must be normalized before policy application. |
| `year` | `year` | Required for review-only grouping. |
| `collection_date_raw` | `Assembly BioSample Collection date` | Counts toward metadata richness. |
| `host` | `Assembly BioSample Host` | Counts toward metadata richness. |
| `isolation_source` | `Assembly BioSample Isolation source` | Counts toward metadata richness. |
| `strain` | `Assembly BioSample Strain` | Used to build a normalized sample token. |
| `isolate` | `Assembly BioSample Isolate` | Preferred over `strain` when specific. |
| `sequencing_tech` | `Assembly Sequencing Tech` | Counts toward metadata richness. |

## Group Construction

Duplicate resolution is a two-stage process:

1. Build candidate duplicate groups using strong identifiers.
2. Rank records within each auto-resolvable group to choose one representative.

### Strong Auto-Grouping Keys

Records belong to the same auto-resolvable duplicate group if they share any non-empty strong identifier below. Grouping should use connected components, so a record linked by any strong key joins the same final group.

| Priority | Grouping key | Auto-collapse allowed | Reason |
| --- | --- | --- | --- |
| 1 | `assembly_accession_root` | Yes | Joins `GCA_` and `GCF_` mirror pairs or the same assembly root across versions. |
| 2 | `biosample_accession` | Yes | Strong sample-level accession. |
| 3 | `sra_run_accession` | Yes | Strong read-level accession after `GC-02`. |
| 4 | `ena_run_accession` | Yes | Strong read-level accession after `GC-02`. |

### Review-Only Candidate Groups

If no strong identifier matches, a record may be placed into a review-only candidate group when all of the following are true:

- same normalized `strain_or_isolate_token`,
- same normalized `country`,
- same `year`,
- same non-empty `bioproject_accession`.

These groups are not auto-collapsed. They should be written to a manual-review queue and left unresolved unless external evidence confirms they are true replicates.

### Guardrails Against False Collapse

The following are never sufficient on their own to auto-group records:

- same country and year only,
- same strain text only,
- same isolate text only,
- same genome size or contig count only,
- generic sample tokens such as `Bordetella pertussis`, `missing`, `unknown`, `not applicable`, or blank values.

When `isolate` is generic or missing, scripts may fall back to `strain` for review grouping, but not for automatic collapse without a strong accession match.

## Representative Selection Within Auto-Resolvable Groups

Within a duplicate group built from strong identifiers, compare records in the exact order below. The first rule that separates two records decides the winner. If a rule is tied, continue to the next rule.

### Rule Order

1. Prefer `raw_reads_available = true` over `false`.
2. Prefer higher `assembly_level_rank`.
3. Prefer lower `n_contigs`.
4. Prefer higher `contig_n50`.
5. Prefer higher `metadata_completeness_score`.
6. Prefer higher `source_database_rank`.
7. Prefer earlier `assembly_release_date`.
8. Prefer lexicographically smaller `assembly_accession` as the final deterministic tie-break.

### Rank Definitions

#### `assembly_level_rank`

| Assembly level | Rank |
| --- | --- |
| `Complete Genome` | 4 |
| `Chromosome` | 3 |
| `Scaffold` | 2 |
| `Contig` | 1 |
| anything else or missing | 0 |

#### `source_database_rank`

| Source database | Rank | Interpretation |
| --- | --- | --- |
| `SOURCE_DATABASE_GENBANK` | 2 | Primary submitter record |
| `SOURCE_DATABASE_REFSEQ` | 1 | Repackaged or mirrored record |
| anything else or missing | 0 | Unknown or unmapped source |

This means a `GCA_` record beats the corresponding `GCF_` mirror if earlier rules remain tied.

#### `metadata_completeness_score`

The score is the count of non-missing values across:

- `biosample_accession`
- `bioproject_accession`
- `country`
- `year`
- `collection_date_raw`
- `host`
- `isolation_source`
- `strain_or_isolate_token`
- `sequencing_tech`

Treat the following as missing after trimming and lowercasing:

- empty string
- `missing`
- `unknown`
- `not applicable`
- `n/a`

If only one of two compared records has a parsable value for a ranking field, the non-missing value wins that comparator.

## Replicate Handling

The main manifest should retain one representative per auto-resolved duplicate group.

Additional records should be handled as follows:

- mirror or duplicate records dropped by strong-ID grouping receive status `drop_duplicate`.
- confirmed technical or biological replicates receive status `keep_replicate_sensitivity_only`.
- unresolved review groups receive status `manual_review_required` and must not be collapsed automatically.

## Recommended Audit Columns For `GC-03`

`bp_duplicate_resolution.tsv` should include at least these fields:

- `duplicate_group_id`
- `duplicate_group_type`
- `duplicate_evidence_basis`
- `record_decision`
- `kept_assembly_accession`
- `decisive_rule_id`
- `decisive_rule_value`
- `review_required`
- `decision_note`

These columns make it possible to explain every retained or dropped row without re-reading source metadata.

## Worked Examples

### Example 1. Exact GenBank and RefSeq mirror pair

Current `bp_step1` contains the pair below for the same BioSample and assembly root:

| Assembly accession | BioSample | Source database | Assembly level | Contigs | Contig N50 |
| --- | --- | --- | --- | --- | --- |
| `GCA_019974435.1` | `SAMD00328237` | `SOURCE_DATABASE_GENBANK` | `Complete Genome` | `1` | `4130169` |
| `GCF_019974435.1` | `SAMD00328237` | `SOURCE_DATABASE_REFSEQ` | `Complete Genome` | `1` | `4130169` |

Decision:

- auto-group by `assembly_accession_root` and `biosample_accession`
- retain `GCA_019974435.1`
- drop `GCF_019974435.1`

Reason:

- QC and metadata are tied, so rule 6 applies and the primary `GenBank` record beats the repackaged `RefSeq` mirror.

### Example 2. Same token does not justify automatic collapse

Current `bp_step1` also contains the following USA 1939 records for token `B203`:

| Assembly accession | BioSample | BioProject | Source database | Assembly level | Contigs |
| --- | --- | --- | --- | --- | --- |
| `GCA_001199415.1` / `GCF_001199415.1` | `SAMEA751178` | `PRJEB2274` | GenBank / RefSeq | `Scaffold` | `574` |
| `GCA_001687385.1` / `GCF_001687385.1` | `SAMN03877216` | `PRJNA279196` | GenBank / RefSeq | `Complete Genome` | `1` |

Decision:

- collapse each mirror pair separately by strong identifiers
- do not auto-collapse the two BioSamples into one sample

Reason:

- the isolate token is similar, but the stable accessions and BioProjects disagree
- choosing the complete genome automatically would risk deleting a distinct historical record

### Example 3. Generic isolate text must not merge distinct BioSamples

For Czech Republic 2012, multiple records share isolate text `Bordetella pertussis` within project `PRJEB26966`, but they belong to distinct BioSamples such as `SAMEA4693572`, `SAMEA4693573`, `SAMEA4693574`, and `SAMEA4693575`.

Decision:

- collapse each `GCA_` and `GCF_` mirror pair separately
- keep the four BioSamples as separate retained candidates unless later external evidence proves replication

Reason:

- `Bordetella pertussis` is a generic token and is explicitly blocked from driving duplicate collapse

### Example 4. Raw-read-linked record wins once `GC-02` completes

If a future duplicate group contains two otherwise tied records and only one has linked SRA or ENA runs:

| Assembly accession | BioSample | Raw reads available | Assembly level | Contigs |
| --- | --- | --- | --- | --- |
| `record_A` | same | `true` | `Scaffold` | `54` |
| `record_B` | same | `false` | `Scaffold` | `54` |

Decision:

- retain `record_A`

Reason:

- rule 1 is decisive and raw-read support outranks later QC and source-database rules

## Implementation Notes For Future Scripts

- Build duplicate groups only after country normalization and raw-read-link recovery are available.
- Use connected components for strong-ID grouping so a record linked by BioSample and another linked by SRA are not split into separate groups.
- Store the first decisive rule, not just the final retained accession.
- Leave unresolved review groups visible rather than silently collapsing them.
