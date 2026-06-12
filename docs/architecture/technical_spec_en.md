# Pertussis Gene Technical Specification

## 1. Scope

This document is a code-centric consolidated technical specification. It is intended to replace the need to jump across many small Markdown notes by describing the active code boundaries, the design of each module, why the design was chosen, where each stage writes its outputs, and how those outputs flow into the manuscript and figure pipeline.

## 2. Design Principles

- The repository is organized around a single source of truth and stage-owned responsibilities. `outputs/workflow/manifest/manifest.tsv` is the cross-stage sample manifest, while `modules/step2_typing/outputs/bp_genotype_manifest.tsv` is the standardized typing-layer manifest that feeds the PRN structural layer; step-specific logic stays in `modules/*/bin/`, and `workflow/bin/` plus `workflow/lib/` are reserved for orchestration and shared helpers.

- Active execution paths are explicitly separated from audit-only material. Audit-only scripts, legacy manuscript tables, and run snapshots are removed from the release branch; recover them from git history if a historical comparison is needed.

- Phylogeny and ASR use dual validation. M4 centers on Gubbins + IQ-TREE, with optional ClonalFrameML and RAxML cross-checks; M5 emits both Fitch and PastML outputs so that evolutionary claims do not depend on a single method.

- Assembly-level calling and read-backed validation are intentionally separated. Step4 first classifies `prn` mechanisms from assemblies, then validates selected events with `ismapper + panISa`, because IS insertions and rearrangements are often fragmented or misresolved in assembly-only data.

- The manuscript and figure layer does not re-implement the core analyses. `manuscript/scripts/` stages frozen TSVs for submission, and `manuscript/figures/scripts/` focuses on rendering.

## 3. End-to-End Execution

| Stage | Summary | Main outputs |
| --- | --- | --- |
| M0 | Build the unified manifest, run the readiness checks, and snapshot versions. | `outputs/workflow/manifest/`, `outputs/workflow/checkpoints/`, `outputs/workflow/versions.txt` |
| M1/M2 | Completeness checks, read-tracing and download planning, assembly QC, and missingness modeling. | `outputs/workflow/reads_plan/`, `outputs/workflow/assembly_qc/`, `outputs/workflow/missingness_model/` |
| M3 | Build the bootstrap SNP alignment with contig-mode Snippy. | `outputs/workflow/snippy_ctg/`, `outputs/workflow/phylo/core.full.aln` |
| M4 | Pre-filter by missingness, run Gubbins, mask recombination, build the ML tree, optionally cross-check with CFML/RAxML. | `outputs/workflow/phylo/` |
| M5 | Re-root on the reference, run dual-track ASR, and package independent origin events. | `outputs/workflow/asr/`, `outputs/workflow/asr_sensitivity/` |
| M6 / Step4 | Validate `prn` mechanism calls with reads and derive hotspot summaries. | `modules/step4_prn_validation/outputs/` |
| M7 / PH + Step6 | Build the public-health master table, join it with genomic summaries, then fit ecological and transmission models. | `modules/public_health/outputs/`, `modules/step6_epi_transmission/outputs/`, `outputs/bp_step6_v2_fixed/` |
| Manuscript / Figures | Stage frozen supplementary tables, figure-data extracts, and final figures. | `manuscript/figure_data/`, `manuscript/supplementary/`, `manuscript/figures/outputs/main/`, `manuscript/figures/outputs/extended_data/` |

## 4. Output Map

- `outputs/workflow/manifest/`: cohort single source of truth, reads linkage, build report.
- `outputs/workflow/checkpoints/`: readiness JSON reports and consolidated readiness report.
- `outputs/workflow/reads_plan/`: read availability and download planning.
- `outputs/workflow/assembly_qc/`: assembly QC statistics and summaries.
- `outputs/workflow/missingness_model/`: interpretablity and selection-bias diagnostics.
- `outputs/workflow/snippy_ctg/`: Snippy plan, batch ledgers, completion summaries.
- `outputs/workflow/phylo/`: core alignment, pre-Gubbins filtering, recombination masks, ML trees, tree comparison.
- `outputs/workflow/asr/`: rooted tree, tip states, Fitch states, PastML states, origin events, event subtrees.
- `outputs/workflow/asr_sensitivity/`: support-threshold and composition-filtered robustness runs.
- `modules/step4_prn_validation/outputs/`: PRN mechanism calls, read validation, IS evidence, TSD, hotspot outputs.
- `modules/public_health/outputs/`: cleaned WHO/WUENIC/AMU tables and country-year master table.
- `modules/step6_epi_transmission/outputs/`: joined analytical table, model tables, cross-validation, transmission outputs.
- `outputs/bp_step4_with_ci/` and `outputs/bp_step6_v2_fixed/`: compatibility output locations still consumed by the current figure layer.
- `manuscript/figure_data/`: figure-ready frozen TSV extracts.
- `manuscript/supplementary/`: canonical supplementary tables and data files.
- `manuscript/figures/outputs/main/`, `manuscript/figures/outputs/extended_data/`: rendered main figures and Supplementary Figure source images.

## 5. File Inventory

### 5.1 Root Orchestration And Workflow Core

- [Snakefile](../Snakefile): Top-level dependency graph that declares what the workflow must materialize while delegating implementation details to `rules/*.smk`. Outputs: default `rule all` targets under `outputs/workflow/`, `modules/step4_prn_validation/outputs/`, and `outputs/workflow/epi/`. Status: active scaffold.
- [config/workflow.yaml](../config/workflow.yaml): Central parameter contract for thresholds, reference paths, and stage inputs so configuration does not drift across shell, Python, and Snakemake. Outputs: consumed by all rules; no direct output file. Status: active.
- [rules/manifest.smk](../rules/manifest.smk): Snakemake rules for manifest building and reads availability so downstream steps share the same cohort boundary. Outputs: `outputs/workflow/manifest/`, `outputs/workflow/checkpoints/reads_availability_report.json`. Status: active.
- [rules/reads_qc.smk](../rules/reads_qc.smk): Encapsulates read-level QC and MultiQC aggregation as declarative workflow targets. Outputs: `outputs/workflow/qc/` and cleaned read staging. Status: active scaffold.
- [rules/assembly_qc.smk](../rules/assembly_qc.smk): Makes assembly QC and missingness diagnostics first-class workflow products. Outputs: `outputs/workflow/qc/assembly_qc_report.tsv`, `outputs/workflow/qc/missingness_model.json`. Status: active.
- [rules/snippy.smk](../rules/snippy.smk): Declares Snippy and core-alignment dependencies so M3 can be scheduled reproducibly. Outputs: `outputs/workflow/phylo/core.full.aln` and related Snippy outputs. Status: active scaffold.
- [rules/recomb_filter.smk](../rules/recomb_filter.smk): Declares recombination filtering as a repeatable stage so the tree chain explicitly depends on non-recombinant sites. Outputs: `outputs/workflow/phylo/recomb_filtered.aln` and Gubbins artifacts. Status: active scaffold.
- [rules/ml_tree.smk](../rules/ml_tree.smk): Produces ML-tree artifacts and comparison outputs so ASR depends on a stable tree interface instead of ad hoc commands. Outputs: `outputs/workflow/phylo/iqtree2/`, `outputs/workflow/phylo/raxmlng/`. Status: active scaffold.
- [rules/asr.smk](../rules/asr.smk): Packages rooting, Fitch, PastML, and origin-event generation in a single dependency chain so event-level conclusions remain auditable. Outputs: `outputs/workflow/asr/`. Status: active scaffold.
- [rules/is_detection.smk](../rules/is_detection.smk): Connects Step4 read-validation products to the root workflow without duplicating Step4 logic, preserving step ownership. Outputs: `modules/step4_prn_validation/outputs/`. Status: active thin wrapper.
- [rules/epi_models.smk](../rules/epi_models.smk): Brings exposure-index building, IPW, and panel models into the workflow graph so the ecological layer is reproducible. Outputs: `outputs/workflow/epi/`. Status: active scaffold.
- [workflow/bin/run_full_workflow.sh](../workflow/bin/run_full_workflow.sh): Canonical shell orchestrator for M0-M5, with resume points and artifact verification after each stage. Outputs: logs in `logs/pipeline/` and all verified stage outputs under `outputs/workflow/`. Status: active primary runner.
- [workflow/bin/m0_foundation.sh](../workflow/bin/m0_foundation.sh): Builds the manifest, runs the readiness checks, and snapshots versions so the rest of the pipeline starts from a frozen readiness checkpoint. Outputs: `outputs/workflow/manifest/`, `outputs/workflow/checkpoints/`, `outputs/workflow/versions.txt`. Status: active.
- [workflow/bin/m1_m2_qc.sh](../workflow/bin/m1_m2_qc.sh): Wraps completeness checks, read planning, assembly QC, and the missingness model into one support stage. Outputs: `outputs/workflow/reads_plan/`, `outputs/workflow/assembly_qc/`, `outputs/workflow/missingness_model/`. Status: active.
- [workflow/bin/m3_snippy.sh](../workflow/bin/m3_snippy.sh): Builds the bootstrap alignment through plan generation, Snippy batch execution, and snippy-core aggregation, supporting staged batching on limited hardware. Outputs: `outputs/workflow/snippy_ctg/`, `outputs/workflow/phylo/core.full.aln`. Status: active.
- [workflow/bin/m4_phylogeny.sh](../workflow/bin/m4_phylogeny.sh): Wraps missingness filtering, Gubbins, recombination masking, IQ-TREE2, and optional CFML/RAxML into a host-friendly phylogeny pipeline that can survive environment variability. Outputs: `outputs/workflow/phylo/`. Status: active.
- [workflow/bin/m5_asr.sh](../workflow/bin/m5_asr.sh): Main M5 wrapper for rooting, Fitch, PastML, origin packaging, and summary generation from a single rooted tree and manifest contract. Outputs: `outputs/workflow/asr/`. Status: active.
- [workflow/bin/m5_asr_sensitivity.sh](../workflow/bin/m5_asr_sensitivity.sh): Runs branch-support and composition-filtered sensitivity scenarios so robustness becomes an auditable artifact rather than a narrative claim. Outputs: `outputs/workflow/asr_sensitivity/`. Status: active.
- [workflow/lib/run_foundation_checks.py](../workflow/lib/run_foundation_checks.py): Consolidates reads availability, vaccine-variable coverage, and validation feasibility into JSON reports so readiness decisions are repeatable and machine-readable. Outputs: `outputs/workflow/checkpoints/*.json`. Status: active.
- [workflow/lib/build_analysis_manifest.py](../workflow/lib/build_analysis_manifest.py): Starts from the Step4 mechanism table as the manuscript cohort and enriches it with Step5/Step1 metadata plus the Step2 standardized typing manifest so the project has a single cohort SSOT spanning the typing layer and the PRN structural layer. Outputs: `outputs/workflow/manifest/manifest.tsv`, `manifest_build_report.json`. Status: active.
- [workflow/lib/aggregate_assembly_qc.py](../workflow/lib/aggregate_assembly_qc.py): Aggregates per-sample assembly QC into a unified table for the root workflow contract. Outputs: `outputs/workflow/qc/` or `outputs/workflow/assembly_qc/`. Status: active.
- [workflow/lib/missingness_model.py](../workflow/lib/missingness_model.py): Models `prn` interpretability and missingness mechanisms so selection bias is explicit and measurable. Outputs: `outputs/workflow/missingness_model/` or `outputs/workflow/qc/missingness_model.json`. Status: active.
- [workflow/lib/filter_alignment_by_missingness.py](../workflow/lib/filter_alignment_by_missingness.py): Filters high-missingness tips before Gubbins so extreme samples do not dominate recombination and tree inference. Outputs: `outputs/workflow/phylo/core.filtered.aln`, `pre_gubbins_missingness.tsv`. Status: active.
- [workflow/lib/mask_recombination.py](../workflow/lib/mask_recombination.py): Masks recombinant sites using the Gubbins GFF so downstream trees are based on non-recombinant sequence. Outputs: `outputs/workflow/phylo/recomb_filtered.aln`, mask summary JSON. Status: active.
- [workflow/lib/extract_iqtree_composition_report.py](../workflow/lib/extract_iqtree_composition_report.py): Parses IQ-TREE composition warnings into a TSV so composition-based sensitivity analysis is scriptable. Outputs: `outputs/workflow/phylo/iqtree2/ml_tree.composition.tsv`. Status: active.
- [workflow/lib/compare_trees.py](../workflow/lib/compare_trees.py): Compares IQ-TREE and RAxML tip sets and consistency so tree robustness is machine-readable. Outputs: `outputs/workflow/phylo/tree_comparison_report.json`. Status: active.
- [workflow/lib/root_tree_on_tip.py](../workflow/lib/root_tree_on_tip.py): Re-roots the ML tree on the reference tip and emits node metadata so Fitch and PastML consume the same rooted topology. Outputs: `outputs/workflow/asr/rooted_ml_tree.reference_rooted.nwk`, `rooted_tree_node_metadata.tsv`. Status: active.
- [workflow/lib/asr_parsimony.py](../workflow/lib/asr_parsimony.py): Implements Fitch parsimony and emits tip states, node states, transitions, and the PastML input table in one pass. Outputs: `outputs/workflow/asr/tip_states.tsv`, `parsimony_states.tsv`, `parsimony_transitions.tsv`, `pastml_input.tsv`. Status: active.
- [workflow/lib/asr_pastml.py](../workflow/lib/asr_pastml.py): Normalizes raw PastML output and compares it against Fitch events so the likelihood track is directly comparable to the parsimony track. Outputs: `outputs/workflow/asr/pastml_states.tsv`, `pastml_origin_events.tsv`, `track_comparison.tsv`. Status: active.
- [workflow/lib/origin_events.py](../workflow/lib/origin_events.py): Scans intact-to-disrupted transitions and packages descendant tip subsets so independent origins become event-level artifacts. Outputs: `outputs/workflow/asr/origin_events.tsv`, `outputs/workflow/asr/event_subtrees/`. Status: active.
- [workflow/lib/prune_tree_by_tips.py](../workflow/lib/prune_tree_by_tips.py): Prunes Newick trees by tip list without requiring `ete3`, enabling composition-filtered sensitivity on the current host. Outputs: pruned tree files for sensitivity runs. Status: active.
- [workflow/lib/build_ap_exposure_index.py](../workflow/lib/build_ap_exposure_index.py): Converts vaccine-program variables into a standardized aP exposure index for downstream models. Outputs: `outputs/workflow/epi/ap_exposure_index.tsv`. Status: active.
- [workflow/lib/ipw_prevalence.py](../workflow/lib/ipw_prevalence.py): Applies inverse-probability weighting to disrupted prevalence so interpretability bias is reflected in the ecological layer. Outputs: `outputs/workflow/epi/ipw_prevalence.tsv`. Status: active.
- [workflow/lib/panel_model.py](../workflow/lib/panel_model.py): Fits country-year panel models linking the exposure index to disrupted prevalence. Outputs: `outputs/workflow/epi/panel_model_results.tsv`, diagnostics. Status: active.
- [workflow/lib/its_feasibility.py](../workflow/lib/its_feasibility.py): Evaluates whether interrupted time-series analysis is feasible before fitting unsupported causal narratives. Outputs: `outputs/workflow/epi/its_feasibility_report.tsv`. Status: active.
- [manuscript/scripts/freeze/generate_supplementary_table_1.py](../manuscript/scripts/freeze/generate_supplementary_table_1.py): Exports Supplementary Table 1 from pipeline outputs so the cohort metadata table is reproducible. Outputs: `manuscript/supplementary/Supplementary_Table_1_genome_metadata.tsv`. Status: active.
- [workflow/bin/find-big-git-objects.sh](../workflow/bin/find-big-git-objects.sh): Maintenance helper for locating oversized Git objects and keeping the repository manageable. Outputs: stdout only. Status: maintenance only.
- [validate_figure1.py](../validate_figure1.py): One-off validator comparing old and fixed `R_e` data after the transition to `outputs/bp_step6_v2_fixed/`. Outputs: `outputs/figure1_validation.png` and console summary. Status: audit/helper, not part of the main pipeline.

### 5.2 Step1 Core

- [modules/step1_ingest/bin/core/01_fetch_ncbi_report.py](../modules/step1_ingest/bin/core/01_fetch_ncbi_report.py): Fetches raw NCBI metadata; separated from cleaning so upstream refreshes remain simple. Outputs: Step1 raw metadata snapshots under `modules/step1_ingest/outputs/`.
- [modules/step1_ingest/bin/core/02_export_ncbi_tsv.py](../modules/step1_ingest/bin/core/02_export_ncbi_tsv.py): Flattens NCBI export into TSV so downstream cleaners consume a stable tabular interface. Outputs: flattened Step1 TSVs.
- [modules/step1_ingest/bin/core/03_clean_metadata_aggregate.py](../modules/step1_ingest/bin/core/03_clean_metadata_aggregate.py): Normalizes country/date/source fields before cohort construction, reducing source heterogeneity early. Outputs: cleaned Step1 metadata tables.
- [modules/step1_ingest/bin/core/04_download_ncbi_genomes.py](../modules/step1_ingest/bin/core/04_download_ncbi_genomes.py): Downloads genome assemblies so metadata and physical FASTA assets are staged together. Outputs: genome FASTAs in Step1 staging and later `pertussis_data/bp_genomes_qc/assemblies/`.
- [modules/step1_ingest/bin/core/05_verify_genome_extract.py](../modules/step1_ingest/bin/core/05_verify_genome_extract.py): Verifies extracted genome assets before they enter QC. Outputs: verification summaries in `modules/step1_ingest/outputs/`.
- [modules/step1_ingest/bin/manifest/06_build_public_manifest.py](../modules/step1_ingest/bin/manifest/06_build_public_manifest.py): Builds the pre-deduplicated public manifest so maximum coverage is preserved before explicit duplicate resolution. Outputs: public manifest tables.
- [modules/step1_ingest/bin/manifest/07_recover_raw_reads.py](../modules/step1_ingest/bin/manifest/07_recover_raw_reads.py): Recovers raw-read links from SRA/ENA so assembly-only and read-backed paths meet at the manifest layer. Outputs: read linkage fields in Step1 manifests.
- [modules/step1_ingest/bin/manifest/08_resolve_duplicates.py](../modules/step1_ingest/bin/manifest/08_resolve_duplicates.py): Resolves duplicate records so downstream models and trees use a single canonical sample per entity. Outputs: deduplicated manifest tables and decision ledgers.
- [modules/step1_ingest/bin/manifest/09_build_analysis_cohorts.py](../modules/step1_ingest/bin/manifest/09_build_analysis_cohorts.py): Builds analysis cohorts, especially country-year cohort C, so analytic boundaries are explicit and reproducible. Outputs: cohort-specific tables such as `bp_cohort_C_country_year.tsv`.
- [modules/step1_ingest/bin/raw_reads/10_build_download_plan.py](../modules/step1_ingest/bin/raw_reads/10_build_download_plan.py): Builds the raw-read download plan so heavy I/O can be batched. Outputs: Step1 download plan tables.
- [modules/step1_ingest/bin/raw_reads/11_split_download_plan.py](../modules/step1_ingest/bin/raw_reads/11_split_download_plan.py): Splits large download plans into shards for distributed or resumable execution. Outputs: shard plans.
- [modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh](../modules/step1_ingest/bin/raw_reads/12_run_shard_to_assembly.sh): Executes one shard from reads through assembly, enabling distributed batch assembly. Outputs: per-shard assemblies and logs in Step1 staging.
- [modules/step1_ingest/bin/raw_reads/13_preflight_env.sh](../modules/step1_ingest/bin/raw_reads/13_preflight_env.sh): Preflight environment checks before heavy download and assembly jobs. Outputs: stdout/logs only.
- [modules/step1_ingest/bin/raw_reads/14_setup_env.sh](../modules/step1_ingest/bin/raw_reads/14_setup_env.sh): Sets up the execution environment across hosts or workspaces. Outputs: environment side effects and logs.
- [modules/step1_ingest/bin/raw_reads/15_fetch_taxon_read_run_catalog.py](../modules/step1_ingest/bin/raw_reads/15_fetch_taxon_read_run_catalog.py): Fetches a taxon-level run catalog to broaden upstream read-link discovery. Outputs: Step1 run catalog tables.
- [modules/step1_ingest/bin/raw_reads/16_build_external_gapfill.py](../modules/step1_ingest/bin/raw_reads/16_build_external_gapfill.py): Builds an external gap-fill plan for unresolved samples, preserving a fallback acquisition path. Outputs: external gap-fill plans.
- [modules/step1_ingest/bin/raw_reads/17_collect_assembled_genomes.py](../modules/step1_ingest/bin/raw_reads/17_collect_assembled_genomes.py): Collects distributed assembly outputs back into a unified staging area. Outputs: collected assembly tables/directories.
- [modules/step1_ingest/bin/raw_reads/18_qc_assembled_genomes.py](../modules/step1_ingest/bin/raw_reads/18_qc_assembled_genomes.py): QC for newly assembled genomes before they are merged into the main database. Outputs: assembly QC staging tables.
- [modules/step1_ingest/bin/raw_reads/19_merge_qc_passed_genomes.py](../modules/step1_ingest/bin/raw_reads/19_merge_qc_passed_genomes.py): Merges QC-passed assemblies into the main assembly corpus. Outputs: merged assembly sets.
- [modules/step1_ingest/bin/raw_reads/20_build_genome_database.py](../modules/step1_ingest/bin/raw_reads/20_build_genome_database.py): Consolidates genomes into `pertussis_data/bp_genomes_qc/assemblies/` so downstream stages target one assembly root. Outputs: `pertussis_data/bp_genomes_qc/assemblies/` and database manifests. Status: active.
- [modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh](../modules/step1_ingest/bin/raw_reads/21_download_missing_assemblies.sh): Downloads missing assemblies as a standalone completion step. Outputs: additional FASTAs in `pertussis_data/bp_genomes_qc/assemblies/` and download logs. Status: active.
- [modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py](../modules/step1_ingest/bin/raw_reads/22_check_genome_completeness.py): Checks genome-database completeness before deeper QC and modeling. Outputs: completeness reports used by M1. Status: active.
- [modules/step1_ingest/bin/raw_reads/23_retry_missing_assemblies.sh](../modules/step1_ingest/bin/raw_reads/23_retry_missing_assemblies.sh): Retries failed assembly downloads so transient transfer errors are separated from permanent gaps. Outputs: retry logs and refreshed assemblies. Status: active.
- [modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py](../modules/step1_ingest/bin/raw_reads/24_trace_reads_availability.py): Traces read availability from the manifest, supporting the readiness checkpoint and later read validation. Outputs: `outputs/workflow/manifest/runs.tsv`, reads availability JSON. Status: active.
- [modules/step1_ingest/bin/raw_reads/25_build_reads_download_plan.py](../modules/step1_ingest/bin/raw_reads/25_build_reads_download_plan.py): Builds the executable read-download plan from the traced linkage layer. Outputs: `outputs/workflow/reads_plan/reads_download_plan.tsv`. Status: active.
- [modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py](../modules/step1_ingest/bin/raw_reads/26_run_assembly_qc.py): Runs assembly QC so downstream Snippy planning and missingness modeling consume explicit quality statistics. Outputs: `outputs/workflow/assembly_qc/assembly_qc_stats.tsv`. Status: active.
- [modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh](../modules/step1_ingest/bin/raw_reads/27_run_snippy_batch.sh): Batch runner for contig-mode Snippy, designed for staged large-scale execution. Outputs: `outputs/workflow/snippy_ctg/` sample runs and histories. Status: active.
- [modules/step1_ingest/bin/raw_reads/28_build_snippy_ctg_plan.py](../modules/step1_ingest/bin/raw_reads/28_build_snippy_ctg_plan.py): Builds the Snippy contig plan so include/exclude decisions are explicit and auditable. Outputs: `outputs/workflow/snippy_ctg/snippy_ctg_plan.tsv`. Status: active.
- [modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh](../modules/step1_ingest/bin/raw_reads/29_run_snippy_core.sh): Aggregates completed Snippy directories into the core alignment, supporting both batch-local and cumulative rebuilds. Outputs: `outputs/workflow/phylo/core.full.aln`, completion summaries. Status: active.

### 5.3 Step2 Assembly Characterization

- [modules/step2_typing/bin/step2_01_qc_filter.py](../modules/step2_typing/bin/step2_01_qc_filter.py): Assembly QC filter that screens low-quality genomes before marker analysis. Outputs: Step2 QC-filtered tables.
- [modules/step2_typing/bin/step2_02_index_genomes.py](../modules/step2_typing/bin/step2_02_index_genomes.py): Indexes genomes for batch BLAST/MLST execution. Outputs: Step2 genome index files.
- [modules/step2_typing/bin/step2_03_run_mlst.py](../modules/step2_typing/bin/step2_03_run_mlst.py): Runs MLST so sequence type becomes a standard downstream annotation field. Outputs: Step2 MLST result tables.
- [modules/step2_typing/bin/step2_04_marker_scan_blast.py](../modules/step2_typing/bin/step2_04_marker_scan_blast.py): BLAST-based marker scan for `prn`, 23S, and related loci. Outputs: Step2 marker hit tables.
- [modules/step2_typing/bin/step2_05_merge_qc_tables.py](../modules/step2_typing/bin/step2_05_merge_qc_tables.py): Merges QC intermediates so MLST and marker layers consume one QC table. Outputs: merged Step2 QC tables.
- [modules/step2_typing/bin/step2_06_merge_mlst.py](../modules/step2_typing/bin/step2_06_merge_mlst.py): Consolidates MLST outputs into a single Step2 table. Outputs: Step2 merged MLST table.
- [modules/step2_typing/bin/step2_07_mlst_summaries.py](../modules/step2_typing/bin/step2_07_mlst_summaries.py): Summarizes MLST distributions into cohort-level descriptive outputs. Outputs: MLST summary tables.
- [modules/step2_typing/bin/step2_08_extract_marker_alleles.py](../modules/step2_typing/bin/step2_08_extract_marker_alleles.py): Extracts marker alleles and hashes them so known and novel variants are tracked consistently. Outputs: allele tables and extracted sequences.
- [modules/step2_typing/bin/step2_09_call_23s_a2047g.py](../modules/step2_typing/bin/step2_09_call_23s_a2047g.py): Specialized caller for 23S A2047G to create an explicit macrolide-resistance marker column. Outputs: A2047G call tables.
- [modules/step2_typing/bin/step2_10_merge_markers.py](../modules/step2_typing/bin/step2_10_merge_markers.py): Merges marker hits into a single table for downstream consumption. Outputs: merged marker tables.
- [modules/step2_typing/bin/step2_11_marker_summaries.py](../modules/step2_typing/bin/step2_11_marker_summaries.py): Builds marker summaries for lighter downstream consumption. Outputs: marker summary tables.
- [modules/step2_typing/bin/step2_12_build_marker_references.py](../modules/step2_typing/bin/step2_12_build_marker_references.py): Builds marker reference assets so scans and allele extraction share the same reference contract. Outputs: Step2 reference assets.
- [modules/step2_typing/bin/step2_13_joint_summaries.py](../modules/step2_typing/bin/step2_13_joint_summaries.py): Joins QC, MLST, and marker outputs into the main Step2 table. Outputs: `modules/step2_typing/outputs/bp_qc_merged_mlst_markers.tsv`. Status: active.
- [modules/step2_typing/bin/step2_14_harmonize_typing.py](../modules/step2_typing/bin/step2_14_harmonize_typing.py): Harmonizes raw marker hashes into canonical `ptxP`/`fim3`/`fhaB2400_5550`/`23S` labels, joins the frozen typing-profile registry, and emits the standardized genotype manifest consumed by the unified manifest and manuscript package. Outputs: `modules/step2_typing/outputs/bp_genotype_manifest.tsv`. Status: active.

### 5.4 Step3 Preliminary PRN And Phylogeny Prep

- [modules/step3_prn_scan/bin/step3_01_extra_summaries.py](../modules/step3_prn_scan/bin/step3_01_extra_summaries.py): Produces extra cohort summaries before deeper phylogeny and mechanism analysis. Outputs: Step3 summary tables.
- [modules/step3_prn_scan/bin/step3_10_prepare_phylogeny_manifest.py](../modules/step3_prn_scan/bin/step3_10_prepare_phylogeny_manifest.py): Prepares the phylogeny manifest used for tree sampling and balancing. Outputs: Step3/Step5 phylogeny manifest tables.
- [modules/step3_prn_scan/bin/step3_20_prn_disruption_scan.py](../modules/step3_prn_scan/bin/step3_20_prn_disruption_scan.py): First-pass `prn` disruption scan to rapidly identify candidate abnormal assemblies. Outputs: preliminary `prn` disruption calls.
- [modules/step3_prn_scan/bin/step3_21_prn_disruption_summaries.py](../modules/step3_prn_scan/bin/step3_21_prn_disruption_summaries.py): Summarizes Step3 disruption calls into cohort-level tables. Outputs: Step3 summary tables.
- [modules/step3_prn_scan/bin/step3_30_prn_trends_tables.py](../modules/step3_prn_scan/bin/step3_30_prn_trends_tables.py): Builds `prn` trend tables as an early time-series interface. Outputs: Step3 trend tables.
- [modules/step3_prn_scan/bin/step3_40_phylo_annotations.py](../modules/step3_prn_scan/bin/step3_40_phylo_annotations.py): Adds annotations to the phylogeny manifest so mechanism and lineage labels are aligned early. Outputs: annotated phylogeny tables.
- [modules/step3_prn_scan/bin/step3_50_prn_breakpoint_evidence.py](../modules/step3_prn_scan/bin/step3_50_prn_breakpoint_evidence.py): Extracts breakpoint evidence from assemblies before Step4 read validation. Outputs: breakpoint evidence tables.
- [modules/step3_prn_scan/bin/step3_51_prn_breakpoint_summaries.py](../modules/step3_prn_scan/bin/step3_51_prn_breakpoint_summaries.py): Summarizes breakpoint evidence across samples. Outputs: breakpoint summary tables.
- [modules/step3_prn_scan/bin/step3_52_extract_prn_gap_sequences.py](../modules/step3_prn_scan/bin/step3_52_extract_prn_gap_sequences.py): Extracts `prn` gap sequences to support Step4 IS scanning and gap interpretation. Outputs: extracted FASTA/gap assets.
- [modules/step3_prn_scan/bin/step3_60_results_digest.py](../modules/step3_prn_scan/bin/step3_60_results_digest.py): Produces a digest of Step3 outputs for quick review and handoff. Outputs: digest tables and reports.

### 5.5 Step4 PRN Mechanism Resolution And Read Validation

- [modules/step4_prn_validation/bin/step4_00_distributed_raw_reads_lib.sh](../modules/step4_prn_validation/bin/step4_00_distributed_raw_reads_lib.sh): Shared shell library for distributed raw-read operations across hosts. Outputs: no direct file; imported by distributed Step4 wrappers.
- [modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh](../modules/step4_prn_validation/bin/step4_00_launch_distributed_raw_reads.sh): Launches distributed read jobs for large-scale Step4 execution. Outputs: distributed job logs and shard state.
- [modules/step4_prn_validation/bin/step4_00_collect_distributed_status.sh](../modules/step4_prn_validation/bin/step4_00_collect_distributed_status.sh): Collects status across distributed jobs for centralized monitoring. Outputs: distributed status summaries.
- [modules/step4_prn_validation/bin/step4_00_sync_distributed_outputs.sh](../modules/step4_prn_validation/bin/step4_00_sync_distributed_outputs.sh): Synchronizes distributed outputs back into the main repository. Outputs: synced Step4 work/output directories.
- [modules/step4_prn_validation/bin/step4_00_rebalance_distributed_shards.py](../modules/step4_prn_validation/bin/step4_00_rebalance_distributed_shards.py): Rebalances shards across workers to reduce skew in distributed execution. Outputs: updated shard plans.
- [modules/step4_prn_validation/bin/step4_01_build_is_reference.py](../modules/step4_prn_validation/bin/step4_01_build_is_reference.py): Builds the IS reference assets used by both assembly scanning and read validation. Outputs: `modules/step4_prn_validation/references/is_elements/`. Status: active.
- [modules/step4_prn_validation/bin/step4_02_scan_prn_mechanisms.py](../modules/step4_prn_validation/bin/step4_02_scan_prn_mechanisms.py): Scans assemblies for `prn` mechanism classes as the high-throughput first pass before read validation. Outputs: `modules/step4_prn_validation/outputs/bp_prn_mechanism_calls.tsv`, related mechanism tables. Status: active.
- [modules/step4_prn_validation/bin/step4_02b_summarize_is_hits.py](../modules/step4_prn_validation/bin/step4_02b_summarize_is_hits.py): Summarizes IS hits so local evidence is elevated into mechanism-level summaries. Outputs: `bp_prn_is_hits.tsv` and summaries.
- [modules/step4_prn_validation/bin/step4_02c_score_prn_calls.py](../modules/step4_prn_validation/bin/step4_02c_score_prn_calls.py): Scores `prn` calls by confidence tier to distinguish robust assembly calls from insufficient evidence. Outputs: confidence-enhanced mechanism call tables.
- [modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables.py](../modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables.py): Original PRN summary-table builder retained for compatibility. Outputs: legacy summary tables. Status: legacy-compatible.
- [modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables_v2.py](../modules/step4_prn_validation/bin/step4_02d_build_prn_summary_tables_v2.py): Current PRN summary builder with the stabilized manuscript-facing logic and interval outputs used by the figure layer. Outputs: `outputs/bp_step4_with_ci/` and Step4 summary tables. Status: active for figure compatibility.
- [modules/step4_prn_validation/bin/step4_03a_build_validation_subset.py](../modules/step4_prn_validation/bin/step4_03a_build_validation_subset.py): Builds the read-validation subset so scarce read resources are focused on the most informative samples. Outputs: validation subset TSVs. Status: active.
- [modules/step4_prn_validation/bin/step4_03b_assess_validation_feasibility.py](../modules/step4_prn_validation/bin/step4_03b_assess_validation_feasibility.py): Validation-feasibility check before running M6 validation. Outputs: validation-feasibility JSON in `outputs/workflow/checkpoints/`. Status: active.
- [modules/step4_prn_validation/bin/step4_03c_prepare_is_reference.py](../modules/step4_prn_validation/bin/step4_03c_prepare_is_reference.py): Prepares the specialized IS reference used by `ismapper/panISa` so read-validation inputs are standardized. Outputs: Step4 validation reference assets. Status: active.
- [modules/step4_prn_validation/bin/step4_03d_build_read_validation_batch.py](../modules/step4_prn_validation/bin/step4_03d_build_read_validation_batch.py): Converts the validation subset into a batch plan consumed by download and validation runners. Outputs: validation batch TSVs. Status: active.
- [modules/step4_prn_validation/bin/step4_03e_run_is_read_validation.sh](../modules/step4_prn_validation/bin/step4_03e_run_is_read_validation.sh): Main shell entry for M6, orchestrating `ismapper`, `panISa`, and downstream parsing as an operational batch runner. Outputs: Step4 work directories, logs, and final validation tables in `modules/step4_prn_validation/outputs/`. Status: active.
- [modules/step4_prn_validation/bin/step4_03_validate_prn_with_reads.py](../modules/step4_prn_validation/bin/step4_03_validate_prn_with_reads.py): Core parser that reconciles `ismapper` and `panISa` evidence into manuscript-facing validation statuses such as `supported`, `supported_concordant`, `unresolved`, and `no_prn_is_signal_detected`. Outputs: `bp_prn_read_validation.tsv`, IS evidence table, TSD table. Status: active.
- [modules/step4_prn_validation/bin/step4_03f_hotspot_test.py](../modules/step4_prn_validation/bin/step4_03f_hotspot_test.py): Tests insertion hotspots so positional clustering becomes a statistical output rather than a visual impression. Outputs: `bp_prn_is_hotspot_results.tsv` and hotspot plots. Status: active.
- [modules/step4_prn_validation/bin/step4_04_summarize_prn_validation.py](../modules/step4_prn_validation/bin/step4_04_summarize_prn_validation.py): Summarizes read-validation results into cohort-level tables so Step4 emits both detail and overview. Outputs: `bp_prn_validation_summary.tsv` and related summaries. Status: active.
- [modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py](../modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py): Builds the missingness model used by the current M1/M2 wrapper. Outputs: `outputs/workflow/missingness_model/`. Status: active dependency.

### 5.6 Step5 Legacy Balanced-Tree Layer

- [modules/step5_phylogeny_asr/bin/step5_01_build_phylogeny_manifest.py](../modules/step5_phylogeny_asr/bin/step5_01_build_phylogeny_manifest.py): Builds the legacy Step5 phylogeny manifest that predates the workflow-native M3/M4/M5 chain. Outputs: `modules/step5_phylogeny_asr/outputs/bp_phylogeny_manifest_balanced.tsv`. Status: legacy but still referenced by manifest building.
- [modules/step5_phylogeny_asr/bin/step5_02_build_global_phylogeny.sh](../modules/step5_phylogeny_asr/bin/step5_02_build_global_phylogeny.sh): Legacy global phylogeny runner retained for historical traceability. Outputs: legacy Step5 phylogeny outputs. Status: legacy.
- [modules/step5_phylogeny_asr/bin/step5_03_reconstruct_prn_states.py](../modules/step5_phylogeny_asr/bin/step5_03_reconstruct_prn_states.py): Legacy ASR state reconstruction on the balanced tree. Outputs: legacy Step5 state tables. Status: legacy.
- [modules/step5_phylogeny_asr/bin/step5_04_count_independent_origins.py](../modules/step5_phylogeny_asr/bin/step5_04_count_independent_origins.py): Legacy independent-origin counter for the balanced-tree analysis. Outputs: `modules/step5_phylogeny_asr/outputs/bp_prn_independent_origins.tsv`. Status: legacy but still informative.
- [modules/step5_phylogeny_asr/bin/step5_05_summarize_clades.py](../modules/step5_phylogeny_asr/bin/step5_05_summarize_clades.py): Summarizes clade-level outputs from the legacy phylogeny. Outputs: legacy clade summary tables. Status: legacy.
- [modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py](../modules/step5_phylogeny_asr/bin/step5_06_build_missingness_model.py): Missingness-model builder still invoked by the current M1/M2 wrapper, retained because it is already validated in the root orchestration path. Outputs: `outputs/workflow/missingness_model/`. Status: active dependency.
- [modules/step5_phylogeny_asr/bin/step5_07_asr_sensitivity.py](../modules/step5_phylogeny_asr/bin/step5_07_asr_sensitivity.py): Legacy ASR sensitivity script from the balanced-tree era. Outputs: legacy Step5 sensitivity outputs. Status: legacy.

### 5.7 Public Health Ingestion

- [modules/public_health/bin/ph_utils.py](../modules/public_health/bin/ph_utils.py): Shared utilities for normalization, Excel parsing, and TSV output in the public-health layer. Outputs: imported by other PH scripts; no direct dataset.
- [modules/public_health/bin/ph_01_build_source_inventory.py](../modules/public_health/bin/ph_01_build_source_inventory.py): Builds the source inventory so every external input is version-tracked. Outputs: source inventory tables in `modules/public_health/outputs/`.
- [modules/public_health/bin/ph_02_normalize_country_names.py](../modules/public_health/bin/ph_02_normalize_country_names.py): Normalizes country names and ISO3 mappings so WHO/WUENIC/AMU/genomic layers can join cleanly. Outputs: `ph_country_name_map.tsv`. Status: active.
- [modules/public_health/bin/ph_03_clean_who_cases.py](../modules/public_health/bin/ph_03_clean_who_cases.py): Cleans WHO case counts into a country-year joinable table. Outputs: cleaned WHO case table. Status: active.
- [modules/public_health/bin/ph_04_clean_wuenic.py](../modules/public_health/bin/ph_04_clean_wuenic.py): Cleans WUENIC coverage data so DTP3/booster variables fit the country-year panel. Outputs: cleaned WUENIC table. Status: active.
- [modules/public_health/bin/ph_05_clean_vaccine_programs.py](../modules/public_health/bin/ph_05_clean_vaccine_programs.py): Cleans vaccine-program metadata including acellular vs whole-cell and `prn_in_vaccine`, making the exposure definition explicit. Outputs: `ph_country_program_metadata.tsv`. Status: active.
- [modules/public_health/bin/ph_06_clean_glass_amu.py](../modules/public_health/bin/ph_06_clean_glass_amu.py): Cleans GLASS AMU data as one antimicrobial-use source. Outputs: cleaned GLASS AMU table. Status: active.
- [modules/public_health/bin/ph_07_clean_esacnet_amu.py](../modules/public_health/bin/ph_07_clean_esacnet_amu.py): Cleans ESAC-Net AMU data to complement or substitute GLASS. Outputs: cleaned ESAC-Net AMU table. Status: active.
- [modules/public_health/bin/ph_08_build_country_year_master.py](../modules/public_health/bin/ph_08_build_country_year_master.py): Joins cases, vaccine programs, and AMU into the master country-year public-health table used by Step6. Outputs: `modules/public_health/outputs/ph_country_year_master.tsv`. Status: active.
- [modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py](../modules/public_health/bin/ph_09_assess_vaccine_variable_coverage.py): Checks whether vaccine variables are sufficient for downstream ecological analysis and now reads the richer formulation curation file by default. Outputs: `outputs/workflow/checkpoints/vaccine_variable_coverage_report.json`. Status: active.

### 5.8 Step6 Ecology And Transmission

- [modules/step6_epi_transmission/bin/step6_01_build_country_year_genomic_summaries.py](../modules/step6_epi_transmission/bin/step6_01_build_country_year_genomic_summaries.py): Aggregates sample-level mechanism, read support, and A2047G calls into country-year genomic summaries. Outputs: `modules/step6_epi_transmission/outputs/bp_country_year_genomic_summary.tsv`. Status: active.
- [modules/step6_epi_transmission/bin/step6_02_join_public_health.py](../modules/step6_epi_transmission/bin/step6_02_join_public_health.py): Joins genomic country-year summaries with the PH master table to create the unified model input table. Outputs: `modules/step6_epi_transmission/outputs/bp_country_year_analysis_input.tsv`. Status: active.
- [modules/step6_epi_transmission/bin/step6_03_fit_primary_models.py](../modules/step6_epi_transmission/bin/step6_03_fit_primary_models.py): Fits the primary ecological model and emits a manuscript-facing effect table while documenting that random effects were not fit under sparse data. Outputs: `bp_country_year_association_models.tsv`, `bp_country_year_model_diagnostics.tsv`. Status: active.
- [modules/step6_epi_transmission/bin/step6_04_run_sensitivity_models.py](../modules/step6_epi_transmission/bin/step6_04_run_sensitivity_models.py): Runs ecological sensitivity models across alternate filters and covariate choices. Outputs: Step6 sensitivity result tables. Status: active.
- [modules/step6_epi_transmission/bin/step6_05_run_amu_exploratory_sensitivity.py](../modules/step6_epi_transmission/bin/step6_05_run_amu_exploratory_sensitivity.py): Runs exploratory AMU sensitivity models so antimicrobial-use signals remain clearly exploratory. Outputs: AMU overlap and ridge-sensitivity tables consumed by manuscript extracts. Status: active.
- [modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers.py](../modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers.py): First-generation `R_e` estimator retained for traceability. Outputs: initial Step6 transmission outputs. Status: legacy-compatible.
- [modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers_v2.py](../modules/step6_epi_transmission/bin/step6_06_estimate_reproduction_numbers_v2.py): Revised `R_e` estimator with cleaned trajectories and quality flags used by the current figure layer. Outputs: `outputs/bp_step6_v2_fixed/bp_country_year_re_trajectories.tsv`, `bp_re_summary_statistics.tsv`. Status: active figure contract.
- [modules/step6_epi_transmission/bin/step6_06_run_re_estimation.sh](../modules/step6_epi_transmission/bin/step6_06_run_re_estimation.sh): Shell wrapper for `R_e` estimation to standardize IO locations and runtime options. Outputs: Step6 or `outputs/bp_step6_v2_fixed/` transmission tables. Status: active.
- [modules/step6_epi_transmission/bin/step6_07_fit_transmission_models.py](../modules/step6_epi_transmission/bin/step6_07_fit_transmission_models.py): Fits transmission models from `R_e` trajectories and covariates. Outputs: Step6 transmission model TSV/JSON outputs. Status: active.
- [modules/step6_epi_transmission/bin/step6_07_run_transmission_models.sh](../modules/step6_epi_transmission/bin/step6_07_run_transmission_models.sh): Shell wrapper for transmission-model fitting with a stable IO contract. Outputs: `modules/step6_epi_transmission/outputs/*.json`, `*.tsv`. Status: active.
- [modules/step6_epi_transmission/bin/step6_08_cross_validation.py](../modules/step6_epi_transmission/bin/step6_08_cross_validation.py): Performs cross-validation so exploratory models retain a basic generalization check. Outputs: Step6 cross-validation summaries. Status: active.
- [modules/step6_epi_transmission/bin/step6_09_mixed_effects_models.py](../modules/step6_epi_transmission/bin/step6_09_mixed_effects_models.py): Fits mixed-effects variants when the data structure permits more hierarchical modeling. Outputs: Step6 mixed-effects result tables. Status: active exploratory.

### 5.9 Manuscript Staging

- [manuscript/scripts/freeze/extract_key_statistics.py](../manuscript/scripts/freeze/extract_key_statistics.py): Extracts manuscript key statistics into JSON/TXT to avoid manual transcription in the text. Outputs: `manuscript/key_statistics.json`, `manuscript/key_statistics.txt`. Status: active.
- [manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py](../manuscript/scripts/freeze/ms_01_build_figure_data_extracts.py): Builds figure-ready TSV extracts and a data dictionary from Step4, Step6, and ASR outputs so the figure layer reads frozen inputs only. Outputs: `manuscript/figure_data/*.tsv`, `manuscript/figure_data_dictionary.md`. Status: active.
- [manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py](../manuscript/scripts/freeze/ms_02_export_workflow_asr_tables.py): Stages workflow-native ASR outputs into supplementary and figure-data locations with the composition-pruned ASR frame promoted to primary. Outputs: `manuscript/supplementary/Supplementary_Table_3_independent_origins.tsv`, `Supplementary_Table_6_ASR_Sensitivity.tsv`, plus Figure 3 extracts. Status: active.
- [manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py](../manuscript/scripts/freeze/ms_03_build_figure3_workflow_tree.py): Converts the composition-pruned rooted ML tree into precomputed segment and node tables so Figure 3 plotting does not depend on dynamic tree parsing inside R. Outputs: `manuscript/figure_data/figure3_workflow_tree_segments.tsv`, `figure3_workflow_tree_nodes.tsv`. Status: active.
- [manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py](../manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py): Orchestrates manuscript diagnostics bundles across ASR, validation evidence, and context audits. Status: active.
- [manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py](../manuscript/scripts/diagnostics/ms_10_build_submission_diagnostics.py): Builds the consolidated ASR, validation, and context ledgers used by the manuscript package. Outputs: `Supplementary_Table_13` through `Supplementary_Table_19`, `Supplementary_Table_21`, and `Supplementary_Table_23` to `Supplementary_Table_24`, plus associated Figure 3 and validation sidecars.
- [manuscript/scripts/diagnostics/ms_12_build_submission_evidence_summary.py](../manuscript/scripts/diagnostics/ms_12_build_submission_evidence_summary.py): Builds the submission evidence summary and lineage/origin collapse tables used by the submission package. Outputs: `Supplementary_Table_25` through `Supplementary_Table_30` and the consolidated validation matrix.
- [manuscript/scripts/review/ms_15_build_selected_country_review_report.py](../manuscript/scripts/review/ms_15_build_selected_country_review_report.py): Builds the selected-country program-history manifest, Stage 1 / Stage 2 selection tables, and evidence grid used by the current selected-country comparison framework. Outputs: `manuscript/figure_data/selected_country/*`. Status: active.
- [manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py](../manuscript/scripts/sidecars/ms_16_build_analysis_upgrade_sidecars.py): Builds the reviewer-facing sidecar tables for year-sensitivity, architecture turnover, origin-burden bridging, and the PRN-locus structural signal specificity audit. Outputs: `manuscript/figure_data/selected_country/*`, `manuscript/supplementary/Supplementary_Table_35_*.tsv` through `Supplementary_Table_38_*.tsv`. Status: active.
- [manuscript/scripts/source_data/ms_14_build_source_data_manifest.py](../manuscript/scripts/source_data/ms_14_build_source_data_manifest.py): Builds the submission-facing source-data manifest for the active Figure 1-5 and Extended Data Fig. 1-10 contract. Outputs: `manuscript/submission_data/source_data/final_source_data_manifest.tsv`. Status: active.
- [manuscript/scripts/source_data/ms_17_build_source_data_workbook.py](../manuscript/scripts/source_data/ms_17_build_source_data_workbook.py): Builds the panel/file manifests plus one Excel workbook per figure from the source-data manifest. Outputs: `manuscript/submission_data/source_data/`. Status: active.

### 5.10 Figures

- [manuscript/figures/bin/render_main.R](../manuscript/figures/bin/render_main.R): Main-figure batch runner for Figures 1-5. Outputs: `manuscript/figures/outputs/main/`. Status: active.
- [manuscript/figures/bin/render_extended_data.R](../manuscript/figures/bin/render_extended_data.R): Extended Data batch runner for Extended Data Fig. 1-10. Outputs: `manuscript/figures/outputs/extended_data/`. Status: active.
- [manuscript/figures/scripts/main/fig01_public_genome_atlas.R](../manuscript/figures/scripts/main/fig01_public_genome_atlas.R): Renders the public-genome atlas, sampling frame, and recoverable-*prn* boundary. Outputs: `manuscript/figures/outputs/main/fig01_public_genome_atlas.*`. Status: active.
- [manuscript/figures/scripts/main/fig02_prn_structural_solution_space.R](../manuscript/figures/scripts/main/fig02_prn_structural_solution_space.R): Renders the proportional *prn* locus atlas and structural-event reuse panels. Outputs: `manuscript/figures/outputs/main/fig02_prn_structural_solution_space.*`. Status: active.
- [manuscript/figures/scripts/main/fig03_repeated_origin_phylogeny.R](../manuscript/figures/scripts/main/fig03_repeated_origin_phylogeny.R): Renders the fan-tree repeated-origin analysis with outer context tracks and ASR sensitivity panels. Outputs: `manuscript/figures/outputs/main/fig03_repeated_origin_phylogeny.*`. Status: active.
- [manuscript/figures/scripts/main/fig04_country_programme_amplification.R](../manuscript/figures/scripts/main/fig04_country_programme_amplification.R): Renders country-programme amplification small multiples and selected-country trajectory contrasts. Outputs: `manuscript/figures/outputs/main/fig04_country_programme_amplification.*`. Status: active.
- [manuscript/figures/scripts/main/fig05_validation_synthesis.R](../manuscript/figures/scripts/main/fig05_validation_synthesis.R): Renders PRN specificity, support-only exposure sensitivity, missingness contrasts, and the identifiability synthesis grid. Outputs: `manuscript/figures/outputs/main/fig05_validation_synthesis.*`. Status: active.
- [manuscript/figures/scripts/lib/data_utils.R](../manuscript/figures/scripts/lib/data_utils.R): Shared data-loading and validation helpers for the figure layer. It now resolves the repository root dynamically, reads frozen manuscript-facing TSV assets, and keeps figure scripts detached from ad hoc working-directory assumptions. Outputs: imported by figure scripts; no direct file. Status: active helper.
- [manuscript/figures/scripts/lib/theme_nature.R](../manuscript/figures/scripts/lib/theme_nature.R): Shared theme, palette, and save helpers so all figures follow a consistent visual contract. Outputs: imported by figure scripts; saves to `manuscript/figures/outputs/main/`. Status: active helper.
- [manuscript/figures/scripts/extended_data/ed09_ecology_sidecar_impl.R](../manuscript/figures/scripts/extended_data/ed09_ecology_sidecar_impl.R): Renders the support-only ecology/programme diagnostic now mapped to Supplementary Figure 9 source assets. Outputs: `manuscript/figures/outputs/extended_data/Extended_Data_Fig_09_Ecology_Sidecar.*`. Status: active.

### 5.11 Tests

- [tests/test_m4_utils.py](../tests/test_m4_utils.py): CLI tests for `mask_recombination.py` and `compare_trees.py`, protecting the most critical pure-Python M4 helpers. Outputs: test-only; no persistent outputs.
- [tests/test_asr_pastml_parser.py](../tests/test_asr_pastml_parser.py): Tests whether `asr_pastml.py` correctly distinguishes strict and compatible origins, protecting the PastML interpretation layer. Outputs: test-only.
- [tests/test_m5_asr_stage1.py](../tests/test_m5_asr_stage1.py): Tests `asr_parsimony.py` and `origin_events.py` on a small fixture tree to protect core M5 logic. Outputs: test-only.
- [tests/test_step4_read_validation.py](../tests/test_step4_read_validation.py): Tests `step4_03_validate_prn_with_reads.py` and `step4_03f_hotspot_test.py`, protecting the central read-validation and hotspot logic in Step4. Outputs: test-only.

## 6. Practical Reading Order

If you want to understand the active production path rather than the full historical code surface, read in this order:

1. `workflow/bin/run_full_workflow.sh`
2. `workflow/lib/build_analysis_manifest.py` and `workflow/lib/run_foundation_checks.py`
3. `modules/step1_ingest/bin/raw_reads/24-29_*`
4. `workflow/bin/m4_phylogeny.sh` and `workflow/bin/m5_asr.sh`
5. `modules/step4_prn_validation/bin/step4_02*` and `step4_03*`
6. `modules/public_health/bin/ph_03-08_*`
7. `modules/step6_epi_transmission/bin/step6_01-07_*`
8. `manuscript/scripts/freeze/ms_01-03_*` and `manuscript/figures/bin/render_main.R`

## 7. Current Technical Debt

- The figure layer still depends on a broad frozen TSV contract under `manuscript/figure_data/`; further consolidation of repeated path helpers and staging manifests would make future renderer migrations smaller.

- The legacy Step5 balanced-tree layer still enriches parts of the manifest chain, so `modules/step5_phylogeny_asr/` cannot be removed wholesale yet.

- The root Snakemake scaffold and the shell-first production path still coexist, which is a controlled but real dual-entry architecture.

- Distributed Step4 scripts and older Step3/Step5 result paths remain because they still matter for auditability and replaying historical analyses.
