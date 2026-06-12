# bp_step3

Step3 is downstream analysis using the merged table produced in `step2_typing`.

## Inputs

- `../step2_typing/outputs/bp_qc_merged_mlst_markers.tsv` (QC + MLST + marker alleles + 23S call)

## Outputs

Written to `outputs/`.

- Step3A (extra summaries): cross-tabs and diversity metrics.
- Step3B (phylogeny): manifest + optional MashTree-based Newick tree.

## Run

### Step3A

```bash
conda run -n ncbi_ds --no-capture-output bash bin/run_step3A.sh
```

### Step3B

```bash
# If mash/mashtree are missing, the runner prints install instructions
conda run -n ncbi_ds --no-capture-output bash bin/run_step3B.sh

# Defaults (can be overridden via env vars):
# - stratified sampling by mlst_st,year with MAX_GENOMES=800
# - set MAX_GENOMES=0 for the full cohort (can be slow)
# Examples:
#   GROUP_COLS=mlst_st,year MAX_GENOMES=1000 PER_GROUP=10 conda run -n ncbi_ds --no-capture-output bash bin/run_step3B.sh
#   MAX_GENOMES=0 conda run -n ncbi_ds --no-capture-output bash bin/run_step3B.sh
```

### Step3C (prn disruption scan)

Detect whether `prn` looks intact vs fragmented (multi-HSP) vs partial by re-BLASTing the `prn` query against each assembly and summarizing query coverage + HSP count.

```bash
conda run -n ncbi_ds --no-capture-output bash bin/run_step3C.sh
```

Key outputs:

- `outputs/bp_prn_disruption_calls.tsv`
- `outputs/bp_prn_call_counts.tsv`
- `outputs/bp_prn_call_by_year.tsv`, `..._by_country.tsv`, `..._by_mlst_st.tsv`, `..._by_year_mlst_st.tsv`
- `outputs/bp_qc_merged_mlst_markers_prn.tsv` (merged table with extra `prn_call` columns)

### Step3D (paper-ready prn trend tables)

Turn `prn_call` into trend tables suitable for plotting/writing: overall by-year trend, by (year×ST) trend with minimum sample-size filtering, Top ST ranking, and cross-tabs vs `23s_A2047G_call` and `ptxP`.

```bash
conda run -n ncbi_ds --no-capture-output bash bin/run_step3D.sh
```

Key outputs:

- `outputs/bp_prn_trend_by_year_clean.tsv`
- `outputs/bp_prn_trend_by_year_mlst_st_min20.tsv` (threshold via `MIN_GROUP_N=...`)
- `outputs/bp_prn_top_sts_rank.tsv` and `outputs/bp_prn_top_sts_top20.tsv`
- `outputs/bp_prn_vs_23s_overall.tsv`, `outputs/bp_prn_vs_23s_by_year.tsv`
- `outputs/bp_prn_vs_ptxP_overall.tsv`, `outputs/bp_prn_vs_ptxP_by_year.tsv`
- `outputs/bp_step3D_prn_trends_summary.txt`

### Step3E (phylogeny annotation table)

Create a metadata table aligned to the MashTree sample IDs so you can color/label the tree in iTOL/ggtree with `prn_call`, `year`, `country`, `mlst_st`, `23S`.

```bash
conda run -n ncbi_ds --no-capture-output bash bin/run_step3E.sh
```

Key output:

- `outputs/bp_phylo_annotations.tsv`

### Step3F (results digest)

Generate an auto-written Markdown digest summarizing the key Step3 findings and a checklist for figures.

```bash
conda run -n ncbi_ds --no-capture-output bash bin/run_step3F_digest.sh
```

Key output:

- `outputs/bp_results_digest.md`

### Step3G (prn breakpoint evidence; hardening disrupted calls)

Re-scan only `prn_call=disrupted_multi_hsp` genomes with a more detailed BLAST output that includes subject contig coordinates. Classify whether the multi-HSP pattern is more consistent with:

- **insertion-like** (same contig, a clear subject-gap much larger than query-gap; often IS-sized),
- **within-contig** (multi-HSP but no clear gap),
- **fragmented-contigs** (HSPs on different contigs; assembly fragmentation plausible).

```bash
conda run -n ncbi_ds --no-capture-output bash bin/run_step3F.sh
```

Key outputs:

- `outputs/bp_prn_breakpoint_evidence.tsv`
- `outputs/bp_prn_bp_category_counts.tsv` and `..._by_year_mlst_st.tsv`
- `outputs/bp_prn_bp_subject_gap_summary.txt` (gap length distribution summary)
- `outputs/bp_prn_bp_subject_gap_counts_insertion_like.tsv`
- `outputs/bp_prn_insertion_gap_plus_flanks.fasta` (gap sequence + flanks; for BLAST vs IS elements)
- `outputs/bp_prn_insertion_gap_plus_flanks.tsv` (metadata for the extracted sequences)
