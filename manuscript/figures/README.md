# Figures Workspace

This directory contains the active static figure pipeline for the Communications Biology submission source package.

## Active Figure Set

Main figures:

- `main/fig01_public_genome_atlas.R`
- `main/fig02_prn_structural_solution_space.R`
- `main/fig03_repeated_origin_phylogeny.R`
- `main/fig04_country_programme_amplification.R`
- `main/fig05_validation_synthesis.R`

Supplementary Figure source scripts:

- `extended_data/ed01_sample_frame_reconciliation.R`
- `extended_data/ed02_year_composition_sensitivity.R`
- `extended_data/ed03_tree_representativeness.R`
- `extended_data/ed04_asr_robustness.R`
- `extended_data/ed05_validation_package_anchors.R`
- `extended_data/ed06_architecture_turnover.R`
- `extended_data/ed07_origin_burden_bridge.R`
- `extended_data/ed08_prn_specificity_audit.R`
- `extended_data/ed09_ecology_sidecar.R`
- `extended_data/ed10_usa_focal_country.R`
- `extended_data/ed11_identifiability_dag.R`
- `extended_data/ed12_country_year_missingness_audit.R`
- `extended_data/ed13_fig01_context_panels.R`

The output filenames still use `Extended_Data_Fig_*` internally because the source renderer predates the Communications Biology retargeting. In the manuscript and final upload package these panels are labelled Supplementary Figures 1-13.

## Rendering

Run the main figures:

```bash
Rscript manuscript/figures/bin/render_main.R
```

Run the Supplementary Figure source set:

```bash
Rscript manuscript/figures/bin/render_extended_data.R
```

## Outputs

- Main figures: `manuscript/figures/outputs/main/`
- Supplementary Figure source images: `manuscript/figures/outputs/extended_data/`

The old compatibility renderers and duplicate output aliases have been removed. `manuscript/text/commsbio_supplementary_information.md` is now the active Supplementary Information source.
