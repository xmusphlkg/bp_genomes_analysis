#!/usr/bin/env Rscript
# Figure 5: Validation, specificity and identifiability synthesis
# Merged from former Fig 5 (specificity/bias bounds) + Fig 6 (synthesis)
# Panels: A = PRN specificity, B = phenotype bridge, C = AMR/lineage overlay,
#          D = missingness bounds, E = evidence grid, F = identifiability ledger

get_script_dir <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) normalizePath(frame$ofile, winslash = "/", mustWork = FALSE) else ""
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) return(dirname(tail(frame_files, 1)))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  normalizePath(".", winslash = "/", mustWork = FALSE)
}
script_dir <- get_script_dir()

source(file.path(script_dir, "..", "lib", "theme_nature.R"))
source(file.path(script_dir, "..", "lib", "data_utils.R"))
source(file.path(script_dir, "..", "lib", "manuscript_rebuild_helpers.R"))

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(stringr)
library(tidyr)

# setting -----------------------------------------------------------

fig5_base_size <- FIGURE_BASE_SIZE
fig5_tag_size  <- FIGURE_TAG_SIZE
fig5_legend_text_size <- FIGURE_LEGEND_TEXT_SIZE
fig5_legend_title_size <- FIGURE_LEGEND_TITLE_SIZE
fig5_annot_size <- FIGURE_ANNOT_SIZE
fig5_delta_limits <- c(-1, 0.75)

# data loading -------------------------------------------------------

controls <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "prn_specificity_negative_control.tsv"), "PRN specificity controls") %>%
  filter(table_scope == "global_overlap_frame") %>%
  mutate(
    n_locus_interpretable = as.numeric(n_locus_interpretable),
    n_signal_positive = as.numeric(n_signal_positive),
    signal_positive_fraction_among_interpretable = as.numeric(signal_positive_fraction_among_interpretable),
    has_locus_frame = !is.na(n_locus_interpretable) & n_locus_interpretable > 0,
    plot_fraction = if_else(has_locus_frame, replace_na(signal_positive_fraction_among_interpretable, 0), 0),
    locus_label = if_else(locus == "prn", "prn", locus_label),
    locus_label = factor(locus_label, levels = rev(c("prn", "FHA/FhaB", "BrkA", "Phg", "BapC", "Fim2", "Fim3", "TcfA", "Vag8", "SphB1"))),
    group = case_when(
      locus == "prn" ~ "prn",
      !has_locus_frame ~ "No frame",
      TRUE ~ "Comparator"
    )
  )

controls_label_xmax <- max(controls$plot_fraction, na.rm = TRUE) * 1.15
controls_label_dat <- controls %>%
  mutate(
    signal_label = case_when(
      !has_locus_frame ~ "no frame",
      plot_fraction == 0 ~ "0%",
      TRUE ~ percent(plot_fraction, accuracy = 0.1)
    ),
    label_x = pmax(
      plot_fraction + controls_label_xmax * 0.012,
      controls_label_xmax * if_else(has_locus_frame, 0.018, 0.035)
    ),
    label_colour = if_else(group == "prn", unname(npg_colors["red"]), FIGURE_MUTED_TEXT)
  )

biology_external <- safe_load(file.path(FIGURE_DATA_DIR, "biology_bridge_external_context.tsv"), "biology bridge external context") %>%
  mutate(
    n_total_or_frame = as.numeric(n_total_or_frame),
    n_prn_producing = as.numeric(n_prn_producing),
    n_prn_deficient = as.numeric(n_prn_deficient),
    n_amr_or_a2047g = as.numeric(n_amr_or_a2047g)
  )

biology_internal <- safe_load(file.path(FIGURE_DATA_DIR, "biology_bridge_internal_sensitivity.tsv"), "biology bridge internal sensitivity") %>%
  mutate(
    n_records = as.numeric(n_records),
    n_prn_interpretable = as.numeric(n_prn_interpretable),
    n_prn_disrupted = as.numeric(n_prn_disrupted),
    prn_disrupted_fraction_among_interpretable = as.numeric(prn_disrupted_fraction_among_interpretable),
    n_gap1043 = as.numeric(n_gap1043),
    n_ptxP3 = as.numeric(n_ptxP3),
    n_fim3_1 = as.numeric(n_fim3_1)
  )

dr <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "selected_country_dr_missingness_summary.tsv"), "DR missingness summary") %>%
  mutate(across(starts_with("delta_"), as.numeric))

evidence <- safe_load(data_paths$rr_evidence_grid, "cross-country evidence grid")


# ===========================================================================
# Panel A: PRN specificity ???outlying recurrent structural signal
# ===========================================================================

pA <- ggplot(controls, aes(y = locus_label)) +
  geom_col(
    data = controls %>% filter(has_locus_frame),
    aes(x = plot_fraction, fill = group),
    width = 0.65, colour = "white", linewidth = 0.18
  ) +
  geom_point(
    data = controls %>% filter(!has_locus_frame),
    aes(x = 0),
    shape = 21, fill = "white", colour = FIGURE_MUTED_TEXT, size = 1.55, stroke = 0.25
  ) +
  geom_text(
    data = controls_label_dat,
    aes(x = label_x, label = signal_label, colour = label_colour),
    hjust = 0, size = fig5_annot_size
  ) +
  scale_colour_identity() +
  scale_x_continuous(labels = percent,
    limits = c(0, controls_label_xmax),
    expand = c(0, 0)) +
  scale_fill_manual(values = c("prn" = unname(npg_colors["red"]), "Comparator" = FIGURE_GREY, "No frame" = "white"), guide = "none") +
  labs(x = "Signal-positive fraction among interpretable loci", y = NULL) +
  theme_nature(base_size = fig5_base_size) +
  theme(axis.text.y = element_text(size = 5.2))

# ===========================================================================
# Panel B: Published phenotype bridge
# ===========================================================================

phenotype_bridge <- biology_external %>%
  filter(
    str_detect(interpretive_role, "phenotype-genotype bridge|direct PRN-expression bridge"),
    !is.na(n_total_or_frame),
    !is.na(n_prn_deficient),
    n_total_or_frame > 0
  ) %>%
  mutate(
    prn_deficient_fraction = n_prn_deficient / n_total_or_frame,
    cohort_label = recode(context_id,
      japan_otsuka_2012_expression_bridge = "Japan\nOtsuka 2012",
      usa_pawloski_2014_expression_bridge = "USA\nPawloski 2014",
      australia_lam_2014_expression_bridge = "Australia\nLam 2014",
      europe_eupert_1998_2015_expression_bridge = "Europe\nBarkoff 2019",
      slovenia_2002_2020_expression = "Slovenia\nBarkoff 2024",
      belgium_2000_2023_expression_bridge = "Belgium\nMartini 2026",
      .default = country_or_region
    ),
    assay_class = case_when(
      str_detect(evidence_type, regex("ELISA", ignore_case = TRUE)) ~ "ELISA",
      str_detect(evidence_type, regex("immunoblot|Western", ignore_case = TRUE)) ~ "Immunoblot",
      TRUE ~ "Published expression"
    ),
    count_label = paste0(n_prn_deficient, "/", n_total_or_frame),
    cohort_label = factor(cohort_label, levels = rev(cohort_label[order(prn_deficient_fraction)]))
  )

pB <- ggplot(phenotype_bridge, aes(prn_deficient_fraction, cohort_label, fill = assay_class)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.18) +
  geom_text(aes(label = count_label), hjust = -0.08, size = fig5_annot_size, colour = FIGURE_INK) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, max(phenotype_bridge$prn_deficient_fraction, na.rm = TRUE) * 1.25),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "ELISA" = unname(npg_colors["green"]),
      "Immunoblot" = unname(npg_colors["blue"]),
      "Published expression" = unname(npg_colors["peach"])
    ),
    name = NULL
  ) +
  labs(x = "PRN-negative in published phenotype cohorts", y = NULL) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.y = element_text(size = 4.8),
    legend.position = "bottom",
    legend.text = element_text(size = fig5_legend_text_size)
  )

# ===========================================================================
# Panel C: AMR and lineage-background overlay
# ===========================================================================

amr_overlay <- biology_internal %>%
  filter(
    analysis_layer == "internal_genome_context",
    (group_variable == "marker_23s_status" & group_value %in% c("MS", "MR_A2047G")) |
      (group_variable == "published_sublineage" & group_value %in% c(
        "United States Sublineage 1", "MR-MT28-PG1", "MR-MT28-PG2", "MR-MT28-others", "MS-MT28"
      ))
  ) %>%
  mutate(
    display_label = recode(group_value,
      MS = "23S susceptible\nmarker",
      MR_A2047G = "23S A2047G\nmarker",
      `United States Sublineage 1` = "US\nSublineage 1",
      `MR-MT28-PG1` = "MR-MT28\nPG1",
      `MR-MT28-PG2` = "MR-MT28\nPG2",
      `MR-MT28-others` = "MR-MT28\nother",
      `MS-MT28` = "MS-MT28",
      .default = group_value
    ),
    context_class = if_else(group_variable == "marker_23s_status", "23S marker", "Published lineage"),
    count_label = paste0(n_prn_disrupted, "/", n_prn_interpretable),
    display_label = factor(display_label, levels = rev(c(
      "23S susceptible\nmarker", "23S A2047G\nmarker", "US\nSublineage 1",
      "MR-MT28\nPG1", "MR-MT28\nPG2", "MR-MT28\nother", "MS-MT28"
    )))
  )

pC <- ggplot(amr_overlay, aes(prn_disrupted_fraction_among_interpretable, display_label, fill = context_class)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.18) +
  geom_text(aes(label = count_label), hjust = -0.08, size = fig5_annot_size, colour = FIGURE_INK) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, 1.28),
    breaks = c(0, 0.3, 0.6, 0.9),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c("23S marker" = unname(npg_colors["peach"]), "Published lineage" = unname(npg_colors["purple"])),
    name = NULL
  ) +
  labs(x = "Genome-defined prn disruption", y = NULL) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.y = element_text(size = 4.8),
    legend.position = "bottom",
    legend.text = element_text(size = fig5_legend_text_size)
  )

# ===========================================================================
# Panel D: Missingness bounds (multiple estimators per country)
# ===========================================================================

dr_long <- dr %>%
  filter(country_iso3 %in% c("USA", "NZL", "JPN", "AUS")) %>%
  transmute(
    country_iso3 = factor(country_iso3, levels = rev(c("USA", "NZL", "JPN", "AUS"))),
    Naive = delta_naive_prevalence,
    `IPW cap20` = delta_ipw_cap20_prevalence,
    `Overlap observed` = delta_overlap_observed_prevalence,
    AIPW = delta_aipw_prevalence
  ) %>%
  pivot_longer(-country_iso3, names_to = "estimator", values_to = "delta") %>%
  mutate(
    delta_plot = squish(delta, range = fig5_delta_limits),
    clipped = !is.na(delta) & delta != delta_plot,
    clip_hjust = if_else(delta < fig5_delta_limits[1], -0.08, 1.08),
    clip_label = if_else(delta < fig5_delta_limits[1], "< -100 pp", "> 75 pp")
  )

dr_range <- dr_long %>%
  group_by(country_iso3) %>%
  summarise(xmin = min(delta_plot, na.rm = TRUE), xmax = max(delta_plot, na.rm = TRUE), .groups = "drop")

pD <- ggplot(dr_long, aes(delta_plot, country_iso3, fill = estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
  geom_segment(data = dr_range, aes(x = xmin, xend = xmax, y = country_iso3, yend = country_iso3),
    inherit.aes = FALSE, colour = FIGURE_RULE_COLOUR, linewidth = 0.44) +
  geom_point(shape = 21, colour = FIGURE_INK, size = 2.15, stroke = 0.16, alpha = 0.92, na.rm = TRUE) +
  geom_text(
    data = dr_long %>% filter(clipped),
    aes(label = clip_label, hjust = clip_hjust),
    colour = FIGURE_MUTED_TEXT, size = 1.65, show.legend = FALSE, na.rm = TRUE
  ) +
  scale_x_continuous(
    labels = label_number(scale = 100, suffix = " pp"),
    limits = fig5_delta_limits,
    breaks = c(-1, -0.5, 0, 0.5),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_fill_manual(values = c(
    "Naive" = unname(npg_colors["red"]),
    "IPW cap20" = unname(npg_colors["blue"]),
    "Overlap observed" = unname(npg_colors["peach"]),
    "AIPW" = unname(npg_colors["green"])
  ), name = NULL) +
  labs(x = "Post-minus-pre change (percentage points)", y = NULL) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = fig5_legend_text_size)
  )

# ===========================================================================
# Panel E: Cross-country evidence grid (from former Fig 6)
# ===========================================================================

country_levels <- c("USA", "NZL", "JPN", "AUS", "GBR", "CHN", "FRA", "BRA")

glyph_dat <- evidence %>%
  filter(country_iso3 %in% country_levels) %>%
  transmute(
    country_iso3,
    Direction = direction_label(prevalence_direction),
    `Repeated origin` = evidence_label(repeated_origin_evidence),
    `Structural reuse` = evidence_label(structural_reuse_evidence),
    Amplification = evidence_label(amplification_pattern),
    Interpretation = evidence_label(final_interpretation_tier)
  ) %>%
  pivot_longer(-country_iso3, names_to = "dimension", values_to = "status") %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = rev(country_levels)),
    dimension = factor(dimension, levels = c("Direction", "Repeated origin", "Structural reuse", "Amplification", "Interpretation")),
    status_class = case_when(
      status %in% c("Upward", "Strong") ~ "Strong",
      status %in% c("Downward", "Bounded", "Limited") ~ "Bounded",
      status %in% c("Context", "None", "No change") ~ "Context",
      TRUE ~ "Uncertain"
    ),
    dimension_label = recode(as.character(dimension),
      Direction = "Direction",
      `Repeated origin` = "Repeated\norigin",
      `Structural reuse` = "Structural\nreuse",
      Amplification = "Amplification",
      Interpretation = "Interpretation"
    )
  )

pE <- ggplot(glyph_dat, aes(dimension, country_iso3)) +
  geom_tile(fill = FIGURE_PANEL_FILL, colour = "white", linewidth = 0.25) +
  geom_point(aes(fill = status_class), shape = 21, size = 3.0, colour = FIGURE_INK, stroke = 0.22) +
  scale_fill_manual(
    values = c("Strong" = unname(npg_colors["green"]), "Bounded" = unname(npg_colors["teal"]),
               "Context" = FIGURE_LIGHT_GREY, "Uncertain" = FIGURE_GREY),
    name = "Evidence"
  ) +
  scale_x_discrete(labels = setNames(unique(glyph_dat$dimension_label), unique(glyph_dat$dimension))) +
  labs(x = NULL, y = NULL) +
  theme_nature_matrix(base_size = fig5_base_size) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 4.9, lineheight = 0.85),
    axis.text.y = element_text(face = "bold", size = 5.5),
    legend.position = "bottom",
    legend.text = element_text(size = fig5_legend_text_size)
  )

# ===========================================================================
# Panel F: Identifiability ledger (from former Fig 6)
# ===========================================================================

ledger <- tibble::tribble(
  ~object, ~status, ~x,
  "Recoverable-locus prevalence", "Identified", 1,
  "Product-aware exposure index", "Bounded", 2,
  "Reported-case scenario model", "Bounded", 2,
  "Country-epoch contrasts", "Bounded", 2,
  "Repeated origin in analysed trees", "Bounded", 2,
  "Circulation-wide prevalence", "Needs new data", 3,
  "Global vaccine-programme causal effect", "Needs new data", 3,
  "Direct PRN antigen expression", "Needs new data", 3
) %>%
  mutate(
    status = factor(status, levels = c("Identified", "Bounded", "Needs new data")),
    object = factor(object, levels = rev(object))
  )

pF <- ggplot(ledger, aes(x, object)) +
  geom_segment(aes(x = 1, xend = x, yend = object), colour = FIGURE_RULE_COLOUR, linewidth = 0.42) +
  geom_point(aes(fill = status), shape = 21, size = 3.0, colour = FIGURE_INK, stroke = 0.22) +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("Identified", "Bounded", "Needs new data"),
    limits = c(0.72, 3.28),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = identifiability_colors, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.x = element_text(face = "bold", size = 5.5),
    axis.text.y = element_text(size = 5)
  )

# ===========================================================================
# Assembly and save
# ===========================================================================

fig5_layout <- "
AABB
CCDD
EEFF
"

fig5 <- pA + pB + pC + pD + pE + free(pF) +
  plot_layout(
    design = fig5_layout,
    widths = c(1, 1, 1, 1),
    heights = c(0.90, 0.90, 1.0)
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = fig5_tag_size, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_nc_pdf(fig5, "fig05_validation_synthesis.pdf", height = NC_MAX_HEIGHT)
save_nc_png(fig5, "fig05_validation_synthesis.png", height = NC_MAX_HEIGHT)
