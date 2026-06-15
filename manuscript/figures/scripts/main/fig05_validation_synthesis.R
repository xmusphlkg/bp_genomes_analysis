#!/usr/bin/env Rscript
# Figure 5: Validation, specificity and phenotype synthesis
# Panels: A = frame transition, B = PRN specificity,
#          C = external phenotype bridge, D = validation metrics,
#          E = genome-call phenotype tiers
# Context-only AMR/lineage, missingness and country-glyph diagnostics remain in
# supplementary/source-data outputs rather than the main synthesis figure.

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
library(ggforce)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(stringr)
library(tidyr)

# setting -----------------------------------------------------------

fig5_base_size <- FIGURE_BASE_SIZE
fig5_tag_size  <- FIGURE_TAG_SIZE
fig5_annot_size <- FIGURE_ANNOT_SIZE

# data loading -------------------------------------------------------

flow_summary <- safe_load(file.path(FIGURE_DATA_DIR, "fig05_evidence_chain_summary.tsv"), "Fig. 5 evidence chain summary") %>%
  mutate(
    n = as.numeric(n),
    parent_n = suppressWarnings(as.numeric(parent_n))
  )

specificity_locus_order <- c("Fim3", "Fim2", "SphB1", "Vag8", "TcfA", "BapC", "Phg", "BrkA", "FHA/FhaB", "prn")

controls <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "prn_specificity_negative_control.tsv"), "PRN specificity controls") %>%
  filter(table_scope == "global_overlap_frame") %>%
  mutate(
    n_locus_interpretable = as.numeric(n_locus_interpretable),
    n_signal_positive = as.numeric(n_signal_positive),
    signal_positive_fraction_among_interpretable = as.numeric(signal_positive_fraction_among_interpretable),
    has_locus_frame = !is.na(n_locus_interpretable) & n_locus_interpretable > 0,
    plot_fraction = if_else(has_locus_frame, replace_na(signal_positive_fraction_among_interpretable, 0), 0),
    locus_label = if_else(locus == "prn", "prn", locus_label),
    locus_label = factor(locus_label, levels = specificity_locus_order),
    locus_rank = match(as.character(locus_label), specificity_locus_order),
    group = case_when(
      locus == "prn" ~ "prn",
      !has_locus_frame ~ "No frame",
      TRUE ~ "Comparator"
    )
  )

controls_label_xmax <- max(controls$plot_fraction, na.rm = TRUE) * 1.18
controls_label_dat <- controls %>%
  mutate(
    signal_label = case_when(
      !has_locus_frame ~ "no frame",
      plot_fraction == 0 ~ "0%",
      TRUE ~ percent(plot_fraction, accuracy = 0.1)
    ),
    label_x = pmax(
      plot_fraction + controls_label_xmax * 0.013,
      controls_label_xmax * if_else(has_locus_frame, 0.022, 0.040)
    ),
    label_colour = if_else(group == "prn", unname(npg_colors["red"]), FIGURE_MUTED_TEXT)
  )

max_comparator_fraction <- controls %>%
  filter(group == "Comparator", has_locus_frame) %>%
  summarise(max_fraction = max(plot_fraction, na.rm = TRUE), .groups = "drop") %>%
  pull(max_fraction)

biology_external <- safe_load(file.path(FIGURE_DATA_DIR, "biology_bridge_external_context.tsv"), "biology bridge external context") %>%
  mutate(
    n_total_or_frame = as.numeric(n_total_or_frame),
    n_prn_producing = as.numeric(n_prn_producing),
    n_prn_deficient = as.numeric(n_prn_deficient),
    n_amr_or_a2047g = as.numeric(n_amr_or_a2047g)
  )

phenotype_tiers <- safe_load(
  file.path(BASE_DIR, "manuscript", "supplementary", "Supplementary_Table_10_Event_Class_Phenotype_Evidence_Tiers.tsv"),
  "phenotype evidence tiers"
) %>%
  mutate(sample_count = as.numeric(sample_count))

validation_summary <- safe_load(
  file.path(FIGURE_DATA_DIR, "caller_validation_sensitivity_summary.tsv"),
  "caller validation sensitivity summary"
) %>%
  mutate(
    n_compared_rows = suppressWarnings(as.numeric(n_compared_rows)),
    n_supporting_rows = suppressWarnings(as.numeric(n_supporting_rows)),
    fraction = suppressWarnings(as.numeric(fraction))
  )

# ===========================================================================
# Panel A: Frame transition
# ===========================================================================

flow_nodes <- tibble::tribble(
  ~label, ~x, ~y, ~fill,
  "Retained archive\n2,406", 1.00, 2.40, FIGURE_PANEL_FILL,
  "Interpretable prn\n1,325", 2.10, 2.95, unname(npg_colors["teal"]),
  "Missing/uncertain\n1,081", 2.10, 1.55, FIGURE_LIGHT_GREY,
  "Resolved disrupted\n577", 3.30, 3.25, unname(npg_colors["red"]),
  "Intact boundary\n748", 3.30, 2.65, FIGURE_MID_GREY,
  "Bounds\n1,081", 3.30, 1.55, FIGURE_LIGHT_GREY,
  "Tier 1a bridge\n544", 4.55, 3.25, unname(npg_colors["green"]),
  "Other tiers\n33", 4.55, 2.65, unname(npg_colors["teal"])
) %>%
  mutate(
    xmin = x - 0.47,
    xmax = x + 0.47,
    ymin = y - 0.25,
    ymax = y + 0.25
  )

flow_edges <- tibble::tribble(
  ~x, ~y, ~xend, ~yend,
  1.47, 2.40, 1.63, 2.95,
  1.47, 2.40, 1.63, 1.55,
  2.57, 2.95, 2.83, 3.25,
  2.57, 2.95, 2.83, 2.65,
  2.57, 1.55, 2.83, 1.55,
  3.77, 3.25, 4.08, 3.25,
  3.77, 3.25, 4.08, 2.65
)

pA <- ggplot() +
  geom_segment(
    data = flow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(2.4, "pt"), type = "closed"),
    linewidth = 0.34,
    colour = FIGURE_RULE_COLOUR,
    lineend = "round"
  ) +
  geom_rect(
    data = flow_nodes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "white",
    linewidth = 0.28
  ) +
  geom_text(
    data = flow_nodes,
    aes(x = x, y = y, label = label),
    size = 1.58,
    lineheight = 0.82,
    fontface = "bold",
    colour = FIGURE_TEXT_COLOUR
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.44, 5.10), ylim = c(1.06, 3.66), clip = "off") +
  labs(x = NULL, y = NULL, tag = "a") +
  theme_void(base_size = fig5_base_size) +
  theme(
    plot.margin = margin(3, 3, 3, 3)
  )

# ===========================================================================
# Panel B: PRN specificity as an outlying recurrent structural signal
# ===========================================================================

pB <- ggplot(controls, aes(y = locus_rank)) +
  annotate(
    "rect",
    xmin = 0, xmax = max_comparator_fraction,
    ymin = -Inf, ymax = Inf,
    fill = FIGURE_PANEL_FILL,
    colour = NA
  ) +
  geom_rect(
    data = controls %>% filter(has_locus_frame),
    aes(xmin = 0, xmax = plot_fraction, ymin = locus_rank - 0.30, ymax = locus_rank + 0.30, fill = group),
    colour = "white", linewidth = 0.18
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
  geom_vline(xintercept = max_comparator_fraction, linewidth = 0.26, linetype = "dashed", colour = FIGURE_MUTED_TEXT) +
  scale_colour_identity() +
  scale_x_continuous(labels = percent,
    limits = c(0, controls_label_xmax),
    expand = c(0, 0)) +
  scale_y_continuous(
    breaks = seq_along(specificity_locus_order),
    labels = specificity_locus_order,
    expand = expansion(add = 0.45)
  ) +
  scale_fill_manual(values = c("prn" = unname(npg_colors["red"]), "Comparator" = FIGURE_GREY, "No frame" = "white"), guide = "none") +
  labs(
    x = "Signal-positive fraction among interpretable loci",
    y = NULL,
    tag = "b"
  ) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.y = element_text(size = 5.2)
  )

# ===========================================================================
# Panel C/E: Published phenotype bridge and present genome-call tiers
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

pC_external <- ggplot(phenotype_bridge, aes(prn_deficient_fraction, cohort_label, fill = assay_class)) +
  geom_segment(aes(x = 0, xend = prn_deficient_fraction, yend = cohort_label),
    colour = FIGURE_RULE_COLOUR, linewidth = 0.55) +
  geom_point(shape = 21, colour = FIGURE_INK, size = 2.25, stroke = 0.22) +
  geom_text(aes(label = count_label), hjust = -0.08, size = fig5_annot_size, colour = FIGURE_INK) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, max(phenotype_bridge$prn_deficient_fraction, na.rm = TRUE) * 1.42),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "ELISA" = unname(npg_colors["green"]),
      "Immunoblot" = unname(npg_colors["blue"]),
      "Published expression" = unname(npg_colors["peach"])
    ),
    guide = "none"
  ) +
  labs(
    x = "PRN-negative in published phenotype cohorts",
    y = NULL,
    tag = "c"
  ) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.y = element_text(size = 4.8)
  )

tier_levels <- c(
  "Tier 1a junction-class phenotype bridge",
  "Tier 1b lesion-class bridge only",
  "Tier 2 genome-disruption plausible",
  "Tier 3 genome-only disruption",
  "Boundary only, not expression-proven"
)

route_family_order <- c("Intact boundary", "Other disruption", "Rearrangement", "IS481 insertion")
tier_axis_order <- c("Tier 1a", "Tier 1b", "Tier 2", "Tier 3", "Intact")

route_tier <- phenotype_tiers %>%
  mutate(
    phenotype_evidence_tier = factor(phenotype_evidence_tier, levels = tier_levels),
    tier_short = recode(as.character(phenotype_evidence_tier),
      `Tier 1a junction-class phenotype bridge` = "Tier 1a",
      `Tier 1b lesion-class bridge only` = "Tier 1b",
      `Tier 2 genome-disruption plausible` = "Tier 2",
      `Tier 3 genome-only disruption` = "Tier 3",
      `Boundary only, not expression-proven` = "Intact",
      .default = as.character(phenotype_evidence_tier)
    ),
    route_family = case_when(
      phenotype_evidence_tier == "Boundary only, not expression-proven" ~ "Intact boundary",
      str_detect(event_subcategory, regex("IS481", ignore_case = TRUE)) ~ "IS481 insertion",
      str_detect(event_subcategory, regex("Inversion|rearrangement", ignore_case = TRUE)) ~ "Rearrangement",
      TRUE ~ "Other disruption"
    )
  ) %>%
  group_by(route_family, tier_short) %>%
  summarise(
    event_classes = n(),
    genomes = sum(sample_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tier_short = factor(tier_short, levels = tier_axis_order),
    route_family = factor(route_family, levels = route_family_order),
    x_rank = match(as.character(tier_short), tier_axis_order),
    y_rank = match(as.character(route_family), route_family_order),
    label_colour = if_else(genomes >= 100, "white", FIGURE_INK),
    label = comma(genomes)
  )

tier_fill <- c(
  "Tier 1a" = unname(npg_colors["green"]),
  "Tier 1b" = unname(npg_colors["teal"]),
  "Tier 2" = unname(npg_colors["blue"]),
  "Tier 3" = unname(npg_colors["peach"]),
  "Intact" = FIGURE_LIGHT_GREY
)

route_grid <- tidyr::crossing(
  x_rank = seq_along(tier_axis_order),
  y_rank = seq_along(route_family_order)
)

pC_tiers <- ggplot() +
  geom_tile(
    data = route_grid,
    aes(x = x_rank, y = y_rank),
    fill = FIGURE_PANEL_FILL,
    colour = "white",
    linewidth = 0.28,
    width = 0.92,
    height = 0.82
  ) +
  geom_point(
    data = route_tier,
    aes(x = x_rank, y = y_rank, size = genomes, fill = tier_short),
    shape = 21,
    colour = FIGURE_INK,
    stroke = 0.22,
    alpha = 0.96
  ) +
  geom_text(
    data = route_tier,
    aes(x = x_rank, y = y_rank, label = label, colour = label_colour),
    size = 1.45,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_size_area(max_size = 8.2, guide = "none") +
  scale_fill_manual(values = tier_fill, guide = "none") +
	scale_x_continuous(
	    breaks = seq_along(tier_axis_order),
	    labels = c("T1a\nbridge", "T1b\nlesion", "T2\nplaus.", "T3\ngenome", "Intact\nboundary"),
	    limits = c(0.45, length(tier_axis_order) + 0.55),
	    expand = c(0, 0)
	  ) +
  scale_y_continuous(
    breaks = seq_along(route_family_order),
    labels = route_family_order,
    limits = c(0.45, length(route_family_order) + 0.55),
    expand = c(0, 0)
  ) +
	  labs(x = "Phenotype-evidence tier", y = "Route family", tag = "e") +
  theme_nature(base_size = fig5_base_size) +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(size = 4.35, lineheight = 0.84),
    axis.text.y = element_text(size = 4.55, face = "bold"),
    axis.title = element_text(size = 4.8),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 3, 2, 3)
  )

pC <- pC_external / pC_tiers +
  plot_layout(heights = c(1.18, 0.78))

# ===========================================================================
# Panel D: Validation dashboard
# ===========================================================================

validation_key <- tibble::tribble(
  ~evidence_layer, ~comparison_or_risk, ~display_label, ~display_class, ~invert_fraction,
  "PRN status", "repo PRN status versus published PRN status or phenotype annotation", "Published PRN status", "External concordance", FALSE,
  "broad mechanism", "repo broad mechanism versus published broad mechanism among PRN-negative concordant rows", "Broad mechanism", "Crosswalk bounded", FALSE,
  "manual/read junction evidence tiers", "event-class confidence tiers from read-backed TSD, public long-read or hybrid, assembly/rule support or unresolved validation", "Read/TSD event classes", "Read-backed support", FALSE,
  "full caller threshold grid", "HSP identity 85-95%, locus coverage 90-99% and relaxed/default/strict IS-support profiles", "Gap1043 dominant\nacross caller grid", "Threshold robustness", FALSE,
  "intact controls", "read-validation audit of intact control rows", "False signal in\nintact controls", "Negative control", FALSE
)

validation_metrics <- validation_summary %>%
  inner_join(validation_key, by = c("evidence_layer", "comparison_or_risk")) %>%
  mutate(display_order = match(display_label, validation_key$display_label)) %>%
  arrange(display_order) %>%
  mutate(
    plot_fraction = if_else(invert_fraction, 1 - fraction, fraction),
    denominator_label = case_when(
      evidence_layer == "full caller threshold grid" ~ paste0(n_supporting_rows, "/", n_compared_rows, " grids"),
      evidence_layer == "intact controls" ~ paste0(n_supporting_rows, "/", n_compared_rows, " false"),
      TRUE ~ paste0(n_supporting_rows, "/", n_compared_rows)
    ),
    value_label = case_when(
      evidence_layer == "intact controls" ~ "0/7\nfalse signals",
      TRUE ~ paste0(percent(plot_fraction, accuracy = 0.1), "\n", denominator_label)
    ),
	    display_label = recode(display_label,
	      `Published PRN status` = "PRN status",
	      `Read/TSD event classes` = "Read/TSD event classes",
	      `Gap1043 dominant\nacross caller grid` = "Caller grid",
	      `False signal in\nintact controls` = "Intact controls"
	    ),
    display_class = factor(
      display_class,
      levels = c("External concordance", "Crosswalk bounded", "Read-backed support", "Threshold robustness", "Negative control")
    ),
    display_label = factor(display_label, levels = rev(display_label)),
    label_x = case_when(
      plot_fraction >= 0.92 ~ plot_fraction - 0.025,
      plot_fraction <= 0.08 ~ plot_fraction + 0.035,
      TRUE ~ plot_fraction + 0.035
    ),
    label_hjust = if_else(plot_fraction >= 0.92, 1, 0)
  )

validation_colours <- c(
  "External concordance" = unname(npg_colors["green"]),
  "Crosswalk bounded" = unname(npg_colors["peach"]),
  "Read-backed support" = unname(npg_colors["teal"]),
  "Threshold robustness" = unname(npg_colors["blue"]),
  "Negative control" = FIGURE_GREY
)

pD <- ggplot(validation_metrics, aes(plot_fraction, display_label)) +
  geom_segment(
    aes(x = 0, xend = plot_fraction, yend = display_label),
    colour = FIGURE_RULE_COLOUR,
    linewidth = 0.48
  ) +
  geom_point(
    aes(fill = display_class),
    shape = 21,
    size = 2.45,
    colour = FIGURE_INK,
    stroke = 0.20,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, label = value_label, hjust = label_hjust),
    size = 1.58,
    lineheight = 0.84,
    colour = FIGURE_TEXT_COLOUR
  ) +
  scale_fill_manual(
    values = validation_colours,
    guide = "none"
  ) +
  scale_x_continuous(labels = percent, limits = c(0, 1.10), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(
    x = "Fraction within validation layer",
    y = NULL,
    tag = "d"
  ) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.y = element_text(size = 4.8, lineheight = 0.86),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

# ===========================================================================
# Assembly and save
# ===========================================================================

fig5_layout <- "
AD
BC
"

fig5 <- pA + pB + pC + pD +
  plot_layout(
    design = fig5_layout,
    widths = c(1, 1),
    heights = c(0.88, 1.22)
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = fig5_tag_size, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_nc_pdf(fig5, "fig05_validation_synthesis.pdf", height = NC_MAX_HEIGHT * 0.84)
save_nc_png(fig5, "fig05_validation_synthesis.png", height = NC_MAX_HEIGHT * 0.84)
