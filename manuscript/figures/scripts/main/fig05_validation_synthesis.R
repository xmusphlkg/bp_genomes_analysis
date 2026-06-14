#!/usr/bin/env Rscript
# Figure 5: Validation, specificity and phenotype synthesis
# Panels: A = denominator/evidence flow, B = PRN specificity,
#          C = phenotype bridge, D = validation evidence stack
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
# Panel A: Denominator and evidence flow
# ===========================================================================

flow_count <- setNames(flow_summary$n, flow_summary$stage_id)

sankey_paths <- tibble::tribble(
  ~archive_frame, ~endpoint_frame, ~genome_call, ~evidence_boundary, ~flow_group, ~value,
  "Retained\n2,406", "Interpretable\n1,325", "Resolved\n577", "Tier 1a\n544", "Tier 1a", flow_count[["tier1a_phenotype_bridge"]],
  "Retained\n2,406", "Interpretable\n1,325", "Resolved\n577", "Tier 1b-3\n33", "Tier 1b-3", flow_count[["other_disrupted_phenotype_tiers"]],
  "Retained\n2,406", "Interpretable\n1,325", "Intact boundary\n748", "Boundary\n748", "Boundary", flow_count[["intact_boundary"]],
  "Retained\n2,406", "Missing/uncertain\n1,081", "Missing/uncertain\n1,081", "Missingness\nbounds\n1,081", "Missing/uncertain", flow_count[["noninterpretable_uncertain"]]
) %>%
  mutate(value = as.numeric(value))

sankey_dat <- ggforce::gather_set_data(
  sankey_paths,
  c("archive_frame", "endpoint_frame", "genome_call", "evidence_boundary"),
  id_name = "flow_id"
) %>%
  mutate(
    x = factor(
      x,
      levels = c("archive_frame", "endpoint_frame", "genome_call", "evidence_boundary")
    )
)

sankey_flow_colours <- c(
  "Tier 1a" = alpha(unname(npg_colors["green"]), 0.70),
  "Tier 1b-3" = alpha(unname(npg_colors["teal"]), 0.50),
  "Boundary" = alpha(FIGURE_MID_GREY, 0.42),
  "Missing/uncertain" = alpha(FIGURE_GREY, 0.42)
)

pA <- ggplot(sankey_dat, aes(x = x, id = flow_id, split = y, value = value)) +
  geom_parallel_sets(
    aes(fill = flow_group),
    alpha = 0.86,
    axis.width = 0.13,
    sep = 0.035,
    strength = 0.56,
    colour = NA
  ) +
  geom_parallel_sets_axes(
    axis.width = 0.13,
    fill = FIGURE_PANEL_FILL,
    colour = "white",
    linewidth = 0.25
  ) +
  geom_parallel_sets_labels(
    angle = 0,
    size = 1.55,
    lineheight = 0.82,
    colour = FIGURE_INK,
    fontface = "bold"
  ) +
  scale_fill_manual(values = sankey_flow_colours, guide = "none") +
  scale_x_discrete(
    labels = c(
      archive_frame = "Archive frame",
      endpoint_frame = "Endpoint frame",
      genome_call = "Genome call",
      evidence_boundary = "Evidence boundary"
    )
  ) +
  labs(x = NULL, y = NULL, tag = "a") +
  theme_void(base_size = fig5_base_size) +
  theme(
    axis.text.x = element_text(size = 5.1, face = "bold", colour = FIGURE_TEXT_COLOUR),
    plot.margin = margin(2, 2, 2, 2)
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
# Panel C: Published phenotype bridge and present genome-call tiers
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
    labels = c("Tier 1a\nbridge", "Tier 1b\nlesion", "Tier 2\nplausible", "Tier 3\ngenome-only", "Intact\nboundary"),
    limits = c(0.45, length(tier_axis_order) + 0.55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_along(route_family_order),
    labels = route_family_order,
    limits = c(0.45, length(route_family_order) + 0.55),
    expand = c(0, 0)
  ) +
  labs(x = "Present genome-call phenotype tier", y = "Route family") +
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
# Panel D: Validation evidence stack
# ===========================================================================

validation_key <- tibble::tribble(
  ~evidence_layer, ~comparison_or_risk, ~display_label, ~display_class, ~invert_fraction,
  "PRN status", "repo PRN status versus published PRN status or phenotype annotation", "Published PRN status", "External concordance", FALSE,
  "broad mechanism", "repo broad mechanism versus published broad mechanism among PRN-negative concordant rows", "Broad mechanism", "Crosswalk bounded", FALSE,
  "manual/read junction evidence tiers", "event-class confidence tiers from read-backed TSD, public long-read or hybrid, assembly/rule support or unresolved validation", "Read/TSD event classes", "Read-backed support", FALSE,
  "full caller threshold grid", "HSP identity 85-95%, locus coverage 90-99% and relaxed/default/strict IS-support profiles", "Gap1043 dominant\nacross caller grid", "Threshold robustness", FALSE,
  "intact controls", "read-validation audit of intact control rows", "False signal in\nintact controls", "Negative control", FALSE
)

validation_stack <- validation_summary %>%
  inner_join(validation_key, by = c("evidence_layer", "comparison_or_risk")) %>%
  mutate(
    plot_fraction = if_else(invert_fraction, 1 - fraction, fraction),
    denominator_label = case_when(
      evidence_layer == "full caller threshold grid" ~ paste0(n_supporting_rows, "/", n_compared_rows, " grids"),
      evidence_layer == "intact controls" ~ paste0(n_supporting_rows, "/", n_compared_rows, " false"),
      TRUE ~ paste0(n_supporting_rows, "/", n_compared_rows)
    ),
    value_label = case_when(
      evidence_layer == "intact controls" ~ denominator_label,
      TRUE ~ paste0(percent(plot_fraction, accuracy = 0.1), "\n", denominator_label)
    ),
    label_x = if_else(plot_fraction >= 0.85, plot_fraction - 0.04, plot_fraction + 0.04),
    label_hjust = if_else(plot_fraction >= 0.85, 1, 0),
    display_label = factor(display_label, levels = rev(validation_key$display_label)),
    display_class = factor(
      display_class,
      levels = c("External concordance", "Crosswalk bounded", "Read-backed support", "Threshold robustness", "Negative control")
    )
  )

validation_colours <- c(
  "External concordance" = unname(npg_colors["green"]),
  "Crosswalk bounded" = unname(npg_colors["peach"]),
  "Read-backed support" = unname(npg_colors["teal"]),
  "Threshold robustness" = unname(npg_colors["blue"]),
  "Negative control" = FIGURE_GREY
)

pD <- ggplot(validation_stack, aes(plot_fraction, display_label)) +
  geom_segment(
    aes(x = 0, xend = plot_fraction, yend = display_label),
    colour = FIGURE_RULE_COLOUR,
    linewidth = 0.55
  ) +
  geom_point(
    aes(fill = display_class),
    shape = 21,
    colour = FIGURE_INK,
    size = 2.65,
    stroke = 0.22,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = label_x,
      label = value_label,
      hjust = label_hjust
    ),
    size = 1.58,
    lineheight = 0.84,
    colour = FIGURE_INK
  ) +
  geom_vline(xintercept = 0, colour = FIGURE_INK, linewidth = 0.28) +
  scale_fill_manual(
    values = validation_colours,
    guide = "none"
  ) +
  scale_x_continuous(
    labels = percent,
    limits = c(0, 1.12),
    breaks = c(0, 0.5, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Observed validation fraction",
    y = NULL,
    tag = "d"
  ) +
  theme_nature(base_size = fig5_base_size) +
  theme(
    axis.text.y = element_text(size = 4.8, lineheight = 0.86),
    panel.grid.major.y = element_blank()
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
