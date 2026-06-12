#!/usr/bin/env Rscript
# Figure 2: constrained structural solution space for prn disruption
# Deep rebuild: information-dense panels inspired by Nature Genetics / Cell style
# Layout: 5 panels (A-E) focused on the structural-reuse claim.

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
library(cowplot)

# setting -----------------------------------------------------------

# Unified typography (consistent with fig01)
fig2_base_size <- FIGURE_BASE_SIZE
fig2_tag_size  <- FIGURE_TAG_SIZE
fig2_legend_text_size <- FIGURE_LEGEND_TEXT_SIZE
fig2_legend_title_size <- FIGURE_LEGEND_TITLE_SIZE
fig2_repel_size <- 1.75
fig2_annot_size <- FIGURE_ANNOT_SIZE    # in-plot annotation text
fig2_top_event_limit <- 6

fig2_event_label <- function(prn_event_id, mechanism_call = NA_character_) {
  mechanism_call <- rep_len(as.character(mechanism_call), length(prn_event_id))
  base <- event_label(prn_event_id)
  case_when(
    str_detect(prn_event_id, "gap1042") & str_detect(mechanism_call, "is481") ~ "IS481 gap1042",
    str_detect(prn_event_id, "gap1042") ~ "Other gap1042",
    TRUE ~ base
  )
}

fig2_architecture_colors <- c(
  architecture_colors,
  "IS481 gap1042" = unname(npg_colors["orange"]),
  "Other gap1042" = unname(npg_colors["green"])
)

# data loading -------------------------------------------------------

reuse <- safe_load(data_paths$rr_structure_reuse, "selected-country structure reuse") %>%
  mutate(
    country_sample_count = as.numeric(country_sample_count),
    sample_count = as.numeric(sample_count),
    country_total_disrupted = as.numeric(country_total_disrupted),
    event = fig2_event_label(prn_event_id, prn_mechanism_call),
    mechanism = mechanism_label(prn_mechanism_call)
  )

concentration <- safe_load(file.path(FIGURE_DATA_DIR, "structural_event_concentration.tsv"), "structural concentration") %>%
  mutate(across(c(dominant_event_share, top3_share, effective_number,
                  null_dominant_event_share_mean, null_top3_share_mean, null_effective_number_mean), as.numeric))

study_weight <- safe_load(file.path(FIGURE_DATA_DIR, "study_dependence", "structure_reuse_block_reweighted.tsv"), "study reweighted structure reuse") %>%
  mutate(
    dominant_event_share = as.numeric(dominant_event_share),
    top3_share = as.numeric(top3_share),
    effective_number = as.numeric(effective_number)
  )

mechanism_summary <- safe_load(data_paths$prn_mechanisms, "PRN mechanism summary") %>%
  mutate(
    sample_count = as.numeric(sample_count),
    event_count = as.numeric(event_count),
    country_count = as.numeric(country_count),
    year_min = as.numeric(year_min),
    year_max = as.numeric(year_max)
  )

junction_confidence <- safe_load(file.path(FIGURE_DATA_DIR, "prn_junction_confidence_matrix.tsv"), "PRN junction confidence matrix") %>%
  mutate(
    sample_count = as.numeric(sample_count),
    country_count = as.numeric(country_count),
    year_min = as.numeric(year_min),
    year_max = as.numeric(year_max),
    breakpoint_left = as.numeric(breakpoint_left),
    breakpoint_right = as.numeric(breakpoint_right),
    supporting_read_count = as.numeric(supporting_read_count),
    supporting_validation_rows = as.numeric(supporting_validation_rows),
    tsd_supported_validation_rows = as.numeric(tsd_supported_validation_rows),
    confidence_label = case_when(
      str_detect(confidence_tier, "tier_1") ~ "TSD+reads",
      str_detect(confidence_tier, "tier_2") ~ "Long-read",
      str_detect(confidence_tier, "tier_3") ~ "Assembly",
      str_detect(confidence_tier, "tier_4") ~ "Unresolved",
      TRUE ~ "Not audited"
    ),
    confidence_label = factor(
      confidence_label,
      levels = c("TSD+reads", "Long-read", "Assembly", "Unresolved", "Not audited")
    )
  )

structural_grammar <- safe_load(file.path(FIGURE_DATA_DIR, "structural_grammar_evidence.tsv"), "structural grammar evidence") %>%
  mutate(
    dominant_event_share = as.numeric(dominant_event_share),
    dominant_event_collapsed_burden = as.numeric(dominant_event_collapsed_burden)
  )

target_site <- safe_load(file.path(FIGURE_DATA_DIR, "is481_target_site_accessibility.tsv"), "IS481 target-site accessibility") %>%
  mutate(
    reference_start_1based = as.numeric(reference_start_1based),
    reference_end_1based = as.numeric(reference_end_1based),
    locus_length_bp = as.numeric(locus_length_bp),
    observed_gap1043_breakpoint_left = as.numeric(observed_gap1043_breakpoint_left),
    observed_gap1043_breakpoint_right = as.numeric(observed_gap1043_breakpoint_right),
    distance_from_observed_breakpoint_to_nearest_exact_target_bp =
      as.numeric(distance_from_observed_breakpoint_to_nearest_exact_target_bp),
    exact_ACTAGG_or_reverse_complement_count =
      as.numeric(exact_ACTAGG_or_reverse_complement_count)
  )

event_evidence <- safe_load(file.path(FIGURE_DATA_DIR, "prn_event_evidence_manifest.tsv"), "PRN event evidence manifest") %>%
  mutate(
    sample_count = as.numeric(sample_count),
    country_count = as.numeric(country_count),
    year_min = as.numeric(year_min),
    year_max = as.numeric(year_max),
    insertion_subject_gap_bp = as.numeric(insertion_subject_gap_bp),
    read_locus_start = as.numeric(read_locus_start),
    read_locus_end = as.numeric(read_locus_end),
    breakpoint_left = as.numeric(breakpoint_left),
    breakpoint_right = as.numeric(breakpoint_right),
    supporting_read_count = as.numeric(supporting_read_count),
    supporting_validation_rows = as.numeric(supporting_validation_rows),
    event = fig2_event_label(prn_event_id, mechanism_call),
    mechanism = mechanism_label(mechanism_call)
  )

within_origin <- safe_load(file.path(FIGURE_DATA_DIR, "within_origin_structural_concentration.tsv"), "within-origin concentration") %>%
  mutate(
    n_disrupted_descendant_tips = as.numeric(n_disrupted_descendant_tips),
    n_unique_events = as.numeric(n_unique_events),
    dominant_event_count = as.numeric(dominant_event_count),
    dominant_event_share = as.numeric(dominant_event_share),
    country_count = as.numeric(country_count)
  )

confidence_shapes <- c(
  "TSD+reads" = 21,
  "Long-read" = 23,
  "Assembly" = 22,
  "Unresolved" = 4,
  "Not audited" = 1
)

confidence_colours <- c(
  "TSD+reads" = unname(npg_colors["green"]),
  "Long-read" = unname(npg_colors["blue"]),
  "Assembly" = FIGURE_DARK_GREY,
  "Unresolved" = unname(npg_colors["peach"]),
  "Not audited" = FIGURE_MID_GREY
)

fig2_null_model_colors <- c(
  "Equal-event" = figure_discrete_at(2),
  "Accessibility-weighted" = figure_discrete_at(5)
)

fig2_study_weight_colors <- c(
  "Naive" = figure_discrete_at(1),
  "Study-weighted" = figure_discrete_at(2),
  "Drop largest block" = figure_discrete_at(5)
)

gap1043_conf <- junction_confidence %>%
  filter(prn_event_id == "prn_evt_coding_disrupted_is481__is481__gap1043") %>%
  slice_head(n = 1)

target_callout <- target_site %>%
  filter(locus == "prn") %>%
  slice_head(n = 1) %>%
  mutate(
    junction_x = ((observed_gap1043_breakpoint_left + observed_gap1043_breakpoint_right) / 2) -
      reference_start_1based + 1,
    callout_label = paste0(
      "ACTAGG/TSD\n",
      distance_from_observed_breakpoint_to_nearest_exact_target_bp,
      " bp; ",
      if_else(nrow(gap1043_conf) > 0, as.character(gap1043_conf$supporting_read_count[[1]]), "NA"),
      " reads"
    )
  )


# ===========================================================================
# Panel A: prn reference feature and recurrent junction evidence
# ===========================================================================

prn_locus <- target_site %>%
  filter(locus == "prn") %>%
  slice_head(n = 1) %>%
  transmute(
    xmin = 1,
    xmax = locus_length_bp,
    y = 4.34,
    label = paste0("BP1054 / prn; pertactin autotransporter (", comma(locus_length_bp), " bp)"),
    reference_start_1based,
    reference_end_1based
  )

prn_reference_start <- prn_locus$reference_start_1based[[1]]
prn_locus_length <- prn_locus$xmax[[1]]

event_axis <- tibble::tribble(
  ~prn_event_id, ~event_class, ~track_y,
  "prn_evt_coding_disrupted_is481__is481__gap1043", "IS481 gap1043", 3.28,
  "prn_evt_rearrangement__within_contig__cov58", "Rearr. cov58", 2.67,
  "prn_evt_rearrangement__within_contig__cov91", "Rearr. cov91", 2.06,
  "prn_evt_other_disruption__insertion_like__gap1042", "Other insertion-like", 1.45
)

event_tracks <- event_evidence %>%
  inner_join(event_axis, by = "prn_event_id") %>%
  mutate(
    junction_x = case_when(
      breakpoint_coordinate_basis == "read_reference" & !is.na(breakpoint_left) & !is.na(breakpoint_right) ~
        ((breakpoint_left + breakpoint_right) / 2) - prn_reference_start + 1,
      !is.na(read_locus_start) & !is.na(read_locus_end) ~
        ((read_locus_start + read_locus_end) / 2) - prn_reference_start + 1,
      TRUE ~ NA_real_
    ),
    junction_x = pmin(pmax(junction_x, 1), prn_locus_length),
    event_colour = unname(fig2_architecture_colors[event_class]),
    track_label = paste0(event_class, " (n=", comma(sample_count), ")")
  ) %>%
  filter(!is.na(junction_x))

junction_guides <- event_tracks %>%
  distinct(junction_x) %>%
  arrange(junction_x) %>%
  mutate(junction_label = paste0(comma(round(junction_x)), " bp"))

gene_arrow <- tibble::tibble(
  x = c(1, prn_locus_length - 85, prn_locus_length, prn_locus_length - 85, 1),
  y = c(4.18, 4.18, 4.34, 4.50, 4.50)
)

axis_tracks <- event_tracks %>%
  arrange(desc(track_y)) %>%
  transmute(track_y, track_label)

target_callout_a <- target_callout %>%
  mutate(callout_label = paste0("ACTAGG TSD\n", distance_from_observed_breakpoint_to_nearest_exact_target_bp, " bp; ",
                                if_else(nrow(gap1043_conf) > 0, as.character(gap1043_conf$supporting_read_count[[1]]), "NA"),
                                " reads"))

pA <- ggplot() +
  geom_segment(
    data = event_tracks,
    aes(x = 1, xend = prn_locus_length, y = track_y, yend = track_y),
    colour = FIGURE_PANEL_FILL, linewidth = 0.28
  ) +
  geom_segment(
    data = junction_guides,
    aes(x = junction_x, xend = junction_x, y = 1.18, yend = 4.08),
    linetype = "22", linewidth = 0.18, colour = FIGURE_MID_GREY
  ) +
  geom_text(
    data = junction_guides,
    aes(x = junction_x, y = 4.04, label = junction_label),
    size = 1.50, colour = FIGURE_MUTED_TEXT, vjust = 1.1
  ) +
  geom_polygon(
    data = gene_arrow,
    aes(x = x, y = y),
    fill = alpha(unname(npg_colors["teal"]), 0.22),
    colour = FIGURE_INK, linewidth = 0.15
  ) +
  geom_segment(
    data = target_callout_a,
    aes(x = junction_x, xend = junction_x, y = 4.50, yend = 4.68),
    inherit.aes = FALSE, linewidth = 0.16, colour = FIGURE_MUTED_TEXT
  ) +
  geom_label(
    data = target_callout_a,
    aes(x = 2010, y = 4.68, label = callout_label),
    inherit.aes = FALSE, size = 1.48, lineheight = 0.78,
    linewidth = 0.12, label.padding = unit(1.0, "pt"),
    fill = "white", colour = FIGURE_INK
  ) +
  geom_segment(
    data = event_tracks,
    aes(x = junction_x, xend = junction_x, y = track_y - 0.16, yend = track_y + 0.16, colour = event_colour),
    linewidth = 0.42, alpha = 0.92
  ) +
  geom_point(
    data = event_tracks,
    aes(x = junction_x, y = track_y, fill = event_class, size = sample_count),
    shape = 21, stroke = 0.2, colour = FIGURE_INK
  ) +
  geom_text(
    data = prn_locus,
    aes(x = 85, y = 4.34, label = label),
    inherit.aes = FALSE, hjust = 0, size = 1.52, colour = FIGURE_INK
  ) +
  scale_colour_identity() +
  scale_size_area(max_size = 3.7, guide = "none") +
  scale_fill_manual(values = fig2_architecture_colors, na.value = FIGURE_GREY, guide = "none") +
  scale_y_continuous(
    NULL,
    breaks = axis_tracks$track_y,
    labels = axis_tracks$track_label,
    limits = c(1.12, 4.82),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    "prn coding sequence coordinate (bp)",
    limits = c(1, prn_locus_length),
    breaks = c(1, 500, 1000, 1500, 2000, 2500, prn_locus_length),
    labels = comma,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.0, hjust = 1, lineheight = 0.86),
    axis.text.x = element_text(size = 5, colour = FIGURE_MUTED_TEXT),
    axis.title.x = element_text(size = 5.4, colour = FIGURE_MUTED_TEXT, margin = margin(t = 1.5)),
    panel.grid.major.x = element_line(colour = FIGURE_PANEL_FILL, linewidth = 0.14)
  )

# ===========================================================================
# Panel B: Global event burden (horizontal bar) + evidence tier annotation
# ===========================================================================

global_arch <- reuse %>%
  filter(!str_detect(event, "Insufficient"), !is.na(sample_count)) %>%
  distinct(prn_event_id, event, mechanism, sample_count) %>%
  group_by(event) %>%
  summarise(
    mechanism = first(mechanism),
    sample_count = sum(sample_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sample_count)) %>%
  slice_head(n = fig2_top_event_limit) %>%
  mutate(event = factor(event, levels = rev(unique(event))))

# Merge junction-confidence tier from the sidecar
event_tier <- junction_confidence %>%
  mutate(event = fig2_event_label(prn_event_id, mechanism_call)) %>%
  filter(event %in% levels(global_arch$event)) %>%
  group_by(event) %>%
  summarise(
    confidence_label = first(confidence_label),
    confidence_tier = first(confidence_tier),
    n_countries = max(country_count, na.rm = TRUE),
    read_support = sum(supporting_read_count, na.rm = TRUE),
    tsd_rows = sum(tsd_supported_validation_rows, na.rm = TRUE),
    .groups = "drop"
  )

global_arch <- global_arch %>%
  left_join(event_tier %>% select(event, confidence_label, n_countries), by = "event")

pB <- ggplot(global_arch, aes(sample_count, event)) +
  geom_col(aes(fill = event), width = 0.68, colour = "white", linewidth = 0.18) +
  geom_text(aes(label = paste0(comma(sample_count), "; ", n_countries, " ctry")),
            hjust = -0.02, size = fig2_annot_size) +
  # Junction-confidence tier markers on left
  geom_point(aes(x = -12, shape = confidence_label, colour = confidence_label),
             size = 1.65, fill = "white", stroke = 0.32) +
  scale_x_continuous(labels = comma,
                     limits = c(-20, max(global_arch$sample_count, na.rm = TRUE) * 1.78),
                     expand = c(0, 0)) +
  scale_fill_manual(values = fig2_architecture_colors, na.value = FIGURE_GREY, guide = "none") +
  scale_shape_manual(values = confidence_shapes, name = "Junction tier", drop = TRUE) +
  scale_colour_manual(values = confidence_colours, name = "Junction tier", drop = TRUE) +
  guides(
    shape = guide_legend(override.aes = list(fill = "white")),
    fill = "none"
  ) +
  labs(x = "Resolved disrupted genomes", y = NULL) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    legend.position = "none",
    legend.title = element_text(size = fig2_legend_title_size),
    legend.text = element_text(size = fig2_legend_text_size)
  )


# ===========================================================================
# Panel C: Country-by-architecture reuse matrix (bubble matrix)
# ===========================================================================

top_events <- global_arch %>% pull(event) %>% as.character() %>% unique()
country_order <- c("USA", "NZL", "JPN", "AUS", "GBR", "CHN", "FRA", "BRA", "CZE")

composition <- reuse %>%
  filter(country_iso3 %in% country_order, event %in% top_events) %>%
  group_by(country_iso3, event) %>%
  summarise(n = sum(country_sample_count, na.rm = TRUE), .groups = "drop") %>%
  group_by(country_iso3) %>%
  mutate(share = n / sum(n, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(country_iso3 = country_order, event = top_events, fill = list(n = 0, share = 0)) %>%
  mutate(
    country_iso3 = factor(country_iso3, levels = country_order),
    event = factor(event, levels = rev(top_events))
  )

country_totals_c <- reuse %>%
  filter(country_iso3 %in% country_order, !str_detect(event, "Insufficient")) %>%
  group_by(country_iso3) %>%
  summarise(total_disrupted = sum(country_sample_count, na.rm = TRUE), .groups = "drop") %>%
  mutate(country_iso3 = factor(country_iso3, levels = country_order))

country_axis_labels <- setNames(country_order, country_order)
country_axis_labels[as.character(country_totals_c$country_iso3)] <-
  paste0(as.character(country_totals_c$country_iso3), "\n", comma(country_totals_c$total_disrupted))

pC <- ggplot(composition, aes(country_iso3, event)) +
  geom_point(aes(size = n, fill = event), shape = 21, colour = FIGURE_INK, stroke = 0.14, alpha = 0.92) +
  scale_x_discrete(labels = country_axis_labels) +
  scale_size_area(max_size = 5.8, breaks = c(1, 10, 50, 200), labels = comma, name = "Genomes") +
  scale_fill_manual(values = fig2_architecture_colors, na.value = FIGURE_GREY, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    axis.text.x = element_text(face = "bold", size = 5.1, lineheight = 0.82),
    legend.position = "bottom",
    legend.title = element_text(size = fig2_legend_title_size),
    legend.text = element_text(size = fig2_legend_text_size)
  )

# ===========================================================================
# Panel D: Concentration vs null models
# ===========================================================================

null_base <- concentration %>%
  filter(scope == "overall", mechanism_group == "all",
         null_model %in% c("equal_probability_multinomial_over_observed_unique_events",
                           "mutational_accessibility_weighted_multinomial")) %>%
  select(null_model, dominant_event_share, top3_share,
         null_dominant_event_share_mean, null_top3_share_mean) %>%
  distinct()

observed_metric <- null_base %>%
  slice_head(n = 1) %>%
  select(dominant_event_share, top3_share) %>%
  pivot_longer(everything(), names_to = "metric_raw", values_to = "value") %>%
  mutate(
    metric = recode(metric_raw,
      dominant_event_share = "Dominant event",
      top3_share = "Top-three events"
    ),
    metric = factor(metric, levels = c("Dominant event", "Top-three events"))
  )

null_metric <- null_base %>%
  transmute(
    null_model = recode(null_model,
      equal_probability_multinomial_over_observed_unique_events = "Equal-event",
      mutational_accessibility_weighted_multinomial = "Accessibility-weighted"
    ),
    `Dominant event` = null_dominant_event_share_mean,
    `Top-three events` = null_top3_share_mean
  ) %>%
  pivot_longer(c(`Dominant event`, `Top-three events`), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("Dominant event", "Top-three events")),
    null_model = factor(null_model, levels = c("Equal-event", "Accessibility-weighted"))
  )

pD <- ggplot() +
  geom_segment(
    data = observed_metric,
    aes(x = 0, xend = value, y = metric, yend = metric),
    colour = FIGURE_RULE_COLOUR, linewidth = 0.55
  ) +
  geom_point(
    data = null_metric,
    aes(value, metric, shape = null_model, colour = null_model),
    fill = "white", size = 2.3, stroke = 0.42
  ) +
  geom_point(
    data = observed_metric,
    aes(value, metric),
    shape = 21, fill = figure_discrete_at(1), colour = FIGURE_INK, size = 2.7, stroke = 0.16
  ) +
  geom_text(
    data = observed_metric,
    aes(value, metric, label = percent(value, accuracy = 1)),
    hjust = -0.22, size = fig2_annot_size, colour = FIGURE_INK
  ) +
  scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0.02, 0.12))) +
  scale_shape_manual(values = c("Equal-event" = 21, "Accessibility-weighted" = 24), name = "Null model") +
  scale_colour_manual(values = fig2_null_model_colors, name = "Null model") +
  guides(colour = guide_legend(override.aes = list(fill = "white"))) +
  labs(x = "Share of resolved disrupted genomes", y = NULL) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = fig2_legend_text_size),
    legend.title = element_text(size = fig2_legend_title_size),
    axis.text.y = element_text(size = 5, lineheight = 0.85)
  )

# ===========================================================================
# Panel E: Study-block stress test (before/after comparison)
# ===========================================================================

weight_metric <- study_weight %>%
  filter(scope == "overall", mechanism_group == "all",
         row_type %in% c("current_naive_reference", "study_block_equalized", "drop_largest_block_naive")) %>%
  mutate(row_type = recode(row_type,
    current_naive_reference = "Naive",
    study_block_equalized = "Study-weighted",
    drop_largest_block_naive = "Drop largest block"
  )) %>%
  arrange(match(row_type, c("Naive", "Study-weighted", "Drop largest block"))) %>%
  distinct(row_type, .keep_all = TRUE) %>%
  select(row_type, dominant_event_share, top3_share) %>%
  pivot_longer(-row_type, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric,
      dominant_event_share = "Dominant event",
      top3_share = "Top-three events"
    ),
    metric = factor(metric, levels = c("Dominant event", "Top-three events")),
    row_type = factor(row_type, levels = c("Naive", "Study-weighted", "Drop largest block")),
    label_vjust = case_when(
      metric == "Top-three events" & row_type == "Naive" ~ 1.55,
      metric == "Top-three events" & row_type == "Drop largest block" ~ -1.15,
      TRUE ~ -1.1
    ),
    label_hjust = case_when(
      metric == "Top-three events" & row_type == "Naive" ~ 1.15,
      metric == "Top-three events" & row_type == "Drop largest block" ~ 0.35,
      TRUE ~ 0.5
    )
  )

study_weight_note <- structural_grammar %>%
  filter(evidence_layer == "study_block_stress_test",
         collapse_or_weighting_rule == "study_block_equalized") %>%
  slice_head(n = 1) %>%
  mutate(
    note = paste0("study-weighted\ntop: ", event_label(dominant_event_id))
  )

pE <- ggplot(weight_metric, aes(value, metric, fill = row_type)) +
  geom_line(aes(group = metric), colour = FIGURE_RULE_COLOUR, linewidth = 0.5) +
  geom_point(shape = 21, colour = FIGURE_INK, size = 2.8, stroke = 0.16) +
  geom_text(aes(label = percent(value, accuracy = 1), vjust = label_vjust, hjust = label_hjust), size = fig2_annot_size) +
  geom_label(
    data = study_weight_note,
    aes(x = 0.30, y = "Top-three events", label = note),
    inherit.aes = FALSE, size = 1.55, lineheight = 0.78,
    linewidth = 0.12, label.padding = unit(1.0, "pt"),
    fill = "white", colour = FIGURE_MUTED_TEXT
  ) +
  scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0.02, 0.06))) +
  scale_fill_manual(
    values = fig2_study_weight_colors,
    name = NULL
  ) +
  labs(x = "Share of resolved disrupted genomes", y = NULL) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = fig2_legend_text_size)
  )


# ===========================================================================
# Panel F: Within-origin structural concentration (phylogenetic independence)
# Shows that even within independent phylogenetic origins, the same event dominates
# ===========================================================================

within_origin_plot <- within_origin %>%
  filter(n_disrupted_descendant_tips >= 2) %>%
  mutate(
    dominant_event = event_label(dominant_prn_event_id),
    origin_label = paste0("Origin ", str_extract(origin_id, "\\d+"))
  ) %>%
  arrange(desc(n_disrupted_descendant_tips)) %>%
  slice_head(n = 10) %>%
  mutate(
    origin_label = factor(origin_label, levels = rev(unique(origin_label))),
    label_x = if_else(dominant_event_share > 0.9,
                      pmax(dominant_event_share - 0.04, 0),
                      pmin(dominant_event_share + 0.045, 1.12)),
    label_hjust = if_else(dominant_event_share > 0.9, 1, 0),
    event_label_text = paste0(n_unique_events, " event", if_else(n_unique_events == 1, "", "s"))
  )

pF <- ggplot(within_origin_plot, aes(dominant_event_share, origin_label)) +
  geom_segment(aes(x = 0, xend = dominant_event_share, yend = origin_label),
               colour = FIGURE_RULE_COLOUR, linewidth = 0.4) +
  geom_point(aes(size = n_disrupted_descendant_tips, fill = dominant_event),
             shape = 21, colour = FIGURE_INK, stroke = 0.18) +
  geom_text(aes(x = label_x, label = event_label_text, hjust = label_hjust),
            size = fig2_annot_size, colour = FIGURE_MUTED_TEXT) +
  scale_x_continuous(labels = percent, limits = c(0, 1.15), expand = c(0, 0)) +
  scale_size_area(max_size = 5, breaks = c(2, 5, 15, 38), name = "Tips") +
  scale_fill_manual(values = fig2_architecture_colors, na.value = FIGURE_GREY, guide = "none") +
  labs(x = "Dominant-event share within origin", y = NULL) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = fig2_legend_title_size),
    legend.text = element_text(size = fig2_legend_text_size)
  )

# ===========================================================================
# Panel G: Event temporal span + validation (lollipop with evidence markers)
# ===========================================================================

event_catalog <- junction_confidence %>%
  filter(sample_count > 0, !str_detect(prn_event_id, "insufficient")) %>%
  mutate(
    event = fig2_event_label(prn_event_id, mechanism_call),
    mechanism = mechanism_label(mechanism_call),
    count_label_x = pmin(year_max + 0.35, 2025.6),
    read_marker_x = 2029.0,
    count_label = paste0(comma(sample_count), "; ", country_count, " ctry")
  ) %>%
  arrange(desc(sample_count), event) %>%
  slice_head(n = fig2_top_event_limit) %>%
  mutate(event = factor(event, levels = rev(unique(event))))

pG <- ggplot(event_catalog, aes(y = event)) +
  geom_segment(aes(x = year_min, xend = year_max, yend = event),
               colour = FIGURE_RULE_COLOUR, linewidth = 0.8) +
  geom_point(aes(x = year_min), shape = "|", colour = FIGURE_MUTED_TEXT, size = 2) +
  geom_point(aes(x = year_max, size = sample_count, fill = event),
             shape = 21, colour = FIGURE_INK, stroke = 0.18) +
  geom_point(
    aes(x = read_marker_x, shape = confidence_label),
    size = 1.35, fill = "white", colour = FIGURE_INK, stroke = 0.32
  ) +
  geom_text(aes(x = count_label_x, label = count_label),
            hjust = -0.05, size = fig2_annot_size) +
  scale_x_continuous(limits = c(2008, 2030.0), breaks = seq(2010, 2025, 5),
                     expand = expansion(mult = c(0.02, 0.08))) +
  scale_size_area(max_size = 4.8, guide = "none") +
  scale_fill_manual(values = fig2_architecture_colors, na.value = FIGURE_GREY, guide = "none") +
  scale_shape_manual(values = confidence_shapes, name = "Junction tier", drop = TRUE) +
  guides(
    shape = guide_legend(nrow = 1, title.position = "left", override.aes = list(fill = "white")),
    colour = "none"
  ) +
  labs(x = "Observed calendar span; point size = genomes", y = NULL) +
  theme_nature(base_size = fig2_base_size) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = fig2_legend_title_size),
    legend.text = element_text(size = fig2_legend_text_size)
  )

# ===========================================================================
# Panel H: Mechanism-class composition (stacked proportional bar)
# ===========================================================================

mech_comp <- mechanism_summary %>%
  mutate(
    call_class = case_when(
      prn_mechanism_call == "intact" ~ "Intact",
      prn_mechanism_call == "insufficient_data" ~ "Insufficient",
      str_detect(prn_mechanism_call, "uncertain_fragmented") ~ "Uncertain",
      str_detect(prn_mechanism_call, "is481") ~ "IS481",
      str_detect(prn_mechanism_call, "inversion|rearrangement") ~ "Rearrangement",
      TRUE ~ "Other disrupted"
    ),
    call_class = factor(
      call_class,
      levels = c("Intact", "IS481", "Rearrangement", "Other disrupted", "Uncertain", "Insufficient")
    )
  ) %>%
  group_by(call_class) %>%
  summarise(n = sum(sample_count, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    frac = n / sum(n),
    label = if_else(frac > 0.04, paste0(comma(n), "\n", percent(frac, accuracy = 0.1)), "")
  )

pH <- ggplot(mech_comp, aes(x = "", y = n, fill = call_class)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = fig2_annot_size, lineheight = 0.82, colour = FIGURE_TEXT_COLOUR) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c(
    "Intact" = unname(npg_colors["blue"]),
    "IS481" = unname(npg_colors["red"]),
    "Rearrangement" = unname(npg_colors["teal"]),
    "Other disrupted" = unname(npg_colors["green"]),
    "Uncertain" = FIGURE_GREY,
    "Insufficient" = FIGURE_LIGHT_GREY
  ), name = NULL) +
  labs(x = NULL, y = "Genomes") +
  theme_nature(base_size = fig2_base_size) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = fig2_legend_text_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ===========================================================================
# Assembly and save
# ===========================================================================

fig2_layout <- "
AABB
CCDD
CCEE
"

fig2 <- pA + pG + pC + pD + pE +
  plot_layout(
    design = fig2_layout,
    widths = c(1, 1, 1, 1),
    heights = c(1.05, 0.5, 0.5)
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = fig2_tag_size, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.title.position = "plot",
    plot.margin = margin(3, 3, 3, 3)
  )

save_nc_pdf(fig2, "fig02_prn_structural_solution_space.pdf", height = 4)
save_nc_png(fig2, "fig02_prn_structural_solution_space.png", height = 4)
