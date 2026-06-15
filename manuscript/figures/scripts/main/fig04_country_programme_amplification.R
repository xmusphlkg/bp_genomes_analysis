#!/usr/bin/env Rscript
# Figure 4: archive context and sampling heterogeneity
# Restructured: unified typography, no panel subtitles, Nature Portfolio compliant

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

fig4_base_size <- FIGURE_BASE_SIZE
fig4_tag_size  <- FIGURE_TAG_SIZE

focus_countries <- c("USA", "NZL", "JPN", "AUS", "GBR", "CHN", "FRA", "BRA")
focus_year_min <- 1980

fig4_event_colors <- c(
  "Tree package" = figure_discrete_at(1),
  "Archive detection" = figure_discrete_at(4)
)

# data loading -------------------------------------------------------

history <- safe_load(data_paths$rr_program_history_manifest, "program history") %>%
  mutate(
    start_year = as.numeric(start_year),
    end_year = as.numeric(end_year),
    epoch_class = short_epoch_type(epoch_type, prn_in_formulation)
  ) %>%
  filter(country_iso3 %in% focus_countries)

prevalence <- safe_load(data_paths$rr_epoch_prevalence, "epoch prevalence") %>%
  mutate(
    start_year = as.numeric(start_year),
    end_year = as.numeric(end_year),
    n_prn_interpretable = as.numeric(n_prn_interpretable),
    naive_prevalence = as.numeric(naive_prevalence),
    ipw_prevalence = as.numeric(ipw_prevalence),
    epoch_class = short_epoch_type(epoch_type, prn_in_formulation),
    midpoint = (start_year + end_year) / 2
  ) %>%
  filter(country_iso3 %in% focus_countries)

bounds <- safe_load(data_paths$rr_epoch_bounds, "epoch missingness bounds") %>%
  mutate(
    bound_lower_prevalence = as.numeric(bound_lower_prevalence),
    bound_upper_prevalence = as.numeric(bound_upper_prevalence),
    bound_width = as.numeric(bound_width)
  )

relative <- safe_load(data_paths$rr_relative_year_plot_data, "relative-year plot data") %>%
  mutate(
    relative_year = as.numeric(relative_year),
    event_year = as.numeric(event_year),
    ipw_prevalence = as.numeric(ipw_prevalence),
    n_genomes_prn_interpretable = as.numeric(n_genomes_prn_interpretable),
    n_origin_clades_active = as.numeric(n_origin_clades_active),
    event_type = recode(event_type, first_local_origin = "Tree package", first_prn_detection = "Archive detection")
  )

pooled_relative <- safe_load(file.path(FIGURE_DATA_DIR, "figure4_event_centered_pooled.tsv"), "pooled event-centered data") %>%
  mutate(
    relative_year = as.numeric(relative_year),
    mean_value = as.numeric(mean_value),
    median_value = as.numeric(median_value),
    n_countries = as.numeric(n_countries),
    event_type = recode(event_type, first_local_origin = "Tree package", first_prn_detection = "Archive detection")
  )

prn_year <- safe_load(data_paths$prn_country_year, "country-year PRN summary") %>%
  mutate(
    year = as.numeric(year),
    n_genomes_total = as.numeric(n_genomes_total),
    n_genomes_prn_interpretable = as.numeric(n_genomes_prn_interpretable),
    n_prn_disrupted = as.numeric(n_prn_disrupted),
    frac_prn_disrupted = as.numeric(frac_prn_disrupted),
    frac_prn_disrupted = if_else(
      is.na(frac_prn_disrupted) & n_genomes_prn_interpretable > 0,
      n_prn_disrupted / n_genomes_prn_interpretable,
      frac_prn_disrupted
    )
  ) %>%
  filter(country_iso3 %in% focus_countries, between(year, focus_year_min, 2026))

turnover <- safe_load(file.path(FIGURE_DATA_DIR, "selected_country", "country_epoch_architecture_turnover_summary.tsv"), "architecture turnover summary") %>%
  mutate(
    next_dominant_share = as.numeric(next_dominant_share),
    architecture_total_variation_distance = as.numeric(architecture_total_variation_distance),
    next_event = event_label(next_dominant_event_id)
  )

# Panel A: Country programme timelines (aligned tracks) -------------------

event_years <- relative %>%
  filter(country_iso3 %in% focus_countries, relative_year == 0, !is.na(event_year), event_year >= focus_year_min) %>%
  distinct(country_iso3, event_type, event_year)

country_track_levels <- rev(focus_countries)
country_track_lookup <- tibble::tibble(
  country_iso3 = country_track_levels,
  country_y = seq_along(country_track_levels)
)

rail_history <- history %>%
  left_join(country_track_lookup, by = "country_iso3") %>%
  mutate(
    rail_y = country_y + 0.20,
    start_plot = pmax(start_year - 0.50, focus_year_min - 0.50),
    end_plot = pmin(end_year + 0.52, 2026.50)
  )

density_max <- max(log10(prn_year$n_genomes_total + 1), na.rm = TRUE)
if (!is.finite(density_max) || density_max <= 0) density_max <- 1
genome_density_breaks <- c(1, 10, 50, 100, 500)
genome_density_breaks <- genome_density_breaks[log10(genome_density_breaks + 1) <= density_max]
genome_density_scaled_breaks <- log10(genome_density_breaks + 1) / density_max

timeline_dat <- prn_year %>%
  left_join(country_track_lookup, by = "country_iso3") %>%
  mutate(
    density_scaled = log10(n_genomes_total + 1) / density_max,
    frac_prn_disrupted = pmin(pmax(frac_prn_disrupted, 0), 1),
    frac_prn_disrupted = if_else(n_genomes_prn_interpretable > 0, frac_prn_disrupted, NA_real_),
    signal_y = country_y - 0.20
  )

timeline_density_dat <- timeline_dat %>%
  filter(n_genomes_total > 0, !is.na(country_y))

timeline_point_dat <- timeline_dat %>%
  filter(n_genomes_prn_interpretable > 0, !is.na(frac_prn_disrupted), !is.na(country_y))

event_years <- event_years %>%
  left_join(country_track_lookup, by = "country_iso3") %>%
  mutate(event_y = country_y - 0.20 + if_else(event_type == "Tree package", 0.16, -0.16))

programme_backbone <- tibble::tibble(
  country_iso3 = country_track_levels,
  country_y = seq_along(country_track_levels),
  rail_y = country_y + 0.20,
  start_plot = focus_year_min - 0.50,
  end_plot = 2026.50
)

pA <- ggplot() +
     geom_tile(data = timeline_density_dat,
               aes(x = year, y = signal_y, alpha = density_scaled),
               fill = unname(npg_colors["teal"]), width = 1.00, height = 0.34, colour = NA) +
     geom_segment(data = programme_backbone,
                  aes(x = start_plot, xend = end_plot, y = rail_y, yend = rail_y),
                  colour = FIGURE_PANEL_FILL, linewidth = 2.7, lineend = "butt") +
     geom_segment(data = rail_history,
                  aes(x = start_plot, xend = end_plot, y = rail_y, yend = rail_y, colour = epoch_class),
                  linewidth = 2.1, lineend = "butt", alpha = 0.72) +
     scale_colour_manual(values = programme_epoch_colors,
                         drop = TRUE,
                         name = "Archive context\nannotation",
                         guide = guide_legend(ncol = 1, override.aes = list(linewidth = 3.0))) +
     geom_point(data = timeline_point_dat,
                aes(x = year, y = signal_y, size = n_genomes_prn_interpretable, fill = frac_prn_disrupted),
                shape = 21, colour = FIGURE_INK, stroke = 0.16, alpha = 0.96) +
     geom_point(data = event_years,
                aes(x = event_year, y = event_y, shape = event_type),
                colour = FIGURE_INK, fill = "white", size = 1.75, stroke = 0.48) +
     scale_fill_gradientn(colours = red_seq,
                          limits = c(0, 1),
                          breaks = c(0, 0.5, 1),
                          labels = percent_format(accuracy = 1),
                          name = "Observed\ndisrupted\nfraction") +
     scale_alpha_continuous(range = c(0.10, 0.48),
                            breaks = genome_density_scaled_breaks,
                            labels = comma(genome_density_breaks),
                            name = "Genomes\nper year",
                            guide = guide_legend(override.aes = list(fill = unname(npg_colors["teal"]), colour = NA, size = 3))) +
     scale_size_area(max_size = 3.4, guide = "none") +
     scale_shape_manual(values = c("Tree package" = 24, "Archive detection" = 21),
                        name = "Marker",
                        guide = guide_legend(ncol = 1)) +
     scale_x_continuous(breaks = seq(1980, 2020, 5), expand = c(0, 0)) +
     scale_y_continuous(
          breaks = seq_along(country_track_levels),
          labels = country_track_levels,
          expand = expansion(add = c(0.48, 0.48))
     ) +
     labs(x = "Calendar year", y = NULL) +
     coord_cartesian(xlim = c(focus_year_min - 0.50, 2026.50), clip = "on") +
     theme_nature(base_size = fig4_base_size) +
     theme(
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.byrow = F,
          axis.text.y = element_text(face = "bold", size = 5.4),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
     )

# ===========================================================================
# Panel B: Epoch IPW prevalence with missingness bounds
# ===========================================================================

country_order_b <- c("USA", "NZL", "JPN", "AUS", "GBR", "CHN")
country_axis_b <- tibble::tibble(
  country_iso3 = country_order_b,
  country_y = length(country_order_b) - seq_along(country_order_b) + 1
)

epoch_bounds_dat <- prevalence %>%
  left_join(
    bounds %>% select(country_iso3, epoch_id, bound_lower_prevalence, bound_upper_prevalence),
    by = c("country_iso3", "epoch_id")
  ) %>%
  filter(country_iso3 %in% country_order_b, !is.na(ipw_prevalence)) %>%
  mutate(
    country_y = length(country_order_b) - match(country_iso3, country_order_b) + 1,
    epoch_midpoint = (start_year + end_year) / 2
  ) %>%
  arrange(country_iso3, epoch_midpoint) %>%
  group_by(country_iso3) %>%
  mutate(
    epoch_rank = row_number(),
    epoch_n = n(),
    epoch_offset = if_else(epoch_n == 1, 0, (epoch_rank - (epoch_n + 1) / 2) * 0.18),
    country_y_plot = country_y + epoch_offset
  ) %>%
  ungroup()

pB <- ggplot(epoch_bounds_dat) +
  geom_segment(data = country_axis_b,
    aes(x = 0, xend = 1, y = country_y, yend = country_y),
    inherit.aes = FALSE, colour = FIGURE_PANEL_FILL, linewidth = 0.25) +
  geom_segment(aes(x = bound_lower_prevalence, xend = bound_upper_prevalence,
                   y = country_y_plot, yend = country_y_plot, colour = epoch_class),
    linewidth = 0.46, alpha = 0.88, lineend = "round", na.rm = TRUE) +
  geom_point(aes(x = ipw_prevalence, y = country_y_plot, size = n_prn_interpretable, fill = epoch_class),
    shape = 21, colour = FIGURE_INK, stroke = 0.18, alpha = 0.96, na.rm = TRUE) +
  annotate("text", x = 0.03, y = 0.72, label = "bounds",
    hjust = 0, size = 1.55, colour = FIGURE_MUTED_TEXT, fontface = "italic") +
  scale_x_continuous(labels = percent, limits = c(0, 1), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(breaks = seq_along(country_order_b), labels = rev(country_order_b),
                     expand = expansion(add = c(0.38, 0.38))) +
  scale_size_area(max_size = 3.8, guide = "none") +
  scale_colour_manual(values = programme_epoch_colors, guide = "none") +
  scale_fill_manual(values = programme_epoch_colors, guide = "none") +
  labs(
    x = "Archive-frame disrupted fraction sensitivity estimate",
    y = NULL
  ) +
  theme_nature(base_size = fig4_base_size) +
  theme(
    axis.text.y = element_text(face = "bold", size = 5.5),
    legend.position = "none"
  )

# ===========================================================================
# Panel C: Event-centered amplification (pooled + country traces)
# ===========================================================================

rel_dat <- relative %>%
  filter(country_iso3 %in% c("USA", "NZL", "JPN", "AUS", "CHN", "FRA"), !is.na(ipw_prevalence)) %>%
  mutate(country_iso3 = factor(country_iso3, levels = c("USA", "NZL", "JPN", "AUS", "CHN", "FRA")))

pooled_dat <- pooled_relative %>%
  filter(metric_name == "ipw_prevalence", !is.na(median_value))

pooled_label_dat <- pooled_dat %>%
  group_by(event_type) %>%
  filter(relative_year == max(relative_year, na.rm = TRUE)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
	  mutate(
	    label_x = 3.33,
	    label_hjust = 0,
	    label_y = pmax(median_value, 0.055),
	    event_label = recode(event_type, `Tree package` = "Tree pkg", `Archive detection` = "Detection")
	  )

pC <- ggplot() +
  annotate("rect", xmin = -0.45, xmax = 0.45, ymin = -Inf, ymax = Inf,
           fill = FIGURE_PANEL_FILL, colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.26, colour = FIGURE_MID_GREY) +
  geom_line(data = rel_dat,
    aes(relative_year, ipw_prevalence, colour = event_type, group = interaction(country_iso3, event_type)),
    linewidth = 0.32, alpha = 0.20, na.rm = TRUE) +
  geom_line(data = pooled_dat,
    aes(relative_year, median_value, colour = event_type),
    linewidth = 0.84, na.rm = TRUE) +
  geom_point(data = pooled_dat,
    aes(relative_year, median_value, fill = event_type),
    shape = 21, size = 2.15, colour = FIGURE_INK, stroke = 0.18, alpha = 0.95, na.rm = TRUE) +
	  geom_text(data = pooled_label_dat,
	    aes(label_x, label_y, label = event_label, colour = event_type, hjust = label_hjust),
	    size = 1.65, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = fig4_event_colors, guide = "none") +
  scale_fill_manual(values = fig4_event_colors, guide = "none") +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0.02, 0.08))) +
  scale_x_continuous(breaks = -3:3) +
  labs(
    x = "Years from aligned marker",
    y = "Median archive-frame disrupted fraction"
  ) +
  coord_cartesian(xlim = c(-3, 3.62), clip = "off") +
  theme_nature(base_size = fig4_base_size) +
  theme(
    legend.position = "none",
    plot.margin = margin(3, 8, 3, 3)
  )

# ===========================================================================
# Panel D: Architecture turnover between epochs
# ===========================================================================

country_order_d <- c("USA", "NZL", "JPN", "AUS")
country_axis_d <- tibble::tibble(
  country_iso3 = country_order_d,
  country_y = length(country_order_d) - seq_along(country_order_d) + 1
)

turnover_dat <- turnover %>%
  filter(comparison_type == "primary_selected_country_pair", country_iso3 %in% c("USA", "NZL", "JPN", "AUS")) %>%
  mutate(
    country_y = length(country_order_d) - match(country_iso3, country_order_d) + 1,
    next_event_plot = recode(next_event, fragmented_contigs = "Fragmented contigs"),
    next_event_plot = factor(next_event_plot, levels = c("Fragmented contigs", "IS481 gap1043")),
    next_event_label = recode(as.character(next_event_plot),
      `IS481 gap1043` = "IS481 1043",
      `Fragmented contigs` = "fragmented contigs",
      .default = as.character(next_event_plot)
    ),
    direct_label = paste0(next_event_label, "\n", percent(next_dominant_share, accuracy = 1)),
    label_x = if_else(architecture_total_variation_distance > 0.82,
                      architecture_total_variation_distance - 0.060,
                      architecture_total_variation_distance + 0.060),
    label_hjust = if_else(architecture_total_variation_distance > 0.82, 1, 0)
  )

fig4_turnover_colors <- c(
  architecture_colors,
  "Fragmented contigs" = figure_discrete_at(6)
)

pD <- ggplot(turnover_dat, aes(architecture_total_variation_distance, country_y)) +
  geom_vline(xintercept = c(0, 0.5, 1), linewidth = 0.18, colour = FIGURE_PANEL_FILL) +
  geom_segment(data = country_axis_d,
    aes(x = 0, xend = 1, y = country_y, yend = country_y),
    inherit.aes = FALSE, colour = FIGURE_PANEL_FILL, linewidth = 0.25) +
  geom_segment(aes(x = 0, xend = architecture_total_variation_distance, yend = country_y),
    colour = FIGURE_RULE_COLOUR, linewidth = 0.46) +
  geom_point(aes(size = next_dominant_share, fill = next_event_plot),
    shape = 21, colour = FIGURE_INK, stroke = 0.18) +
  geom_text(aes(x = label_x, label = direct_label, hjust = label_hjust),
            size = 1.42, lineheight = 0.82, colour = FIGURE_MUTED_TEXT) +
  annotate("text", x = 0.03, y = 0.55, label = "same", hjust = 0, size = 1.55, colour = FIGURE_MUTED_TEXT) +
  annotate("text", x = 0.97, y = 0.55, label = "replaced", hjust = 1, size = 1.55, colour = FIGURE_MUTED_TEXT) +
  scale_x_continuous(labels = percent, limits = c(0, 1), breaks = c(0, 0.5, 1), expand = expansion(mult = c(0.01, 0.06))) +
  scale_y_continuous(breaks = seq_along(country_order_d), labels = rev(country_order_d),
                     expand = expansion(add = c(0.55, 0.38))) +
  scale_size_area(max_size = 5.2, labels = percent, name = "Dominant share", guide = "none") +
  scale_fill_manual(values = fig4_turnover_colors, na.value = FIGURE_GREY, guide = "none") +
  labs(
    x = "Architecture turnover (0 = same, 100% = replaced)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_nature(base_size = fig4_base_size) +
  theme(
    axis.text.y = element_text(face = "bold", size = 5.5),
    legend.position = "none"
  )

# ===========================================================================
# Assembly and save
# ===========================================================================

fig4_layout <- "
AAA
BCD
"

fig4 <- pA + pB + pC + free(pD) +
  plot_layout(design = fig4_layout, heights = c(1.15, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = fig4_tag_size, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1),
    plot.margin = margin(3, 3, 3, 3)
  )

save_nc_pdf(fig4, "fig04_archive_context_amplification.pdf", height = NC_MAX_HEIGHT)
save_nc_png(fig4, "fig04_archive_context_amplification.png", height = NC_MAX_HEIGHT)
