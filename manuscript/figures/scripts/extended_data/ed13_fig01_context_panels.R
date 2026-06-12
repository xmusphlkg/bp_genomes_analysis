#!/usr/bin/env Rscript
# Extended Data Fig. 13: Figure 1 context panels (former Figure 1 panels D, E, and G)

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
library(ggrepel)
library(patchwork)
library(readr)
library(scales)

fig13_base_size <- FIGURE_BASE_SIZE
fig13_tag_size <- FIGURE_TAG_SIZE
fig13_legend_text_size <- FIGURE_LEGEND_TEXT_SIZE
fig13_legend_title_size <- FIGURE_LEGEND_TITLE_SIZE
fig13_repel_size <- 1.75

prn_year <- safe_load(data_paths$prn_country_year, "Figure 1 PRN country-year summary") %>%
     mutate(
          across(c(year, n_genomes_total, n_genomes_prn_interpretable, n_prn_intact,
                   n_prn_disrupted, n_prn_uncertain_fragmented, n_prn_insufficient), as.numeric)
     )

focus_year_min <- 1980
prn_year_focus <- prn_year %>%
     filter(!is.na(year), year >= focus_year_min)

flow <- safe_load(
     file.path(BASE_DIR, "manuscript", "supplementary", "Supplementary_Table_7_Cohort_Flow_and_Tree_Selection.tsv"),
     "cohort flow"
) %>%
     mutate(n_rows = as.numeric(n_rows))

stage_n <- function(id) {
     value <- flow %>% filter(stage_id == id) %>% pull(n_rows)
     if (length(value) == 0 || is.na(value[1])) NA_real_ else value[1]
}

country_totals <- prn_year_focus %>%
     group_by(country_iso3) %>%
     summarise(
          n_retained = sum(n_genomes_total, na.rm = TRUE),
          n_interpretable = sum(n_genomes_prn_interpretable, na.rm = TRUE),
          n_disrupted = sum(n_prn_disrupted, na.rm = TRUE),
          n_uncertain = sum(n_prn_uncertain_fragmented + n_prn_insufficient, na.rm = TRUE),
          frac_disrupted = if_else(n_interpretable > 0, n_disrupted / n_interpretable, NA_real_),
          frac_noninterpretable = if_else(n_retained > 0, pmax(n_retained - n_interpretable, 0) / n_retained, NA_real_),
          .groups = "drop"
     )

# panel D -------------------------------------------------------------

pD_dat <- country_totals %>%
     mutate(
          has_interpretable = n_interpretable > 0,
          frac_disrupted_plot = if_else(has_interpretable, frac_disrupted, 0)
     )

pD <- ggplot(pD_dat, aes(frac_noninterpretable, frac_disrupted_plot)) +
     geom_point(
          data = pD_dat %>% filter(has_interpretable),
          aes(size = n_retained, fill = n_interpretable),
          shape = 21,
          colour = FIGURE_INK,
          stroke = 0.18,
          alpha = 0.9
     ) +
     geom_point(
          data = pD_dat %>% filter(!has_interpretable),
          aes(size = n_retained),
          shape = 21,
          fill = FIGURE_LIGHT_GREY,
          colour = FIGURE_MID_GREY,
          stroke = 0.18,
          alpha = 0.95
     ) +
     ggrepel::geom_text_repel(
          data = pD_dat %>% filter(has_interpretable),
          aes(label = country_iso3),
          size = fig13_repel_size,
          min.segment.length = 0,
          segment.size = 0.1,
          segment.colour = FIGURE_MID_GREY,
          box.padding = 0.12,
          point.padding = 0.08,
          max.overlaps = 15,
          colour = FIGURE_MUTED_TEXT
     ) +
     scale_x_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
     scale_y_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
     scale_size_area(
          max_size = 4.8,
          breaks = c(1, 10, 50, 100, 500),
          labels = comma,
          name = "Retained\ngenomes"
     ) +
     scale_fill_gradientn(
          colours = blue_seq,
          trans = scales::log1p_trans(),
          breaks = c(0, 1, 10, 50, 100, 500),
          labels = comma,
          name = "Interpretable\ngenomes"
     ) +
     labs(x = "Non-interpretable fraction", y = "Disrupted fraction") +
     theme_nature(base_size = fig13_base_size) +
     theme(
          legend.position = "bottom",
          legend.box = "vertical",
          legend.title = element_text(size = fig13_legend_title_size),
          legend.text = element_text(size = fig13_legend_text_size)
     ) +
     guides(
          fill = guide_colorbar(title.position = "left", title.hjust = 0.5, barwidth = 5, barheight = 0.45, order = 1),
          size = guide_legend(order = 2, override.aes = list(fill = unname(npg_colors["teal"])))
     )

# panel E -------------------------------------------------------------

layer_dat <- tibble::tribble(
     ~layer, ~stage_id, ~group,
     "Combined manifest", "S01", "Manifest",
     "Interpretable prn", "S02", "Recoverable",
     "Disrupted prn", NA_character_, "Signal",
     "Structurally resolved", NA_character_, "Signal",
     "Validation set", "S05", "Audit",
     "Primary tree", "S11", "Audit"
) %>%
     mutate(
          n = case_when(
               layer == "Disrupted prn" ~ 617,
               layer == "Structurally resolved" ~ 577,
               TRUE ~ vapply(stage_id, stage_n, numeric(1))
          ),
          fraction = n / stage_n("S01"),
          layer = factor(layer, levels = layer)
     )

pE <- ggplot(layer_dat, aes(layer, fraction, group = 1)) +
     geom_line(colour = FIGURE_RULE_COLOUR, linewidth = 0.55) +
     geom_point(aes(fill = group, size = n), shape = 21, colour = FIGURE_INK, stroke = 0.22) +
     geom_text(aes(label = paste0(comma(n), "\n", percent(fraction, accuracy = 1))), vjust = -0.8, size = 1.9, lineheight = 0.84) +
     scale_y_continuous(labels = percent, breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
     coord_cartesian(ylim = c(0, 1.28), clip = "off") +
     scale_size_area(max_size = 5.6, guide = "none") +
     scale_fill_manual(values = c("Manifest" = FIGURE_DARK_GREY, "Recoverable" = unname(npg_colors["blue"]), "Signal" = unname(npg_colors["red"]), "Audit" = unname(npg_colors["green"])), guide = "none") +
     labs(x = NULL, y = "Fraction of manifest") +
     theme_nature(base_size = fig13_base_size) +
     theme(
          axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5),
          plot.margin = margin(10, 3, 3, 3)
     )

# panel G: Non-interpretable composition ---------------------------------

noninterp_dat <- prn_year_focus %>%
     group_by(country_iso3) %>%
     summarise(
          Interpretable = sum(n_genomes_prn_interpretable, na.rm = TRUE),
          Fragmented = sum(n_prn_uncertain_fragmented, na.rm = TRUE),
          Insufficient = sum(n_prn_insufficient, na.rm = TRUE),
          .groups = "drop"
     ) %>%
     mutate(
          Total = Interpretable + Fragmented + Insufficient,
          Other = pmax(0, Total - Interpretable - Fragmented - Insufficient)
     ) %>%
     filter(Total > 0) %>%
     select(country_iso3, Interpretable, Fragmented, Insufficient) %>%
     pivot_longer(cols = -country_iso3, names_to = "category", values_to = "n") %>%
     group_by(category) %>%
     summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
     mutate(
          category = factor(category, levels = c("Interpretable", "Fragmented", "Insufficient")),
          frac = n / sum(n),
          label = paste0(comma(n), "\n(", percent(frac, accuracy = 0.1), ")")
     )

pG <- ggplot(noninterp_dat, aes(x = category, y = n, fill = category)) +
     geom_col(width = 0.65, colour = FIGURE_INK, linewidth = 0.18) +
     geom_text(aes(label = label), vjust = -0.3, size = 1.9, lineheight = 0.84) +
     scale_y_continuous(
          labels = comma,
          expand = expansion(mult = c(0, 0.18))
     ) +
     scale_fill_manual(
          values = c(
               "Interpretable" = unname(npg_colors["blue"]),
               "Fragmented" = unname(npg_colors["peach"]),
               "Insufficient" = unname(npg_colors["grey"])
          ),
          guide = "none"
     ) +
     labs(x = NULL, y = "Genomes") +
     theme_nature(base_size = fig13_base_size) +
     theme(
          axis.text.x = element_text(size = 5.5),
          plot.margin = margin(3, 3, 3, 3)
     )

# save -------------------------------------------------------------

ed13 <- free(pD) + pE + pG +
     plot_layout(
          ncol = 3,
          widths = c(1.10, 0.95, 0.95),
          guides = "keep"
     ) +
     plot_annotation(tag_levels = "A") &
     theme(
          plot.tag = element_text(face = "bold", size = fig13_tag_size, colour = FIGURE_TEXT_COLOUR),
          plot.tag.position = c(0, 1),
          plot.title.position = "plot",
          plot.margin = margin(3, 3, 3, 3)
     )

save_ed_pdf(ed13, "Extended_Data_Fig_13_Figure_1_Context_Panels.pdf", height = NC_MAX_HEIGHT / 2)
save_ed_png(ed13, "Extended_Data_Fig_13_Figure_1_Context_Panels.png", height = NC_MAX_HEIGHT / 2)
