#!/usr/bin/env Rscript
# Figure 1: public-genome atlas and recoverable-locus signal boundary

# function --------------------------------------------------------

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
library(maps)
library(patchwork)
library(readr)
library(scales)
library(stringr)
library(tidyr)

# setting -----------------------------------------------------------

# Unified typography for all panels (Nature standard: 7 pt base)
fig1_base_size <- FIGURE_BASE_SIZE
fig1_tag_size  <- FIGURE_TAG_SIZE        # bold panel labels (A-D)
fig1_legend_text_size <- FIGURE_LEGEND_TEXT_SIZE # legend item labels
fig1_legend_title_size <- FIGURE_LEGEND_TITLE_SIZE
fig1_repel_size <- 1.75    # country ISO3 labels in scatter panel (D)

map_bi_dim <- 4

make_bi_palette <- function(dim = 4) {
     low <- grDevices::col2rgb(FIGURE_PANEL_FILL)
     x_high <- grDevices::col2rgb(unname(npg_colors["teal"]))
     y_high <- grDevices::col2rgb(unname(npg_colors["red"]))
     both_high <- grDevices::col2rgb(unname(npg_colors["blue"]))
     grid <- expand.grid(x = seq_len(dim), y = seq_len(dim))
     cols <- apply(grid, 1, function(row) {
          tx <- if (dim == 1) 0 else (row[["x"]] - 1) / (dim - 1)
          ty <- if (dim == 1) 0 else (row[["y"]] - 1) / (dim - 1)
          rgb <- (1 - tx) * (1 - ty) * low +
               tx * (1 - ty) * x_high +
               (1 - tx) * ty * y_high +
               tx * ty * both_high
          grDevices::rgb(rgb[1], rgb[2], rgb[3], maxColorValue = 255)
     })
     stats::setNames(cols, paste(grid$x, grid$y, sep = "-"))
}

bi_palette <- make_bi_palette(map_bi_dim)

map_bi_x_labels <- c("<10", "<50", "<100", "100+")
map_bi_y_labels <- c("0%", "<5%", "<10%", "10%+")
map_ocean_fill <- FIGURE_MAP_FILL

# Keep the full map region names so subregions like Hong Kong can be matched explicitly.
world_region_lookup <- c(
     Argentina = "ARG",
     Australia = "AUS",
     Brazil = "BRA",
     Canada = "CAN",
     China = "CHN",
     `Czech Republic` = "CZE",
     Denmark = "DNK",
     Spain = "ESP",
     Finland = "FIN",
     France = "FRA",
     UK = "GBR",
     Guatemala = "GTM",
     Haiti = "HTI",
     India = "IND",
     Italy = "ITA",
     Japan = "JPN",
     Kenya = "KEN",
     Mexico = "MEX",
     Netherlands = "NLD",
     Norway = "NOR",
     `New Zealand` = "NZL",
     Poland = "POL",
     Russia = "RUS",
     Sweden = "SWE",
     Taiwan = "TWN",
     `South Africa` = "ZAF",
     USA = "USA"
)

fortify_world_map <- function(map_obj) {
     piece_id <- cumsum(is.na(map_obj$x) & is.na(map_obj$y)) + 1
     tibble::tibble(
          long = map_obj$x,
          lat = map_obj$y,
          group = piece_id,
          order = seq_along(map_obj$x),
          region_raw = map_obj$names[piece_id]
     ) %>%
          filter(!is.na(long), !is.na(lat)) %>%
          mutate(
               region_root = sub(":.*$", "", region_raw),
               country_iso3 = dplyr::recode(region_root, !!!world_region_lookup, .default = NA_character_),
               country_iso3 = if_else(startsWith(region_raw, "China:Hong Kong"), "HKG", country_iso3)
          )
}

prn_year <- safe_load(data_paths$prn_country_year, "Figure 1 PRN country-year summary") %>%
     mutate(
          across(c(year, n_genomes_total, n_genomes_prn_interpretable, n_prn_intact,
                   n_prn_disrupted, n_prn_uncertain_fragmented, n_prn_insufficient), as.numeric)
     )

focus_year_min <- 1980
prn_year_focus <- prn_year %>%
     filter(!is.na(year), year >= focus_year_min)

score <- safe_load(data_paths$rr_selection_scorecard, "selection scorecard") %>%
     mutate(
          n_retained_genomes = as.numeric(n_retained_genomes),
          n_prn_interpretable = as.numeric(n_prn_interpretable),
          interpretability_fraction = as.numeric(interpretability_fraction),
          selection_state = factor(
               selection_state,
               levels = c("primary_and_triangulated", "primary_only", "context_only")
          )
     )

stage_n <- function(id) {
     value <- flow %>% filter(stage_id == id) %>% pull(n_rows)
     if (length(value) == 0 || is.na(value[1])) NA_real_ else value[1]
}

country_xy <- tibble::tribble(
     ~country_iso3, ~lon, ~lat,
     "ARG", -64, -34, "AUS", 134, -25, "BRA", -52, -10, "CAN", -106, 56,
     "CHN", 104, 35, "CZE", 15, 49.8, "DNK", 10, 56, "ESP", -4, 40,
     "FIN", 26, 64, "FRA", 2, 46, "GBR", -2, 54, "IRL", -8, 53,
     "ITA", 12, 43, "JPN", 138, 37, "KEN", 37, 0, "MEX", -102, 23,
     "NLD", 5, 52, "NOR", 10, 62, "NZL", 174, -41, "SWE", 15, 62,
     "TWN", 121, 24, "USA", -98, 39, "ZAF", 24, -29
)

# panel A -------------------------------------------------------------

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
     ) %>%
     left_join(country_xy, by = "country_iso3")

country_bi_classes <- country_totals %>%
     mutate(
          # Use all countries with at least one interpretable genome so the filled map is not
          # restricted to the higher-support subset.
          map_retained = if_else(!is.na(frac_disrupted), n_retained, NA_real_),
          map_disrupted = if_else(!is.na(frac_disrupted), frac_disrupted, NA_real_)
     ) %>%
     filter(!is.na(map_retained), !is.na(map_disrupted)) %>%
     mutate(
          retained_bin = case_when(
               map_retained < 10 ~ 1L,
               map_retained < 50 ~ 2L,
               map_retained < 100 ~ 3L,
               TRUE ~ 4L
          ),
          disrupted_bin = case_when(
               map_disrupted == 0 ~ 1L,
               map_disrupted < 0.05 ~ 2L,
               map_disrupted < 0.10 ~ 3L,
               TRUE ~ 4L
          ),
          bi_class = paste(retained_bin, disrupted_bin, sep = "-")
     ) %>%
     select(country_iso3, bi_class)

world_map <- map("world", fill = TRUE, plot = FALSE)
world <- fortify_world_map(world_map)
world_fill <- world %>%
     left_join(country_bi_classes, by = "country_iso3") %>%
     filter(!is.na(bi_class))

map_bi_legend <- expand.grid(x = seq_len(map_bi_dim), y = seq_len(map_bi_dim)) %>%
     mutate(bi_class = paste(x, y, sep = "-")) %>%
     ggplot(aes(x, y, fill = bi_class)) +
     geom_tile() +
     scale_fill_manual(values = bi_palette, guide = "none") +
     scale_x_continuous(
          breaks = seq_len(map_bi_dim),
          labels = map_bi_x_labels,
          expand = c(0, 0)
     ) +
     scale_y_continuous(
          breaks = seq_len(map_bi_dim),
          labels = map_bi_y_labels,
          expand = c(0, 0)
     ) +
     coord_equal(expand = FALSE, clip = "off") +
     labs(x = "Retained\ngenomes", y = "Disrupted %") +
     theme_bw(base_size = 5) +
     theme(
          plot.background = element_rect(fill = map_ocean_fill, colour = NA),
          plot.margin = margin(0, 0, 0, 0),
          axis.title = element_text(size = 4.5),
          axis.text.x = element_text(size = 3.5, margin = margin(t = 1)),
          axis.text.y = element_text(size = 3.5, margin = margin(r = 1)),
          axis.ticks = element_blank()
     )

pA_map <- ggplot() +
     geom_polygon(
          data = world,
          aes(long, lat, group = group),
          fill = "white",
          colour = FIGURE_RULE_COLOUR,
          linewidth = 0.10
     ) +
     geom_polygon(
          data = world_fill,
          aes(long, lat, group = group, fill = bi_class),
          colour = "white",
          linewidth = 0.10,
          alpha = 0.96
     ) +
     # Rectangle indicating Europe inset extent
     annotate("rect", xmin = -12, xmax = 35, ymin = 35, ymax = 72,
              fill = NA, colour = FIGURE_MUTED_TEXT, linewidth = 0.25, linetype = "dashed") +
     scale_fill_manual(values = bi_palette, na.value = "white", guide = "none") +
     coord_quickmap(xlim = c(-170, 180), ylim = c(-55, 72), expand = FALSE) +
     labs(x = NULL, y = NULL) +
     theme_nature_map(base_size = fig1_base_size) +
     theme(
          legend.position = "none",
          panel.background = element_rect(fill = map_ocean_fill, colour = NA)
     )

# Europe inset
pA_europe <- ggplot() +
     geom_polygon(
          data = world,
          aes(long, lat, group = group),
          fill = "white",
          colour = FIGURE_RULE_COLOUR,
          linewidth = 0.12
     ) +
     geom_polygon(
          data = world_fill,
          aes(long, lat, group = group, fill = bi_class),
          colour = "white",
          linewidth = 0.12,
          alpha = 0.96
     ) +
     scale_fill_manual(values = bi_palette, na.value = "white", guide = "none") +
     coord_quickmap(xlim = c(-12, 35), ylim = c(35, 72), expand = FALSE) +
     labs(x = NULL, y = NULL) +
     theme_nature_map(base_size = fig1_base_size) +
     theme(
          legend.position = "none",
          panel.border = element_rect(colour = FIGURE_MUTED_TEXT, fill = NA, linewidth = 0.35),
          plot.background = element_rect(fill = map_ocean_fill, colour = NA),
          panel.background = element_rect(fill = map_ocean_fill, colour = NA),
          plot.margin = margin(0, 0, 0, 0)
     )

library(cowplot)

pA <- ggdraw(pA_map) +
     draw_plot(
          map_bi_legend,
          x = -0.08,
          y = 0.04,
          width  = 0.45,
          height = 0.45
     ) +
     draw_plot(
          pA_europe,
          x = 0.28,
          y = 0.2,
          width  = 0.37,
          height = 0.46
     )

# panel B -------------------------------------------------------------

country_levels <- country_totals %>%
     arrange(desc(n_retained)) %>%
     pull(country_iso3)

heatmap_years <- seq(focus_year_min, max(prn_year_focus$year, na.rm = TRUE), by = 1)
heatmap_dat <- prn_year_focus %>%
     filter(country_iso3 %in% country_levels, !is.na(year)) %>%
     mutate(
          country_iso3 = factor(country_iso3, levels = rev(country_levels)),
          year = as.integer(year)
     ) %>%
     complete(
          country_iso3,
          year = heatmap_years,
          fill = list(n_genomes_total = 0)
     ) %>%
     mutate(country_iso3 = factor(country_iso3, levels = rev(country_levels)))

pB <- ggplot(heatmap_dat, aes(year, country_iso3, fill = n_genomes_total)) +
     geom_tile(width = 0.9, height = 0.74) +
     scale_fill_gradientn(
          colours = blue_seq,
          trans = scales::log1p_trans(),
          breaks = c(0, 1, 10, 50, 100, 500),
          labels = comma,
          name = "Genomes"
     ) +
     scale_x_continuous(
          limits = c(focus_year_min - 0.45, max(heatmap_years) + 0.45),
          breaks = seq(1980, 2020, 10),
          expand = c(0, 0)
     ) +
     labs(x = "Collection year", y = NULL) +
     theme_nature(base_size = fig1_base_size) +
     theme(
          axis.text.y = element_text(size = 4.8),
          legend.position = "bottom",
          legend.title = element_text(size = fig1_legend_title_size),
          legend.text = element_text(size = fig1_legend_text_size)
     ) +
     guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5, barwidth = 5, barheight = 0.45))

# panel C -------------------------------------------------------------

annual <- prn_year_focus %>%
     filter(!is.na(year)) %>%
     group_by(year) %>%
     summarise(
          n_total = sum(n_genomes_total, na.rm = TRUE),
          n_interpretable = sum(n_genomes_prn_interpretable, na.rm = TRUE),
          n_disrupted = sum(n_prn_disrupted, na.rm = TRUE),
          frac = if_else(n_interpretable > 0, n_disrupted / n_interpretable, NA_real_),
          .groups = "drop"
     ) %>%
     filter(n_interpretable >= 5) %>%
     rowwise() %>%
     mutate(
          ci = list(binom.test(n_disrupted, n_interpretable)$conf.int),
          lo = ci[[1]],
          hi = ci[[2]]
     ) %>%
     ungroup()

pC <- ggplot(annual, aes(year, frac)) +
     geom_ribbon(aes(ymin = lo, ymax = hi), fill = alpha(npg_colors["red"], 0.13), colour = NA) +
     geom_line(colour = npg_colors["red"], linewidth = 0.42) +
     geom_point(aes(size = n_interpretable), shape = 21, fill = npg_colors["red"], colour = "white", stroke = 0.18) +
     # annotate("text", x = focus_year_min + 1, y = 0.05, label = "Pooled; not sampling-adjusted",
     #          hjust = 0, size = 1.7, colour = FIGURE_MUTED_TEXT, fontface = "italic") +
     scale_y_continuous(labels = percent,
                        limits = c(0, 1),
                        expand = expansion(mult = c(0.02, 0.05))) +
     scale_size_area(max_size = 4.8,
                     limits = c(5, 200),
                     breaks = c(5, 10, 50, 200), name = "Interpretable genomes") +
     labs(x = "Collection year", y = "Disrupted fraction") +
     theme_nature(base_size = fig1_base_size)

# panel D -------------------------------------------------------------

score_top_labels <- score %>%
     filter(country_iso3 %in% c("NLD", "BRA", "CZE")) %>%
     mutate(
          label_x = case_when(
               country_iso3 == "NLD" ~ 3.0,
               country_iso3 == "BRA" ~ 6.3,
               country_iso3 == "CZE" ~ 11.5,
               TRUE ~ n_prn_interpretable
          ),
          label_y = case_when(
               country_iso3 == "BRA" ~ 1.065,
               TRUE ~ 1.025
          )
     )

pF <- score %>%
     ggplot(aes(n_prn_interpretable, interpretability_fraction, label = country_iso3)) +
     geom_point(
          aes(fill = selection_state),
          shape = 21,
          size = 2.35,
          colour = FIGURE_INK,
          stroke = 0.18,
          alpha = 0.95
     ) +
     ggrepel::geom_text_repel(
          data = score %>% filter(!country_iso3 %in% c("NLD", "BRA", "CZE")),
          size = fig1_repel_size,
          min.segment.length = 0,
          segment.size = 0.1,
          segment.colour = FIGURE_MID_GREY,
          box.padding = 0.12,
          point.padding = 0.08,
          max.overlaps = Inf,
          colour = FIGURE_MUTED_TEXT
     ) +
     geom_text(
          data = score_top_labels,
          aes(x = label_x, y = label_y, label = country_iso3),
          inherit.aes = FALSE,
          size = fig1_repel_size,
          colour = FIGURE_MUTED_TEXT
     ) +
     scale_x_continuous(
          trans = "log10",
          breaks = c(1, 5, 20, 100, 1000),
          labels = comma,
          expand = expansion(mult = c(0.05, 0.08))
     ) +
     scale_y_continuous(
          labels = percent,
          limits = c(0, 1.08),
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          expand = expansion(mult = c(0.02, 0.05))
     ) +
     scale_fill_manual(
          values = c(
               primary_and_triangulated = unname(npg_colors["green"]),
               primary_only = unname(npg_colors["peach"]),
               context_only = FIGURE_GREY
          ),
          breaks = c("primary_and_triangulated", "primary_only", "context_only"),
          labels = c("Primary + triangulated", "Primary only", "Context only"),
          guide = guide_legend(nrow = 1, byrow = TRUE)
     ) +
     labs(x = "Interpretable genomes", y = "Interpretability fraction") +
     theme_nature(base_size = fig1_base_size) +
     theme(
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = fig1_legend_text_size)
     )

# save -------------------------------------------------------------

fig1_layout <- "
AAB
CDB
"

fig1 <- free(pA) + free(pB) + pC + pF +
     plot_layout(
          design = fig1_layout,
          widths = c(1, 1, 1),
          heights = c(1.1, 0.90)
     ) +
     plot_annotation(tag_levels = "a") &
     theme(
          plot.tag = element_text(face = "bold", size = fig1_tag_size, colour = FIGURE_TEXT_COLOUR),
          plot.tag.position = c(0, 1),
          plot.title.position = "plot",
          plot.margin = margin(3, 3, 3, 3)
     )

save_nc_pdf(fig1, "fig01_public_genome_atlas.pdf", height = 3.9)
save_nc_png(fig1, "fig01_public_genome_atlas.png", height = 3.9)
