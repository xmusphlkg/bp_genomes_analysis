#!/usr/bin/env Rscript
# Extended Data Fig. 12: country-year interpretability and study-block audit

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

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)

audit <- safe_load(
  file.path(BASE_DIR, "manuscript", "figure_data", "epidemiology_revision_country_year_audit.tsv"),
  "country-year interpretability and block audit"
) %>%
  mutate(
    year = as.integer(year),
    n_total_records = as.numeric(n_total_records),
    interpretability_fraction = as.numeric(interpretability_fraction),
    dominant_base_block_share_total = as.numeric(dominant_base_block_share_total),
    dominant_base_block_share_interpretable = as.numeric(dominant_base_block_share_interpretable),
    tile_alpha = pmin(1, 0.28 + log10(pmax(n_total_records, 1)) / 1.5)
  )

country_levels <- audit %>%
  group_by(country_iso3, country_name) %>%
  summarise(total_records = sum(n_total_records, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_records), country_iso3) %>%
  pull(country_iso3)

audit <- audit %>%
  mutate(country_iso3 = factor(country_iso3, levels = rev(country_levels)))

pA <- ggplot(audit, aes(year, country_iso3)) +
  geom_tile(aes(fill = interpretability_fraction, alpha = tile_alpha), colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = teal_seq,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = percent,
    name = "Interpretability",
    guide = guide_colorbar(title.position = "top", barwidth = grid::unit(48, "pt"), barheight = grid::unit(3, "pt"))
  ) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = c(1940, 1960, 1980, 2000, 2020)) +
  labs(x = NULL, y = NULL) +
  theme_nature_matrix() +
  theme(
    axis.text.y = element_text(size = 6.4, face = "bold")
  )

pB <- ggplot(audit, aes(year, country_iso3)) +
  geom_tile(aes(fill = dominant_base_block_share_total, alpha = tile_alpha), colour = "white", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = orange_seq,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = percent,
    name = "Dominant study block share",
    guide = guide_colorbar(title.position = "top", barwidth = grid::unit(48, "pt"), barheight = grid::unit(3, "pt"))
  ) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = c(1940, 1960, 1980, 2000, 2020)) +
  labs(x = "Collection year", y = NULL) +
  theme_nature_matrix() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

fig <- pA + pB +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = FIGURE_TAG_SIZE, colour = FIGURE_TEXT_COLOUR),
    plot.tag.position = c(0, 1)
  )

save_ed_pdf(fig, "Extended_Data_Fig_12_Country_Year_Missingness_Audit.pdf", height = 6.0)
save_ed_png(fig, "Extended_Data_Fig_12_Country_Year_Missingness_Audit.png", height = 6.0)
