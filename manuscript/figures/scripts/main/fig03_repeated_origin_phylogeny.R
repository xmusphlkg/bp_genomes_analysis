#!/usr/bin/env Rscript
# Figure 3: repeated-acquisition phylogeny reconstructed from frozen tree tables
# Layout: AAB / AAC / AAD (Nature double-column, 3-row composite)
# A = optimized fan phylogeny with annotation rings
# B = ASR Fitch package count by scenario (range plot)
# C = Stochastic mapping package count (range plot)
# D = Event-specific minimum tree-level package burden for recurrent architectures

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

library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(patchwork)
library(readr)
library(scales)
library(stringr)
library(tidyr)

# setting -----------------------------------------------------------

fig3_base_size <- FIGURE_BASE_SIZE
fig3_tag_size  <- FIGURE_TAG_SIZE

scenario_label <- function(x) {
     dplyr::case_when(
          x %in% c("primary") ~ "Primary",
          x %in% c("rooting", "rooting sensitivity") ~ "Rooting sensitivity",
          x %in% c("support threshold") ~ "Support threshold",
          x %in% c("composition filtered") ~ "Composition-filtered",
          x %in% c("resampling study block balanced") ~ "Study-block balanced",
          x %in% c("resampling country balanced") ~ "Country-balanced",
          x %in% c("resampling time balanced") ~ "Time-balanced stress",
          TRUE ~ stringr::str_to_sentence(x)
     )
}

scenario_levels <- c(
     "Time-balanced stress",
     "Country-balanced",
     "Study-block balanced",
     "Composition-filtered",
     "Support threshold",
     "Rooting sensitivity",
     "Primary"
)

# data loading --------------------------------------------------------

tip_calls <- safe_load(file.path(FIGURE_DATA_DIR, "fig02_prn_mechanism_calls.tsv"), "tip-level PRN calls") %>%
     select(
          assembly_accession,
          call_prn_mechanism = prn_mechanism_call,
          call_prn_event = prn_event_id,
          prn_call_confidence,
          read_validation_status,
          read_validation_support
     )

nodes <- safe_load(data_paths$fig3_tree_nodes, "workflow tree nodes") %>%
     left_join(tip_calls, by = c("tip_label" = "assembly_accession")) %>%
     mutate(
          x_cladogram = as.numeric(x_cladogram),
          x_branch_length = as.numeric(x_branch_length),
          is_fitch_origin = as_logical_flag(is_fitch_origin),
          is_reference = as_logical_flag(is_reference),
          year = as.numeric(year),
          call_prn_mechanism = coalesce(call_prn_mechanism, prn_mechanism_call),
          tip_state = state_label(if_else(is_reference, "Reference", coalesce(observed_prn_state, call_prn_mechanism))),
          mechanism_class = mechanism_label(call_prn_mechanism),
          architecture = if_else(tip_state == "Disrupted", event_label(call_prn_event), NA_character_),
          label = if_else(node_type == "tip", tip_label, node_id),
          country_focus = if_else(country_iso3 %in% names(country_role_colors), country_iso3, "Other"),
          year_bin = case_when(
               is_reference ~ "Reference",
               is.na(year) ~ "Year unresolved",
               year < 1990 ~ "<1990",
               year < 2000 ~ "1990s",
               year < 2010 ~ "2000s",
               year < 2020 ~ "2010s",
               TRUE ~ "2020+"
          ),
          read_track = case_when(
               str_detect(read_validation_support, "read|clip|is") & !read_validation_support %in% c("not_evaluated", "") ~ "Read-supported",
               read_validation_status %in% c("not_run", "", NA_character_) ~ "Not run",
               TRUE ~ "Read-linked unresolved"
          ),
          confidence_track = case_when(
               prn_call_confidence == "rule_high" ~ "High-confidence call",
               prn_call_confidence == "rule_medium" ~ "Medium-confidence call",
               prn_call_confidence == "rule_low" ~ "Low/fragmented call",
               TRUE ~ "No call confidence"
          )
     )

asr <- safe_load(file.path(FIGURE_DATA_DIR, "asr_scenario_registry.tsv"), "ASR scenario registry") %>%
     mutate(
          fitch_origin_events = as.numeric(fitch_origin_events),
          largest_disrupted_clade_share = as.numeric(largest_disrupted_clade_share),
          disrupted_tip_count = as.numeric(disrupted_tip_count),
          scenario_class = scenario_label(str_replace_all(scenario_class, "_", " "))
     )

simmap <- safe_load(file.path(FIGURE_DATA_DIR, "asr_stochastic_mapping_summary.tsv"), "stochastic mapping summary") %>%
     mutate(
          stochastic_origin_count_median = as.numeric(stochastic_origin_count_median),
          stochastic_origin_count_lower_95 = as.numeric(stochastic_origin_count_lower_95),
          stochastic_origin_count_upper_95 = as.numeric(stochastic_origin_count_upper_95),
          scenario_class = scenario_label(str_replace_all(scenario_class, "_", " "))
     )

event_specific <- safe_load(file.path(FIGURE_DATA_DIR, "event_specific_acquisition_summary.tsv"), "event-specific acquisition summary") %>%
     mutate(
          sample_count = as.numeric(sample_count),
          n_country_year_cells = as.numeric(n_country_year_cells),
          n_mlst_st = as.numeric(n_mlst_st),
          acquisition_package_count = as.numeric(acquisition_package_count),
          non_singleton_package_count = as.numeric(non_singleton_package_count),
          largest_package_disrupted_tips = as.numeric(largest_package_disrupted_tips),
          rank_by_genome_burden = as.numeric(rank_by_genome_burden)
     )

within_origin <- safe_load(file.path(FIGURE_DATA_DIR, "within_origin_structural_concentration.tsv"), "within-origin structural concentration") %>%
     mutate(
          n_disrupted_descendant_tips = as.numeric(n_disrupted_descendant_tips),
          n_unique_events = as.numeric(n_unique_events),
          dominant_event_count = as.numeric(dominant_event_count),
          dominant_event_share = as.numeric(dominant_event_share),
          top3_event_share = as.numeric(top3_event_share),
          country_count = as.numeric(country_count),
          mlst_st_count = as.numeric(mlst_st_count),
          year_min = as.numeric(year_min),
          year_max = as.numeric(year_max)
     )

# panel A: optimized fan phylogeny with annotation rings ----------------

build_phylo <- function(node_table) {
     tips <- node_table %>% filter(node_type == "tip") %>% arrange(y_order, label)
     internals <- node_table %>% filter(node_type != "tip") %>% arrange(node_id)
     node_map <- c(
          setNames(seq_len(nrow(tips)), tips$node_id),
          setNames(nrow(tips) + seq_len(nrow(internals)), internals$node_id)
     )
     edge_df <- node_table %>%
          filter(!is.na(parent_node_id), nzchar(parent_node_id), node_id %in% names(node_map), parent_node_id %in% names(node_map)) %>%
          left_join(
               node_table %>% select(parent_node_id_lookup = node_id, parent_x_cladogram = x_cladogram),
               by = c("parent_node_id" = "parent_node_id_lookup")
          ) %>%
          transmute(
               parent = unname(node_map[parent_node_id]),
               child = unname(node_map[node_id]),
               edge_length = pmax(replace_na(x_cladogram - parent_x_cladogram, 0), 1e-8)
          )
     phy <- list(
          edge = as.matrix(edge_df[, c("parent", "child")]),
          edge.length = edge_df$edge_length,
          tip.label = tips$tip_label,
          Nnode = nrow(internals),
          node.label = internals$node_id
     )
     class(phy) <- "phylo"
     phy <- ape::reorder.phylo(phy)
     tip_depths <- ape::node.depth.edgelength(phy)[seq_along(phy$tip.label)]
     max_tip_depth <- max(tip_depths, na.rm = TRUE)
     if (is.finite(max_tip_depth) && max_tip_depth > 0) {
          phy$edge.length <- phy$edge.length / max_tip_depth
     }
     phy
}

phy <- build_phylo(nodes)
tree_meta <- nodes %>%
     transmute(
          label = if_else(node_type == "tip", tip_label, node_id),
          node_type,
          tip_label,
          is_fitch_origin,
          is_reference,
          country_iso3,
          year,
          year_bin,
          read_track,
          confidence_track,
          mechanism_class = factor(mechanism_class, levels = names(mechanism_colors)),
          architecture = factor(architecture, levels = names(architecture_colors)),
          tip_state = factor(tip_state, levels = names(prn_colors)),
          country_focus
     )

# panel A: fan phylogeny with annotation rings -------------------------

p_tree_base <- ggtree(
     phy,
     # Topology-only fan: encoding SNP branch lengths detaches dense annotation
     # rings from many tips at final figure size.
     branch.length = "none",
     layout = "fan",
     open.angle = 8,
     size = 0.13,
     colour = FIGURE_GREY
)

tip_track_meta <- tree_meta %>%
     filter(node_type == "tip") %>%
     mutate(
          architecture_ring = case_when(
               as.character(tip_state) == "Insufficient data" ~ "Unresolved",
               is.na(as.character(architecture)) ~ "No disrupted event",
               as.character(architecture) == "IS481 gap1043" ~ "IS481 gap1043",
               str_detect(as.character(architecture), "^IS481") ~ "Other IS481",
               str_detect(as.character(architecture), "^Rearr") ~ "Rearrangement",
               TRUE ~ "Other disruption"
          ),
          country_ring = case_when(
               country_focus %in% c("USA", "CHN", "NZL", "JPN") ~ country_focus,
               country_focus %in% c("AUS", "GBR", "FRA", "BRA", "CZE", "FIN") ~ "Other selected",
               TRUE ~ "Other"
          ),
          year_ring = replace_na(year_bin, "Year unresolved")
     ) %>%
     arrange(match(label, phy$tip.label))

# Keep only the non-redundant main-text context rings. Tip points already encode
# PRN state, and the tree topology carries lineage context.
tip_tracks <- as.data.frame(tip_track_meta %>%
                                 select(
                                      Architecture = architecture_ring,
                                      Country = country_ring,
                                      Year = year_ring
                                 ))
rownames(tip_tracks) <- tip_track_meta$label

track_values <- sort(unique(unlist(tip_tracks, use.names = FALSE)))
track_palette <- c(
     prn_colors,
     architecture_colors,
     country_role_colors,
     "Other" = FIGURE_LIGHT_GREY,
     "Other selected" = FIGURE_MID_GREY,
     "No disrupted event" = FIGURE_PANEL_FILL,
     "Unresolved" = FIGURE_LIGHT_GREY,
     "Other IS481" = unname(npg_colors["orange"]),
     "Rearrangement" = unname(npg_colors["blue"]),
     "Other disruption" = unname(npg_colors["green"]),
     "Reference" = FIGURE_INK,
     "<1990" = FIGURE_LIGHT_GREY,
     "1990s" = FIGURE_GREY,
     "2000s" = FIGURE_MID_GREY,
     "2010s" = unname(npg_colors["teal"]),
     "2020+" = unname(npg_colors["blue"]),
     "Year unresolved" = FIGURE_PANEL_FILL
)
missing_track_values <- setdiff(track_values, names(track_palette))
if (length(missing_track_values) > 0) {
     track_palette <- c(track_palette, setNames(rep(FIGURE_LIGHT_GREY, length(missing_track_values)), missing_track_values))
}

pA_tree <- p_tree_base %<+% tree_meta +
     geom_tippoint(aes(colour = tip_state), size = 0.82, alpha = 0.96, na.rm = TRUE) +
     geom_point2(
          aes(subset = is_fitch_origin, shape = "Fitch package"),
          size = 1.75,
          stroke = 0.42,
          fill = "white",
          colour = unname(npg_colors["red"]),
          na.rm = TRUE
     ) +
     scale_colour_manual(values = prn_colors, na.value = FIGURE_GREY, name = "Tip state") +
     scale_shape_manual(values = c("Fitch package" = 21), name = NULL) +
     guides(
          colour = guide_legend(order = 1, override.aes = list(size = 2)),
          shape = guide_legend(
               order = 2,
               override.aes = list(fill = "white", colour = unname(npg_colors["red"]), size = 2.4, stroke = 0.45)
          )
     )

pA_body <- ggtree::gheatmap(
     pA_tree,
     tip_tracks,
     offset = 0.06,
     width = 0.32,
     color = NA,
     colnames = FALSE,
     font.size = 1.6
) +
     scale_fill_manual(values = track_palette, na.value = FIGURE_PANEL_FILL, guide = "none") +
     theme_nature_tree(base_size = fig3_base_size) +
     theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

panel_a_key_rows <- bind_rows(
     tibble::tibble(
          track = "State",
          label = c("intact", "disrupt.", "insuff.", "ref.", "package"),
          fill = c(unname(prn_colors["Intact"]), unname(prn_colors["Disrupted"]), unname(prn_colors["Insufficient data"]), unname(prn_colors["Reference"]), "white"),
          stroke = c(rep(FIGURE_INK, 4), unname(npg_colors["red"])),
          shape = c(rep(22, 4), 21),
          y = 4
     ),
     tibble::tibble(
          track = "Arch.",
          label = c("intact/no event", "unres.", "IS481 1043", "other IS481", "rearr.", "other disrupt."),
          fill = c(
               FIGURE_PANEL_FILL,
               FIGURE_LIGHT_GREY,
               unname(architecture_colors["IS481 gap1043"]),
               unname(npg_colors["orange"]),
               unname(npg_colors["blue"]),
               unname(npg_colors["green"])
          ),
          stroke = FIGURE_INK,
          shape = 22,
          y = 3
     ),
     tibble::tibble(
          track = "Country",
          label = c("USA", "CHN", "NZL", "JPN", "other selected", "Other"),
          fill = c(unname(country_role_colors[c("USA", "CHN", "NZL", "JPN")]), FIGURE_MID_GREY, FIGURE_LIGHT_GREY),
          stroke = FIGURE_INK,
          shape = 22,
          y = 2
     ),
     tibble::tibble(
          track = "Year",
          label = c("<1990", "1990s", "2000s", "2010s", "2020+", "unres.", "ref."),
          fill = c(FIGURE_LIGHT_GREY, FIGURE_GREY, FIGURE_MID_GREY, unname(npg_colors["teal"]), unname(npg_colors["blue"]), FIGURE_PANEL_FILL, FIGURE_INK),
          stroke = FIGURE_INK,
          shape = 22,
          y = 1
     )
) %>%
     group_by(track) %>%
     mutate(x = seq(1.45, 11.25, length.out = n())) %>%
     ungroup()

panel_a_track_labels <- panel_a_key_rows %>%
     distinct(track, y)

pA_key <- ggplot(panel_a_key_rows, aes(x = x, y = y)) +
     geom_text(
          data = panel_a_track_labels,
          aes(x = 0.05, y = y, label = track),
          inherit.aes = FALSE,
          hjust = 0,
          size = 1.75,
          fontface = "bold",
          colour = FIGURE_TEXT_COLOUR
     ) +
     geom_point(aes(fill = fill, colour = stroke, shape = shape), size = 1.9, stroke = 0.22) +
     geom_text(aes(x = x + 0.17, label = label), hjust = 0, size = 1.55, colour = FIGURE_TEXT_COLOUR) +
     scale_fill_identity() +
     scale_colour_identity() +
     scale_shape_identity() +
     coord_cartesian(xlim = c(0, 12.45), ylim = c(0.45, 5.05), clip = "off") +
     theme_void(base_size = fig3_base_size) +
     theme(plot.margin = margin(0, 0, 0, 0))

pA_composite <- pA_body / pA_key + plot_layout(heights = c(1, 0.26))
pA <- wrap_elements(full = pA_composite)

# panel B: ASR Fitch package count by scenario -------------------------

asr_summary <- asr %>%
     filter(!is.na(fitch_origin_events)) %>%
     group_by(scenario_class) %>%
     summarise(
          min_origin = min(fitch_origin_events, na.rm = TRUE),
          median_origin = median(fitch_origin_events, na.rm = TRUE),
          max_origin = max(fitch_origin_events, na.rm = TRUE),
          n_scenarios = n(),
          .groups = "drop"
     ) %>%
     mutate(scenario_class = factor(scenario_class, levels = scenario_levels))

pB <- ggplot(asr_summary, aes(median_origin, scenario_class)) +
     geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
     geom_segment(aes(x = min_origin, xend = max_origin, yend = scenario_class), colour = FIGURE_GREY, linewidth = 0.68) +
     geom_point(size = 2.5, fill = npg_colors["red"], shape = 21, colour = FIGURE_INK, stroke = 0.18) +
     scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
     labs(x = "Fitch package count", y = NULL) +
     theme_nature(base_size = fig3_base_size) +
     theme(axis.text.y = element_text(size = 4.9))

# panel C: Stochastic mapping package count ----------------------------

simmap_summary <- simmap %>%
     group_by(scenario_class) %>%
     summarise(
          lower = min(stochastic_origin_count_lower_95, na.rm = TRUE),
          median = median(stochastic_origin_count_median, na.rm = TRUE),
          upper = max(stochastic_origin_count_upper_95, na.rm = TRUE),
          .groups = "drop"
     ) %>%
     mutate(scenario_class = factor(scenario_class, levels = scenario_levels))

pC <- ggplot(simmap_summary, aes(median, scenario_class)) +
     geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
     geom_segment(aes(x = lower, xend = upper, yend = scenario_class), colour = FIGURE_GREY, linewidth = 0.68) +
     geom_point(size = 2.5, fill = npg_colors["purple"], shape = 21, colour = FIGURE_INK, stroke = 0.18) +
     labs(x = "Stochastic-map package count (95% CI)", y = NULL) +
     theme_nature(base_size = fig3_base_size) +
     theme(axis.text.y = element_text(size = 4.9))

# panel D: Event-specific acquisition packages ------------------------

event_package_plot <- event_specific %>%
     filter(rank_by_genome_burden <= 3) %>%
     arrange(desc(sample_count)) %>%
     mutate(
          event_display = .env$event_label(prn_event_id),
          event_display = factor(event_display, levels = rev(unique(event_display))),
          package_label = paste0(
               acquisition_package_count, " pkg; ",
               n_country_year_cells, " cells\n",
               comma(sample_count), " genomes; ",
               n_mlst_st, " ST"
          ),
          label_x = acquisition_package_count + 0.48
     )

pD <- ggplot(event_package_plot, aes(acquisition_package_count, event_display)) +
     geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.25, colour = FIGURE_MID_GREY) +
     geom_segment(aes(x = 0, xend = acquisition_package_count, yend = event_display),
                  colour = FIGURE_RULE_COLOUR, linewidth = 0.58) +
     geom_point(aes(size = sample_count, fill = event_display),
                shape = 21, colour = FIGURE_INK, stroke = 0.20, alpha = 0.96) +
     geom_text(aes(x = label_x, label = package_label),
               hjust = 0, size = 1.45, lineheight = 0.82, colour = FIGURE_MUTED_TEXT) +
     scale_x_continuous(limits = c(0, max(event_package_plot$label_x, na.rm = TRUE) + 4.3),
                        breaks = 0:9,
                        expand = expansion(mult = c(0.02, 0.02))) +
     scale_size_area(max_size = 4.2, guide = "none") +
     scale_fill_manual(values = architecture_colors, guide = "none") +
     labs(x = "Minimum tree-level\nacquisition packages", y = NULL) +
     coord_cartesian(clip = "off") +
     theme_nature(base_size = fig3_base_size) +
     theme(
          axis.text.y = element_text(size = 5.2, face = "bold"),
          plot.margin = margin(3, 18, 3, 3)
     )

# assembly: multi-panel evidence layout --------------------------------

fig3_layout <- "
AAB
AAC
AAD
"

fig3 <- free(pA) + pB + pC + free(pD) +
     plot_layout(
          design = fig3_layout,
          widths = c(1, 1, 0.78),
          heights = c(0.9, 0.9, 0.9)
     ) +
     plot_annotation(tag_levels = "a") &
     theme(
          plot.tag = element_text(face = "bold", size = fig3_tag_size, colour = FIGURE_TEXT_COLOUR),
          plot.tag.position = c(0, 1),
          plot.title.position = "plot",
          plot.margin = margin(3, 3, 3, 3)
     )

save_nc_pdf(fig3, "fig03_repeated_origin_phylogeny.pdf", height = 5.5)
save_nc_png(fig3, "fig03_repeated_origin_phylogeny.png", height = 5.5)
