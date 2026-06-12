# Shared helpers for the manuscript visual rebuild.

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(tidyr)
library(patchwork)

as_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

pct_label <- function(x, accuracy = 1) {
  scales::percent(x, accuracy = accuracy)
}

short_epoch_type <- function(epoch_type, prn_in_formulation = NA_character_) {
  out <- dplyr::case_when(
    epoch_type %in% c("wP_only", "wp_only") ~ "wP only",
    epoch_type %in% c("aP_with_PRN", "routine_ap_prn_positive") ~ "aP + PRN",
    epoch_type %in% c("aP_without_PRN", "routine_ap_prn_negative") ~ "aP - PRN",
    epoch_type %in% c("mixed_brand_heterogeneous") ~ "Mixed PRN",
    epoch_type %in% c("transition_mixed") ~ "Transition",
    prn_in_formulation %in% c("yes", "prn_positive") ~ "aP + PRN",
    prn_in_formulation %in% c("no", "prn_negative") ~ "aP - PRN",
    prn_in_formulation %in% c("mixed") ~ "Mixed PRN",
    TRUE ~ "Uncertain"
  )
  factor(out, levels = names(programme_epoch_colors))
}

state_label <- function(x) {
  dplyr::case_when(
    x %in% c("disrupted", "Disrupted", "coding_disrupted") ~ "Disrupted",
    x %in% c("intact", "Intact") ~ "Intact",
    x %in% c("Reference", "reference") ~ "Reference",
    TRUE ~ "Insufficient data"
  )
}

mechanism_label <- function(x) {
  dplyr::case_when(
    str_detect(x, "is481") ~ "IS481 insertion",
    str_detect(x, "inversion|rearrangement") ~ "Rearrangement",
    str_detect(x, "insufficient") ~ "Insufficient data",
    TRUE ~ "Other disruptions"
  )
}

event_label <- function(x) {
  dplyr::case_when(
    str_detect(x, "gap1043") ~ "IS481 gap1043",
    str_detect(x, "gap1042") ~ "gap1042",
    str_detect(x, "gap1045") ~ "gap1045",
    str_detect(x, "gap1041") & str_detect(x, "is481") ~ "IS481 gap1041",
    str_detect(x, "gap1041") ~ "Other gap1041",
    str_detect(x, "gap1040") ~ "Other gap1040",
    str_detect(x, "gap1044") ~ "Other gap1044",
    str_detect(x, "gap54") ~ "gap54",
    str_detect(x, "gap204") ~ "gap204",
    str_detect(x, "cov58") ~ "Rearr. cov58",
    str_detect(x, "cov91") ~ "Rearr. cov91",
    str_detect(x, "cov94") ~ "Rearr. cov94",
    str_detect(x, "insufficient") ~ "Insufficient",
    TRUE ~ str_replace_all(x, "^prn_evt_|coding_disrupted_", "")
  )
}

direction_label <- function(x) {
  dplyr::case_when(
    x %in% c("upward", "increase", "increase_with_overlapping_bounds") ~ "Upward",
    x %in% c("downward", "decrease", "decrease_with_overlapping_bounds") ~ "Downward",
    x %in% c("no_change") ~ "No change",
    TRUE ~ "Uncertain"
  )
}

evidence_label <- function(x) {
  dplyr::case_when(
    x %in% c("strong", "clear_post-origin_amplification", "clear_post-detection_amplification",
             "strongly selection-compatible") ~ "Strong",
    x %in% c("present", "bounded", "compatible but bounded", "limited_or_nonpersistent_post-detection_change",
             "biologically interesting but uncertain") ~ "Bounded",
    x %in% c("limited", "event_level_only", "epidemiologically estimable but not triangulated") ~ "Limited",
    x %in% c("none", "not_observed", "context only") ~ "None",
    TRUE ~ "Uncertain"
  )
}

ordered_countries <- function(countries) {
  preferred <- c("USA", "NZL", "JPN", "AUS", "GBR", "CHN", "FRA", "BRA", "CZE", "FIN")
  present <- preferred[preferred %in% countries]
  c(present, sort(setdiff(countries, present)))
}

panel_title <- function(title) {
  ggtitle(title) +
    theme(plot.title = element_text(face = "bold", size = 7.4, margin = margin(b = 3)))
}

compact_legend_bottom <- theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(5.5, "pt")
)
