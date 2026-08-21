# ==============================================================================
# figure_07_temperature_categorical_ors_by_svn.R
# FIGURE 8 - Categorical admission-temperature odds ratios by SVN group
#            KCH (Malawi) and SMCH (Zimbabwe), 2022-2025
#
# Shows the odds of neonatal death in each temperature category versus the
# normal range (36.5-37.5 degC), for the overall cohort and each SVN group.
# Demonstrates that hypothermia raises risk in every phenotype, while only
# Term-AGA turns upward (U-shape) at hyperthermia.
#
# Source data: outputs/<pop>/dep_table09_temp_cat_or_overall.csv
#              outputs/<pop>/dep_table09_temp_cat_or_by_svn.csv
#              (produced by 36_temperature_categorical_OR.R)
# NOTE: ASCII characters only.
# ==============================================================================

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

cat("\n===== Temperature categorical OR by SVN group =====\n")

ov <- read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_overall.csv"),
               show_col_types = FALSE) %>%
  transmute(series = "Overall", temp_cat, OR)

bg <- read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_by_svn.csv"),
               show_col_types = FALSE) %>%
  transmute(series = group, temp_cat, OR)

cat_levels <- c("1_SevereHypothermia", "2_Hypothermia", "3_LowNormal",
                "4_Normal (Reference)", "5_Hyperthermia")
cat_labels <- c("Severe\nhypothermia", "Hypothermia", "Mild\nhypothermia",
                "Normal\n(reference)", "Hyperthermia")

plot_df <- bind_rows(ov, bg) %>%
  filter(temp_cat %in% cat_levels) %>%
  mutate(
    temp_f = factor(temp_cat, levels = cat_levels, labels = cat_labels),
    series = factor(series, levels = c("Overall", SVN_LEVELS))
  )

series_cols <- c("Overall" = "#222222", SVN_COLOURS)

fig8 <- ggplot(plot_df, aes(x = temp_f, y = OR, colour = series, group = series)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_colour_manual(values = series_cols, name = NULL) +
  labs(
    title    = "Temperature-mortality by SVN group",
    subtitle = "Odds of neonatal death vs normal range (36.5-37.5 degC); log scale",
    x        = NULL,
    y        = "Odds ratio (log scale)",
    caption  = svn_caption("Reference: normal range. Source: script 36 (categorical temperature OR).")
  ) +
  theme_svn(legend_position = "right")

save_figure(fig8, "figure_07_temperature_categorical_ors_by_svn", width = 9, height = 6)

# ------------------------------------------------------------------------------
# COMPANION VALUES TABLE: OR (95% CI) and Wald p vs the normal range, per group
# Written so the figure's underlying numbers and test are reproducible from R.
# ------------------------------------------------------------------------------
plain_lab <- c("1_SevereHypothermia" = "Severe hypothermia",
               "2_Hypothermia"        = "Hypothermia",
               "3_LowNormal"          = "Mild hypothermia",
               "4_Normal (Reference)" = "Normal (reference)",
               "5_Hyperthermia"       = "Hyperthermia")

tbl8 <- bind_rows(
  read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_overall.csv"), show_col_types = FALSE) %>%
    mutate(group = "Overall"),
  read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_by_svn.csv"), show_col_types = FALSE)
) %>%
  filter(temp_cat %in% names(plain_lab)) %>%
  transmute(
    group                = factor(group, levels = c("Overall", SVN_LEVELS)),
    temperature_category = unname(plain_lab[temp_cat]),
    OR                   = round(OR, 3),
    ci_lower             = round(ci_lower, 3),
    ci_upper             = round(ci_upper, 3),
    p_value_vs_normal    = round(p_value, 4)
  ) %>%
  arrange(group)

write_csv(tbl8, file.path(OUTPUT_DIR, "figure_07_temperature_by_svn_table.csv"))
cat("  Saved: figure_07_temperature_by_svn_table.csv\n")
cat("Temperature categorical OR by SVN: complete.\n")
