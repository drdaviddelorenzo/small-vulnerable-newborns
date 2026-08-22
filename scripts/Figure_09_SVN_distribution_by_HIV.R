# ==============================================================================
# Purpose: Figure 9 -- SVN phenotype distribution by maternal HIV status,
# combined sites. HIV-exposed infants shift away from Term-AGA towards the
# SGA phenotypes (chi-square 21.3, p<0.001), consistent with an effect via
# fetal growth rather than direct mortality.
#
# Reads: dep_figure09_svn_distribution_by_hiv.csv
#   (written by dep_Figure09_hiv_distribution_by_svn.R -- run that first)
# ==============================================================================

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

cat("\n===== SVN distribution by maternal HIV status =====\n")

plot_df <- read_csv(file.path(OUTPUT_DIR, "dep_figure09_svn_distribution_by_hiv.csv"),
                    show_col_types = FALSE) %>%
  filter(facility == "Combined",
         hiv_status %in% c("HIV Negative", "HIV Positive")) %>%
  mutate(
    svn_f = svn_factor(svn_category),
    hiv   = factor(ifelse(hiv_status == "HIV Positive", "HIV-exposed", "HIV-unexposed"),
                   levels = c("HIV-unexposed", "HIV-exposed"))
  )

fig9 <- ggplot(plot_df, aes(x = svn_f, y = pct_within_hiv_group, fill = hiv)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", pct_within_hiv_group)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("HIV-unexposed" = "#74ADD1", "HIV-exposed" = "#B2182B"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "HIV-exposed infants shift toward SGA phenotypes",
    subtitle = "SVN phenotype distribution within HIV groups (combined sites, tested mothers)",
    x        = NULL,
    y        = "% of HIV group",
    caption  = svn_caption("Source: script 37 (HIV distribution by SVN).")
  ) +
  theme_svn(legend_position = "top") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_figure(fig9, "figure_09_svn_distribution_by_hiv", width = 8, height = 6)

# ------------------------------------------------------------------------------
# COMPANION VALUES TABLE: counts and % within each HIV group, plus the
# chi-square test of association (from dep_Figure09_hiv_distribution_by_svn.R).
# ------------------------------------------------------------------------------
tbl9 <- plot_df %>%
  transmute(
    svn_category         = as.character(svn_f),
    hiv_group            = as.character(hiv),
    n                    = n,
    pct_within_hiv_group = pct_within_hiv_group
  ) %>%
  arrange(hiv_group, svn_category)
write_csv(tbl9, file.path(OUTPUT_DIR, "figure_09_hiv_distribution_table.csv"))

test9 <- read_csv(file.path(OUTPUT_DIR, "dep_figure09_chisq_results.csv"), show_col_types = FALSE) %>%
  filter(comparison == "HIV status ~ svn_category")
write_csv(test9, file.path(OUTPUT_DIR, "figure_09_hiv_distribution_chisq_test.csv"))

cat("  Saved: figure_09_hiv_distribution_table.csv + figure_09_hiv_distribution_chisq_test.csv\n")
cat("SVN distribution by maternal HIV status: complete.\n")
