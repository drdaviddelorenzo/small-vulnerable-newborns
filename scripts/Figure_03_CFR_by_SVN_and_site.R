# ==============================================================================
# figure_03_cfr_by_svn_and_site.R
# FIGURE 3 (manuscript) - Case Fatality Rate by SVN phenotype and site,
#   computed on the ANALYSIS COHORT (n=30,449; KCH 10,517, SMCH 19,932),
#   the same population as manuscript Table 5.
#
# Unlike fig5A (which uses the larger all-matched-admissions denominator and so
# differs slightly from Table 5), this reads the cohort-basis CFRs with Wilson
# CIs straight from dep_figure03_cfr_comparison.csv (classification == "SVN").
#
# Source: outputs/<population>/nnu/dep_figure03_cfr_comparison.csv
# ASCII-only output (DSH requirement).
# ==============================================================================

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

cat("\n===== CFR by SVN phenotype and site (analysis cohort) =====\n")

# ------------------------------------------------------------------------------
# LOAD DATA - cohort-basis CFRs with Wilson CIs (already 0-100 %)
# ------------------------------------------------------------------------------
cfr <- read_csv(
  file.path(OUTPUT_DIR_NNU, "dep_figure03_cfr_comparison.csv"),
  show_col_types = FALSE
)

cfr_fac <- cfr %>%
  filter(classification == "SVN", category %in% SVN_CATS) %>%
  mutate(
    svn_cat      = svn_factor(category),
    facility_lbl = recode(facility, !!!FACILITY_LABELS),
    # point CFR recomputed at full precision so labels match Table 5 rounding,
    # expressed PER 1,000 admissions (consistent with the rest of the manuscript).
    # Stored cfr_pct / ci_lower / ci_upper are on the 0-100 % scale, so x10 -> per 1,000.
    cfr   = 1000 * as.numeric(n_deaths) / as.numeric(n_total),
    ci_lo = 10 * as.numeric(ci_lower),
    ci_hi = 10 * as.numeric(ci_upper)
  )

# ------------------------------------------------------------------------------
# FIGURE 3 (reuses fig5A styling)
# ------------------------------------------------------------------------------
fig3 <- ggplot(cfr_fac, aes(
  x = svn_cat, y = cfr,
  fill = facility_lbl
)) +
  geom_col(
    position = position_dodge(width = 0.72), width = 0.68
  ) +
  geom_errorbar(
    aes(ymin = ci_lo, ymax = ci_hi),
    position = position_dodge(width = 0.72),
    width = 0.28, colour = "grey35", linewidth = 0.65
  ) +
  geom_text(
    aes(y = ci_hi + 8, label = sprintf("%.0f", cfr)),
    position = position_dodge(width = 0.72),
    vjust = 0, size = 2.8, colour = "grey20"
  ) +
  scale_fill_manual(
    # Blue/red palette matching Figure 9 (HIV distribution).
    values = c(
      "KCH (Malawi)"    = "#74ADD1",
      "SMCH (Zimbabwe)" = "#B2182B"
    ),
    name = "Facility"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title    = "Case fatality rate by SVN phenotype and site",
    subtitle = "Error bars: 95% Wilson confidence intervals.",
    x        = NULL,
    y        = "Case fatality rate (per 1,000 admissions)",
    caption  = NULL
  ) +
  theme_svn(legend_position = "right") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_figure(fig3, "figure_03_cfr_by_svn_and_site", width = 8, height = 4.8)

# ------------------------------------------------------------------------------
# VERIFY rendered values against manuscript Table 5
# ------------------------------------------------------------------------------
chk <- cfr_fac %>%
  arrange(facility, svn_cat) %>%
  transmute(facility, category = as.character(svn_cat),
            label = sprintf("%.0f (%.0f-%.0f)", cfr, ci_lo, ci_hi))
cat("\nRendered CFR per 1,000 (95% CI) by site and phenotype:\n")
print(chk, n = nrow(chk))

cat("CFR by SVN phenotype and site: complete.\n")
