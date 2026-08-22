# ==============================================================================
# Purpose: generates the case-fatality data that Figure_03_CFR_by_SVN_and_site.R
# reads.
#
# Definitions: case fatality (deaths per known outcome, %) with Wilson 95%
# CI, by site, for both classification systems (birth-weight LBW categories
# and SVN phenotypes), on the analysis cohort (valid BW + GA). Deaths use the
# shared outcome_bin (infants brought in dead counted as deaths). Respects
# POPULATION ("all" default, "USS" for the ultrasound subset).
#
# Output: dep_figure03_cfr_comparison.csv
#   columns: facility, category, n_total, n_known, n_deaths, n_alive,
#            cfr_pct, ci_lower, ci_upper, classification ("LBW" | "SVN")
# ==============================================================================
# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")
suppressPackageStartupMessages({library(dplyr)})
cat("\n=== 05b_cfr_comparison.R ===\n")

data <- process_nnu_data(load_nnu_data())$svn_filtered %>%
    mutate(bw_fac  = factor(bw_cat,       levels = BW_LEVELS),
           svn_fac = factor(svn_category, levels = SVN_LEVELS))

cfr_group <- function(df, gcol, label) {
    df %>%
        group_by(facility, category = .data[[gcol]]) %>%
        summarise(n_total  = n(),
                  n_known  = sum(!is.na(outcome_bin)),
                  n_deaths = sum(outcome_bin == 1, na.rm = TRUE),
                  n_alive  = sum(outcome_bin == 0, na.rm = TRUE),
                  cfr_pct  = round(n_deaths / n_known * 100, 2),
                  ci_lower = round(mapply(function(d, t) prop_ci_wilson(d, t)[["lower"]] * 100,
                                          n_deaths, n_known), 2),
                  ci_upper = round(mapply(function(d, t) prop_ci_wilson(d, t)[["upper"]] * 100,
                                          n_deaths, n_known), 2),
                  classification = label, .groups = "drop") %>%
        filter(!is.na(category))
}

cfr_all <- bind_rows(cfr_group(data, "bw_fac",  "LBW"),
                     cfr_group(data, "svn_fac", "SVN"))
write.csv(cfr_all, file.path(OUTPUT_DIR_NNU, "dep_figure03_cfr_comparison.csv"), row.names = FALSE)

cat(sprintf("Cohort: %d; deaths (outcome_bin==1): %d\n",
            nrow(data), sum(data$outcome_bin == 1, na.rm = TRUE)))
print(as.data.frame(cfr_all), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(OUTPUT_DIR_NNU, "dep_figure03_cfr_comparison.csv")))
cat("05b_cfr_comparison.R complete.\n")
