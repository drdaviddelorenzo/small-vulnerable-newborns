# ==============================================================================
# Table_11_cfr_by_resuscitation.R
# Manuscript TABLE 11 - Case fatality by resuscitation requirement at birth.
#
# Fixes the reviewer query "why does this only add to 17,385?": the previous
# table showed only the four KNOWN resuscitation categories and hid the large
# not-recorded group in a footnote. This version:
#   - runs on the 30,449 analysis cohort (svn_filtered, four SVN phenotypes,
#     valid outcome) - the SAME population as every other table;
#   - shows "Not recorded" as its own row, so every admission is accounted for
#     and the rows sum to 30,449 / 5,348 deaths;
#   - reports n, deaths and CFR per 1,000 admissions (whole numbers) with 95%
#     Wilson confidence intervals.
#
# Resuscitation coding follows Table_09 (variable `resus`), but keeps blank /
# missing / sentinel values as an explicit "Not recorded" level instead of
# folding them into "None".
#
# Output (outputs/): table_10_cfr_by_resuscitation.csv
# ASCII only (DSH requirement).
# ==============================================================================
# --- Locate shared setup after the July-2026 SCRIPTS reorganisation -----------
# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")
OUT2 <- OUTPUT_DIR   # 4DSH: alias so ported OUT2 references resolve to outputs/<pop>/tables
suppressPackageStartupMessages({library(dplyr); library(stringr)})
cat("\n=== Table_11_cfr_by_resuscitation.R ===\n")

# Wilson 95% CI for a proportion, returned on the per-1,000 scale.
wilson1000 <- function(x, n){
    if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
    p <- x / n; z <- 1.96; d <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / d
    hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
    c((ctr - hw) * 1000, (ctr + hw) * 1000)
}

PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")

# ------------------------------------------------------------------------------
# 30,449 analysis cohort with resuscitation categorised (Not recorded kept)
# ------------------------------------------------------------------------------
nnu <- process_nnu_data(load_nnu_data())
if (!"resus" %in% names(nnu$svn_filtered))
    stop("Column 'resus' not found in svn_filtered.")

cohort <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(
        resus_clean = toupper(trimws(as.character(resus))),
        resus_cat = case_when(
            is.na(resus) | resus_clean %in% c("", "NA", "NR", "U", "UNK", "UNKNOWN",
                                              "-6", "-7", "-8", "-9") ~ "Not recorded",
            resus_clean %in% c("NONE", "NO", "N", "FALSE", "0")       ~ "None",
            grepl("CPR|CARDIO", resus_clean)                          ~ "Cardiopulmonary resuscitation",
            grepl("BVM|BAG.*VALVE|MASK", resus_clean)                 ~ "Bag-valve-mask ventilation",
            grepl("STIM|SUCT|O2|OXYGEN", resus_clean)                 ~ "Stimulation / suction / oxygen",
            TRUE                                                       ~ "Not recorded"),
        resus_cat = factor(resus_cat,
            levels = c("None", "Stimulation / suction / oxygen",
                       "Bag-valve-mask ventilation", "Cardiopulmonary resuscitation",
                       "Not recorded")))

# ------------------------------------------------------------------------------
# n, deaths, CFR per 1,000 (95% Wilson CI) by category + Total
# ------------------------------------------------------------------------------
by_cat <- cohort %>%
    group_by(resus_cat, .drop = FALSE) %>%
    summarise(n = n(), deaths = sum(outcome_bin == 1), .groups = "drop")

tot <- data.frame(resus_cat = "Total", n = nrow(cohort),
                  deaths = sum(cohort$outcome_bin == 1), stringsAsFactors = FALSE)

tab <- bind_rows(as.data.frame(by_cat), tot) %>%
    rowwise() %>%
    mutate(cfr_per1000 = ifelse(n > 0, round(1000 * deaths / n), NA_real_),
           ci_lo = round(wilson1000(deaths, n)[1]),
           ci_hi = round(wilson1000(deaths, n)[2]),
           cfr_ci = ifelse(is.na(cfr_per1000), "-",
                           sprintf("%d (%d-%d)", cfr_per1000, ci_lo, ci_hi))) %>%
    ungroup() %>%
    transmute(`Resuscitation at birth` = as.character(resus_cat),
              `Admissions, n` = n,
              `Deaths, n` = deaths,
              `CFR per 1,000 (95% CI)` = cfr_ci)

write.csv(tab, file.path(OUT2, "table_10_cfr_by_resuscitation.csv"), row.names = FALSE)
cat("Saved: table_10_cfr_by_resuscitation.csv\n")
print(tab, row.names = FALSE)

# Reconciliation check
stopifnot(sum(by_cat$n) == nrow(cohort),
          sum(by_cat$deaths) == sum(cohort$outcome_bin == 1))
cat(sprintf("Reconciles: admissions %d, deaths %d\n",
            nrow(cohort), sum(cohort$outcome_bin == 1)))
cat("Table_11_cfr_by_resuscitation.R complete.\n")
