# 37_hiv_distribution_by_svn.R
# HIV Prevalence and Distribution by SVN and LBW Category
#
# Purpose:
#   Addresses meeting question (9 June 2026, Hannah Calder, Rachel Poffley,
#   David de Lorenzo): "what is the distribution of HIV exposure across SVN
#   and LBW categories, and is maternal HIV associated with SVN category?"
#
#   This is distinct from Scripts 32/34 which show CFR stratified by HIV
#   status. Script 37 asks the REVERSE question: do HIV-positive mothers
#   have babies with a different SVN/LBW distribution?
#
# Analyses:
#   1. HIV positivity rate within each SVN category (among tested)
#   2. HIV positivity rate within each LBW category (among tested)
#   3. SVN distribution within HIV+ vs HIV- groups (column percentages)
#   4. LBW distribution within HIV+ vs HIV- groups
#   5. Chi-squared test: association between SVN category and HIV status
#   6. Chi-squared test: association between LBW category and HIV status
#   7. CFR by HIV status within each SVN group (confirmation of Script 32/34)
#
# Output files (in OUTPUT_DIR):
#   dep_figure09_hiv_prevalence_by_svn.csv    — HIV+ rate within each SVN group
#   dep_figure09_hiv_prevalence_by_lbw.csv    — HIV+ rate within each LBW group
#   dep_figure09_svn_distribution_by_hiv.csv  — SVN distribution within HIV groups
#   dep_figure09_lbw_distribution_by_hiv.csv  — LBW distribution within HIV groups
#   dep_figure09_cfr_by_hiv_svn.csv           — CFR by HIV x SVN (pooled sites)
#   dep_figure09_chisq_results.csv            — Chi-squared test results
#   37_hiv_svn_distribution.xlsx       — Combined Excel workbook

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

library(dplyr)
library(tidyr)
library(openxlsx)

cat("=== Script 37: HIV Distribution by SVN and LBW Category ===\n\n")

# ============================================================
# 1. LOAD AND PREPARE DATA
# ============================================================

cat("Loading NNU data...\n")
raw <- load_nnu_data()
nnu <- process_nnu_data(raw)
base <- nnu$svn_filtered %>%
  filter(!is.na(outcome_bin))

cat(sprintf("Base population: %d  (KCH: %d, SMCH: %d)\n",
            nrow(base),
            sum(base$facility == "KCH"),
            sum(base$facility == "SMCH")))

# ============================================================
# 2. CLASSIFY HIV STATUS
# ============================================================

# HIV classification: Positive = tested and Reactive
# Negative = tested and Non-Reactive
# Unknown = tested but result missing/other
# Not Tested = mathivtest == FALSE or NA

classify_hiv <- function(mathivtest_col, hivtestresult_col) {
  tested  <- toupper(trimws(as.character(mathivtest_col))) %in% c("TRUE","T","YES","Y","1")
  not_tested <- toupper(trimws(as.character(mathivtest_col))) %in% c("FALSE","F","NO","N","0")
  result  <- toupper(trimws(as.character(hivtestresult_col)))

  dplyr::case_when(
    tested  & result == "R"                          ~ "HIV Positive",
    tested  & result == "NR"                         ~ "HIV Negative",
    tested  & !result %in% c("R","NR")               ~ "Tested (unknown result)",
    not_tested                                        ~ "Not Tested",
    TRUE                                              ~ "Missing"
  )
}

analysis_df <- base %>%
  mutate(
    hiv_status = classify_hiv(mathivtest, hivtestresult),
    hiv_binary = dplyr::case_when(
      hiv_status == "HIV Positive" ~ "HIV Positive",
      hiv_status == "HIV Negative" ~ "HIV Negative",
      TRUE                         ~ NA_character_
    )
  )

cat("\nHIV status distribution:\n")
print(table(analysis_df$hiv_status, useNA = "ifany"))
cat("\nHIV status by facility:\n")
print(table(analysis_df$facility, analysis_df$hiv_status, useNA = "ifany"))

# ============================================================
# 3. HIV PREVALENCE WITHIN SVN CATEGORIES
# ============================================================

cat("\n=== HIV prevalence by SVN group ===\n")

hiv_prev_svn <- analysis_df %>%
  filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
  group_by(facility, svn_category) %>%
  summarise(
    n_tested    = n(),
    n_positive  = sum(hiv_status == "HIV Positive"),
    n_negative  = sum(hiv_status == "HIV Negative"),
    pct_positive = round(n_positive / n_tested * 100, 1),
    .groups = "drop"
  ) %>%
  bind_rows(
    # Combined sites
    analysis_df %>%
      filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
      group_by(svn_category) %>%
      summarise(
        facility    = "Combined",
        n_tested    = n(),
        n_positive  = sum(hiv_status == "HIV Positive"),
        n_negative  = sum(hiv_status == "HIV Negative"),
        pct_positive = round(n_positive / n_tested * 100, 1),
        .groups = "drop"
      )
  ) %>%
  arrange(facility, svn_category)

print(hiv_prev_svn)

# ============================================================
# 4. HIV PREVALENCE WITHIN LBW CATEGORIES
# ============================================================

cat("\n=== HIV prevalence by LBW group ===\n")

hiv_prev_lbw <- analysis_df %>%
  filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
  group_by(facility, bw_cat) %>%
  summarise(
    n_tested    = n(),
    n_positive  = sum(hiv_status == "HIV Positive"),
    n_negative  = sum(hiv_status == "HIV Negative"),
    pct_positive = round(n_positive / n_tested * 100, 1),
    .groups = "drop"
  ) %>%
  bind_rows(
    analysis_df %>%
      filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
      group_by(bw_cat) %>%
      summarise(
        facility    = "Combined",
        n_tested    = n(),
        n_positive  = sum(hiv_status == "HIV Positive"),
        n_negative  = sum(hiv_status == "HIV Negative"),
        pct_positive = round(n_positive / n_tested * 100, 1),
        .groups = "drop"
      )
  ) %>%
  arrange(facility, bw_cat)

print(hiv_prev_lbw)

# ============================================================
# 5. SVN DISTRIBUTION WITHIN HIV GROUPS
# ============================================================

cat("\n=== SVN distribution within HIV groups ===\n")

svn_dist_by_hiv <- analysis_df %>%
  filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
  group_by(facility, hiv_status, svn_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(facility, hiv_status) %>%
  mutate(
    total_in_hiv_group = sum(n),
    pct_within_hiv_group = round(n / total_in_hiv_group * 100, 1)
  ) %>%
  ungroup() %>%
  bind_rows(
    analysis_df %>%
      filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
      group_by(hiv_status, svn_category) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(hiv_status) %>%
      mutate(
        facility = "Combined",
        total_in_hiv_group = sum(n),
        pct_within_hiv_group = round(n / total_in_hiv_group * 100, 1)
      ) %>%
      ungroup()
  )

print(svn_dist_by_hiv)

# ============================================================
# 6. LBW DISTRIBUTION WITHIN HIV GROUPS
# ============================================================

lbw_dist_by_hiv <- analysis_df %>%
  filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
  group_by(facility, hiv_status, bw_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(facility, hiv_status) %>%
  mutate(
    total_in_hiv_group = sum(n),
    pct_within_hiv_group = round(n / total_in_hiv_group * 100, 1)
  ) %>%
  ungroup() %>%
  bind_rows(
    analysis_df %>%
      filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
      group_by(hiv_status, bw_cat) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(hiv_status) %>%
      mutate(
        facility = "Combined",
        total_in_hiv_group = sum(n),
        pct_within_hiv_group = round(n / total_in_hiv_group * 100, 1)
      ) %>%
      ungroup()
  )

# ============================================================
# 7. CFR BY HIV STATUS AND SVN GROUP
# ============================================================

cat("\n=== CFR by HIV status and SVN group ===\n")

cfr_hiv_svn <- analysis_df %>%
  filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
  group_by(facility, hiv_status, svn_category) %>%
  summarise(
    n         = n(),
    n_deaths  = sum(outcome_bin, na.rm = TRUE),
    cfr_pct   = round(n_deaths / n * 100, 1),
    .groups = "drop"
  ) %>%
  bind_rows(
    analysis_df %>%
      filter(hiv_status %in% c("HIV Positive","HIV Negative")) %>%
      group_by(hiv_status, svn_category) %>%
      summarise(
        facility  = "Combined",
        n         = n(),
        n_deaths  = sum(outcome_bin, na.rm = TRUE),
        cfr_pct   = round(n_deaths / n * 100, 1),
        .groups = "drop"
      )
  )

print(cfr_hiv_svn)

# ============================================================
# 8. CHI-SQUARED TESTS
# ============================================================

cat("\n=== Chi-squared tests ===\n")

run_chisq <- function(df, group_var, strata_var, fac_label) {
  # Guard: need at least 1 row and both HIV groups present
  if (nrow(df) == 0) {
    cat(sprintf("  Skipping chi-sq (%s ~ %s, %s): no rows\n",
                strata_var, group_var, fac_label))
    return(NULL)
  }
  ct <- table(df[[strata_var]], df[[group_var]])
  # Drop any levels with zero marginal sum (e.g. facility with no HIV results)
  ct <- ct[rowSums(ct) > 0, colSums(ct) > 0, drop = FALSE]
  if (any(dim(ct) < 2)) {
    cat(sprintf("  Skipping chi-sq (%s ~ %s, %s): insufficient variation\n",
                strata_var, group_var, fac_label))
    return(tibble(
      comparison = paste("HIV status ~", group_var),
      facility   = fac_label,
      n          = nrow(df),
      chi_sq     = NA_real_,
      df         = NA_real_,
      p_value    = NA_real_,
      note       = "Skipped: no HIV result data at this site"
    ))
  }
  cs <- suppressWarnings(chisq.test(ct))
  tibble(
    comparison = paste("HIV status ~", group_var),
    facility   = fac_label,
    n          = sum(ct),
    chi_sq     = round(cs$statistic, 3),
    df         = cs$parameter,
    p_value    = round(cs$p.value, 6),
    note       = ifelse(any(cs$expected < 5),
                        "Some expected counts < 5 (interpret with caution)",
                        "All expected counts >= 5")
  )
}

tested_only <- analysis_df %>%
  filter(hiv_status %in% c("HIV Positive","HIV Negative"))

tested_only_kch  <- tested_only %>% filter(facility == "KCH")
tested_only_smch <- tested_only %>% filter(facility == "SMCH")

cat(sprintf("  Tested (HIV pos/neg): Combined=%d, KCH=%d, SMCH=%d\n",
            nrow(tested_only), nrow(tested_only_kch), nrow(tested_only_smch)))

chisq_results <- bind_rows(
  run_chisq(tested_only,      "svn_category", "hiv_status", "Combined"),
  run_chisq(tested_only_kch,  "svn_category", "hiv_status", "KCH"),
  run_chisq(tested_only_smch, "svn_category", "hiv_status", "SMCH"),
  run_chisq(tested_only,      "bw_cat",        "hiv_status", "Combined"),
  run_chisq(tested_only_kch,  "bw_cat",        "hiv_status", "KCH"),
  run_chisq(tested_only_smch, "bw_cat",        "hiv_status", "SMCH")
)

print(chisq_results)

# ============================================================
# 9. SAVE OUTPUTS
# ============================================================

cat("\n=== Saving outputs ===\n")

write.csv(hiv_prev_svn,    file.path(OUTPUT_DIR, "dep_figure09_hiv_prevalence_by_svn.csv"),   row.names = FALSE)
write.csv(hiv_prev_lbw,    file.path(OUTPUT_DIR, "dep_figure09_hiv_prevalence_by_lbw.csv"),   row.names = FALSE)
write.csv(svn_dist_by_hiv, file.path(OUTPUT_DIR, "dep_figure09_svn_distribution_by_hiv.csv"), row.names = FALSE)
write.csv(lbw_dist_by_hiv, file.path(OUTPUT_DIR, "dep_figure09_lbw_distribution_by_hiv.csv"), row.names = FALSE)
write.csv(cfr_hiv_svn,     file.path(OUTPUT_DIR, "dep_figure09_cfr_by_hiv_svn.csv"),          row.names = FALSE)
write.csv(chisq_results,   file.path(OUTPUT_DIR, "dep_figure09_chisq_results.csv"),            row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "HIV prev by SVN");    writeData(wb, "HIV prev by SVN",    hiv_prev_svn)
addWorksheet(wb, "HIV prev by LBW");    writeData(wb, "HIV prev by LBW",    hiv_prev_lbw)
addWorksheet(wb, "SVN dist by HIV");    writeData(wb, "SVN dist by HIV",    svn_dist_by_hiv)
addWorksheet(wb, "LBW dist by HIV");    writeData(wb, "LBW dist by HIV",    lbw_dist_by_hiv)
addWorksheet(wb, "CFR HIV x SVN");      writeData(wb, "CFR HIV x SVN",      cfr_hiv_svn)
addWorksheet(wb, "Chi-sq tests");       writeData(wb, "Chi-sq tests",       chisq_results)
saveWorkbook(wb, file.path(OUTPUT_DIR, "37_hiv_svn_distribution.xlsx"), overwrite = TRUE)

cat(sprintf("\nOutputs saved to: %s\n", OUTPUT_DIR))
cat("dep_figure09_hiv_prevalence_by_svn.csv\n")
cat("dep_figure09_hiv_prevalence_by_lbw.csv\n")
cat("dep_figure09_svn_distribution_by_hiv.csv\n")
cat("dep_figure09_lbw_distribution_by_hiv.csv\n")
cat("dep_figure09_cfr_by_hiv_svn.csv\n")
cat("dep_figure09_chisq_results.csv\n")
cat("37_hiv_svn_distribution.xlsx\n")
cat("\nScript 37 complete.\n")
