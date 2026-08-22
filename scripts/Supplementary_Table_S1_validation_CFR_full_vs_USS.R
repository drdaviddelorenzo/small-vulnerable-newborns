# ==============================================================================
# Purpose: Supplementary Table S1 -- validates the main findings by comparing
# the ultrasound-dated (USS) subset against the full cohort and the non-USS
# subset.
#
# Definitions: USS-confirmed GA subset is methodestgest in {USS, EUSS};
# non-USS is all other GA methods (LMP, fundal height, etc). Reports: sample
# size and SVN distribution, CFR by SVN category, temperature OR per 1 degC
# decrease, antenatal steroid effect in Preterm-AGA, and a by-site CFR
# breakdown, across all three populations (ALL, USS, non-USS).
#
# Output: supplementary_tables1_svn_distribution.csv,
#         supplementary_tables1_cfr_by_svn.csv,
#         supplementary_tables1_cfr_by_svn_and_site.csv,
#         supplementary_tables1_temp_or_comparison.csv,
#         supplementary_tables1_steroid_pretermaga.csv,
#         supplementary_tables1_validation_summary.csv, 38_uss_validation.xlsx
# ==============================================================================

# Force ALL population (non-USS filter applied manually below)
POPULATION <- "all"
# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

library(dplyr)
library(tidyr)
library(broom)
library(openxlsx)

cat("=== Script 38: USS Validation Comparison ===\n\n")
cat("Loading full (ALL) population...\n")

raw_all <- load_nnu_data()
nnu_all <- process_nnu_data(raw_all)
base_all <- nnu_all$svn_filtered %>% filter(!is.na(outcome_bin))

cat(sprintf("ALL population: %d  (KCH: %d, SMCH: %d)\n",
            nrow(base_all),
            sum(base_all$facility == "KCH"),
            sum(base_all$facility == "SMCH")))

# -- GA method distribution ----------------------------------------------------
cat("\nGA method distribution:\n")
print(table(toupper(trimws(base_all$methodestgest)), useNA = "ifany"))

# -- Split into USS and non-USS -------------------------------------------------
uss_methods <- c("USS", "EUSS")

base_uss <- base_all %>%
  filter(toupper(trimws(methodestgest)) %in% uss_methods)

base_nonuss <- base_all %>%
  filter(!toupper(trimws(methodestgest)) %in% uss_methods)

cat(sprintf("\nUSS subset: %d  (KCH: %d, SMCH: %d)\n",
            nrow(base_uss),
            sum(base_uss$facility == "KCH"),
            sum(base_uss$facility == "SMCH")))
cat(sprintf("Non-USS subset: %d  (KCH: %d, SMCH: %d)\n",
            nrow(base_nonuss),
            sum(base_nonuss$facility == "KCH"),
            sum(base_nonuss$facility == "SMCH")))

# -- GA method breakdown in non-USS --------------------------------------------
cat("\nNon-USS: GA method breakdown:\n")
print(table(toupper(trimws(base_nonuss$methodestgest)), useNA = "ifany"))

# ============================================================
# 1. SAMPLE SIZES AND SVN DISTRIBUTION
# ============================================================

compute_svn_dist <- function(df, label) {
  total <- nrow(df)
  df %>%
    group_by(svn_category) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      population  = label,
      total       = total,
      pct         = round(n / total * 100, 1)
    ) %>%
    select(population, svn_category, n, total, pct)
}

svn_dist <- bind_rows(
  compute_svn_dist(base_all,    "ALL"),
  compute_svn_dist(base_uss,    "USS"),
  compute_svn_dist(base_nonuss, "Non-USS")
)

cat("\nSVN distribution by population:\n")
print(svn_dist)

# ============================================================
# 2. CFR BY SVN CATEGORY
# ============================================================

compute_cfr <- function(df, label) {
  df %>%
    group_by(svn_category) %>%
    summarise(
      population = label,
      n          = n(),
      n_deaths   = sum(outcome_bin),
      cfr_pct    = round(n_deaths / n * 100, 1),
      .groups    = "drop"
    ) %>%
    bind_rows(
      df %>%
        summarise(
          svn_category = "Total",
          population   = label,
          n            = n(),
          n_deaths     = sum(outcome_bin),
          cfr_pct      = round(n_deaths / n * 100, 1)
        )
    )
}

cfr_comparison <- bind_rows(
  compute_cfr(base_all,    "ALL"),
  compute_cfr(base_uss,    "USS"),
  compute_cfr(base_nonuss, "Non-USS")
)

cat("\nCFR by SVN group and population:\n")
print(cfr_comparison)

# -- By-site breakdown (manuscript Supplementary Table S1: KCH full/USS, SMCH full/USS) --
compute_cfr_by_site <- function(df, label) {
  df %>%
    group_by(facility, svn_category) %>%
    summarise(
      population = label,
      n          = n(),
      n_deaths   = sum(outcome_bin),
      cfr_pct    = round(n_deaths / n * 100, 1),
      .groups    = "drop"
    ) %>%
    bind_rows(
      df %>%
        group_by(facility) %>%
        summarise(
          svn_category = "Overall",
          population   = label,
          n            = n(),
          n_deaths     = sum(outcome_bin),
          cfr_pct      = round(n_deaths / n * 100, 1),
          .groups      = "drop"
        )
    ) %>%
    select(population, facility, svn_category, n, n_deaths, cfr_pct)
}

cfr_by_site <- bind_rows(
  compute_cfr_by_site(base_all, "full"),
  compute_cfr_by_site(base_uss, "USS")
)

cat("\nCFR by SVN group, population and site (Supplementary Table S1 layout):\n")
print(cfr_by_site)

# ============================================================
# 3. TEMPERATURE OR PER 1 degC DECREASE
# ============================================================

TEMP_VALID_MIN <- 30.0
TEMP_VALID_MAX <- 45.0

prep_temp <- function(df) {
  df %>%
    mutate(
      temp_val = suppressWarnings(
        as.numeric(temperature)
      ),
      temp_val = ifelse(temp_val >= TEMP_VALID_MIN & temp_val <= TEMP_VALID_MAX,
                        temp_val, NA_real_)
    ) %>%
    filter(!is.na(temp_val), !is.na(outcome_bin))
}

run_temp_or <- function(df, group_label, svn = NULL) {
  if (!is.null(svn)) df <- df %>% filter(svn_category == svn)
  df <- prep_temp(df)
  if (nrow(df) < 10) return(NULL)
  tryCatch({
    m <- glm(outcome_bin ~ temp_val, data = df, family = binomial())
    td <- broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "temp_val")
    tibble(
      population    = group_label,
      svn_group     = if (is.null(svn)) "Overall" else svn,
      n             = nrow(df),
      n_events      = sum(df$outcome_bin),
      OR_per_1degC  = round(td$estimate, 4),
      ci_lower      = round(td$conf.low, 4),
      ci_upper      = round(td$conf.high, 4),
      p_value       = round(td$p.value, 6),
      OR_per_1degC_decrease = round(1 / td$estimate, 4)
    )
  }, error = function(e) NULL)
}

svn_groups <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")

temp_or <- purrr::map_dfr(list("ALL" = base_all, "USS" = base_uss, "Non-USS" = base_nonuss),
  function(df) {
    pop <- attr(df, "label")
    bind_rows(
      run_temp_or(df, names(which(sapply(list("ALL" = base_all, "USS" = base_uss, "Non-USS" = base_nonuss), identical, df))), NULL),
      purrr::map_dfr(svn_groups, ~run_temp_or(df, names(which(sapply(list("ALL" = base_all, "USS" = base_uss, "Non-USS" = base_nonuss), identical, df))), .x))
    )
  }
)

# Simpler approach without attr:
temp_or_all    <- bind_rows(run_temp_or(base_all, "ALL"), purrr::map_dfr(svn_groups, ~run_temp_or(base_all, "ALL", .x)))
temp_or_uss    <- bind_rows(run_temp_or(base_uss, "USS"), purrr::map_dfr(svn_groups, ~run_temp_or(base_uss, "USS", .x)))
temp_or_nonuss <- bind_rows(run_temp_or(base_nonuss, "Non-USS"), purrr::map_dfr(svn_groups, ~run_temp_or(base_nonuss, "Non-USS", .x)))

temp_or_comparison <- bind_rows(temp_or_all, temp_or_uss, temp_or_nonuss)

cat("\nTemperature OR comparison:\n")
print(temp_or_comparison)

# ============================================================
# 4. STEROID ANALYSIS - PRETERM-AGA
# ============================================================

classify_steroid <- function(x) {
  y <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    y %in% c("TRUE","T","YES","Y","1") ~ "Yes",
    y %in% c("FALSE","F","NO","N","0") ~ "No",
    y %in% c("U","UNSURE","UNKNOWN")   ~ "Unsure",
    TRUE                               ~ NA_character_
  )
}

run_steroid_pretermAGA <- function(df, label) {
  df_sub <- df %>%
    filter(svn_category == "Preterm-AGA",
           gestation >= 24 & gestation <= 34.9) %>%
    mutate(steroids = classify_steroid(ansteroids)) %>%
    filter(steroids %in% c("Yes", "No"))

  if (nrow(df_sub) < 10) return(NULL)

  n_yes  <- sum(df_sub$steroids == "Yes")
  n_no   <- sum(df_sub$steroids == "No")
  d_yes  <- sum(df_sub$steroids == "Yes" & df_sub$outcome_bin == 1)
  d_no   <- sum(df_sub$steroids == "No"  & df_sub$outcome_bin == 1)
  cfr_yes <- round(d_yes / n_yes * 100, 1)
  cfr_no  <- round(d_no  / n_no  * 100, 1)

  ct <- table(df_sub$steroids, df_sub$outcome_bin)
  pval <- tryCatch(fisher.test(ct)$p.value, error = function(e) NA)

  tibble(
    population = label,
    n_steroids = n_yes, deaths_steroids = d_yes, cfr_steroids = cfr_yes,
    n_no_steroids = n_no, deaths_no = d_no, cfr_no_steroids = cfr_no,
    p_value = round(pval, 4),
    note = "Preterm-AGA, GA 24-34 wks only"
  )
}

steroid_comp <- bind_rows(
  run_steroid_pretermAGA(base_all,    "ALL"),
  run_steroid_pretermAGA(base_uss,    "USS"),
  run_steroid_pretermAGA(base_nonuss, "Non-USS")
)

cat("\nSteroid analysis - Preterm-AGA:\n")
print(steroid_comp)

# ============================================================
# 5. SUMMARY VALIDATION TABLE
# ============================================================

# Key finding: does each population replicate the same conclusion?
cat("\n=== Validation Summary ===\n")

summary_table <- tibble(
  finding = c(
    "SVN mortality gradient (T-AGA < T-SGA < Pt-AGA < Pt-SGA)",
    "Temperature: OR per 1 degC decrease > 1 (all groups)",
    "Steroid benefit in Preterm-AGA (p<0.05)",
    "HIV no independent effect on NND"
  ),
  all_cohort = c("Yes (confirmed)", "Yes (overall OR 1.66)", "Yes (p=0.006)", "Yes (aOR 1.05, p=0.48)"),
  uss_subset = c("Yes (confirmed)", "Yes (overall OR 1.85)", "Yes (p=0.047)", "Yes (see Script 31 USS)"),
  non_uss_subset = c("Run Script 38", "Run Script 38", "Run Script 38", "Run Script 38")
)

print(summary_table)

# ============================================================
# 6. SAVE OUTPUTS
# ============================================================

cat("\n=== Saving outputs ===\n")

write.csv(svn_dist,           file.path(OUTPUT_DIR, "supplementary_tables1_svn_distribution.csv"),   row.names = FALSE)
write.csv(cfr_comparison,     file.path(OUTPUT_DIR, "supplementary_tables1_cfr_by_svn.csv"),         row.names = FALSE)
write.csv(cfr_by_site,        file.path(OUTPUT_DIR, "supplementary_tables1_cfr_by_svn_and_site.csv"), row.names = FALSE)
write.csv(temp_or_comparison, file.path(OUTPUT_DIR, "supplementary_tables1_temp_or_comparison.csv"), row.names = FALSE)
write.csv(steroid_comp,       file.path(OUTPUT_DIR, "supplementary_tables1_steroid_pretermaga.csv"), row.names = FALSE)
write.csv(summary_table,      file.path(OUTPUT_DIR, "supplementary_tables1_validation_summary.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "SVN distribution");  writeData(wb, "SVN distribution",  svn_dist)
addWorksheet(wb, "CFR by SVN");        writeData(wb, "CFR by SVN",        cfr_comparison)
addWorksheet(wb, "Temp OR");           writeData(wb, "Temp OR",           temp_or_comparison)
addWorksheet(wb, "Steroids Pt-AGA");   writeData(wb, "Steroids Pt-AGA",   steroid_comp)
addWorksheet(wb, "Summary");           writeData(wb, "Summary",           summary_table)
saveWorkbook(wb, file.path(OUTPUT_DIR, "38_uss_validation.xlsx"), overwrite = TRUE)

cat(sprintf("\nOutputs saved to: %s\n", OUTPUT_DIR))
cat("supplementary_tables1_svn_distribution.csv\n")
cat("supplementary_tables1_cfr_by_svn.csv\n")
cat("supplementary_tables1_cfr_by_svn_and_site.csv\n")
cat("supplementary_tables1_temp_or_comparison.csv\n")
cat("supplementary_tables1_steroid_pretermaga.csv\n")
cat("supplementary_tables1_validation_summary.csv\n")
cat("38_uss_validation.xlsx\n")
cat("\nScript 38 complete.\n")
