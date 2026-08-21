# 36_temperature_categorical_OR.R
# Categorical Temperature Analysis — Odds Ratios with Normal Range as Reference
#
# Purpose:
#   Extends Script 35 (continuous) and Script 18 (raw counts) by computing
#   formal odds ratios for each temperature category using 36.5–37.5°C as
#   the reference group (clinical normal range). Covers both hypothermic and
#   hyperthermic ends of the spectrum. Analysis by SVN group and overall.
#
# Background:
#   Raised at team meeting 9 June 2026 (Hannah Calder, Rachel Poffley,
#   David de Lorenzo). The existing continuous analysis (Script 35) models
#   temperature linearly; this script shows the categorical OR structure
#   to capture any U-shaped relationship, particularly at the hyperthermic
#   end (>37.5°C).
#
# Temperature categories (matching Script 18 intervals):
#   1. Severe hypothermia:  < 32.0°C
#   2. Hypothermia:         32.0 – 35.9°C
#   3. Low-normal:          36.0 – 36.4°C
#   4. Normal (REFERENCE):  36.5 – 37.5°C
#   5. Hyperthermia:        > 37.5°C
#
# Output files (in OUTPUT_DIR):
#   dep_table09_temp_cat_or_overall.csv    — ORs vs normal, overall
#   dep_table09_temp_cat_or_by_svn.csv     — ORs vs normal, by SVN group
#   dep_table09_temp_cat_or_by_lbw.csv     — ORs vs normal, by LBW category
#   dep_table09_temp_cat_sample_sizes.csv  — N per temperature category per group
#   36_temperature_categorical.xlsx  — Combined Excel workbook

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

cat("=== Script 36: Categorical Temperature OR Analysis ===\n\n")

# ============================================================
# 1. TEMPERATURE CLASSIFICATION
# ============================================================

TEMP_VALID_MIN <- 30.0
TEMP_VALID_MAX <- 45.0

classify_temp_cat <- function(temp) {
  dplyr::case_when(
    is.na(temp)                        ~ NA_character_,
    temp < TEMP_VALID_MIN              ~ NA_character_,
    temp > TEMP_VALID_MAX              ~ NA_character_,
    temp <  32.0                       ~ "1_SevereHypothermia",
    temp >= 32.0 & temp <  36.0       ~ "2_Hypothermia",
    temp >= 36.0 & temp <  36.5       ~ "3_LowNormal",
    temp >= 36.5 & temp <= 37.5       ~ "4_Normal",
    temp >  37.5                       ~ "5_Hyperthermia",
    TRUE                               ~ NA_character_
  )
}

TEMP_LABELS <- c(
  "1_SevereHypothermia" = "Severe hypothermia (<32.0°C)",
  "2_Hypothermia"       = "Hypothermia (32.0–35.9°C)",
  "3_LowNormal"         = "Low-normal (36.0–36.4°C)",
  "4_Normal"            = "Normal (36.5–37.5°C) [Reference]",
  "5_Hyperthermia"      = "Hyperthermia (>37.5°C)"
)

# ============================================================
# 2. LOAD AND PREPARE ANALYSIS DATASET
# ============================================================

cat("Loading NNU data...\n")
raw  <- load_nnu_data()
nnu  <- process_nnu_data(raw)
base <- nnu$svn_filtered %>%
  filter(!is.na(outcome_bin))

cat(sprintf("Base population: %d  (KCH: %d, SMCH: %d)\n",
            nrow(base),
            sum(base$facility == "KCH"),
            sum(base$facility == "SMCH")))

cat("Preparing temperature analysis dataset...\n")

analysis_df <- base %>%
  mutate(
    temp_val = suppressWarnings(
      as.numeric(temperature)
    ),
    temp_val = ifelse(temp_val >= TEMP_VALID_MIN & temp_val <= TEMP_VALID_MAX,
                      temp_val, NA_real_),
    temp_cat = classify_temp_cat(temp_val),
    temp_cat = factor(temp_cat,
                      levels = c("4_Normal",
                                 "1_SevereHypothermia",
                                 "2_Hypothermia",
                                 "3_LowNormal",
                                 "5_Hyperthermia"))
  ) %>%
  filter(
    !is.na(temp_cat),
    !is.na(outcome_bin)
  )

cat(sprintf("Analysis N with valid temperature and outcome: %d\n", nrow(analysis_df)))
cat(sprintf("Temperature category distribution:\n"))
print(table(analysis_df$temp_cat, useNA = "ifany"))

# ============================================================
# 3. HELPER: RUN CATEGORICAL LOGISTIC REGRESSION
# ============================================================

run_cat_logistic <- function(df, group_label) {
  if (nrow(df) < 10) {
    return(tibble(
      group      = group_label,
      temp_cat   = NA_character_,
      n_total    = nrow(df),
      n_events   = sum(df$outcome_bin, na.rm = TRUE),
      OR         = NA_real_,
      ci_lower   = NA_real_,
      ci_upper   = NA_real_,
      p_value    = NA_real_,
      note       = "Insufficient data"
    ))
  }

  # Check reference group has enough events
  ref_n <- sum(df$temp_cat == "4_Normal", na.rm = TRUE)
  if (ref_n < 5) {
    return(tibble(
      group      = group_label,
      temp_cat   = NA_character_,
      n_total    = nrow(df),
      n_events   = sum(df$outcome_bin, na.rm = TRUE),
      OR         = NA_real_,
      ci_lower   = NA_real_,
      ci_upper   = NA_real_,
      p_value    = NA_real_,
      note       = "Reference group too small"
    ))
  }

  tryCatch({
    model <- glm(outcome_bin ~ temp_cat,
                 data   = df,
                 family = binomial(link = "logit"))

    tidy_out <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        group    = group_label,
        temp_cat = sub("temp_cat", "", term),
        n_total  = nrow(df),
        n_events = sum(df$outcome_bin, na.rm = TRUE),
        OR       = round(estimate, 4),
        ci_lower = round(conf.low, 4),
        ci_upper = round(conf.high, 4),
        p_value  = round(p.value, 6),
        note     = paste0("OR vs Normal (36.5-37.5°C), n=", ref_n, " in reference group")
      ) %>%
      select(group, temp_cat, n_total, n_events, OR, ci_lower, ci_upper, p_value, note)

    # Add reference row
    ref_row <- tibble(
      group    = group_label,
      temp_cat = "4_Normal (Reference)",
      n_total  = nrow(df),
      n_events = sum(df$outcome_bin, na.rm = TRUE),
      OR       = 1.00,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value  = NA_real_,
      note     = paste0("Reference group, n=", ref_n)
    )

    bind_rows(ref_row, tidy_out)

  }, error = function(e) {
    tibble(
      group    = group_label,
      temp_cat = NA_character_,
      n_total  = nrow(df),
      n_events = sum(df$outcome_bin, na.rm = TRUE),
      OR       = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value  = NA_real_,
      note     = paste0("Model error: ", e$message)
    )
  })
}

# ============================================================
# 4. OVERALL ANALYSIS
# ============================================================

cat("\n=== Running overall analysis ===\n")

or_overall <- run_cat_logistic(analysis_df, "Overall")
cat("Overall results:\n")
print(or_overall)

# ============================================================
# 5. BY SVN GROUP
# ============================================================

cat("\n=== Running by SVN group ===\n")

svn_groups <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")

or_svn <- purrr::map_dfr(svn_groups, function(grp) {
  df_grp <- analysis_df %>% filter(svn_category == grp)
  cat(sprintf("  %s: n=%d\n", grp, nrow(df_grp)))
  run_cat_logistic(df_grp, grp)
})

# ============================================================
# 6. BY LBW CATEGORY
# ============================================================

cat("\n=== Running by LBW category ===\n")

lbw_groups <- c("ELBW", "VLBW", "LBW", "NBW", "HBW")

or_lbw <- purrr::map_dfr(lbw_groups, function(cat_name) {
  df_grp <- analysis_df %>% filter(bw_cat == cat_name)
  cat(sprintf("  %s: n=%d\n", cat_name, nrow(df_grp)))
  run_cat_logistic(df_grp, cat_name)
})

# ============================================================
# 7. SAMPLE SIZE TABLE
# ============================================================

cat("\n=== Generating sample size table ===\n")

sample_sizes <- analysis_df %>%
  group_by(svn_category, temp_cat, .drop = FALSE) %>%
  summarise(
    n         = n(),
    n_deaths  = sum(outcome_bin, na.rm = TRUE),
    cfr_pct   = round(n_deaths / n * 100, 1),
    .groups   = "drop"
  ) %>%
  mutate(
    temp_label = TEMP_LABELS[as.character(temp_cat)]
  )

# ============================================================
# 8. SAVE OUTPUTS
# ============================================================

cat("\n=== Saving outputs ===\n")

write.csv(or_overall, file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_overall.csv"), row.names = FALSE)
write.csv(or_svn,     file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_by_svn.csv"),  row.names = FALSE)
write.csv(or_lbw,     file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_by_lbw.csv"),  row.names = FALSE)
write.csv(sample_sizes, file.path(OUTPUT_DIR, "dep_table09_temp_cat_sample_sizes.csv"), row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Overall OR");     writeData(wb, "Overall OR",     or_overall)
addWorksheet(wb, "By SVN Group");   writeData(wb, "By SVN Group",   or_svn)
addWorksheet(wb, "By LBW Cat");     writeData(wb, "By LBW Cat",     or_lbw)
addWorksheet(wb, "Sample Sizes");   writeData(wb, "Sample Sizes",   sample_sizes)
saveWorkbook(wb, file.path(OUTPUT_DIR, "36_temperature_categorical.xlsx"), overwrite = TRUE)

cat(sprintf("\nOutputs saved to: %s\n", OUTPUT_DIR))
cat("dep_table09_temp_cat_or_overall.csv\n")
cat("dep_table09_temp_cat_or_by_svn.csv\n")
cat("dep_table09_temp_cat_or_by_lbw.csv\n")
cat("dep_table09_temp_cat_sample_sizes.csv\n")
cat("36_temperature_categorical.xlsx\n")
cat("\nScript 36 complete.\n")
