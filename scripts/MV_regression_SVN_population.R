# ==============================================================================
# 40_svn_multivariable_regression.R
# NEW ANALYSIS (July 2026): Independent predictors of death WITHIN the SVN
# population.
#
# Background (co-author decisions, 7 July 2026):
#   1a  Population is restricted to the SVN phenotypes only (Term-SGA,
#       Preterm-AGA, Preterm-SGA). Term-AGA (the non-vulnerable group) is
#       EXCLUDED. Birth-weight category is NOT used as a predictor (it is
#       collinear with the phenotype, which is defined by weight + gestation).
#       Gestation band is also dropped for the same reason - prematurity is
#       carried by the phenotype term.
#   2b  The birth-weight mortality gradient (ELBW vs NBW) is NOT modelled here.
#       It stays descriptive (existing Table 5 / Figure 3). See script 10/12.
#   3a  We fit a MULTIVARIABLE model (all predictors together) and also report
#       the univariable estimate for each predictor. The covariate list is
#       printed and saved so the manuscript can state exactly what the model is
#       adjusted for.
#   4a  The hyperthermia / Term-AGA finding is kept but sourced from a
#       WHOLE-COHORT analysis - see companion script 41.
#
# Vulnerability axis in this model: svn_cat3 (phenotype), reference = Term-SGA
# (the least-severe SVN phenotype; Term-AGA is not in the population).
#
# Outputs (in ../SCRIPTS2-JULY2026/outputs):
#   40_svn_population_counts.csv       - denominator + events (verify 17,956)
#   40_svn_univariable.csv             - univariable OR per predictor
#   40_svn_multivariable.csv           - multivariable model (overall)
#   40_svn_multivariable.xlsx          - univariable + multivariable + covariates
#   40_svn_multivariable_missingness.csv
#
# NOTE: ASCII characters only (DSH requirement).
# ==============================================================================

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

cat("\n=== MV_regression_SVN_population.R ===\n")

OUT2 <- OUTPUT_DIR   # alias so ported OUT2 references resolve to output/<pop>/tables

# ==============================================================================
# LOAD DATA AND RESTRICT TO THE SVN POPULATION
# ==============================================================================
cat("Loading NNU data...\n")
raw  <- load_nnu_data()
nnu  <- process_nnu_data(raw)

SVN3 <- c("Term-SGA", "Preterm-AGA", "Preterm-SGA")   # Term-AGA EXCLUDED

data <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin)) %>%
    filter(svn_category %in% SVN3) %>%
    mutate(
        svn_cat3 = factor(as.character(svn_category), levels = SVN3),
        facility = factor(as.character(facility), levels = c("SMCH", "KCH"))
    )

# Verify the denominator (expect total 17,956; KCH 6,761; SMCH 11,195)
counts <- data %>%
    group_by(facility, svn_cat3) %>%
    summarise(n = n(), deaths = sum(outcome_bin == 1), .groups = "drop")
totals <- data %>%
    summarise(facility = "TOTAL", svn_cat3 = "ALL SVN",
              n = n(), deaths = sum(outcome_bin == 1))
site_tot <- data %>%
    group_by(facility) %>%
    summarise(svn_cat3 = "ALL SVN", n = n(),
              deaths = sum(outcome_bin == 1), .groups = "drop")
count_tbl <- bind_rows(counts, site_tot, totals)
write.csv(count_tbl, file.path(OUT2, "40_svn_population_counts.csv"),
          row.names = FALSE)

cat(sprintf("SVN population: %d  (KCH %d, SMCH %d); deaths %d\n",
            nrow(data),
            sum(data$facility == "KCH"), sum(data$facility == "SMCH"),
            sum(data$outcome_bin == 1)))
cat("  (Expect total 17,956; KCH 6,761; SMCH 11,195 - verify against 40_svn_population_counts.csv)\n")

# ==============================================================================
# DERIVE RISK-FACTOR VARIABLES (Version A categorisations; identical coding to
# script 31 so results are comparable to the previous pipeline). Birth-weight
# and gestation-band predictors are intentionally NOT derived here.
# ==============================================================================
data <- data %>%
    mutate(
        # HIV
        hiv_status = case_when(
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y")  &
                toupper(trimws(as.character(hivtestresult))) == "R"  ~ "HIV Positive",
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y")  &
                toupper(trimws(as.character(hivtestresult))) == "NR" ~ "HIV Negative",
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y")  ~ "HIV Tested (unknown result)",
            toupper(trimws(as.character(mathivtest))) %in% c("FALSE","N") ~ "Not Tested",
            TRUE ~ NA_character_
        ),
        hiv_status = factor(hiv_status,
                            levels = c("HIV Negative","HIV Positive",
                                       "HIV Tested (unknown result)","Not Tested")),

        # Parity (Version A)
        par_num = suppressWarnings(as.numeric(par)),
        parity = factor(case_when(
            is.na(par_num) ~ NA_character_,
            par_num == 1   ~ "Primiparous (1st delivery)",
            par_num <= 4   ~ "Multiparous (2nd-4th)",
            TRUE           ~ "Grand Multipara (5th+)"
        ), levels = c("Primiparous (1st delivery)",
                      "Multiparous (2nd-4th)","Grand Multipara (5th+)")),

        # Maternal age (reference 20-24 years)
        matage_num = suppressWarnings(as.numeric(matageyrs)),
        matage_cat = case_when(
            is.na(matage_num)  ~ NA_character_,
            matage_num < 16    ~ "<16 years",
            matage_num < 20    ~ "16-19 years",
            matage_num < 25    ~ "20-24 years",
            matage_num < 35    ~ "25-34 years",
            TRUE               ~ ">=35 years"
        ),
        matage_cat = factor(matage_cat,
                            levels = c("20-24 years","<16 years","16-19 years",
                                       "25-34 years",">=35 years")),

        # ANC (Version A: WHO threshold 4 visits)
        anc_num = suppressWarnings(as.numeric(antenatalcare)),
        anc = factor(case_when(
            is.na(anc_num) ~ NA_character_,
            anc_num == 0   ~ "None (0 visits)",
            anc_num <= 3   ~ "Inadequate (1-3 visits)",
            TRUE           ~ "Adequate (>=4 visits)"
        ), levels = c("Adequate (>=4 visits)","None (0 visits)","Inadequate (1-3 visits)")),

        # Type of pregnancy
        typebirth_cat = case_when(
            is.na(typebirth) | trimws(typebirth) == "" ~ NA_character_,
            toupper(trimws(typebirth)) %in% c("S","SINGLETON") ~ "Singleton",
            toupper(trimws(typebirth)) %in% c("TW1","TW2","TWIN1","TWIN2","TWIN") ~ "Twin",
            toupper(trimws(typebirth)) %in% c("TR1","TR2","TR3","TRIPLET") ~ "Triplet",
            TRUE ~ NA_character_
        ),
        typebirth_cat = factor(typebirth_cat, levels = c("Singleton","Twin","Triplet")),

        # Preeclampsia (Version A: PRECLAM only)
        preeclampsia = factor(
            ifelse(is.na(pregconditions), NA,
                   ifelse(grepl("PRECLAM", toupper(as.character(pregconditions))),
                          "Preeclampsia/Eclampsia", "No Preeclampsia")),
            levels = c("No Preeclampsia","Preeclampsia/Eclampsia")),

        # Mode of delivery (Version A)
        modedelivery_num = suppressWarnings(as.numeric(modedelivery)),
        mode_del = factor(case_when(
            is.na(modedelivery_num) ~ NA_character_,
            modedelivery_num == 1   ~ "SVD",
            modedelivery_num %in% c(2,3,6,7) ~ "Assisted Vaginal",
            modedelivery_num == 4   ~ "ELCS",
            modedelivery_num %in% c(5,8) ~ "EMCS",
            TRUE ~ NA_character_
        ), levels = c("SVD","Assisted Vaginal","ELCS","EMCS")),

        # Resuscitation
        resus_clean = toupper(trimws(as.character(resus))),
        resus_cat = case_when(
            is.na(resus) | resus_clean %in% c("","NA","NONE","NO","FALSE") ~ "None",
            grepl("CPR|CARDIO", resus_clean)       ~ "CPR",
            grepl("BVM|BAG.*VALVE|MASK", resus_clean) ~ "BVM",
            grepl("STIM|SUCT|O2|OXYGEN", resus_clean) ~ "Stimulation/Suction/O2",
            TRUE ~ "Resuscitation (type unknown)"
        ),
        resus_cat = factor(resus_cat,
                           levels = c("None","Stimulation/Suction/O2","BVM","CPR",
                                      "Resuscitation (type unknown)")),

        # Inborn / Outborn
        inout_cat = case_when(
            !is.na(inorout) ~ case_when(
                toupper(trimws(as.character(inorout))) %in% c("TRUE","INBORN","1","IN","YES") ~ "Inborn",
                toupper(trimws(as.character(inorout))) %in% c("FALSE","OUTBORN","0","OUT","NO") ~ "Outborn",
                TRUE ~ NA_character_),
            !is.na(admittedfrom) ~ case_when(
                toupper(trimws(as.character(admittedfrom))) %in% c("LW","LABOUR WARD","LW ") ~ "Inborn",
                TRUE ~ "Outborn"),
            TRUE ~ NA_character_
        ),
        inout_cat = factor(inout_cat, levels = c("Inborn","Outborn")),

        # Apgar 5 min (Version A: 3 groups)
        apgar5_num = suppressWarnings(as.numeric(apgar5)),
        apgar5_num = ifelse(apgar5_num >= 0 & apgar5_num <= 10, apgar5_num, NA_real_),
        apgar5_cat = factor(case_when(
            is.na(apgar5_num)  ~ NA_character_,
            apgar5_num <= 3    ~ "Severely depressed (0-3)",
            apgar5_num <= 6    ~ "Moderately depressed (4-6)",
            TRUE               ~ "Normal (7-10)"
        ), levels = c("Normal (7-10)","Moderately depressed (4-6)","Severely depressed (0-3)")),

        # Sex
        sex_cat = case_when(
            toupper(trimws(as.character(gender_std))) %in% c("M","MALE","BOY")    ~ "Male",
            toupper(trimws(as.character(gender_std))) %in% c("F","FEMALE","GIRL") ~ "Female",
            TRUE ~ NA_character_
        ),
        sex_cat = factor(sex_cat, levels = c("Male","Female")),

        # Temperature (categorical; reference Normal 36.5-37.5 degC)
        temp_val = suppressWarnings(as.numeric(temperature)),
        temp_val = ifelse(temp_val >= 30 & temp_val <= 45, temp_val, NA_real_),
        temp_cat = case_when(
            is.na(temp_val)    ~ NA_character_,
            temp_val <  32.0   ~ "Severe Hypothermia (<32 degC)",
            temp_val <  36.0   ~ "Moderate Hypothermia (32-35.9 degC)",
            temp_val <  36.5   ~ "Mild Hypothermia (36.0-36.4 degC)",
            temp_val <= 37.5   ~ "Normal (36.5-37.5 degC)",
            TRUE               ~ "Hyperthermia (>37.5 degC)"
        ),
        temp_cat = factor(temp_cat,
                          levels = c("Normal (36.5-37.5 degC)",
                                     "Severe Hypothermia (<32 degC)",
                                     "Moderate Hypothermia (32-35.9 degC)",
                                     "Mild Hypothermia (36.0-36.4 degC)",
                                     "Hyperthermia (>37.5 degC)"))
    )

# ==============================================================================
# COVARIATE SET (printed and saved so the manuscript can state adjustments)
# svn_cat3 = SVN phenotype (ref Term-SGA); facility adjusts for site.
# NO birth-weight and NO gestation-band term (carried by phenotype).
# ==============================================================================
# Pre-eclampsia is EXCLUDED: within the SVN population it is degenerate
# (univariable returns "Insufficient data" - too few informative cases once the
# non-vulnerable Term-AGA are removed), so it cannot be estimated.
selected_vars <- c("svn_cat3", "facility", "hiv_status", "parity", "matage_cat",
                    "anc", "typebirth_cat", "mode_del",
                    "resus_cat", "inout_cat", "apgar5_cat", "sex_cat", "temp_cat")

# "Not recorded" is kept as an explicit level ONLY for variables with
# substantive missingness. For near-complete variables a "Not recorded" level
# is a tiny, perfectly separated cell that yields nonsensical odds ratios
# (e.g. HIV 6 missing -> OR ~19000). Those few rows are handled by complete-case
# instead. Threshold: keep the level when >1% is missing.
NOTREC_VARS <- c("parity", "matage_cat", "anc", "apgar5_cat")

cat("\nCovariates in the multivariable model (adjustment set):\n")
cat(paste0("  - ", selected_vars, collapse = "\n"), "\n")
writeLines(selected_vars, file.path(OUT2, "40_covariate_list.txt"))

# ==============================================================================
# UNIVARIABLE ESTIMATES (one predictor at a time)
# ==============================================================================
run_univar <- function(df, predictor_var, outcome_var = "outcome_bin") {
    df2 <- df %>%
        select(all_of(c(outcome_var, predictor_var))) %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[predictor_var]]))
    if (nrow(df2) < 10 || length(unique(df2[[outcome_var]])) < 2 ||
        !is.factor(df2[[predictor_var]]) ||
        length(levels(droplevels(df2[[predictor_var]]))) < 2) {
        return(tibble(predictor = predictor_var, term = NA_character_,
                      n_total = nrow(df2), n_events = sum(df2[[outcome_var]] == 1),
                      OR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                      p_value = NA_real_, note = "Insufficient data"))
    }
    fmla <- as.formula(paste(outcome_var, "~", predictor_var))
    fit  <- glm(fmla, data = df2, family = binomial(link = "logit"))
    res  <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    ref_level <- levels(df2[[predictor_var]])[1]
    ref_row <- tibble(term = paste0(predictor_var, ref_level),
                      estimate = 1.0, conf.low = NA_real_,
                      conf.high = NA_real_, p.value = NA_real_)
    res %>%
        filter(term != "(Intercept)") %>%
        bind_rows(ref_row) %>%
        mutate(predictor = predictor_var, n_total = nrow(df2),
               n_events = sum(df2[[outcome_var]] == 1),
               OR = round(estimate, 3), ci_lower = round(conf.low, 3),
               ci_upper = round(conf.high, 3), p_value = round(p.value, 4),
               note = ifelse(term == paste0(predictor_var, ref_level), "Reference", "")) %>%
        select(predictor, term, n_total, n_events, OR, ci_lower, ci_upper, p_value, note)
}

cat("\nRunning univariable regressions...\n")
uni <- bind_rows(lapply(selected_vars, function(v) run_univar(data, v)))
write.csv(uni, file.path(OUT2, "40_svn_univariable.csv"), row.names = FALSE)
cat("  Saved: 40_svn_univariable.csv\n")

# ==============================================================================
# MULTIVARIABLE MODEL (all predictors together; "Not recorded" kept as a level)
# ==============================================================================
miss_rows <- lapply(selected_vars, function(v) {
    n_miss <- sum(is.na(data[[v]]))
    data.frame(variable = v, n_total = nrow(data), n_missing = n_miss,
               pct_missing = round(n_miss / nrow(data) * 100, 1),
               stringsAsFactors = FALSE)
})
miss_summary <- do.call(rbind, miss_rows)
write.csv(miss_summary, file.path(OUT2, "40_svn_multivariable_missingness.csv"),
          row.names = FALSE)

add_not_recorded <- function(x) {
    if (is.factor(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "Not recorded" }
    x
}
# Only the high-missingness variables get a "Not recorded" level; the remaining
# selected variables keep their NAs and are handled by complete-case in the
# model fit (run_multivar drops rows with any remaining NA).
data <- data %>% mutate(across(all_of(NOTREC_VARS), add_not_recorded))

n_cc <- sum(complete.cases(data[, selected_vars, drop = FALSE]))
cat(sprintf("Analysis sample after missing-data handling: %d / %d (%.1f%%)\n",
            n_cc, nrow(data), n_cc / nrow(data) * 100))
cat("  ('Not recorded' kept for: ", paste(NOTREC_VARS, collapse = ", "),
    "; other variables complete-case.)\n", sep = "")

run_multivar <- function(df, vars, outcome_var = "outcome_bin", label = "Overall") {
    df2 <- df %>%
        select(all_of(c(outcome_var, vars))) %>%
        filter(!is.na(.data[[outcome_var]])) %>%
        tidyr::drop_na()
    cat(sprintf("  %s: %d complete cases, %d events\n",
                label, nrow(df2), sum(df2[[outcome_var]] == 1)))
    if (nrow(df2) < 20 || length(unique(df2[[outcome_var]])) < 2) return(NULL)

    active_vars <- vars[vapply(vars, function(v) {
        col <- df2[[v]]
        if (is.factor(col)) length(levels(droplevels(col))) >= 2
        else length(unique(col[!is.na(col)])) >= 2
    }, logical(1))]
    dropped <- setdiff(vars, active_vars)
    if (length(dropped) > 0)
        cat(sprintf("  Dropped (single level): %s\n", paste(dropped, collapse = ", ")))
    if (length(active_vars) == 0) return(NULL)

    fmla <- as.formula(paste(outcome_var, "~", paste(active_vars, collapse = " + ")))
    fit  <- glm(fmla, data = df2, family = binomial(link = "logit"))
    # Profiled CIs preferred; Wald fallback if profiling fails (separation in a
    # near-empty 'Not recorded' cell).
    tidied <- tryCatch(
        broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE),
        error = function(e) {
            cat(sprintf("  Profiled CI failed (%s); using Wald CIs.\n", e$message))
            co <- summary(fit)$coefficients
            ci <- suppressMessages(confint.default(fit))
            tibble::tibble(term = rownames(co), estimate = exp(co[, "Estimate"]),
                           std.error = co[, "Std. Error"], statistic = co[, "z value"],
                           p.value = co[, "Pr(>|z|)"],
                           conf.low = exp(ci[, 1]), conf.high = exp(ci[, 2]))
        })
    tidied %>%
        filter(term != "(Intercept)") %>%
        mutate(model_group = label, n_total = nrow(df2),
               n_events = sum(df2[[outcome_var]] == 1),
               OR = round(estimate, 3), ci_lower = round(conf.low, 3),
               ci_upper = round(conf.high, 3), p_value = round(p.value, 4),
               OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, ci_lower, ci_upper)) %>%
        select(model_group, term, n_total, n_events,
               OR, ci_lower, ci_upper, p_value, OR_CI)
}

cat("\nFitting multivariable model (SVN population, overall)...\n")
m1 <- run_multivar(data, selected_vars, label = "SVN population (overall)")

if (!is.null(m1)) {
    write.csv(m1, file.path(OUT2, "40_svn_multivariable.csv"), row.names = FALSE)
    cat("  Saved: 40_svn_multivariable.csv\n")

    wb <- createWorkbook()
    addWorksheet(wb, "Multivariable")
    writeData(wb, "Multivariable", m1)
    addWorksheet(wb, "Univariable")
    writeData(wb, "Univariable", uni)
    addWorksheet(wb, "Covariates")
    writeData(wb, "Covariates", data.frame(covariate = selected_vars))
    addWorksheet(wb, "Missingness")
    writeData(wb, "Missingness", miss_summary)
    addWorksheet(wb, "Population_counts")
    writeData(wb, "Population_counts", count_tbl)
    saveWorkbook(wb, file.path(OUT2, "40_svn_multivariable.xlsx"), overwrite = TRUE)
    cat("  Saved: 40_svn_multivariable.xlsx\n")

    cat("\nKey multivariable ORs (SVN population):\n")
    print(m1 %>% select(term, OR_CI, p_value))
}

cat("\n40_svn_multivariable_regression.R complete.\n")
