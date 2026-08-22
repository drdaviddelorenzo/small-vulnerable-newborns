# ==============================================================================
# Purpose: Table 11 -- case fatality with vs without antenatal corticosteroids
# among eligible preterm neonates, restricted to the clinically eligible
# window (24+0 to 34+6 weeks gestation).
#
# Definitions: the `ansteroids` variable is grouped into Steroids Used (Y) /
# Not Used (N) / Unsure (U) / Missing. Reports CFR (deaths/n) for "Steroids
# used" vs "not used", by facility and by LBW/SGA category, with Fisher's
# exact test p-values.
#
# Output: table_11_full_*.csv, table_11_ga24to34wk_*.csv, table_11_cfr_by_steroids.csv
# ==============================================================================

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")

# Steroids Groups
STEROIDS_GROUPS <- c(
    "Steroids Used (Y)",
    "Steroids Not Used (N)",
    "Unsure (U)",
    "Missing"
)

# Normalize outcomes
normalise_outcome <- function(x) {
    y <- tolower(trimws(as.character(x)))
    y <- gsub("[^a-z0-9]+", " ", y)
    y <- gsub("\\s+", " ", y)

    dplyr::case_when(
        grepl("\\bdc\\b|discharged|\\bdcr\\b|\\bdpc\\b|\\bdckmc\\b|kangaroo mother care", y) ~ "Alive",
        grepl("\\babs\\b|abscond|\\bdama\\b|against medical advice", y) ~ "Alive",
        grepl("\\btrh\\b|\\btro\\b|transfer", y) ~ "Alive",
        grepl("\\bnnd\\b|neonatal death|\\bdda\\b|died during admission", y) ~ "Neonatal Death",
        grepl("\\bbid\\b|brought in dead", y) ~ "Neonatal Death",
        grepl("\\bstb\\b|stillbirth|still born", y) ~ "Stillbirth",
        y == "" ~ NA_character_,
        TRUE ~ "Unknown/Other"
    )
}

# Normalize Antenatal Steroids
# Values: Y, N, U, or Missing
normalise_steroids <- function(x) {
    y <- toupper(trimws(as.character(x)))
    dplyr::case_when(
        y == "Y" ~ "Steroids Used (Y)",
        y == "N" ~ "Steroids Not Used (N)",
        y == "U" ~ "Unsure (U)",
        is.na(y) | y == "" | y == "NA" ~ "Missing",
        TRUE ~ "Missing" # Treat unexpected values as missing or check them
    )
}

# ===== LOAD AND PROCESS DATA =====

cat("\n=== Loading NNU Data ===\n")
all_admissions <- load_nnu_data()

# Extract antenatal steroids column
all_admissions <- all_admissions %>%
    mutate(steroids_raw = ansteroids)

# Process
cat("Processing variables...\n")
analysis_df <- all_admissions %>%
    mutate(
        bw_kg = bw_kg_canonical(.),
        gestation_weeks = suppressWarnings(as.numeric(gestation)),

        # Outcome
        outcome_norm = normalise_outcome(neotreeoutcome),

        # Steroids Categorization
        steroids_cat = normalise_steroids(steroids_raw),

        # Validations
        is_valid_bw = is.finite(bw_kg) & bw_kg >= BW_MIN_KG & bw_kg <= BW_MAX_KG,
        is_invalid_ga = is.na(gestation_weeks) | gestation_weeks < GA_MIN | gestation_weeks > GA_MAX,
        sga_status = dplyr::if_else(gestation_weeks >= GA_MIN & !is.na(gestation_weeks), classify_sga(bw_kg, gestation_weeks, if ("gender" %in% names(.)) gender else NA_character_), NA_character_),
        sga_cat_temp = sga_category_from_gestational_age(sga_status, gestation_weeks),
        is_valid_sga = is_valid_bw & !is.na(sga_cat_temp),

        # Categorizations including "Discarded"
        bw_category = if_else(is_valid_bw, bw_category_from_kg(bw_kg), "Discarded"),
        lbw_explanation = if_else(is_valid_bw, "Included in analysis", "Missing or invalid birth weight"),
        sga_explanation = case_when(
            is_valid_sga ~ "Included in analysis",
            !is_valid_bw & is_invalid_ga ~ "Missing both BW and GA",
            !is_valid_bw ~ "Missing BW",
            is_invalid_ga ~ "Missing GA",
            TRUE ~ "Missing GA"
        ),
        sga_category = if_else(is_valid_sga, sga_cat_temp, sga_explanation)
    )

# Filter for valid data (Valid Outcome) - Birthweight mappings managed securely within explicit matrices natively gracefully natively accurately explicitly
valid_df <- analysis_df %>%
    filter(outcome_norm %in% c("Alive", "Neonatal Death"))

cat("Total Analysis Population (Valid Outcome):", nrow(valid_df), "\n")

# ===== ANALYSIS FUNCTION =====

run_q21_analysis <- function(valid_df, analysis_df, file_prefix, label) {
    cat(sprintf("\n=== Running Analysis: %s ===\n", label))
    cat("Total Analysis Population (Valid Outcome):", nrow(valid_df), "\n")

    # ----- CALC STATISTICS -----
    cat("Calculating Statistics...\n")

    # Helper for group stats
    calc_group_stats <- function(data, group_name, filter_expr, group_cols, group_levels) {
        category_col <- group_cols[1]
        stats <- data %>%
            filter(eval(parse(text = filter_expr))) %>%
            filter(!is.na(.data[[category_col]])) %>%
            group_by(facility, category = .data[[category_col]]) %>%
            summarise(
                NND = sum(outcome_norm == "Neonatal Death"),
                DC = sum(outcome_norm == "Alive"),
                .groups = "drop"
            ) %>%
            mutate(steroids_group = group_name)

        totals <- data %>%
            filter(eval(parse(text = filter_expr))) %>%
            group_by(facility) %>%
            summarise(
                NND = sum(outcome_norm == "Neonatal Death"),
                DC = sum(outcome_norm == "Alive"),
                .groups = "drop"
            ) %>%
            mutate(category = "Total", steroids_group = group_name)
        bind_rows(stats, totals)
    }

    # Helper for wide outcome
    calc_wide_outcome <- function(data, group_cols, group_levels) {
        defs <- list(
            list(name = "Steroids Used (Y)", expr = "steroids_cat == 'Steroids Used (Y)'"),
            list(name = "Steroids Not Used (N)", expr = "steroids_cat == 'Steroids Not Used (N)'"),
            list(name = "Unsure (U)", expr = "steroids_cat == 'Unsure (U)'"),
            list(name = "Missing", expr = "steroids_cat == 'Missing'")
        )
        all_stats <- list()
        for (d in defs) {
            all_stats[[d$name]] <- calc_group_stats(data, d$name, d$expr, group_cols, group_levels)
        }
        summary_long <- bind_rows(all_stats)
        complete_grid <- expand_grid(
            facility = unique(data$facility),
            category = c(group_levels, "Discarded", "Missing BW", "Missing GA", "Missing both BW and GA", "Total"),
            steroids_group = STEROIDS_GROUPS
        )
        summary_complete <- complete_grid %>%
            left_join(summary_long, by = c("facility", "category", "steroids_group")) %>%
            replace_na(list(NND = 0, DC = 0))

        summary_wide <- summary_complete %>%
            pivot_wider(
                names_from = category,
                values_from = c(NND, DC),
                names_glue = "{category}_{.value}"
            )

        col_order <- c("facility", "steroids_group")
        for (cat in c(group_levels, "Discarded", "Missing BW", "Missing GA", "Missing both BW and GA", "Total")) {
            col_name_nnd <- paste0(cat, "_NND")
            col_name_dc <- paste0(cat, "_DC")
            if (col_name_nnd %in% names(summary_wide)) {
                col_order <- c(col_order, col_name_nnd, col_name_dc)
            }
        }
        summary_wide %>% select(any_of(col_order)) %>% arrange(facility, match(steroids_group, STEROIDS_GROUPS))
    }

    # LBW Categories
    lbw_levels <- c("ELBW", "VLBW", "LBW", "NBW", "HBW")
    lbw_wide <- calc_wide_outcome(valid_df, c("bw_category", "lbw_explanation"), lbw_levels) %>%
        mutate(classification = "Birthweight")

    # SGA Categories
    sga_levels <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")
    sga_wide <- calc_wide_outcome(valid_df, c("sga_category", "sga_explanation"), sga_levels) %>%
        mutate(classification = "SGA")

    # ----- SAVE RESULTS -----
    cat("Saving Results...\n")

    bw_note <- sprintf("is.finite(bw_kg) & bw_kg >= %.1f kg & bw_kg <= %.1f kg", BW_MIN_KG, BW_MAX_KG)
    ga_note <- sprintf("!is.na(gestation) & gestation >= %d wks & gestation <= %d wks", GA_MIN, GA_MAX)
    sga_note <- paste0(bw_note, " AND ", ga_note, " AND !is.na(sga_category)")
    outcome_note <- "normalise_outcome(neotreeoutcome) %in% c(\"Alive\", \"Neonatal Death\")"
    
    n_total <- nrow(analysis_df)
    sample_size_summary <- tibble(
        description = c(
            "Total admissions",
            "Analysis population (valid outcome)",
            "Excluded: missing/invalid outcome",
            "Analysis population - Steroids Used (Y)",
            "Analysis population - Steroids Not Used (N)",
            "Analysis population - Unsure (U)",
            "Analysis population - Missing steroids record",
            "Analysis population (KCH) - Steroids Used (Y)",
            "Analysis population (KCH) - Steroids Not Used (N)",
            "Analysis population (KCH) - Unsure (U)",
            "Analysis population (KCH) - Missing steroids record",
            "Analysis population (SMCH) - Steroids Used (Y)",
            "Analysis population (SMCH) - Steroids Not Used (N)",
            "Analysis population (SMCH) - Unsure (U)",
            "Analysis population (SMCH) - Missing steroids record",
            "Analysis population with valid birthweight (LBW sample)",
            "Discarded for LBW (missing/invalid birthweight)",
            "Analysis population with valid SGA classification",
            "Discarded for SGA (Missing BW)",
            "Discarded for SGA (Missing GA)",
            "Discarded for SGA (Missing both BW and GA)"
        ),
        n = c(
            n_total,
            nrow(valid_df),
            sum(!analysis_df$outcome_norm %in% c("Alive", "Neonatal Death")),
            sum(valid_df$steroids_cat == "Steroids Used (Y)"),
            sum(valid_df$steroids_cat == "Steroids Not Used (N)"),
            sum(valid_df$steroids_cat == "Unsure (U)"),
            sum(valid_df$steroids_cat == "Missing"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Steroids Used (Y)"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Steroids Not Used (N)"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Unsure (U)"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Missing"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Steroids Used (Y)"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Steroids Not Used (N)"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Unsure (U)"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Missing"),
            sum(valid_df$is_valid_bw),
            sum(!valid_df$is_valid_bw),
            sum(valid_df$is_valid_sga),
            sum(valid_df$sga_explanation == "Missing BW"),
            sum(valid_df$sga_explanation == "Missing GA"),
            sum(valid_df$sga_explanation == "Missing both BW and GA")
        ),
        pct_of_total = round(100 * n / n_total, 1),
        Filter_notes = c(
            NA_character_, outcome_note, paste0("NOT (", outcome_note, ")"),
            paste0(outcome_note, " AND steroids == 'Y'"), paste0(outcome_note, " AND steroids == 'N'"),
            paste0(outcome_note, " AND steroids == 'U'"), paste0(outcome_note, " AND steroids is NA/blank/other"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids == 'Y'"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids == 'N'"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids == 'U'"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids is Missing"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids == 'Y'"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids == 'N'"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids == 'U'"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids is Missing"),
            paste0(outcome_note, " AND ", bw_note),
            paste0(outcome_note, " AND NOT (", bw_note, ")"),
            paste0(outcome_note, " AND ", sga_note),
            "Invalid/missing birth weight", "Invalid/missing gestational age", "Both BW and GA invalid/missing"
        )
    )

    write_csv(sample_size_summary, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_00_sample_size_summary.csv")))
    write_csv(lbw_wide, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_01_steroids_outcomes_LBW_wide.csv")))
    write_csv(sga_wide, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_02_steroids_outcomes_SGA_wide.csv")))

    # ----- CFR TESTS -----
    cat("Running CFR Chi-Square and Pairwise Tests...\n")
    
    calc_cfr_tests <- function(data, category_col, valid_categories, group_col, group_levels, group_label) {
        d <- data %>% filter(.data[[group_col]] %in% group_levels, .data[[category_col]] %in% valid_categories)
        facilities <- unique(d$facility)
        results <- list()
        for (fac in facilities) {
            for (cat in valid_categories) {
                sub <- d %>% filter(facility == fac, .data[[category_col]] == cat)
                ct_df <- sub %>% group_by(group = .data[[group_col]]) %>%
                    summarise(NND = sum(outcome_norm == "Neonatal Death"), DC  = sum(outcome_norm == "Alive"), .groups = "drop") %>%
                    filter((NND + DC) > 0)
                if (nrow(ct_df) < 2) next
                ct_mat <- as.matrix(ct_df[, c("NND", "DC")]); rownames(ct_mat) <- ct_df$group
                cfr_vals <- ct_df %>% mutate(total = NND + DC, CFR = round(100 * NND / total, 1))
                chisq_res <- tryCatch(suppressWarnings(chisq.test(ct_mat, correct = FALSE)), error = function(e) NULL)
                overall_p <- if (!is.null(chisq_res)) chisq_res$p.value else NA_real_
                overall_stat <- if (!is.null(chisq_res)) chisq_res$statistic else NA_real_
                overall_df   <- if (!is.null(chisq_res)) chisq_res$parameter else NA_real_
                expected_warn <- if (!is.null(chisq_res)) any(chisq_res$expected < 5) else NA
                grp_names <- rownames(ct_mat); pairs <- combn(grp_names, 2, simplify = FALSE); n_pairs <- length(pairs)
                pairwise_rows <- map_dfr(pairs, function(pair) {
                    g1 <- pair[1]; g2 <- pair[2]; sub_mat <- ct_mat[c(g1, g2), , drop = FALSE]
                    fish_p <- tryCatch(fisher.test(sub_mat)$p.value, error = function(e) NA_real_)
                    tibble(facility = fac, category = cat, comparison = paste0(g1, " vs. ", g2), group1 = g1, group2 = g2,
                           CFR_group1_pct = cfr_vals$CFR[cfr_vals$group == g1], n_group1 = cfr_vals$total[cfr_vals$group == g1], NND_group1 = cfr_vals$NND[cfr_vals$group == g1],
                           CFR_group2_pct = cfr_vals$CFR[cfr_vals$group == g2], n_group2 = cfr_vals$total[cfr_vals$group == g2], NND_group2 = cfr_vals$NND[cfr_vals$group == g2],
                           fisher_p_raw = fish_p, n_pairs_tested = n_pairs)
                })
                pairwise_rows <- pairwise_rows %>% mutate(
                    fisher_p_bonferroni = pmin(fisher_p_raw * n_pairs_tested, 1),
                    fisher_p_bonferroni = round(fisher_p_bonferroni, 6),
                    fisher_p_raw = round(fisher_p_raw, 6),
                    chisq_statistic = round(as.numeric(overall_stat), 4),
                    chisq_df = as.integer(overall_df),
                    chisq_p_value = round(overall_p, 6),
                    low_expected_count_warning = expected_warn,
                    exposure_variable = group_label,
                    classification = if (category_col == "bw_category") "Birthweight" else "SGA"
                )
                results[[paste(fac, cat, sep = "_")]] <- pairwise_rows
            }
        }
        bind_rows(results)
    }

    testable_steroids <- c("Steroids Used (Y)", "Steroids Not Used (N)", "Unsure (U)")
    cfr_tests_lbw <- calc_cfr_tests(valid_df, "bw_category", lbw_levels, "steroids_cat", testable_steroids, "Antenatal Steroids")
    cfr_tests_sga <- calc_cfr_tests(valid_df, "sga_category", sga_levels, "steroids_cat", testable_steroids, "Antenatal Steroids")
    cfr_tests_all <- bind_rows(cfr_tests_lbw, cfr_tests_sga) %>% arrange(facility, classification, category, comparison)

    write_csv(cfr_tests_all, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_03_steroids_CFR_chisquare_pairwise.csv")))
    cat(sprintf("Saved: 03_steroids_CFR_chisquare_pairwise.csv (Rows: %d)\n", nrow(cfr_tests_all)))
    invisible(cfr_tests_all)
}

# ===== EXECUTE ANALYSES =====

# 1. Full Sample Analysis
run_q21_analysis(valid_df, analysis_df, "table_11_full", "Full Sample")

# 2. Gestation Subsample Analysis (24 to 34 weeks)
cat("\n=== Creating Gestation Subsample (24 to 34 weeks) ===\n")
sub_analysis_df <- analysis_df %>%
    filter(gestation_weeks >= 24 & gestation_weeks < 35)

sub_valid_df <- sub_analysis_df %>%
    filter(outcome_norm %in% c("Alive", "Neonatal Death"))

res_sub <- run_q21_analysis(sub_valid_df, sub_analysis_df, "table_11_ga24to34wk",
                            "Gestation Subsample (24 to 34 weeks)")

# ===== MANUSCRIPT TABLE 11 =====
# CFR expressed as deaths per 1,000 admissions, with deaths/n shown, No vs Yes
# steroids, for the preterm SVN phenotypes in the eligible GA 24-34+6 wk window.
#   p (No vs Yes) = Fisher exact on the two shown groups (matches a 2-column table)
#   p (3-group)   = 3-group (N/Y/U) chi-square, reported for reference
fmtp    <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
comma   <- function(x) formatC(x, format = "d", big.mark = ",")
cfr1000 <- function(nnd, n) round(1000 * nnd / n)
cell    <- function(nnd, n) sprintf("%d (%s/%s)", cfr1000(nnd, n), comma(nnd), comma(n))

table13 <- res_sub %>%
    filter(classification == "SGA",
           comparison == "Steroids Not Used (N) vs. Steroids Used (Y)",
           category %in% c("Preterm-AGA", "Preterm-SGA")) %>%
    mutate(.fac = factor(facility, levels = c("KCH", "SMCH")),
           .cat = factor(category, levels = c("Preterm-AGA", "Preterm-SGA"))) %>%
    arrange(.fac, .cat) %>%
    transmute(
        Group = paste(facility, category),
        `No steroids, CFR per 1,000 (deaths/n)` = cell(NND_group1, n_group1),
        `Steroids, CFR per 1,000 (deaths/n)`    = cell(NND_group2, n_group2),
        `p (No vs Yes, Fisher)`  = fmtp(fisher_p_raw),
        `p (3-group chi-square)` = fmtp(chisq_p_value))

write_csv(table13, file.path(OUTPUT_DIR_NNU, "table_11_cfr_by_steroids.csv"))
cat("\n=== Manuscript Table 11: CFR (per 1,000) by antenatal steroids, eligible preterm ===\n")
print(as.data.frame(table13), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(OUTPUT_DIR_NNU, "table_11_cfr_by_steroids.csv")))

cat("\nAll analyses complete!\n")
