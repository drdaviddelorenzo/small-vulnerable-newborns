# ==============================================================================
# run_all.R
# Orchestrates the full small-vulnerable-newborns (SVN) analysis pipeline.
#
# Usage:
#   Rscript run_all.R
#
# What it does:
#   1. Checks that input/ contains the six required (confidential) data files
#      -- see input/README.md -- and stops with a clear message if not.
#   2. Installs any R packages used by the pipeline that are not already
#      installed.
#   3. Runs every script in scripts/ once for POPULATION="all" and once for
#      POPULATION="USS", each as its own fresh R process, in dependency
#      order (a small number of Figure_*/helper_* scripts read a CSV written
#      by an earlier Table_*/dep_* script in the same population pass).
#   4. A script that errors is logged and skipped -- it does NOT stop the
#      rest of the run. A full pass/fail summary is printed at the end.
#
# Outputs land in output/all/{tables,graphs}/ and output/USS/{tables,graphs}/
# (script code writes to a lowercase output/uss/ folder for the USS pass --
# see 00a_setup_shared.R -- this is expected, not a bug).
# ==============================================================================

SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)

# ------------------------------------------------------------------------------
# 1. Check input/ for the required data files
# ------------------------------------------------------------------------------
REQUIRED_INPUT_FILES <- c(
    "mwi-master-data-20220101-to-20251231.csv",
    "mwi-matched-nnu-data-20220101-to-20251231.csv",
    "mwi-maternal-data-20220101-to-20251231.csv",
    "zim-master-data-20220101-to-20251231.csv",
    "zim-matched-nnu-data-20220101-to-20251231.csv",
    "zim-maternal-data-20220101-to-20251231.csv"
)
missing_inputs <- REQUIRED_INPUT_FILES[
    !file.exists(file.path("input", REQUIRED_INPUT_FILES))
]
if (length(missing_inputs) > 0) {
    cat("\n==================================================================\n")
    cat("Cannot run the pipeline -- input/ is missing required data file(s):\n\n")
    for (f in missing_inputs) cat("  -", f, "\n")
    cat("\nThese are confidential patient-level extracts and are NOT included\n")
    cat("in this repository. See input/README.md for what is expected and\n")
    cat("where to obtain them.\n")
    cat("==================================================================\n\n")
    quit(status = 1, save = "no")
}
cat("input/ check passed -- all 6 required data files found.\n\n")

# ------------------------------------------------------------------------------
# 2. Ensure required R packages are installed
# ------------------------------------------------------------------------------
REQUIRED_PACKAGES <- c(
    "dplyr", "readr", "tidyr", "ggplot2", "scales", "purrr",
    "openxlsx", "broom", "stringr", "lubridate"
)
have_pkg <- vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
missing_pkgs <- REQUIRED_PACKAGES[!have_pkg]
if (length(missing_pkgs) > 0) {
    cat("Installing missing R package(s):", paste(missing_pkgs, collapse = ", "), "\n")
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
    still_missing <- missing_pkgs[!vapply(missing_pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(still_missing) > 0) {
        cat("\nERROR: could not install: ", paste(still_missing, collapse = ", "),
            "\nInstall these manually and re-run.\n", sep = "")
        quit(status = 1, save = "no")
    }
}
cat("Package check passed.\n\n")

# ------------------------------------------------------------------------------
# 3. Script list, in dependency order
#    (Layer 0 = reads only raw data; Layer 1 = reads a Layer-0 script's CSV
#    output from the SAME population pass, so must run after it.)
# ------------------------------------------------------------------------------
SCRIPTS <- c(
    # --- Layer 0: independent (raw data only) ---
    "helper_Table01_SVN_by_site_and_stillbirth_OR.R",
    "Table_01a_maternal_characteristics.R",
    "Table_01b_stillbirth_rate_and_OR.R",
    "Figure_01_cohort_flowchart.R",
    "Figure_02_stillbirth_OR_by_SVN.R",
    "Figure_02_SVN_distribution_births_vs_admissions.R",
    "Table_02_NNU_characteristics.R",
    "Table_03_04_distribution_and_CFR.R",
    "dep_Figure03_CFR_comparison.R",
    "Table_05_timing_of_death.R",
    "Figure_05_timing_of_death_by_SVN.R",
    "Table_06_length_of_stay.R",
    "helper_Table03_06_kruskal_wallis_cfr_los.R",
    "Table_07a_causes_of_death.R",
    "Table_07bc_causes_of_death_by_site.R",
    "Table_08_univariable_predictors_AllSVN.R",
    "Table_09_temperature_categorical_ORs.R",
    "dep_Table09_temperature_categorical_OR.R",
    "helper_Table09_temperature_svn_interaction.R",
    "Table_10_cfr_by_resuscitation.R",
    "Table_11_CFR_by_antenatal_steroids.R",
    "Table_12_CFR_by_sex.R",
    "MV_regression_SVN_population.R",
    "MV_regression_SVN_population_continuous_temperature.R",
    "MV_regression_whole_cohort.R",
    "MV_regression_whole_cohort_continuous_temperature.R",
    "dep_Figure09_hiv_distribution_by_svn.R",
    "Appendix_Table_Sx_univariable_full_vs_USS_vs_nonUSS.R",
    "Supplementary_Table_S1_validation_CFR_full_vs_USS.R",
    # --- Layer 1: depends on a Layer-0 script's output CSV ---
    "Figure_03_CFR_by_SVN_and_site.R",          # <- dep_Figure03_CFR_comparison.R
    "Figure_06_univariable_ORs_AllSVN.R",       # <- Table_08_univariable_predictors_AllSVN.R
    "Figure_07_temperature_categorical_ORs_by_SVN.R", # <- dep_Table09_temperature_categorical_OR.R
    "Figure_09_SVN_distribution_by_HIV.R",      # <- dep_Figure09_hiv_distribution_by_svn.R
    "helper_Table09_temperature_AllSVN_column.R" # <- dep_Table09_temperature_categorical_OR.R
)

missing_scripts <- SCRIPTS[!file.exists(file.path("scripts", SCRIPTS))]
if (length(missing_scripts) > 0) {
    cat("ERROR: run_all.R references scripts that do not exist in scripts/:\n")
    for (s in missing_scripts) cat("  -", s, "\n")
    quit(status = 1, save = "no")
}

# ------------------------------------------------------------------------------
# 4. Run every script, for both populations, each as its own R process
# ------------------------------------------------------------------------------
results <- list()

for (population in c("all", "USS")) {
    cat(sprintf("\n##################################################################\n"))
    cat(sprintf("### POPULATION = %s\n", population))
    cat(sprintf("##################################################################\n\n"))

    for (script in SCRIPTS) {
        cat(sprintf("--> Running %s  [POPULATION=%s]...\n", script, population))
        script_path <- normalizePath(file.path("scripts", script))
        log_file <- tempfile(fileext = ".log")

        exit_status <- system2(
            "Rscript",
            args = shQuote(script_path),
            env = paste0("POPULATION=", population),
            stdout = log_file,
            stderr = log_file
        )

        ok <- identical(exit_status, 0L)
        results[[paste(population, script)]] <- list(
            population = population, script = script, ok = ok
        )

        if (ok) {
            cat("    OK\n")
        } else {
            cat(sprintf("    FAILED (exit status %s) -- log below:\n", exit_status))
            log_lines <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) character())
            cat(paste("   |", log_lines, collapse = "\n"), "\n")
        }
        unlink(log_file)
    }
}

# ------------------------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------------------------
cat("\n##################################################################\n")
cat("### RUN SUMMARY\n")
cat("##################################################################\n")

ok_flags <- vapply(results, function(r) r$ok, logical(1))
cat(sprintf("%d / %d script runs succeeded.\n", sum(ok_flags), length(ok_flags)))

if (any(!ok_flags)) {
    cat("\nFAILED runs:\n")
    for (r in results) {
        if (!r$ok) cat(sprintf("  - [%s] %s\n", r$population, r$script))
    }
    cat("\nSee the log output above for each failure's error message.\n")
} else {
    cat("\nAll scripts completed successfully for both populations.\n")
}
cat("\nOutputs are in output/all/ and output/uss/ (tables/ and graphs/ subfolders).\n")
