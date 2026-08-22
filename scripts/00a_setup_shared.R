# ==============================================================================
# Purpose: shared setup sourced by every script in this folder -- data
# loading, SVN/SGA classification, birthweight handling, and plotting theme.
#
# Definitions: POPULATION is "all" (default, every newborn with valid GA) or
# "USS" (only newborns with GA estimated by ultrasound scan); set it before
# sourcing this file, or leave unset for "all". Study period 1 January 2022
# to 31 December 2025; facilities KCH (Malawi) and SMCH (Zimbabwe).
#
# Expects, relative to this file's own folder:
#   ../input/   the six source CSVs (not committed to the repo)
#   ../output/  all/{tables,graphs}/ and uss/{tables,graphs}/
#
# Data files (see input/README.md for the full list):
#   NNU, matched only:   zim-matched-nnu-data-*.csv, mwi-matched-nnu-data-*.csv
#   NNU, all admissions: zim-master-data-*.csv, mwi-master-data-*.csv
#                        (used only by Figure_01_cohort_flowchart.R, for the
#                        pre-match denominator)
#   Maternal:            zim-maternal-data-*.csv, mwi-maternal-data-*.csv
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ggplot2)
    library(scales)
    library(forcats)
    library(purrr)
    library(openxlsx)
    library(broom)
    library(stringr)
    library(lubridate)
})

# ==============================================================================
# POPULATION CONFIG
# ==============================================================================
if (!exists("POPULATION")) {
    .pp <- Sys.getenv("POPULATION"); POPULATION <- if (nzchar(.pp)) .pp else "all"
}
stopifnot(POPULATION %in% c("all", "USS"))
POP_FOLDER <- if (POPULATION == "USS") "uss" else "all"   # output subfolder name

cat(sprintf("=== POPULATION: %s (outputs -> outputs/%s/) ===\n", POPULATION, POP_FOLDER))

# ==============================================================================
# PATHS (all relative to this script's own directory, so the repository can
# be moved or copied anywhere and everything still resolves correctly)
# ==============================================================================
# SCRIPT_DIR is set by the small header at the top of every script in this
# folder, BEFORE it sources this file. Fall back to getwd() if sourced some
# other way (e.g. interactively with setwd() already done).
SCRIPT_DIR <- if (exists("SCRIPT_DIR")) SCRIPT_DIR else getwd()

# Data folder: sibling of scripts/, i.e. small-vulnerable-newborns/input/
DATA_DIR <- if (exists("DATA_DIR")) DATA_DIR else
    file.path(SCRIPT_DIR, "..", "input")
if (!dir.exists(DATA_DIR))
    stop("Cannot find the data folder at ", DATA_DIR,
         ". Expected input/ next to scripts/ -- see input/README.md for the ",
         "exact files required.")

# --- Data source files ---
SMCH_NNU_FILE       <- file.path(DATA_DIR, "zim-matched-nnu-data-20220101-to-20251231.csv")
KCH_NNU_FILE        <- file.path(DATA_DIR, "mwi-matched-nnu-data-20220101-to-20251231.csv")
SMCH_MATERNAL_FILE  <- file.path(DATA_DIR, "zim-maternal-data-20220101-to-20251231.csv")
KCH_MATERNAL_FILE   <- file.path(DATA_DIR, "mwi-maternal-data-20220101-to-20251231.csv")
SMCH_NNU_MASTER_FILE <- file.path(DATA_DIR, "zim-master-data-20220101-to-20251231.csv")
KCH_NNU_MASTER_FILE  <- file.path(DATA_DIR, "mwi-master-data-20220101-to-20251231.csv")

for (.f in c(SMCH_NNU_FILE, KCH_NNU_FILE, SMCH_MATERNAL_FILE, KCH_MATERNAL_FILE)) {
    if (!file.exists(.f)) stop("Required data file not found: ", .f)
}
rm(.f)

# --- Output directories --------------------------------------------------
# output/<all|uss>/tables/ and output/<all|uss>/graphs/. OUTPUT_DIR_MATERNAL
# and OUTPUT_DIR_NNU are aliases for OUTPUT_DIR (some scripts use one name,
# some the other; all three resolve to the same folder).
OUTPUT_ROOT <- file.path(SCRIPT_DIR, "..", "output")
OUTPUT_DIR_POP      <- file.path(OUTPUT_ROOT, POP_FOLDER)
OUTPUT_DIR          <- file.path(OUTPUT_DIR_POP, "tables")
OUTPUT_DIR_MATERNAL <- OUTPUT_DIR
OUTPUT_DIR_NNU      <- OUTPUT_DIR
OUTPUT_DIR_GRAPHS   <- file.path(OUTPUT_DIR_POP, "graphs")
for (.d in c(OUTPUT_DIR, OUTPUT_DIR_GRAPHS)) {
    if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}
rm(.d)

# ==============================================================================
# STUDY PERIOD AND FACILITY FILTERS
# ==============================================================================
STUDY_START <- as.Date("2022-01-01")
STUDY_END   <- as.Date("2025-12-31")

VALID_FACILITIES <- c("KCH", "SMCH")

# ==============================================================================
# VALIDITY THRESHOLDS
# ==============================================================================
BW_MIN_KG  <- 0.3    # 300 g
BW_MAX_KG  <- 7.0    # 7000 g
GA_MIN_NNU <- 24     # minimum GA for NNU SGA analysis (weeks)
GA_MIN_MAT <- 22     # minimum GA for maternal analysis (weeks)
GA_MAX     <- 42

# ==============================================================================
# INTERGROWTH-21st 10TH PERCENTILE REFERENCE (grams)
# Weeks 14-23: Fetal EFW; Weeks 24-42: Neonatal birthweight
# ==============================================================================
PCT10_UNISEX <- c(
    "14"=76,  "15"=98,  "16"=125, "17"=158, "18"=198, "19"=244,
    "20"=298, "21"=359, "22"=429, "23"=508,
    "24"=560, "25"=640, "26"=740, "27"=840, "28"=940,  "29"=1060,
    "30"=1190,"31"=1340,"32"=1500,
    "33"=1580,"34"=1780,"35"=2000,"36"=2220,"37"=2440, "38"=2650,
    "39"=2850,"40"=3010,"41"=3150,"42"=3260
)
PCT10_BOYS <- c(
    "14"=76,  "15"=98,  "16"=125, "17"=158, "18"=198, "19"=244,
    "20"=298, "21"=359, "22"=429, "23"=508,
    "24"=580, "25"=660, "26"=760, "27"=870, "28"=980,  "29"=1110,
    "30"=1250,"31"=1410,"32"=1580,
    "33"=1600,"34"=1810,"35"=2030,"36"=2250,"37"=2480, "38"=2700,
    "39"=2900,"40"=3070,"41"=3210,"42"=3320
)
PCT10_GIRLS <- c(
    "14"=76,  "15"=98,  "16"=125, "17"=158, "18"=198, "19"=244,
    "20"=298, "21"=359, "22"=429, "23"=508,
    "24"=540, "25"=620, "26"=720, "27"=810, "28"=910,  "29"=1030,
    "30"=1150,"31"=1300,"32"=1460,
    "33"=1540,"34"=1750,"35"=1960,"36"=2180,"37"=2400, "38"=2610,
    "39"=2800,"40"=2960,"41"=3090,"42"=3190
)

# ==============================================================================
# CATEGORY LEVELS (for ordered factors)
# ==============================================================================
BW_LEVELS  <- c("ELBW","VLBW","LBW","NBW","HBW")
SVN_LEVELS <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
SVN_CATS   <- SVN_LEVELS   # alias used by graph scripts

GA_MIN <- GA_MIN_NNU

SVN_COLOURS <- c(
    "Term-AGA"    = "#2166AC",
    "Term-SGA"    = "#74ADD1",
    "Preterm-AGA" = "#D73027",
    "Preterm-SGA" = "#FDAE61"
)
SVN_COLOURS_LABELLED <- c(
    "Term AGA"    = "#2166AC",
    "Term SGA"    = "#74ADD1",
    "Preterm AGA" = "#D73027",
    "Preterm SGA" = "#FDAE61"
)
FACILITY_COLOURS <- c("KCH" = "#1B7837", "SMCH" = "#762A83")
FACILITY_LABELS  <- c(KCH = "KCH (Malawi)", SMCH = "SMCH (Zimbabwe)")

SVN_LABELS <- c(
    "Term-AGA"    = "Term AGA",
    "Term-SGA"    = "Term SGA",
    "Preterm-AGA" = "Preterm AGA",
    "Preterm-SGA" = "Preterm SGA"
)

svn_factor <- function(x) {
    factor(x, levels = SVN_LEVELS, labels = unname(SVN_LABELS))
}

# ==============================================================================
# HELPER FUNCTIONS - BIRTHWEIGHT
# ==============================================================================
# Canonical birthweight is produced UPSTREAM by the Neotree cleaning pipeline
# as a single column, birthweight_g, in GRAMS (coalesce of the legacy
# era/site columns birthweight, bwt, bwtdis, with a kg->g guard). Read
# birthweight_g directly; the fallback below is for pre-fix inputs only.
BW_CANONICAL_COL <- "birthweight_g"

to_kg <- function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    pos   <- x_num[is.finite(x_num) & x_num > 0]
    if (length(pos) > 0) {
        med <- median(pos, na.rm = TRUE)
        if (is.finite(med) && med > 20 && med <= 7000) return(x_num / 1000)
    }
    x_num
}

get_bw_column <- function(df) {
    if (BW_CANONICAL_COL %in% names(df)) return(as.character(df[[BW_CANONICAL_COL]]))
    result <- rep(NA_character_, nrow(df))
    for (col in c("bwt","birthweight","bwtdis")) {
        if (col %in% names(df))
            result <- ifelse(is.na(result) & !is.na(df[[col]]), df[[col]], result)
    }
    result
}

bw_kg_canonical <- function(df) {
    if (BW_CANONICAL_COL %in% names(df)) {
        return(suppressWarnings(as.numeric(df[[BW_CANONICAL_COL]])) / 1000)
    }
    to_kg(get_bw_column(df))
}

# ==============================================================================
# HELPER FUNCTIONS - SGA / SVN CLASSIFICATION
# ==============================================================================

get_10th_pct <- function(ga_weeks, gender_val = NA_character_) {
    if (is.na(ga_weeks)) return(NA_real_)
    ga  <- round(ga_weeks)
    gac <- as.character(ga)
    if (ga < 14 || ga > 42) return(NA_real_)
    g <- toupper(trimws(as.character(gender_val)))
    pct <- if (!is.na(g) && g %in% c("M","MALE","BOY"))   PCT10_BOYS[gac]
      else if (!is.na(g) && g %in% c("F","FEMALE","GIRL")) PCT10_GIRLS[gac]
      else PCT10_UNISEX[gac]
    if (is.na(pct)) NA_real_ else as.numeric(pct) / 1000
}

classify_sga <- function(bw_kg, ga_weeks, gender_val = NA_character_) {
    pct10 <- mapply(get_10th_pct, ga_weeks, gender_val, USE.NAMES = FALSE)
    ifelse(is.na(bw_kg) | is.na(pct10), NA_character_,
           ifelse(bw_kg < pct10, "SGA", "AGA"))
}

make_svn_category <- function(sga, ga_weeks) {
    dplyr::case_when(
        is.na(sga) | is.na(ga_weeks)          ~ NA_character_,
        ga_weeks >= 37 & sga == "SGA"         ~ "Term-SGA",
        ga_weeks >= 37 & sga == "AGA"         ~ "Term-AGA",
        ga_weeks <  37 & sga == "SGA"         ~ "Preterm-SGA",
        ga_weeks <  37 & sga == "AGA"         ~ "Preterm-AGA",
        TRUE                                   ~ NA_character_
    )
}

bw_category <- function(bw_kg) {
    dplyr::case_when(
        is.na(bw_kg)          ~ NA_character_,
        bw_kg <  1.000        ~ "ELBW",
        bw_kg <  1.500        ~ "VLBW",
        bw_kg <  2.500        ~ "LBW",
        bw_kg <= 4.000        ~ "NBW",
        TRUE                  ~ "HBW"
    )
}

bw_category_from_kg              <- bw_category
sga_category_from_gestational_age <- make_svn_category

# ==============================================================================
# HELPER FUNCTIONS - OUTCOMES
# ==============================================================================

normalise_outcome <- function(x) {
    y <- tolower(trimws(as.character(x)))
    y <- gsub("[^a-z0-9]+", " ", y)
    y <- trimws(gsub("\\s+", " ", y))
    dplyr::case_when(
        grepl("\\bdc\\b|discharged|\\bdcr\\b|\\bdpc\\b|\\bdckmc\\b|kangaroo mother care", y) ~ "Alive",
        grepl("\\babs\\b|abscond|\\bdama\\b|against medical advice", y) ~ "Alive",
        grepl("\\btrh\\b|\\btro\\b|transfer", y) ~ "Alive",
        y == "lb" ~ "Alive",
        grepl("\\bnnd\\b|neonatal death|\\bdda\\b|died during admission", y) ~ "Neonatal Death (in-facility)",
        grepl("\\bbid\\b|brought in dead", y) ~ "Neonatal Death (pre-facility)",
        grepl("\\bennd\\b|early neonatal death", y) ~ "Neonatal Death (in-facility)",
        grepl("\\bstb\\b|\\bsbf\\b|\\bsbm\\b|\\bstbf\\b|\\bstbm\\b|stillbirth|still born", y) ~ "Stillbirth",
        y == "" ~ NA_character_,
        TRUE    ~ "Unknown/Other"
    )
}

binary_outcome <- function(outcome_norm) {
    dplyr::case_when(
        outcome_norm %in% c("Neonatal Death (in-facility)",
                            "Neonatal Death (pre-facility)") ~ 1L,
        outcome_norm == "Alive"                              ~ 0L,
        TRUE                                                  ~ NA_integer_
    )
}

# ==============================================================================
# HELPER FUNCTION - LENGTH OF STAY
# ==============================================================================
MAX_LOS_DAYS <- 180

get_los_days <- function(df) {
    if (all(c("datetimeadmission","datetimedischarge") %in% names(df))) {
        los <- suppressWarnings(as.numeric(difftime(
            lubridate::parse_date_time(df$datetimedischarge,
                                       orders = c("Ymd HMS", "ymd HMS", "Ymd", "ymd"),
                                       quiet = TRUE),
            lubridate::parse_date_time(df$datetimeadmission,
                                       orders = c("Ymd HMS", "ymd HMS", "Ymd", "ymd"),
                                       quiet = TRUE),
            units = "days"
        )))
        los <- ifelse(is.finite(los) & los >= 0 & los <= MAX_LOS_DAYS, los, NA_real_)
        return(los)
    }
    if ("timespent" %in% names(df)) {
        ts <- suppressWarnings(as.numeric(df$timespent))
        los <- ts / 24
        los <- ifelse(is.finite(los) & los >= 0 & los <= MAX_LOS_DAYS, los, NA_real_)
        return(los)
    }
    rep(NA_real_, nrow(df))
}

# ==============================================================================
# HELPER FUNCTION - DATE FILTERING
# ==============================================================================

filter_study_period <- function(df, date_col = "adm_date_parsed") {
    if (!date_col %in% names(df)) {
        fallbacks <- c("dateadmission", "datetimeadmission", "startedat")
        for (fb in fallbacks) {
            if (fb %in% names(df)) {
                date_col <- fb
                break
            }
        }
    }
    if (!date_col %in% names(df)) {
        warning("No date column found for study period filter. Returning all rows.")
        return(df)
    }

    df <- df %>%
        mutate(.adm_date = suppressWarnings(
            as.Date(lubridate::parse_date_time(
                .data[[date_col]],
                orders = c("Ymd HMS", "ymd HMS", "Ymd", "ymd"),
                quiet = TRUE
            ))
        ))

    n_before <- nrow(df)
    df <- df %>%
        filter(!is.na(.adm_date),
               .adm_date >= STUDY_START,
               .adm_date <= STUDY_END)
    n_after <- nrow(df)

    cat(sprintf("  Study period filter (%s to %s): %d -> %d rows (%d excluded)\n",
                STUDY_START, STUDY_END, n_before, n_after, n_before - n_after))
    df
}

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

load_nnu_data <- function(sites = c("KCH","SMCH")) {
    datasets <- list()

    if ("KCH" %in% sites) {
        cat("  Reading KCH NNU data...\n")
        d <- suppressMessages(read_csv(KCH_NNU_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        if (!"facility" %in% names(d)) d <- d %>% mutate(facility = "KCH")
        d <- d %>% filter(facility %in% "KCH")
        datasets[["KCH"]] <- d
    }

    if ("SMCH" %in% sites) {
        cat("  Reading SMCH NNU data...\n")
        d <- suppressMessages(read_csv(SMCH_NNU_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        d <- d %>% filter(facility %in% "SMCH")
        datasets[["SMCH"]] <- d
    }

    combined <- bind_rows(datasets)

    combined <- filter_study_period(combined, date_col = "adm_date_parsed")

    if (POPULATION == "USS") {
        n_before_uss <- nrow(combined)
        combined <- combined %>%
            filter(!is.na(methodestgest) &
                   toupper(trimws(methodestgest)) %in% c("USS", "EUSS"))
        cat(sprintf("  USS filter: %d -> %d rows\n",
                    n_before_uss, nrow(combined)))
    }

    cat(sprintf("  NNU data loaded: %d rows (KCH: %d, SMCH: %d)\n",
                nrow(combined),
                sum(combined$facility == "KCH"),
                sum(combined$facility == "SMCH")))

    combined
}

load_maternal_data <- function(sites = c("KCH","SMCH")) {
    datasets <- list()

    if ("KCH" %in% sites) {
        cat("  Reading KCH maternal data...\n")
        d <- suppressMessages(read_csv(KCH_MATERNAL_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        if (!"facility" %in% names(d)) d <- d %>% mutate(facility = "KCH")
        d <- d %>% filter(facility %in% "KCH")
        datasets[["KCH"]] <- d
    }

    if ("SMCH" %in% sites) {
        cat("  Reading SMCH maternal data...\n")
        d <- suppressMessages(read_csv(SMCH_MATERNAL_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        d <- d %>% filter(facility %in% "SMCH")
        datasets[["SMCH"]] <- d
    }

    combined <- bind_rows(datasets)

    combined <- filter_study_period(combined, date_col = "dateadmission")

    cat(sprintf("  Maternal data loaded: %d rows (KCH: %d, SMCH: %d)\n",
                nrow(combined),
                sum(combined$facility == "KCH"),
                sum(combined$facility == "SMCH")))

    combined
}

# ==============================================================================
# CORE PROCESSING - NNU DATA
# ==============================================================================

process_nnu_data <- function(raw) {
    gender_col <- if ("gender_std" %in% names(raw)) "gender_std"
                  else if ("gender" %in% names(raw)) "gender"
                  else NULL

    processed <- raw %>%
        mutate(
            bw_kg          = bw_kg_canonical(.),
            gestation_weeks= suppressWarnings(as.numeric(gestation)),
            gender_std     = if (!is.null(gender_col)) .data[[gender_col]] else NA_character_,
            outcome_norm   = normalise_outcome(neotreeoutcome),
            outcome_bin    = binary_outcome(outcome_norm)
        )

    lbw <- processed %>%
        filter(is.finite(bw_kg), bw_kg >= BW_MIN_KG, bw_kg <= BW_MAX_KG) %>%
        mutate(bw_cat = bw_category(bw_kg))

    svn <- lbw %>%
        filter(is.finite(gestation_weeks),
               gestation_weeks >= GA_MIN_NNU,
               gestation_weeks <= GA_MAX) %>%
        mutate(
            sga_status   = classify_sga(bw_kg, gestation_weeks, gender_std),
            svn_category = make_svn_category(sga_status, gestation_weeks)
        ) %>%
        filter(!is.na(svn_category))

    analysis <- svn %>%
        filter(!is.na(outcome_bin))

    list(raw = raw, processed = processed,
         lbw_filtered = lbw, svn_filtered = svn, analysis = analysis)
}

# ==============================================================================
# CORE PROCESSING - MATERNAL DATA
# ==============================================================================

process_maternal_data <- function(raw) {
    raw %>%
        mutate(
            bw_kg          = bw_kg_canonical(.),
            gestation_weeks= suppressWarnings(as.numeric(gestation)),
            outcome_norm   = normalise_outcome(neotreeoutcome),
            matageyrs_num  = suppressWarnings(as.numeric(matageyrs))
        )
}

# ==============================================================================
# GGPLOT2 THEME - publication-ready
# ==============================================================================
theme_svn <- function(base_size = 11, legend_position = "right") {
    theme_classic(base_size = base_size) %+replace%
        theme(
            axis.line   = element_line(colour = "grey30", linewidth = 0.5),
            axis.ticks  = element_line(colour = "grey50", linewidth = 0.4),
            axis.text   = element_text(colour = "grey20"),
            axis.title  = element_text(colour = "grey10", face = "bold"),
            panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.position  = legend_position,
            legend.title     = element_text(face = "bold", size = rel(0.9)),
            legend.text      = element_text(size = rel(0.85)),
            legend.key.size  = unit(0.45, "cm"),
            legend.background = element_rect(fill = "white", colour = "grey80"),
            strip.background = element_rect(fill = "grey94", colour = "grey70"),
            strip.text       = element_text(face = "bold", colour = "grey15"),
            plot.title       = element_text(face = "bold", hjust = 0,
                                            size = rel(1.15), colour = "grey10",
                                            margin = margin(b = 4)),
            plot.subtitle    = element_text(hjust = 0, colour = "grey30",
                                            size = rel(0.95), margin = margin(b = 6)),
            plot.caption     = element_text(hjust = 1, colour = "grey50",
                                            size = rel(0.75), margin = margin(t = 6)),
            plot.margin      = margin(10, 14, 8, 10)
        )
}

theme_set(theme_svn())

# ==============================================================================
# FIGURE SAVING
# ==============================================================================

save_figure <- function(plot, name, width = 10, height = 7, dpi = 300,
                        outdir = OUTPUT_DIR_GRAPHS) {
    pdf_file <- file.path(outdir, paste0(name, ".pdf"))
    png_file <- file.path(outdir, paste0(name, ".png"))
    ggsave(pdf_file, plot = plot, width = width, height = height, device = "pdf")
    ggsave(png_file, plot = plot, width = width, height = height, dpi = dpi, device = "png")
    cat(sprintf("  Saved: %s (PDF + PNG) [%.0f x %.0f in]\n", name, width, height))
}

# ==============================================================================
# WILSON 95% CI FOR A PROPORTION
# ==============================================================================

prop_ci_wilson <- function(x, n, conf = 0.95) {
    if (is.na(n) || n == 0) return(c(lower = NA_real_, upper = NA_real_))
    z <- qnorm(1 - (1 - conf) / 2)
    p <- x / n
    d <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / d
    marg <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
    c(lower = max(0, ctr - marg), upper = min(1, ctr + marg))
}

prop_ci <- function(x, n, conf = 0.95) {
    res <- prop_ci_wilson(x, n, conf)
    list(lower = unname(res["lower"]), upper = unname(res["upper"]))
}

add_prop_ci <- function(df, n_col = "deaths", total_col = "total_admissions",
                        rate_col = "mortality_rate") {
    df %>%
        rowwise() %>%
        mutate(
            ci_res   = list(prop_ci_wilson(.data[[n_col]], .data[[total_col]])),
            ci_lower = ci_res[["lower"]] * 100,
            ci_upper = ci_res[["upper"]] * 100
        ) %>%
        ungroup() %>%
        select(-ci_res)
}

# ==============================================================================
# STANDARD CAPTIONS
# ==============================================================================

svn_caption <- function(extra = NULL) {
    base <- paste0(
        "KCH = Kamuzu Central Hospital, Malawi; ",
        "SMCH = Sally Mugabe Children's Hospital, Zimbabwe. 2022-2025.")
    pop_note <- if (POPULATION == "USS") " [USS-confirmed GA subset]" else ""
    full <- paste0(base, pop_note)
    if (!is.null(extra)) paste0(full, "\n", extra) else full
}

# Suppress default Rplots.pdf generation
pdf(file = NULL)
options(device = function() { pdf(file = NULL) })

cat(sprintf("00a_setup_shared.R loaded. Data: %s | Output: %s\n", DATA_DIR, OUTPUT_DIR_POP))
