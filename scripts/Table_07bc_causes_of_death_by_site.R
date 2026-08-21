# ==============================================================================
# Table_08_causes_of_death_by_svn.R
# Manuscript TABLE 8 - Leading causes of neonatal death by SVN phenotype.
#
# Produces the five collapsed cause groups (as used in the manuscript) as
# n (% of deaths in that column) for each SVN phenotype, plus pooled All-SVN
# and Overall columns, with a Total deaths row:
#
#   Rows    : Prematurity / RDS
#             Birth asphyxia / HIE
#             Neonatal sepsis / infection
#             Congenital anomalies (incl. gastroschisis)
#             Other / unspecified
#             Total deaths
#   Columns : Term-AGA | Term-SGA | Preterm-AGA | Preterm-SGA | All-SVN | Overall
#
# Population = svn_filtered (current 30,449 cohort; deaths = outcome_bin == 1;
# total deaths = 5,348), identical to Table 5 / script 45, so the death totals
# reconcile with the rest of the manuscript.
#
# The five-group mapping is taken verbatim from SCRIPTS V3/17_causes_of_death.R
# (map_cause_to_category), extended only to catch three free-text prematurity
# variants and one anomaly variant that pass through normalise_cause as raw
# strings (PREM, "PREMATURITY WITH RD", "PREMATURE (32-36 WEEKS)",
# MYELOMENINGOCELE); these additions are flagged inline below.
#
# Outputs (outputs/):
#   table_07bc_causes_by_svn.csv         - wide grid, cells "n (x.x%)"
#   table_07bc_causes_by_svn_and_site.csv- same grid split KCH / SMCH
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
suppressPackageStartupMessages({library(dplyr); library(stringr); library(tidyr)})
cat("\n=== Table_08_causes_of_death_by_svn.R ===\n")

# ------------------------------------------------------------------------------
# FIVE-GROUP CAUSE MAPPING (verbatim from 17_causes_of_death.R::map_cause_to_category,
# operating on the RAW causedeath code). Additions to catch free-text variants
# are marked with  ## +added .
# ------------------------------------------------------------------------------
map_cause5 <- function(code) {
    u <- toupper(trimws(code))
    dplyr::case_when(
        u %in% c(
            "PREMATURITY/RDS",
            "PRRDS", "PR", "RDNT", "TTN", "LBW", "RDS",
            "PREMRD", "TERMRD", "EXPREM", "VPREM", "EXLBW", "VLBW", "CLEFTRD",
            "PREM", "PREMATURITY WITH RD", "PREMATURE (32-36 WEEKS)"  ## +added free-text prematurity
        ) ~ "Prematurity / RDS",

        u %in% c(
            "BIRTH ASPHYXIA/HIE",
            "BA", "HIE", "ASP", "SHIE", "MA", "LAS", "AP", "MAS", "BBA"
        ) ~ "Birth asphyxia / HIE",

        u %in% c(
            "NEONATAL SEPSIS",
            "NEONATAL SEPSIS (EARLY ONSET)",
            "NEONATAL SEPSIS (LATE ONSET)",
            "NS", "EONS", "LONS", "SEPS", "HYPOGSY",
            "MEN", "PN", "NT", "NEOT", "BRON", "GA", "GE"
        ) ~ "Neonatal sepsis / infection",

        u %in% c(
            "CONGENITAL ANOMALY", "GASTROSCHISIS", "OMPHALOCELE",
            "CONG", "CHD", "BI", "BDN", "HYDR",
            "BO", "EXOM", "G", "IA", "NEC", "J",
            "CA", "GASTRO", "GSCH", "GSCHIS",
            "OMPH", "OM", "OCA", "MD", "MYELO", "CLEFTLIPPALATE",
            "MYELOMENINGOCELE"  ## +added free-text anomaly
        ) ~ "Congenital anomalies (incl. gastroschisis)",

        TRUE ~ "Other / unspecified"
    )
}

CAUSE_GROUPS <- c(
    "Prematurity / RDS",
    "Birth asphyxia / HIE",
    "Neonatal sepsis / infection",
    "Congenital anomalies (incl. gastroschisis)",
    "Other / unspecified"
)
PH <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")

# ------------------------------------------------------------------------------
# LOAD DEATHS ON THE CURRENT COHORT
# ------------------------------------------------------------------------------
nnu <- process_nnu_data(load_nnu_data())
if (!"causedeath" %in% names(nnu$svn_filtered))
    stop("Column 'causedeath' not found in svn_filtered.")

deaths <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), outcome_bin == 1, svn_category %in% PH) %>%
    mutate(
        facility   = as.character(facility),
        phenotype  = factor(as.character(svn_category), levels = PH),
        cod        = str_trim(as.character(causedeath)),
        cod        = ifelse(cod == "" | is.na(cod) | cod == "NA", "Unknown", cod),
        cause_grp  = factor(map_cause5(cod), levels = CAUSE_GROUPS)
    )

cat(sprintf("Total deaths on cohort: %d (Term-AGA %d; All-SVN %d)\n",
            nrow(deaths), sum(deaths$phenotype == "Term-AGA"),
            sum(deaths$phenotype != "Term-AGA")))

# ------------------------------------------------------------------------------
# BUILD ONE COLUMN OF n (%) FOR A GIVEN SUBSET
# ------------------------------------------------------------------------------
col_np <- function(df) {
    tot <- nrow(df)
    counts <- df %>% count(cause_grp, .drop = FALSE) %>%
        mutate(pct = if (tot > 0) 100 * n / tot else rep(0, length(n)),
               np  = sprintf("%d (%.1f%%)", n, pct))
    out <- setNames(counts$np, as.character(counts$cause_grp))[CAUSE_GROUPS]
    c(out, "Total deaths" = sprintf("%d (100%%)", tot))
}

make_grid <- function(df) {
    data.frame(
        Cause         = c(CAUSE_GROUPS, "Total deaths"),
        `Term-AGA`    = col_np(df %>% filter(phenotype == "Term-AGA")),
        `Term-SGA`    = col_np(df %>% filter(phenotype == "Term-SGA")),
        `Preterm-AGA` = col_np(df %>% filter(phenotype == "Preterm-AGA")),
        `Preterm-SGA` = col_np(df %>% filter(phenotype == "Preterm-SGA")),
        `All-SVN`     = col_np(df %>% filter(phenotype != "Term-AGA")),
        `Overall`     = col_np(df),
        check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE
    )
}

grid <- make_grid(deaths)
write.csv(grid, file.path(OUT2, "table_07bc_causes_by_svn.csv"), row.names = FALSE)
cat("Saved: table_07bc_causes_by_svn.csv\n")
print(grid, row.names = FALSE)

# ------------------------------------------------------------------------------
# SAME GRID FOR EACH HOSPITAL (Table 8a = KCH, Table 8b = SMCH), same columns as
# the combined table, so all three tables reconcile on the current cohort.
# ------------------------------------------------------------------------------
grid_kch  <- make_grid(deaths %>% filter(facility == "KCH"))
grid_smch <- make_grid(deaths %>% filter(facility == "SMCH"))

write.csv(grid_kch,  file.path(OUT2, "table_07bc_causes_by_svn_kch.csv"),  row.names = FALSE)
write.csv(grid_smch, file.path(OUT2, "table_07bc_causes_by_svn_smch.csv"), row.names = FALSE)
cat("Saved: table_07bc_causes_by_svn_kch.csv\n")
cat("Saved: table_07bc_causes_by_svn_smch.csv\n")

# Combined stacked version (both sites in one file) retained for convenience
grid_site <- bind_rows(
    grid_kch  %>% mutate(Site = "KCH",  .before = 1),
    grid_smch %>% mutate(Site = "SMCH", .before = 1)
)
write.csv(grid_site, file.path(OUT2, "table_07bc_causes_by_svn_and_site.csv"),
          row.names = FALSE)
cat("Saved: table_07bc_causes_by_svn_and_site.csv\n")

cat("\n--- KCH ---\n");  print(grid_kch,  row.names = FALSE)
cat("\n--- SMCH ---\n"); print(grid_smch, row.names = FALSE)
cat("Table_08_causes_of_death_by_svn.R complete.\n")
