# ==============================================================================
# Purpose: Table 12 -- case fatality by infant sex within each SVN
# phenotype, both sites combined and by site. Reports % and deaths per 1,000
# admissions with Wilson 95% CI.
#
# Definitions: population = svn_filtered. Deaths use the shared outcome_bin
# (infants brought in dead counted as deaths). Respects POPULATION ("all"
# default, "USS" for the ultrasound subset).
#
# Output: table_12_cfr_by_sex.csv
# ==============================================================================
# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")
OUT2 <- OUTPUT_DIR

suppressPackageStartupMessages({library(dplyr)})
cat("\n=== 54_cfr_by_sex.R ===\n")
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")

wilson <- function(x, n){
    if (n == 0) return(c(NA, NA)); p <- x/n; z <- 1.96; d <- 1 + z^2/n
    lo <- (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d
    hi <- (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d
    round(c(1000*lo, 1000*hi), 1)
}

d <- process_nnu_data(load_nnu_data())$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(svn_category = factor(as.character(svn_category), levels = PH),
           facility = factor(as.character(facility), levels = c("KCH","SMCH")),
           sex_cat = case_when(
               toupper(trimws(as.character(gender_std))) %in% c("M","MALE","BOY")    ~ "Male",
               toupper(trimws(as.character(gender_std))) %in% c("F","FEMALE","GIRL") ~ "Female",
               TRUE ~ NA_character_),
           sex_cat = factor(sex_cat, levels = c("Male","Female"))) %>%
    filter(!is.na(sex_cat))

cfr_by_sex <- function(df, site){
    df %>% group_by(svn_category, sex_cat) %>%
        summarise(n = n(), deaths = sum(outcome_bin == 1), .groups = "drop") %>%
        rowwise() %>%
        mutate(site = site, cfr_pct = round(100*deaths/n, 1),
               per1000 = round(1000*deaths/n, 1),
               ci = list(wilson(deaths, n)),
               ci95_per1000 = sprintf("%.1f-%.1f", ci[1], ci[2])) %>%
        ungroup() %>%
        select(site, svn_category, sex_cat, n, deaths, cfr_pct, per1000, ci95_per1000)
}
res <- bind_rows(cfr_by_sex(d, "Combined"),
                 cfr_by_sex(filter(d, facility == "KCH"),  "KCH"),
                 cfr_by_sex(filter(d, facility == "SMCH"), "SMCH")) %>%
    arrange(match(site, c("Combined","KCH","SMCH")), svn_category, sex_cat)
write.csv(res, file.path(OUT2, "table_12_cfr_by_sex.csv"), row.names = FALSE)

cat("CFR by sex within SVN phenotype (Combined = Table 12):\n")
print(as.data.frame(filter(res, site == "Combined")), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(OUT2, "table_12_cfr_by_sex.csv")))
cat("54_cfr_by_sex.R complete.\n")
