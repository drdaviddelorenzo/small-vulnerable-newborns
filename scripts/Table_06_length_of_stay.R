# ==============================================================================
# 51_los_by_category.R  (Table 7 and 7b)
# Length of stay (days) among SURVIVORS by SVN phenotype (Table 7) and by
# birth-weight category (Table 7b), by site and combined: median (IQR), mean,
# and n. Ported to SCRIPTS v4 to run on the CURRENT dataset.
#
# LOS variable = lengthofstay (days); stays > 180 days excluded as implausible.
# Survivors = outcome_bin == 0. Population = svn_filtered (matches scripts 42/48).
# Respects POPULATION ("all" default; set to "USS" in 51u for the USS subset).
#
# Output (outputs/ or outputs/USS/): table_06_length_of_stay_by_svn.csv, table_06_length_of_stay_by_bw.csv
# ASCII only (DSH).
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

suppressPackageStartupMessages({library(dplyr)})
cat("\n=== 51_los_by_category.R ===\n")
MAX_LOS_DAYS <- 180
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")

nnu <- process_nnu_data(load_nnu_data())
surv <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), outcome_bin == 0, svn_category %in% PH) %>%
    mutate(facility = factor(as.character(facility), levels=c("KCH","SMCH")),
           svn_category = factor(as.character(svn_category), levels=PH),
           bw_kg = bw_kg_canonical(.),
           bw_category = factor(ifelse(is.finite(bw_kg)&bw_kg>=BW_MIN_KG&bw_kg<=BW_MAX_KG,
                                       bw_category_from_kg(bw_kg), NA_character_),
                                levels=c("ELBW","VLBW","LBW","NBW","HBW")),
           # LOS is not a stored column in the current joined data; compute it
           # from timestamps (datetimedischarge - datetimeadmission), as in the
           # V3 pipeline (13_length_of_stay.R).
           adm_ts = suppressWarnings(lubridate::parse_date_time(datetimeadmission, orders=c("Ymd HMS","Ymd HM","Ymd"))),
           dis_ts = suppressWarnings(lubridate::parse_date_time(datetimedischarge, orders=c("Ymd HMS","Ymd HM","Ymd"))),
           los = as.numeric(difftime(dis_ts, adm_ts, units="days")),
           los = ifelse(is.finite(los) & los >= 0 & los <= MAX_LOS_DAYS, los, NA_real_))

los_by <- function(df, gcol){
    f <- function(d, site){
        d %>% group_by(category=.data[[gcol]]) %>%
            summarise(n=sum(!is.na(los)),
                      median_los=round(median(los,na.rm=TRUE),1),
                      iqr=sprintf("%.1f-%.1f", quantile(los,.25,na.rm=TRUE), quantile(los,.75,na.rm=TRUE)),
                      mean_los=round(mean(los,na.rm=TRUE),1), .groups="drop") %>%
            mutate(site=site)
    }
    bind_rows(f(filter(df,facility=="KCH"),"KCH"),
              f(filter(df,facility=="SMCH"),"SMCH"),
              f(df,"Combined")) %>%
        transmute(site, category, n, median_los, iqr, mean_los) %>%
        filter(!is.na(category)) %>% arrange(category, site)
}
svn_tab <- los_by(surv, "svn_category")
bw_tab  <- los_by(surv, "bw_category")
write.csv(svn_tab, file.path(OUT2,"table_06_length_of_stay_by_svn.csv"), row.names=FALSE)
write.csv(bw_tab,  file.path(OUT2,"table_06_length_of_stay_by_bw.csv"),  row.names=FALSE)

cat(sprintf("Survivors with valid LOS: %d\n", sum(!is.na(surv$los))))
cat("Median LOS (days) by SVN phenotype:\n"); print(as.data.frame(svn_tab), row.names=FALSE)
cat("\nMedian LOS (days) by birth weight:\n"); print(as.data.frame(bw_tab), row.names=FALSE)
cat("\nSaved table_06_length_of_stay_by_svn.csv, table_06_length_of_stay_by_bw.csv\n")
cat("51_los_by_category.R complete.\n")
