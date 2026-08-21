# ==============================================================================
# fig2_svn_births_vs_nnu_grouped.R  (SCRIPTS v4)  -- SELF-CONTAINED (data + graph)
# Figure 2 (GROUPED): SVN phenotype distribution across THREE series
#   (1) Births (Zimbabwe/SMCH)      - all SMCH births with valid BW+GA
#   (2) NNU admissions (SMCH only)  - like-for-like within-site comparison
#   (3) NNU admissions (both sites) - KCH + SMCH combined
# Requested (David, 9 Jul): "Births.." legend; n AND % on every bar; SMCH-only
# admissions series; the script must derive BOTH the data and the graph.
#
# Data are DERIVED here from the pipeline loaders (no dependence on intermediate
# CSVs), mirroring script 55 (maternal SVN) and script 53 (NNU svn_filtered):
#   Births  : process_maternal_data(load_maternal_data()), valid BW+GA, SMCH
#   NNU     : process_nnu_data(load_nnu_data())$svn_filtered, by facility
# Each series is % within its SVN-classifiable total (sums to 100%).
# ASCII only (DSH).  Run:  Rscript fig2_svn_births_vs_nnu_grouped.R
# ==============================================================================
# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(scales)})
cat("\n=== fig2_svn_births_vs_nnu_grouped.R (self-contained) ===\n")
OUT2G <- OUTPUT_DIR_GRAPHS; dir.create(OUT2G, recursive = TRUE, showWarnings = FALSE)
PH <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")
G_BIRTHS <- "Births (SMCH only)"
G_SMCH   <- "NNU admissions (SMCH only)"
G_BOTH   <- "NNU admissions (both sites)"
GA_MIN_M <- if (exists("GA_MIN_MAT")) GA_MIN_MAT else 22

# ---- (1) DERIVE births SVN distribution (SMCH, all births, valid BW+GA) -------
mat <- process_maternal_data(load_maternal_data()) %>%
    filter(is.finite(bw_kg), bw_kg >= BW_MIN_KG, bw_kg <= BW_MAX_KG,
           is.finite(gestation_weeks), gestation_weeks >= GA_MIN_M, gestation_weeks <= GA_MAX,
           as.character(facility) == "SMCH") %>%
    mutate(svn_category = make_svn_category(classify_sga(bw_kg, gestation_weeks), gestation_weeks))
births <- mat %>% filter(svn_category %in% PH) %>%
    count(svn_category, name = "n") %>% mutate(group = G_BIRTHS)

# ---- (2/3) DERIVE NNU SVN distribution (svn_filtered), by site + combined -----
nnu <- process_nnu_data(load_nnu_data())$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(facility = as.character(facility), svn_category = as.character(svn_category))
nnu_smch <- nnu %>% filter(facility == "SMCH") %>%
    count(svn_category, name = "n") %>% mutate(group = G_SMCH)
nnu_both <- nnu %>% count(svn_category, name = "n") %>% mutate(group = G_BOTH)

# ---- assemble + percentages within each series -------------------------------
plot_df <- bind_rows(births, nnu_smch, nnu_both) %>%
    group_by(group) %>% mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
    mutate(svn_cat = factor(svn_category, levels = PH),
           group   = factor(group, levels = c(G_BIRTHS, G_SMCH, G_BOTH)),
           lab     = paste0(comma(n), "\n(", sprintf("%.1f", pct), "%)"))
write.csv(plot_df %>% select(group, svn_category, n, pct),
          file.path(OUT2G, "fig2_svn_births_vs_nnu_grouped_counts.csv"), row.names = FALSE)

# ---- graph -------------------------------------------------------------------
GROUP_COLS <- setNames(c("#8FBF8F", "#E0A458", "#C39BD3"), c(G_BIRTHS, G_SMCH, G_BOTH))
p <- ggplot(plot_df, aes(svn_cat, pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.78,
             colour = "grey30", linewidth = 0.3) +
    geom_text(aes(label = lab), position = position_dodge(width = 0.85),
              vjust = -0.25, size = 2.7, lineheight = 0.9, colour = "grey15") +
    scale_fill_manual(values = GROUP_COLS, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 65), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "SVN phenotype distribution: births vs neonatal-unit admissions",
         x = NULL, y = "% of newborns") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
ggsave(file.path(OUT2G, "fig2_svn_births_vs_nnu_grouped.png"), p, width = 9.5, height = 5.4, dpi = 300)
cat("  Saved:", file.path(OUT2G, "fig2_svn_births_vs_nnu_grouped.png"), "\n")
print(as.data.frame(plot_df %>% select(group, svn_category, n, pct)))
cat("fig2_svn_births_vs_nnu_grouped.R complete.\n")
