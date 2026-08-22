# ==============================================================================
# Purpose: forest plot of univariable odds ratios for neonatal death within
# the pooled small-vulnerable-newborn (All-SVN) population, plus a companion
# comparison plot against Term-AGA (for the appendix).
#
# Definitions: restricted to the All-SVN group (Term-SGA + Preterm-AGA +
# Preterm-SGA); birth weight and gestation band are excluded (definitional
# components of the phenotype).
#
# Reads: table_08_univariable_svn.csv, table_08_univariable_termaga.csv
#   (written by Table_08_univariable_predictors_AllSVN.R -- run that first)
# Output: figure_06_univariable_ors_allsvn.png/.pdf,
#         figure_06_univariable_ors_svn_vs_termaga.png/.pdf
# ==============================================================================

# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")
OUT2 <- OUTPUT_DIR

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
})

cat("\n=== 43_fig4_svn_univariable_forest.R ===\n")

# OUT2 defined in the setup block above (this folder's 'outputs')

svn_path <- file.path(OUT2, "table_08_univariable_svn.csv")
ref_path <- file.path(OUT2, "table_08_univariable_termaga.csv")
if (!file.exists(svn_path))
    stop("Cannot find table_08_univariable_svn.csv - run Table_08_univariable_predictors_AllSVN.R first.")

# ==============================================================================
# LABELS: map predictor variable + level to readable text and group headers.
# The vector order here is the TOP-TO-BOTTOM order in the plot.
# ==============================================================================
PRED_LABEL <- c(
    facility      = "Site",
    hiv_status    = "Maternal HIV status",
    parity        = "Parity",
    matage_cat    = "Maternal age",
    anc           = "Antenatal care",
    typebirth_cat = "Type of pregnancy",
    mode_del      = "Mode of delivery",
    resus_cat     = "Resuscitation at birth",
    inout_cat     = "Place of birth",
    apgar5_cat    = "Apgar score (5 min)",
    sex_cat       = "Infant sex",
    temp_cat      = "Admission temperature"
)
PRED_ORDER <- names(PRED_LABEL)

# Tidy the level text (strip the wordy parenthetical qualifiers for the axis).
tidy_level <- function(x) {
    x <- gsub("\\s*\\(.*?\\)", "", x)   # drop "(1-3 visits)" etc.
    x <- trimws(x)
    x[x == ""] <- "(ref)"
    x
}

prep <- function(path, group_name) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    df %>%
        filter(predictor %in% PRED_ORDER) %>%
        filter(is.na(note) | note != "Insufficient data") %>%
        mutate(
            group        = group_name,
            pred_label   = PRED_LABEL[predictor],
            level_label  = tidy_level(level),
            is_ref       = !is.na(note) & note == "Reference",
            predictor    = factor(predictor, levels = PRED_ORDER)
        )
}

svn <- prep(svn_path, "All-SVN")

# Row layout: walk the predictors top-to-bottom and reserve one BLANK line for
# each factor header before its levels (reference level first). Giving the
# header its own row means it never collides with a data point or axis label.
svn <- svn %>% arrange(predictor, desc(is_ref))
pos       <- numeric(nrow(svn))   # vertical slot for each level row (1 = top)
hdr_pos   <- numeric(0)           # vertical slot for each factor header
hdr_lab   <- character(0)
cur       <- 0
prev_pred <- NULL
for (i in seq_len(nrow(svn))) {
    p <- as.character(svn$predictor[i])
    if (is.null(prev_pred) || p != prev_pred) {
        cur     <- cur + 1                        # blank line -> header sits here
        hdr_pos <- c(hdr_pos, cur)
        hdr_lab <- c(hdr_lab, svn$pred_label[i])
        prev_pred <- p
    }
    cur      <- cur + 1
    pos[i]   <- cur
}
maxpos       <- cur + 1                            # invert so slot 1 is at the top
svn$row_id   <- maxpos - pos
svn$row_label <- svn$level_label
header_df    <- data.frame(pred_label = hdr_lab, row_id = maxpos - hdr_pos,
                           stringsAsFactors = FALSE)

# Significance flag for colour (non-reference rows only).
svn <- svn %>%
    mutate(sig = case_when(
        is_ref                         ~ "Reference",
        !is.na(p_value) & p_value < 0.05 ~ "p < 0.05",
        TRUE                            ~ "Not significant"
    ))

# X-axis limits: log scale, clamp absurdly wide upper CIs so the plot stays
# readable (flagged points get an arrow at the cap).
XMAX <- 60
XMIN <- 0.1
svn <- svn %>%
    mutate(ci_upper_plot = pmin(ci_upper, XMAX),
           capped        = !is.na(ci_upper) & ci_upper > XMAX)

# header_df (factor headers, one per reserved blank row) was built in the
# layout loop above.

pal <- c("Reference" = "#777777", "p < 0.05" = "#B2182B",
         "Not significant" = "#2166AC")

p_main <- ggplot(svn, aes(x = OR, y = row_id, colour = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper_plot),
                   height = 0.25, na.rm = TRUE) +
    geom_point(aes(shape = is_ref), size = 2.4, na.rm = TRUE) +
    geom_text(data = subset(svn, capped),
              aes(x = XMAX, label = ">"), nudge_y = 0, hjust = 0,
              show.legend = FALSE, size = 3) +
    scale_x_log10(limits = c(XMIN, XMAX),
                  breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32),
                  labels = c("0.1","0.25","0.5","1","2","4","8","16","32")) +
    scale_y_continuous(breaks = svn$row_id, labels = svn$row_label,
                       expand = expansion(add = c(0.8, 1.2))) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16), guide = "none") +
    geom_text(data = header_df, aes(x = XMIN, y = row_id, label = pred_label),
              inherit.aes = FALSE, hjust = 0, fontface = "bold", size = 3.1) +
    labs(x = "Odds ratio for neonatal death (95% CI, log scale)", y = NULL,
         title = "Univariable predictors of neonatal death in the SVN population",
         subtitle = "All-SVN group (Term-SGA, Preterm-AGA, Preterm-SGA); reference levels shown open") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(OUTPUT_DIR_GRAPHS, "figure_06_univariable_ors_allsvn.png"),
       p_main, width = 8.5, height = 10, dpi = 300)
ggsave(file.path(OUTPUT_DIR_GRAPHS, "figure_06_univariable_ors_allsvn.pdf"),
       p_main, width = 8.5, height = 10)
cat("  Saved: figure_06_univariable_ors_allsvn.png / .pdf\n")

# ==============================================================================
# APPENDIX COMPARISON: All-SVN vs Term-AGA, non-reference levels, side by side.
# ==============================================================================
if (file.exists(ref_path)) {
    ref <- prep(ref_path, "Term-AGA")
    both <- bind_rows(svn %>% select(predictor, pred_label, level_label, is_ref,
                                     OR, ci_lower, ci_upper, p_value, group),
                      ref %>% select(predictor, pred_label, level_label, is_ref,
                                     OR, ci_lower, ci_upper, p_value, group)) %>%
        filter(!is_ref, !is.na(OR)) %>%
        arrange(predictor, level_label) %>%
        mutate(lab = paste0(pred_label, ": ", level_label),
               lab = factor(lab, levels = rev(unique(lab))),
               ci_upper_plot = pmin(ci_upper, XMAX),
               group = factor(group, levels = c("All-SVN", "Term-AGA")))

    p_cmp <- ggplot(both, aes(x = OR, y = lab, colour = group)) +
        geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
        geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper_plot),
                       height = 0.3, position = position_dodge(width = 0.6),
                       na.rm = TRUE) +
        geom_point(size = 2.2, position = position_dodge(width = 0.6),
                   na.rm = TRUE) +
        scale_x_log10(limits = c(XMIN, XMAX),
                      breaks = c(0.1, 0.5, 1, 2, 4, 8, 16, 32),
                      labels = c("0.1","0.5","1","2","4","8","16","32")) +
        scale_colour_manual(values = c("All-SVN" = "#B2182B",
                                       "Term-AGA" = "#2166AC"), name = NULL) +
        labs(x = "Odds ratio for neonatal death (95% CI, log scale)", y = NULL,
             title = "Univariable predictors: All-SVN vs Term-AGA") +
        theme_minimal(base_size = 10) +
        theme(panel.grid.minor = element_blank(),
              legend.position = "bottom")

    ggsave(file.path(OUTPUT_DIR_GRAPHS, "figure_06_univariable_ors_svn_vs_termaga.png"),
           p_cmp, width = 9, height = 11, dpi = 300)
    ggsave(file.path(OUTPUT_DIR_GRAPHS, "figure_06_univariable_ors_svn_vs_termaga.pdf"),
           p_cmp, width = 9, height = 11)
    cat("  Saved: figure_06_univariable_ors_svn_vs_termaga.png / .pdf\n")
} else {
    cat("  (table_08_univariable_termaga.csv not found; comparison plot skipped.)\n")
}

cat("\n43_fig4_svn_univariable_forest.R complete.\n")
