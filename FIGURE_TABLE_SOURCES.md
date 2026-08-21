# Figure / table sources

Maps each item in the manuscript to the script that produces it. Every
script runs for both `POPULATION = "all"` and `POPULATION = "USS"`; unless
noted, the manuscript uses the `all` output and the `USS` output is the
sensitivity/validation counterpart.

**Manuscript numbering note:** as of the current working draft
(`SVN_Paper_V16b.docx`), several sections appear twice in the document under
different Appendix/main-text numbering (e.g. the temperature-OR tables,
resuscitation CFR table, steroids CFR table, and the univariable-OR forest
plot appear once as "Table/Figure N" and again as "Appendix Table/Figure N").
This is a known, not-yet-resolved issue in the manuscript text itself (not a
gap in this pipeline) — each duplicated item below is produced **once** by
this pipeline; the manuscript's duplicate labelling needs reconciling
separately.

| Manuscript item | Source script(s) |
|---|---|
| Table 1 (maternal characteristics + stillbirth by phenotype) | `Table_01a_maternal_characteristics.R`, `Table_01b_stillbirth_rate_and_OR.R`, `helper_Table01_SVN_by_site_and_stillbirth_OR.R` |
| Figure 1 (CONSORT-style cohort flowchart) | `Figure_01_cohort_flowchart.R` |
| Figure 2 (stillbirth odds ratio by SVN phenotype) | `Figure_02_stillbirth_OR_by_SVN.R` |
| Figure 2 (SVN distribution: births vs admissions — appears in the sensitivity/validation section under the same "Figure 2" label) | `Figure_02_SVN_distribution_births_vs_admissions.R` |
| Table 2 (NNU population characteristics) | `Table_02_NNU_characteristics.R` |
| Table 3 (CFR by phenotype x birth-weight category) / Table 4 (admissions by phenotype x birth-weight category) | `Table_03_04_distribution_and_CFR.R` |
| Figure 3 (CFR by SVN phenotype and site) | `dep_Figure03_CFR_comparison.R` (data) + `Figure_03_CFR_by_SVN_and_site.R` (chart) |
| Table 5 (timing of neonatal death) | `Table_05_timing_of_death.R` |
| Figure 5 (timing of death chart) | `Figure_05_timing_of_death_by_SVN.R` |
| Table 6 (length of stay) | `Table_06_length_of_stay.R` |
| Table 7a/7b/7c (causes of death, combined / KCH / SMCH) | `Table_07a_causes_of_death.R`, `Table_07bc_causes_of_death_by_site.R` |
| Table 8 (univariable predictors, All-SVN vs Term-AGA) | `Table_08_univariable_predictors_AllSVN.R` |
| Figure 6 / Figure 8 (univariable OR forest plot, All-SVN population) | `Figure_06_univariable_ORs_AllSVN.R` |
| Figure 4b (univariable OR forest plot, USS-SVN sensitivity population) | `Figure_06_univariable_ORs_AllSVN.R` run under `POPULATION="USS"` |
| Table 9a/9b (categorical admission-temperature ORs by phenotype) | `dep_Table09_temperature_categorical_OR.R` (data) + `Table_09_temperature_categorical_ORs.R` + `helper_Table09_temperature_AllSVN_column.R` (adds the pooled All-SVN column) |
| Figure 7 (categorical temperature-OR chart) | `Figure_07_temperature_categorical_ORs_by_SVN.R` |
| **Table 9 — univariable + multivariable regression, SVN population (Term-SGA ref)** | `MV_regression_SVN_population.R` |
| **Table 9 — univariable + multivariable regression, whole cohort (Term-AGA ref)** | `MV_regression_whole_cohort.R` |
| Table 10 (CFR by resuscitation requirement) | `Table_10_cfr_by_resuscitation.R` |
| Table 11 (CFR with/without antenatal corticosteroids) | `Table_11_CFR_by_antenatal_steroids.R` |
| Table 12 (CFR by sex within phenotype) | `Table_12_CFR_by_sex.R` |
| Figure 9 (SVN distribution by HIV exposure) | `dep_Figure09_hiv_distribution_by_svn.R` (data) + `Figure_09_SVN_distribution_by_HIV.R` (chart) |
| Appendix Table Sx (univariable, full vs USS vs non-USS) | `Appendix_Table_Sx_univariable_full_vs_USS_vs_nonUSS.R` |
| Supplementary Table S1 (validation: CFR, full vs USS) | `Supplementary_Table_S1_validation_CFR_full_vs_USS.R` |

## Open item — not yet resolved

The manuscript's HIV section references a **"Supplementary Table S2"**
alongside "Appendix Figure 9" — no script in this repository has been
identified as its source yet. This needs chasing down as part of the
manuscript-numbering cleanup (see the duplication note above); it may turn
out to be produced already under a different name, or it may be a genuine
gap. Flagging here rather than guessing.

## Regression note

`MV_regression_SVN_population.R` and `MV_regression_whole_cohort.R` were
brought into this pipeline from an earlier working-copy stage
(`SCRIPTS v4 - MANUSCRIPT MODIFICATIONS JULY2026/40_svn_multivariable_regression.R`
and `40b_wholecohort_regression_termAGA_ref.R`) because they are the only
scripts that reproduce the manuscript's multivariable regression numbers —
verified against the manuscript text to 2 decimal places on every adjusted
odds ratio checked. `Table_08_univariable_predictors_AllSVN.R` is a
univariable-only table (a separate, later analytical decision) and is not a
substitute for these two scripts.
