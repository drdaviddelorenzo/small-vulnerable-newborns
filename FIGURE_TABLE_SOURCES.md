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
| Table 5 (timing of neonatal death) / Supplementary Table S2 (timing of neonatal death **by site**) | `Table_05_timing_of_death.R` — writes both the by-phenotype breakdown (Table 5) and the by-site breakdown (Supplementary Table S2) in one run; see the note below |
| Figure 5 (timing of death chart) | `Figure_05_timing_of_death_by_SVN.R` |
| Table 6 (length of stay) | `Table_06_length_of_stay.R` |
| Distribution comparisons across SVN phenotype / birth-weight category (case fatality, length of stay) | `helper_Table03_06_kruskal_wallis_cfr_los.R` |
| Table 7a/7b/7c (causes of death, combined / KCH / SMCH) | `Table_07a_causes_of_death.R`, `Table_07bc_causes_of_death_by_site.R` |
| Table 8 (univariable predictors, All-SVN vs Term-AGA) | `Table_08_univariable_predictors_AllSVN.R` |
| Figure 6 / Figure 8 (univariable OR forest plot, All-SVN population) | `Figure_06_univariable_ORs_AllSVN.R` |
| Figure 4b (univariable OR forest plot, USS-SVN sensitivity population) | `Figure_06_univariable_ORs_AllSVN.R` run under `POPULATION="USS"` |
| Table 9a/9b (categorical admission-temperature ORs by phenotype) | `dep_Table09_temperature_categorical_OR.R` (data) + `Table_09_temperature_categorical_ORs.R` + `helper_Table09_temperature_AllSVN_column.R` (adds the pooled All-SVN column) |
| Figure 7 (categorical temperature-OR chart) | `Figure_07_temperature_categorical_ORs_by_SVN.R` |
| Temperature x SVN-phenotype interaction test (admission-temperature paragraph) | `helper_Table09_temperature_svn_interaction.R` |
| **Table 9 — univariable + multivariable regression, SVN population (Term-SGA ref)** | `MV_regression_SVN_population.R` |
| **Table 9 — univariable + multivariable regression, whole cohort (Term-AGA ref)** | `MV_regression_whole_cohort.R` |
| Adjusted continuous (per-1degC) temperature effect, SVN population / whole cohort | `MV_regression_SVN_population_continuous_temperature.R` / `MV_regression_whole_cohort_continuous_temperature.R` — same adjustment set as the two scripts above, with `temp_cat` replaced by the continuous `temp_val` |
| Table 10 (CFR by resuscitation requirement) | `Table_10_cfr_by_resuscitation.R` |
| Table 11 (CFR with/without antenatal corticosteroids) | `Table_11_CFR_by_antenatal_steroids.R` |
| Table 12 (CFR by sex within phenotype) | `Table_12_CFR_by_sex.R` |
| Figure 9 (SVN distribution by HIV exposure) | `dep_Figure09_hiv_distribution_by_svn.R` (data) + `Figure_09_SVN_distribution_by_HIV.R` (chart) |
| Appendix Table Sx (univariable, full vs USS vs non-USS) | `Appendix_Table_Sx_univariable_full_vs_USS_vs_nonUSS.R` |
| Supplementary Table S1 (validation: CFR by SVN phenotype and site, full cohort vs ultrasound-dated subset) | `Supplementary_Table_S1_validation_CFR_full_vs_USS.R` — writes `supplementary_tables1_cfr_by_svn_and_site.csv` (added 2026-08-21; the script previously only produced a both-sites-pooled CFR, not the by-site breakdown the manuscript table actually needs) |

## Resolved: "Supplementary Table S2"

The manuscript's HIV section cites **"Supplementary Table S2"** alongside
"Appendix Figure 9" — but S2's own caption, in the Supplementary material
section, is *"Timing of neonatal death by site"*, which has nothing to do
with HIV. This is a manuscript cross-reference slip, not a gap in this
pipeline: `Table_05_timing_of_death.R` already produces the by-site timing
breakdown (`table_05_timing_by_site.csv` / `_cfr.csv`), and its admission
counts and case-fatality rates match the manuscript's S2 table exactly
(verified against real data 2026-08-21 — every count and CFR matched; e.g.
KCH 0-1 day: 342 deaths, CFR 33 per 1,000). The only cosmetic difference is
the percentage base: the manuscript's S2 shows each cell as a percentage of
that column's *recorded-age* deaths only, while the script's own
console/summary percentages include the "Unknown age" column in the
denominator — the underlying counts are identical either way, this is a
display convention, not a data problem. Whoever reconciles the manuscript
text (see the duplication note above) should either fix the S2 citation in
the HIV paragraph to point at whatever HIV-specific table was actually
intended, or confirm S2 was meant to be about timing all along and drop the
stray HIV cross-reference.

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

## Discussion adjusted-OR citations (numbers reconciliation, 2026-09-08)

Three adjusted odds ratios quoted in the Discussion trace to specific rows of
the multivariable-regression outputs. Recorded here so each cited figure is
checkable against a single CSV cell (all paths relative to
`output/all/tables/`).

| Discussion statement | Value | Script | Output CSV | Row |
|---|---|---|---|---|
| Maternal HIV infection / HIV exposure not independently associated with mortality | aOR 0.96 (95% CI 0.83-1.11) | `MV_regression_SVN_population.R` | `40_svn_multivariable.csv` | `hiv_statusHIV Positive` |
| Outborn status and the hazards of inter-facility transfer | aOR 2.37 (95% CI 2.19-2.56) | `MV_regression_whole_cohort.R` | `40b_wholecohort_multivariable.csv` | `inout_catOutborn` |
| Adjusted odds of death per 1 degC fall in admission temperature | aOR 1.34 (95% CI 1.29-1.38); a 34% rise | `MV_regression_whole_cohort_continuous_temperature.R` | `41b_wholecohort_multivariable_tempcontinuous.csv` | `temp_val` (`OR_per_1degC_decrease` column = 1 / the per-1 degC-rise OR of 0.749) |

The univariable "67% per 1 degC fall (OR 1.67, 95% CI 1.62-1.72)" quoted in
the Results temperature paragraph is the **whole-cohort univariable**
estimate: `41b_wholecohort_univariable_tempcontinuous.csv`, row `temp_val`
(OR 0.599 per 1 degC rise; 1 / 0.599 = 1.669). It is not an adjusted figure;
the SVN-population univariable equivalent is OR 1.56 per 1 degC fall.

Kruskal-Wallis "p<0.001" for case fatality across SVN phenotypes and across
birth-weight categories: `helper_Table03_06_kruskal_wallis_cfr_los.R` ->
`helper_table03_06_kruskal_wallis_cfr_los.csv` (Pearson chi-square also
reported; both p<0.001).

Temperature-by-SVN interaction "chi-square 10.5, df 3, p=0.015":
`helper_Table09_temperature_svn_interaction.R` ->
`helper_table09_temp_svn_interaction.csv` (`lr_chisq` 10.475, `p_value`
0.01493).
