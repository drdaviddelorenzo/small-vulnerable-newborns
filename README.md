# small-vulnerable-newborns

Analysis code for the study **"Small vulnerable newborn outcomes across
Malawi and Zimbabwe: secondary analysis of a prospective neonatal cohort"**
(Neotree, Kamuzu Central Hospital and Sally Mugabe Central Hospital,
2022-2025).

This repository contains the complete R analysis pipeline so that every
table and figure in the manuscript is reproducible. The underlying clinical
data are confidential and are **not** included — see
[DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

## What the study does

Compares the INTERGROWTH-21st small-vulnerable-newborn (SVN) framework
(Term-AGA, Term-SGA, Preterm-AGA, Preterm-SGA) against birth-weight-alone
(low birth weight, LBW) classification for describing outcomes among live
births and neonatal-unit admissions, examines independent predictors of
neonatal death (univariable and multivariable), and tests whether SVN
classification derived from routinely estimated gestational age holds up
when validated against an ultrasound-dated subset.

## Repository layout

```
scripts/                  All R analysis and figure scripts
input/                    Confidential data goes here (not committed) -- see input/README.md
output/                   Pipeline outputs land here (not committed)
run_all.R                 Runs the whole pipeline end to end
FIGURE_TABLE_SOURCES.md   Map of each manuscript table/figure -> source script
DATA_AVAILABILITY.md      Data/code availability statement and access terms
```

## Pipeline

- Every script sources `scripts/00a_setup_shared.R`, which defines the
  INTERGROWTH-21st SGA lookup, SVN classification, plausible-range
  outlier handling, colour palettes, plotting theme, and the six data-file
  paths (all resolved relative to the script's own location, so the repo
  can be cloned anywhere).
- Most scripts read only the raw `input/` data. A small number of
  `Figure_*`/`helper_*` scripts read a CSV written by an earlier
  `Table_*`/`dep_*` script in the same run — `run_all.R` runs everything in
  the correct order automatically.
- Two scripts (`MV_regression_SVN_population.R`,
  `MV_regression_whole_cohort.R`) fit the multivariable logistic regression
  models reported in the manuscript's predictor tables, alongside the
  univariable-only tables produced by the rest of the pipeline.
- `run_all.R` runs the whole pipeline twice: once for the full
  gestational-age-valid cohort (`POPULATION = "all"`) and once for the
  ultrasound-dated subset (`POPULATION = "USS"`).

### Run

```bash
# 1. Place the six required data files in input/ -- see input/README.md
# 2. From the repository root:
Rscript run_all.R
```

`run_all.R` checks `input/` is complete before doing anything else, installs
any missing R packages, and prints a running log plus a final pass/fail
summary. A failure in one script does not stop the rest of the run.

Outputs are written to `output/all/{tables,graphs}/` and
`output/uss/{tables,graphs}/`.

## Key analytical choices

- SGA: birth weight below the sex-specific INTERGROWTH-21st 10th centile
  (unisex table where sex is missing); preterm: <37 weeks gestation.
- The multivariable regression models exclude birth weight and gestational
  age as predictors (they define the SVN phenotype) and pre-eclampsia
  (not estimable within the SVN-only population); missing values for
  variables with substantive missingness are retained as an explicit
  "Not recorded" level, other predictors are analysed on complete cases.
- Antenatal corticosteroid analysis is restricted to the clinically eligible
  window (24+0 to 34+6 weeks gestation).
- Admission temperature is modelled both continuously (per 1 degC) and
  categorically (reference 36.5-37.5 degC).
- All scripts are ASCII-only (a legacy requirement from the UCL Data Safe
  Haven environment the analysis was originally run in).

## Citation

If you use this code, please cite the manuscript (in preparation) and this
repository — see [CITATION.cff](CITATION.cff) and
[DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).
