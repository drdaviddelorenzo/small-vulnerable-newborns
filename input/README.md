# input/

This folder is where you place the confidential, patient-level Neotree data
extracts before running `run_all.R`. **These files are not included in this
repository** — see `DATA_AVAILABILITY.md` at the repo root for how to request
access.

Place exactly these six files directly in this folder (no subfolders):

| Filename | Contents |
|---|---|
| `mwi-matched-nnu-data-20220101-to-20251231.csv` | KCH (Malawi) neonatal-unit admissions, matched to a discharge record |
| `zim-matched-nnu-data-20220101-to-20251231.csv` | SMCH (Zimbabwe) neonatal-unit admissions, matched to a discharge record |
| `mwi-maternal-data-20220101-to-20251231.csv` | KCH (Malawi) maternal delivery/outcome records |
| `zim-maternal-data-20220101-to-20251231.csv` | SMCH (Zimbabwe) maternal delivery/outcome records |
| `mwi-master-data-20220101-to-20251231.csv` | KCH (Malawi) — every admission, matched or not (used only for the CONSORT flowchart denominator) |
| `zim-master-data-20220101-to-20251231.csv` | SMCH (Zimbabwe) — every admission, matched or not (used only for the CONSORT flowchart denominator) |

All six cover the study period **1 January 2022 – 31 December 2025** and the
two study facilities: Kamuzu Central Hospital (KCH, Malawi) and Sally Mugabe
Central Hospital (SMCH, Zimbabwe). Column requirements are documented in
`scripts/00a_setup_shared.R`.

Run `Rscript run_all.R` from the repository root once these files are in
place. `run_all.R` checks for all six by exact filename before doing
anything else, and will tell you exactly which are missing if the check
fails.
