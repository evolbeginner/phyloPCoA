# Changelog
All notable changes to this project will be documented in this file.

## [v0.11.3]
### Fixed
- fixed occassional bugs of no output when root age >= 100

## [v0.11.2]
### Improved
- outpus in adonis/ reorganized

## [v0.11.1]
### Added
- Clearer descriptions into `Readme.md`

## [v0.10.8]
### Added
- `prop2.tbl` as the back-transformed proportion table from `P`
- A fourth PCoA panel using Bray-Curtis distance on `prop2.tbl`; each `pcoa.pdf` page now contains four PCoA plots and one tree

## [v0.10.7]
### Added
- `--pagel_lam_mode auto` to compare `global`, `per_feature`, `hierarchical`, and `none` by AIC
- `--pagel_lam_sim` / `--sim_pagel_lam` to simulate fixed or beta-drawn Pagel lambda values in simulation mode
- `pagel_lam_model_selection.tbl` for AIC-based lambda mode selection

## [v0.10.6]
### Improved
- Pagel's lambda improved (beta distr hierarchical, global, none)

## [v0.10.5]
### Added
- discrete trait simulated by CTMC rather than the liability model
- output adonis

## [v0.10.4]
### Improved
- Moved the script history notes from `phyloPCoA.R` into this changelog
### Added
- Pagel's lambda

## [0.10.3] - 2026-02-04
### Added
- Species exclusion via `--species_exclude`
- Species mismatch btwn the tree and otu_table

## 2026-02-02
### Improved
- Able to read more than one state of traits with `-g grp_file -t trait`

## 2026-01-31
### Added
- `get_grp_info()`

## 2025-11-04
### Improved
- All output is written to `outdir`

## 2025-11-01
### Improved
- Correlation with the true values instead of `Rho`

## 2025-10-31
### Added
- `LDA`, `fdr`, and `BDI`

## 2025-10-24
### Fixed
- `Rho` solved by forcing `v %*% t(v)` where `v` is i.i.d. `U(0.5, 1)`

## 2025-10-22
### Improved
- Reorganized into functions
- `R` improved with `nearPD`
- `pcoa` plots with polygon

## 2025-10-21
### Improved
- Greatly improved overall workflow
- `prop`, `abun`, and `log_prop_geomean` now replace `Q`, `Y`, and `P`

## 2025-07-18
### Improved
- Wrapped code to make it easier to read
- Added `-T 30 -B 8 -p 0.2`

## 2025-07-16
### Improved
- Some small update for XQ

## 2025-07-16
### Improved
- Some more updates for XQ
