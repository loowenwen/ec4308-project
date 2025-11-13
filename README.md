# EC4308 Machine Learning and Economic Forecasting

## Research Question
**Can modern machine-learning methods outperform traditional econometric benchmarks at 1-, 3-, 6-, and 12-month-ahead horizons in predicting U.S. inflation?**

## Repository Overview
| Location | What lives here |
| --- | --- |
| `data/` | Raw FRED-MD file (`2025-09-MD.csv`), cleaned splits (`fred_train/test.csv`), engineered targets, and model outputs (`rf_results/`, `ar(p)_results/`, ensemble RDS files, etc.). |
| `figures/` | Exported plots comparing forecasts across methods. |
| `src/` | All scripts/Rmd notebooks for cleaning, feature engineering, model estimation, and ensembling (see highlights below). |

## Data Pipeline (`src/data_cleaning.R`)
- Source: FRED-MD vintage 2025-09 with monthly data from 1960 onward. CPI (`CPIAUCSL`) is transformed to monthly inflation (log first difference, ×100).
- Missingness: drops very sparse series (`ACOGNO`, `TWEXAFEGSMTHx`, `UMCSENTx`, `ANDENOx`, `VIXCLSx`) and forward-fills only short April-2020 gaps in `CP3Mx` and `COMPAPFFx`. Future partial months (July/August 2025) are removed.
- Transformations: applies FRED-MD codes 1–7 with safer log handling, stores the codes in `data/transform_codes_.csv`, and enforces code 5 for CPI at runtime.
- Stationarity & filtering: trims rows lost to differencing, drops remaining `NA`s, runs ADF/UR tests (trend & no-trend) for documentation, and writes transformed predictors (`predictors_transformed.csv`) plus CPI target (`cpi_target_full.csv`).
- Train/test split: keeps pre-2020 data for model training and reserves the last 100 months (2020–2025) for out-of-sample evaluation (`y_train/test.csv`, `fred_train/test.csv`).

## Feature Engineering (`src/feature_engineering.R`)
- Builds the Coulombe et al. Z-matrix families combining CPI lags, stationary predictor lags (`X_t`), raw levels (`Z_Ht`), MARX blocks (lag-averaged panels), and optional level information.
- Generates 12-lag MARX and moving-average factors (MAF), and prepares PCA-ready data for factor-augmented specifications.
- Saves every engineered design matrix to `src/all_Z_matrices.rds`, which is the single entry point for the modeling scripts.

## Modeling Suite
### Econometric Benchmarks
- `src/randomwalk_model.R`: classic random-walk with drift-style forecasts built from CPI lags.
- `src/ar_model.R` & `data/ar(p)_results/`: rolling AR(p) (p = 1…12) with summaries in `summary_best_model.rds`.
- `src/arima_model.Rmd`: rolling auto-ARIMA with exogenous regressors.

### Machine-Learning Models
- Linear penalized models: `src/lasso.Rmd` (plug-in LASSO via `hdm`), `src/lasso temp.R`, `src/en.R`, `src/ridge.R`, and `src/ridge_plots_clean_fixed/` assets. Each loops over Z-matrices and horizons with 100-step rolling windows.
- Tree-based models: `src/randomforest_model.R` trains (i) direct RF on transformed predictors, (ii) RF with rolling PCA factors, and (iii) RF on enriched Z-matrices; outputs live in `data/rf_results/`.
- Gradient boosting: `src/boosting/Boosting_final.R` plus helper scripts explore PCA-on-the-fly boosting with MARX inputs; forecasts are saved as `boosted_all_final.rds` and per-Z prediction files.

### Ensembling & Comparison
- `src/ensemble.Rmd` ingests every model’s `*.rds`, aligns horizons, and builds:
  - Simple equal-weight averages across ARIMA + penalized + tree models.
  - A Granger-Ramanathan (non-negative, sum-to-one) constrained ensemble via `lsei`.
- Best-in-class RMSE/MAE combinations are written to `src/all_best_models_rmse.rds` for reporting.

## Evaluation Design & Key Findings
- **Rolling window:** all models forecast the final 100 observations (2020-01 to 2025-06) using expanding windows; targets are shifted by h = {1, 3, 6, 12} months.
- **Metrics:** RMSE and MAE, always computed on aligned, horizon-specific CPI inflation.
- **Benchmarks:** AR(p) results (`data/ar(p)_results/summary_best_model.rds`) highlight AR(7) as best 1–6 month performer while AR(2–5) dominate the 12-month MAE front—these serve as reference points for machine-learning gains.

| Horizon | Best model + feature set | RMSE | MAE |
| --- | --- | --- | --- |
| 1 month | Ridge regression on level-based `Z_Ht` | 0.2372 | 0.1724 |
| 3 months | BIC LASSO on `Z_Ht` | 0.2432 | 0.1782 |
| 6 months | BIC LASSO on `Z_Ht` | 0.2529 | 0.1826 |
| 12 months | BIC LASSO on `Z_Ht` | 0.2407 | 0.1722 |

*Source: `src/all_best_models_rmse.rds`. Ridge/elastic net also perform strongly at short horizons, while tree-based models shine once factor + MARX features are added. Both ensembles underperform the best single models but provide robustness checks.*

## Reproducing the Pipeline
1. **Install dependencies (R ≥ 4.3):** `tidyverse`, `lubridate`, `tseries`, `urca`, `glmnet`, `hdm`, `randomForest`, `gbm`, `lsei`, `forecast`, `doParallel`, `HDeconometrics`, etc.
2. **Data prep:** `Rscript src/data_cleaning.R` (writes all CSV artifacts to `data/`).
3. **Feature engineering:** `Rscript -e "source('src/feature_engineering.R')"` (produces `src/all_Z_matrices.rds`, `X_t.rds`, `y_t.rds`).
4. **Model families:** render or run the relevant scripts/Rmds, e.g.
   - `Rscript -e "rmarkdown::render('src/lasso.Rmd')"`
   - `Rscript -e "rmarkdown::render('src/arima_model.Rmd')"`
   - `Rscript src/randomforest_model.R`
   - `Rscript src/boosting/Boosting_final.R`
5. **Ensemble + comparison:** `Rscript -e "rmarkdown::render('src/ensemble.Rmd')"` to refresh `src/all_best_models_rmse.rds` and tables.
6. **Visuals (optional):** rerun plotting chunks inside each Rmd; images land in `figures/`.

> ⚠️ The rolling estimators and boosting scripts are multi-hour jobs unless you parallelize (the LASSO/Ridge scripts already register cores via `doParallel`). Consider trimming horizons or feature sets when testing changes.

## Outputs at a Glance
- **Clean data:** `data/fred_train.csv`, `data/fred_test.csv`, `data/predictors_transformed.csv`, `data/cpi_target_full.csv`.
- **Model artifacts:** `data/rf_results/`, `data/ar(p)_results/`, `src/elasticnet_results_100.rds`, `src/ridge_results_100.rds`, `src/rlasso_results.rds`, boosting prediction files, etc.
- **Top-line comparison:** `src/all_best_models_rmse.rds` (table above) plus `figures/` for visuals ready to drop into slides/papers.
