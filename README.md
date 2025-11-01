# EC4308-PROJECT

## Data Cleaning Summary (src/data_cleaning.R)

This project uses the FRED-MD dataset (CSV: `data/2025-09-MD.csv`). The `src/data_cleaning.R` script performs the following tasks and design decisions:

- Missingness handling:
	- Drops columns with very high NA rates (e.g., `ACOGNO`, `TWEXAFEGSMTHx`, `UMCSENTx`, `ANDENOx`, `VIXCLSx`).
	- Inspects specific problematic dates (April 2020, July/August 2025) and forward-fills short gaps for `CP3Mx` and `COMPAPFFx` only (one-month gap). This is a pragmatic choice and should be justified in analysis notes.

- Transformations (FRED-MD codes 1–7):
	- Implements the standard FRED-MD transforms: levels, differences, logs, log differences, second differences, and 12-month log differences.
	- Safer log handling: log transforms return `NA` for non-positive values to avoid NaNs.
	- The script extracts transformation codes from the first data row (when present) and saves them to `transform_codes.csv` (in the data folder).

- CPI transformation override:
	- Following Prof. Denis recommendation, the script forces `CPIAUCSL` to use transform code `5` (log first difference). The override is applied at runtime and the `transform_codes.csv` file is updated to reflect the change.

- Row trimming & NA removal:
	- Drops initial rows lost to differencing/logging (determined by the maximum required lag across transforms) and then drops any remaining rows with `NA` after transformation.

- Stationarity checks:
	- Runs Augmented Dickey-Fuller tests (via `tseries::adf.test`) and `ur.df` tests with trend; reports variables flagged as non-stationary and saves results to CSVs (`nonstationary_no_trend.csv`, `nonstationary_with_trend.csv`). Note: the script currently uses short lag choices (k = 0); lag selection may be refined.

- Train/test split & outputs:
	- Splits transformed series into training (pre-2020) and testing (2020 onward) sets and writes `fred_train.csv` and `fred_test.csv` to the data folder.
	- Saves additional artifacts for reproducibility: `fred_complete_levels.csv`, `fred_data_transformed.csv`, `transform_codes.csv`, stationarity CSVs.

- Warnings and notes:
	- The script may emit warnings about date parsing if the input date format differs from `%m/%d/%Y`. Consider robust parsing (`lubridate::parse_date_time`) if needed.
	- Log transforms can produce NaNs when applied to non-positive values; the implementation now uses safer checks but some upstream warnings may still appear if raw `log()` is called on such values.