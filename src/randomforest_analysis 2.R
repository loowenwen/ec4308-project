library(tidyverse)

# define path to rf_results folder
rf_path <- "../data/rf_results"

# define the Z matrices and horizons
Z_list <- c("predictors_raw", "predictors_transformed", "Z_Fstationary", "Z_Ht", "Z_Level_Fstationary",
            "Z_X", "Z_X_MARX", "Z_X_MARX_Fstationary", "Z_X_MARX_Level_Fstationary")

horizons <- c(1, 3, 6, 12)

# ---------- FUNCTION TO READ AND EXTRACT METRICS ----------
read_rf_result <- function(Z, h) {
  file_name <- sprintf("%s/%s_h%d.rds", rf_path, Z, h)
  if (!file.exists(file_name)) return(NULL)
  
  # read the individual result
  res <- readRDS(file_name)
  
  # extract metrics safely
  rmse <- res$metrics$RMSE %||% NA
  mae  <- res$metrics$MAE  %||% NA
  preds <- res$results$yhat
  
  list(pred = preds, rmse = rmse, mae = mae)
}

# ---------- BUILD THE BIG MASTER LIST ----------
rf_master <- list()

for (Z in Z_list) {
  cat("Combining results for:", Z, "\n")
  
  # for each Z, gather all horizons into a list
  horizon_list <- list()
  for (h in horizons) {
    horizon_list[[paste0("h", h)]] <- read_rf_result(Z, h)
  }
  
  # store inside main list
  rf_master[[Z]] <- horizon_list
}

# ---------- SAVE FINAL COMBINED LIST ----------
saveRDS(rf_master, file = "../data/rf_results/rf_master_results.rds")

cat("Combined results saved to ../data/rf_results/rf_master_results.rds\n")





# load combined master results
rf_master <- readRDS("../data/rf_results/rf_master_results.rds")

# ---------- BUILD SUMMARY TABLE ----------
summary_tbl <- map_dfr(names(rf_master), function(Z) {
  Z_entry <- rf_master[[Z]]
  
  map_dfr(names(Z_entry), function(h_key) {
    horizon_num <- as.numeric(str_extract(h_key, "\\d+"))
    rmse_val <- Z_entry[[h_key]]$rmse
    
    tibble(
      Z = Z,
      Horizon = horizon_num,
      RMSE = rmse_val
    )
  })
})

# ---------- COMPUTE RANKS ----------
summary_ranked <- summary_tbl %>%
  group_by(Horizon) %>%
  mutate(Rank = rank(RMSE, ties.method = "first")) %>%
  arrange(Horizon, Rank) %>%
  ungroup()

# ---------- VIEW / SAVE ----------
print(summary_ranked)

# save to .rds
saveRDS(summary_ranked, "../data/rf_results/rf_summary_ranked.rds")

cat("Summary table with RMSE ranking saved!\n")
