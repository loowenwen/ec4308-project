library(ggplot2)
library(dplyr)
library(lubridate)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)
library(dplyr)
library(gt)
library(reshape2)

#Load Data
predictors_df <- read_csv("../data/predictors_transformed.csv")
cpi_df <- read_csv("../data/cpi_target_full.csv") %>%
  mutate(sasdate = ymd(sasdate))
data_full <- predictors_df %>%
  inner_join(cpi_df, by = "sasdate")

#rename for convenience
data_full <- data_full %>%
  rename(inflation = CPI_t)

# ---------------------- 1. Summary Stats Of CPI_t ----------------------
summary(cpi_df$CPI_t)


# ---------------------- 2. Time Series Plot (CPI_t, 1960-2025) ----------------------

#time range for rolling window evaluation
start_date <- ymd("2017-03-01")
end_date   <- ymd("2025-06-01")

#plot time series
ggplot(cpi_df, aes(x = sasdate, y = CPI_t)) + 
  geom_rect(
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    fill = "lightgray", alpha = 0.02,   
    inherit.aes = FALSE
  ) + 
  geom_line(color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.4) +
  
  # First and last origin: clean dashed lines
  geom_vline(xintercept = c(start_date, end_date), 
             linetype = "dashed", color = "gray30", linewidth = 0.5) +
  
  
  # Minimal label
  annotate("text", x = ymd("2021-06-01"), y = 1.7, 
           label = "Rolling-Window\nEvaluation Period", 
           color = "black", size = 3.5, fontface = "plain", hjust = 0.5) +
  
  scale_x_date(
    limits = c(ymd("1960-01-01"), ymd("2025-06-01") + months(12)),     
    breaks = seq(ymd("1960-01-01"), ymd("2025-01-01"), by = "5 years"), 
    labels = function(x) ifelse(year(x) == 2025, "2025 (Jun)", year(x)),  
    expand = c(0, 0)  # No padding
  ) +
  
  labs(
    title = "U.S. Monthly Inflation Rate (Jan 1960-June 2025)",
    subtitle = "Target Variable: 100 × Δlog(CPIAUCSL)",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "gray70")
  )

# Save plot
ggsave("time_series_plot.png", width = 11, height = 6, dpi = 300, bg = "white")



# ---------------------- 3. Correlation with Target ------------
# Compute correlations with inflation target
corr_with_target <-  data_full %>%
  select(-sasdate) %>%
  summarise(across(-inflation, ~ cor(., data_full$inflation, use = "complete.obs"))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "correlation") %>%
  arrange(desc(abs(correlation)))

# Top correlated predictors
head(corr_with_target, 20)

top_corr <- corr_with_target %>%
  arrange(desc(abs(correlation))) %>%
  head(20) %>%
  mutate(
    correlation = round(correlation, 3),
    Rank = row_number()
  ) %>%
  select(Rank, variable, correlation)

# Identify highest positive & negative correlations 
max_pos_idx <- which.max(top_corr$correlation)
max_neg_idx <- which.min(top_corr$correlation)

# Define alternating row colors 
row_bg <- rep(c("white", "grey95"), length.out = nrow(top_corr))

# Override background for extremes
row_bg[max_pos_idx] <- "#c6dbef" # light blue
row_bg[max_neg_idx] <- "#fcbba1" # light red

#  Define a custom table theme 
table_theme <- ttheme_default(
  core = list(
    fg_params = list(cex = 0.9),
    bg_params = list(fill = row_bg)
  ),
  colhead = list(fg_params = list(fontface = "bold", cex = 1))
)

# Build tableGrob 
corr_table <- tableGrob(top_corr, rows = NULL, theme = table_theme)

# Bold the extreme correlation values 
# Column 3 = correlation
core_fg_idx <- which(corr_table$layout$name == "core-fg")
# grobs[[1]] contains the matrix of cells
for (i in c(max_pos_idx, max_neg_idx)) {
  corr_table$grobs[[core_fg_idx]]$children[[1]]$gp[i, 3] <- gpar(fontface = "bold")
}

# Add title and subtitle 
title <- textGrob(
  "Top 20 Predictors Correlated with CPI Growth",
  gp = gpar(fontsize = 14, fontface = "bold")
)
subtitle <- textGrob(
  "Highest positive (blue) & negative (red) correlations highlighted",
  gp = gpar(fontsize = 10)
)

# Combine title + subtitle + table 
full_table <- grid.arrange(title, subtitle, corr_table, ncol = 1, heights = c(0.15, 0.1, 2))

# Save as PNG 
png("../figures/top20_correlations_colored.png", width = 900, height = 500)
grid.draw(full_table)
dev.off()

cat("✅ Saved PNG table to '../figures/top20_correlations_colored.png'\n")




# ---------------------- 4. Correlation Matrix -----------------
predictor_mat <- predictors_df %>% select(-sasdate)
corr_matrix <- cor(predictor_mat, use = "pairwise.complete.obs")

# Extract upper triangle only
corr_upper <- corr_matrix
corr_upper[lower.tri(corr_upper, diag = TRUE)] <- NA

# Get strong correlations
threshold <- 0.8
high_corr_pairs <- which(abs(corr_upper) > threshold, arr.ind = TRUE)

# Convert to data frame for inspection
high_corr_df <- data.frame(
  var1 = rownames(corr_upper)[high_corr_pairs[,1]],
  var2 = colnames(corr_upper)[high_corr_pairs[,2]],
  correlation = corr_upper[high_corr_pairs]
) %>%
  arrange(desc(abs(correlation)))

# Preview
head(high_corr_df, 20)


corr_melt <- melt(corr_matrix, varnames = c("Var1", "Var2"), value.name = "correlation")

# Remove self-correlations and duplicates (upper triangle only)
corr_melt <- corr_melt %>%
  filter(Var1 != Var2) %>%
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(..1, ..2)), collapse = "_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

# Select top N strongest correlations 
top_n <- 20
top_corr <- corr_melt %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:top_n)

#smmller matrix for heatmap
vars_to_plot <- unique(c(top_corr$Var1, top_corr$Var2))
corr_subset <- corr_matrix[vars_to_plot, vars_to_plot]

corr_subset_melt <- melt(corr_subset)
colnames(corr_subset_melt) <- c("Var1", "Var2", "correlation")  # <--- rename here

#plot-
heatmap_plot <- ggplot(corr_subset_melt, aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#1a9850", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Top 20 Strongest Predictor Correlations",
    subtitle = "Absolute correlation between predictor pairs",
    fill = "Correlation"
  )

#save
png("../figures/top20_predictor_corr.png", width = 1000, height = 800)
print(heatmap_plot)
dev.off()
library(scales)


#######SCREE plots (decide the pc number)############
set.seed("1234")
pca <- prcomp(X_full, center = TRUE, scale. = TRUE) 
var_explained <- pca$sdev^2 / sum(pca$sdev^2) 
cum_var_explained <- cumsum(var_explained) 
variance_table <- data.frame( PC = paste0("PC", 1:length(var_explained)), Variance_Explained = var_explained, Cumulative_Variance = cum_var_explained ) 
print(variance_table)

png("../figures/scree_plot_simple.png",   # change path if you like
    width = 1000, height = 600, res = 140)

plot(pca, type = "l", main = "Scree Plot of PCA Factors (On Stationary Xs)", npcs = 125)
abline(v = 32, col = "blue", lty = 2, lwd = 2)
text(x = 32 + 2, y = max(pca$sdev^2) * 0.8,
     labels = "PC32 (80% Cumulative Variance)", col = "blue", font = 0.4, adj = c(0, 0.5), cex = 1.2)
dev.off()

# Add 80% cumulative line
cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
pc80 <- which(cum_var >= 0.8)[1]
abline(v = pc80, col = "purple", lty = 2)
text(pc80, max(pca$sdev^2)*0.8, paste("PC", pc80, "\n≥80%"), col = "purple", adj = c(-0.2, 1))







