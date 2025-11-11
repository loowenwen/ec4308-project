library(ggplot2)
library(dplyr)
library(lubridate)
library(RColorBrewer)
library(caret)

#######Time Series Plot (target variable over time)############
cpi <- read_csv("../data/cpi_target_full.csv") %>%
  mutate(sasdate = ymd(sasdate))

start_date <- ymd("2017-03-01")
end_date   <- ymd("2025-06-01")

ggplot(cpi, aes(x = sasdate, y = CPI_t)) + 
  geom_rect(
    aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
    fill = "lightgray", alpha = 0.02,   # Subtle red, very transparent
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


# Save
ggsave("time_series_plot.png", width = 11, height = 6, dpi = 300, bg = "white")




#######Correlation Heatmap############


predictors_df <- fred_data_values %>% select(-any_of("CPIAUCSL"))
predictors_cov <- predictors_df %>%
  select(-sasdate) %>%            
  select(where(is.numeric)) %>%      
  na.omit()

# -----------------------------
# 1.Look at correlation of all predictors with target -> then pick top 30 -> plot
# -----------------------------

abs_corr <- apply(predictors_cov, 2,
                 function(x) abs(cor(x, cpi$CPI_t,
                                     use = "pairwise.complete.obs")))
top30_names <- names(sort(abs_corr, decreasing = TRUE))[1:30]

predictors_top30 <- predictors_clean %>% select(all_of(top30_names))
cormat <- cor(predictors_top30, use = "pairwise.complete.obs")
cormat <- round(cormat, 2)
melted <- melt(cormat, na.rm = TRUE)

#plot
p <- ggplot(melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#d73027",    # red
    mid = "white",
    high = "#1a9850",   # green (or use "#4575b4" for blue)
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  geom_text(aes(label = value), color = "black", size = 2.3, fontface = "plain") +
  labs(
    title = "Correlation Heatmap: Top 30 Predictors Highly Correlated With CPI_t",
    subtitle = "FRED-MD transformed series, 1961M01–2025M06",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, color = "gray50", hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.6, "cm")
  ) +
  coord_fixed()
print(p)

ggsave("correlation_heatmap.png", width = 20, height = 10, dpi = 500, bg = "white")



# -----------------------------
# 2.Look at correlation within predictors only -> pairs (|r| > 0.95 but < 1) -> extract the unique variables -> do heatmap on them
# -----------------------------
corr_matrix <- cor(predictors_cov, use = "pairwise.complete.obs")
upper_tri <- corr_matrix[upper.tri(corr_matrix)]
high_corr_indices <- which(abs(corr_matrix) > 0.95 & abs(corr_matrix) < 1, arr.ind = TRUE)
high_corr_indices <- high_corr_indices[high_corr_indices[,1] < high_corr_indices[,2], ]

# extract variable names and correlation values
high_pairs <- data.frame(
  Var1 = rownames(corr_matrix)[high_corr_indices[,1]],
  Var2 = colnames(corr_matrix)[high_corr_indices[,2]],
  Correlation = corr_matrix[high_corr_indices]
)
# Subset correlation matrix to variables involved in high-correlation pairs
vars_to_plot <- unique(c(high_pairs$Var1, high_pairs$Var2))
corr_subset <- corr_matrix[vars_to_plot, vars_to_plot]

# Melt for ggplot
melted_corr <- melt(corr_subset, na.rm = TRUE)

# -----------------------------
# Optional: reorder by hierarchical clustering
# -----------------------------
hc <- hclust(dist(corr_subset), method = "complete")
ord <- rownames(corr_subset)[hc$order]
melted_corr$Var1 <- factor(melted_corr$Var1, levels = ord)
melted_corr$Var2 <- factor(melted_corr$Var2, levels = rev(ord))

#plot
ggplot(melted_corr, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#d73027", mid = "white", high = "#1a9850",
    midpoint = 0, limits = c(-1, 1), name = "Correlation"
  ) +
  geom_text(aes(label = round(value, 2)), color = "black", size = 2.5) +
  labs(
    title = "Heatmap of Highly Correlated Predictors",
    subtitle = paste0("Pairs with |r| > 0.95 (excluding self-correlations);\n", length(vars_to_plot), " predictors involved"),
    x = "Predictor 1", y = "Predictor 2"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(color = "gray50", hjust = 0.5, size = 9.5),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed()

ggsave("highlycorrelated_predictors_heatmap.png", width = 20, height = 10, dpi = 500, bg = "white")
