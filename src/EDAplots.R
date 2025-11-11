library(ggplot2)
library(dplyr)
library(lubridate)
library(RColorBrewer)

#######Time Series Plot (target variable over time)############
cpi <- read_csv("../data/cpi_target_full.csv") %>%
  mutate(sasdate = ymd(sasdate))

start_date <- ymd("2017-03-01")
end_date   <- ymd("2025-06-01")

ggplot(cpi, aes(x = sasdate, y = CPI_t)) +
  geom_line(color = "darkblue", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.4) +
  
  
  # First and last origin: clean dashed lines
  geom_vline(xintercept = c(start_date, end_date), 
             linetype = "dashed", color = "gray30", linewidth = 0.5) +
  
  # Minimal label
  annotate("text", x = ymd("2021-06-01"), y = 1.7, 
           label = "Rolling-Window\nEvaluation Period", 
           color = "black", size = 3.2, fontface = "plain", hjust = 0.5) +
  
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  
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
ggsave("EDA_inflation_rolling_window.png", width = 11, height = 6, dpi = 300, bg = "white")







