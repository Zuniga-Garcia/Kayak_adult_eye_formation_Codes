library(tidyverse)
library(ggpubr)
library(scales)

# upload data
data <- read_csv("rdf_sin_retinin.csv")

# Convert genotype to a factor (be sure that"control" be the first level ) 
data$genotipo <- as.factor(data$genotipo)
data$genotipo <- fct_relevel(data$genotipo, "control")

# Colors
genotype_levels <- levels(data$genotipo)
color_palette <- setNames(
  c("goldenrod2", "#08306B", "#2171B5", "#4292C6", "#6BAED6")[1:length(genotype_levels)],
  genotype_levels
)

# Function to find the first peak
find_first_peak <- function(r, g_r) {
  # Fiend the first maximum 
  smoothed <- predict(loess(g_r ~ r, span = 0.2))
  peaks <- which(diff(sign(diff(smoothed))) == -2) + 1
  if (length(peaks) == 0) return(list(position = NA, height = NA))
  return(list(position = r[peaks[1]], height = g_r[peaks[1]]))
}

# Compute the statistics
results <- data %>%
  group_by(genotipo, experimento) %>%
  group_modify(~ {
    # Order for radial distance
    sorted_data <- .x %>% arrange(r_um)
    r <- sorted_data$r_um
    g_r <- sorted_data$g_r
    g_r_rand <- sorted_data$g_r_rand
    
    # Find peaks
    real_peak <- find_first_peak(r, g_r)
    rand_peak <- find_first_peak(r, g_r_rand)
    
    # Compute adjust metrics
    r_squared <- 1 - sum((g_r - g_r_rand)^2)/sum((g_r - mean(g_r))^2)
    
    tibble(
      peak_pos = real_peak$position,
      peak_height = real_peak$height,
      rand_peak_pos = rand_peak$position,
      rand_peak_height = rand_peak$height,
      r_squared = r_squared,
      correlation = cor(g_r, g_r_rand, use = "complete.obs")
    )
  }) %>%
  ungroup()

# Compute agregate statistics
stats <- results %>%
  group_by(genotipo) %>%
  summarise(
    peak_pos_mean = mean(peak_pos, na.rm = TRUE),
    peak_height_mean = mean(peak_height, na.rm = TRUE),
    peak_pos_se = sd(peak_pos, na.rm = TRUE)/sqrt(n()),
    peak_height_se = sd(peak_height, na.rm = TRUE)/sqrt(n())
  )

# Graph of mean curves
curve_plot <- data %>%
  group_by(genotipo, r_um) %>%
  summarise(
    mean_g_r = mean(g_r),
    sd_g_r = sd(g_r),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = r_um)) +
  geom_ribbon(aes(ymin = mean_g_r - sd_g_r, 
                  ymax = mean_g_r + sd_g_r, 
                  fill = genotipo), 
              alpha = 0.2) +
  geom_line(aes(y = mean_g_r, color = genotipo), linewidth = 1) +
  geom_point(data = stats, 
             aes(x = peak_pos_mean, 
                 y = peak_height_mean, 
                 color = genotipo),
             size = 4) +
  geom_line(data = data %>% 
              group_by(r_um) %>% 
              summarise(g_r_rand = mean(g_r_rand)),
            aes(y = g_r_rand), 
            color = "black", 
            linetype = "dashed", 
            linewidth = 1) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Corneal nipples RDF",
       x = "μm", 
       y = "g(r)") +
  theme_minimal()

# Graph adjust
r2_plot <- ggplot(results, aes(x = genotipo, y = r_squared, color = genotipo)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "red") +
  geom_hline(yintercept = 0.9, color = "red", linetype = "dashed") +
  scale_color_manual(values = color_palette) +
  labs(title = "Coefficient of Determination (R²)", y = "R²") +
  theme_minimal() +
  theme(legend.position = "none")

# Graph peak positions
peak_pos_plot <- ggplot(results, aes(x = genotipo, color = genotipo)) +
  geom_jitter(aes(y = peak_pos), width = 0.2, size = 3, alpha = 0.7) +
  stat_summary(aes(y = peak_pos), fun = mean, geom = "point", shape = 18, size = 5, color = "red") +
  scale_color_manual(values = color_palette) +
  labs(title = "First peack position", y = "μm") +
  theme_minimal() +
  theme(legend.position = "none")

# Conbine graphs
final_plot <- ggarrange(
  curve_plot,
  ggarrange(r2_plot, peak_pos_plot, ncol = 2),
  nrow = 2,
  heights = c(2, 1)
)

# Final graph
print(final_plot)

# Save results
ggsave("rdf_analysis_final.png", final_plot, width = 12, height = 8, dpi = 300)
write_csv(results, "resultados_peaks_rdf_final.csv")
