rm(list = ls())

library(ggplot2)
library(afex)
library(tidyverse)
library(simr)

# ---- create the simulation function ----
simfun_rmanova <- function(n, f){

  sd_error <- 1
  effect_mean <- f * sd_error

  numerosity_label <- c("8", "12", "16", "20")
  numerosity_level <- c(1, 2, 3, 4) # assume equal distance
  numerosity_factor <- numerosity_level / sd(numerosity_level) * effect_mean

  # main effect by distractor
  distractor_label <- c("non_dis", "dis")
  distractor_level <- c(0, -1) # distractor supress the gamma burst
  distractor_factor <- distractor_level / sd(distractor_level) * effect_mean

  # interaction of numerosity and distractor
  # interaction_level <- c(0, 0, 0, 0, # non_dis
  #                         0, 0, 0, 0) # dis
  # interaction_factor <- interaction_level / sd(interaction_level) * effect_mean

  # create data
  subject_intercept <- rnorm(n, mean = 100)

  design <- expand.grid(
    id = factor(seq(1:n), level = seq(1:n)),
    numerosity = factor(numerosity_label, levels = c(8, 12, 16, 20)),
    trial = factor(distractor_label, level = distractor_label)
  )

  trial_mean <- subject_intercept[as.integer(design$id)] +
    numerosity_factor[as.numeric(design$numerosity)] +
    distractor_factor[as.numeric(design$trial)] # +
    # interaction_factor[(as.numeric(design$trial) - 1) * 4 +
    #                      as.numeric(design$numerosity)]

  # add trial error
  design$gamma <- rnorm(length(trial_mean), mean = trial_mean, sd = sd_error) 

  aov_gamma <- aov_ez(
    id = "id",
    dv = "gamma",
    within = c("numerosity", "trial"),
    data = design
  )

  # plot.new()
  # p <- ggplot(design, aes(x = numerosity,
  #               y = gamma,
  #               color = trial,
  #               group = trial)) +

  #   stat_summary(fun = mean, geom = "line") +
  #   stat_summary(fun = mean, geom = "point", size = 3) +
  #   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1)

  anova(aov_gamma)["numerosity", "Pr(>F)"] < 0.05
  # print(anova(aov_gamma))
  # print(p)
}

n_list <- seq(20, 100, by = 10)
f_list <- seq(0.1, 0.4, length.out = 10)
results <- expand.grid(
  n = n_list,
  f = f_list
)
results$pwr <- NA

format_elapsed <- function(seconds) {
  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  seconds <- seconds %% 60
  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

update_progress <- function(step, total_steps, start_time, width = 30) {
  progress <- step / total_steps
  filled <- round(width * progress)
  bar <- paste0(strrep("=", filled), strrep(" ", width - filled))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  cat(sprintf(
    "\r[%s] %3.0f%% (%d/%d) elapsed: %s",
    bar,
    progress * 100,
    step,
    total_steps,
    format_elapsed(elapsed)
  ))
  flush.console()
}

total_steps <- length(n_list) * length(f_list)
current_step <- 0
start_time <- Sys.time()

for (n in n_list){
  for (f in f_list){
    data <- replicate(10000, simfun_rmanova(n, f))
    idx <- which(results$n == n & results$f == f)
    results$pwr[idx] <- mean(data)

    current_step <- current_step + 1
    update_progress(current_step, total_steps, start_time)
  }
}

results$se <- sqrt(results$pwr * (1 - results$pwr) / 10000)

cat("\n")

plot_data <- results |>
  filter(n >= 20, n <= 80) |>
  mutate(
    n_factor = factor(n, levels = n_list[n_list >= 20 & n_list <= 80]),
    f_factor = factor(round(f, 2), levels = round(f_list, 2)),
    pwr_blue = if_else(pwr >= 0.8, pwr, NA_real_),
    pwr_label = sprintf("%.2f", pwr)
  )

pwr_heatmap <- ggplot(plot_data, aes(x = n_factor, y = f_factor)) +
  geom_tile(aes(fill = pwr_blue), color = "white", linewidth = 0.4) +
  geom_text(
    data = filter(plot_data, pwr < 0.8),
    aes(label = pwr_label),
    color = "grey20",
    size = 5
  ) +
  geom_text(
    data = filter(plot_data, pwr >= 0.8),
    aes(label = pwr_label),
    color = "white",
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    name = "Power",
    low = "#B5D7E8",
    high = "#06407B",
    limits = c(0.8, 1),
    breaks = seq(0.8, 1, by = 0.05),
    na.value = "grey85"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barheight = grid::unit(4, "cm"),
      barwidth = grid::unit(0.7, "cm")
    )
  ) +
  labs(
    x = "Sample size (N)",
    y = "Cohen's f",
    title = "Power Heatmap for Numerosity Effect",
    caption = "Grey cells indicate power < 0.80; blue cells indicate power >= 0.80."
  ) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.caption = element_text(size = 18, hjust = 0.5)
  )

print(pwr_heatmap)
ggsave("power_heatmap.png", pwr_heatmap, width = 10, height = 8, dpi = 300)
