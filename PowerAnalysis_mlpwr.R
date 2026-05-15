rm(list = ls())

library(mlpwr)
library(pwr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(afex)
library(tidyverse)

set.seed(123)

i <- 0
f_list <- seq(0.1, 0.4, length.out = 10)
res_list <- vector("list", length(f_list))
for (f in f_list){
  i <- i + 1

  # create the simulation function
  simfun_rmanova <- function(N){

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
    interaction_level <- c(0, 0, 0, 0, # non_dis
                          1, 2, 3, 4) # dis
    interaction_factor <- interaction_level / sd(interaction_level) * effect_mean

    # create data
    subject_intercept <- rnorm(N, mean = 100)

    design <- expand.grid(
      id = factor(seq(1:N), level = seq(1:N)),
      numerosity = factor(numerosity_label, levels = c(8, 12, 16, 20)),
      trial = factor(distractor_label, level = distractor_label)
    )
    trial_mean <- subject_intercept[as.integer(design$id)] + 
                numerosity_factor[as.numeric(design$numerosity)] + 
                distractor_factor[as.numeric(design$trial)] #+
                #interaction_factor[(as.numeric(design$trial)-1)*4 + as.numeric(design$numerosity)]
    design$gamma <- rnorm(length(trial_mean), mean = trial_mean, sd = sd_error) # add trial error
    # design$gamma <- trial_mean

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

    print(anova(aov_gamma)$'Pr(>F)'[1] < 0.05)
    # print(anova(aov_gamma))
    # print(p)
  }


  res_list[[i]] <- find.design(
    simfun = simfun_rmanova,
    boundaries = c(10, 200),
    power = 0.8
  )
  
}

res <- data.frame(
  f = f_list,
  N = sapply(res_list, function(x) x$final$design$N)
)
plot.new()
res_plot <- ggplot(res, aes(x = f_list, y = N)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = res$N[res$f == 0.2], color = "#D55E00") +
  labs(
    x = "Cohen's f",
    y = "Required N",
    title = "Power = 0.8"
  ) +
  theme_minimal(base_size = 30) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(res_plot)

# ggsave("pwr0.8.png", res_plot)



tar <- res_list[[which(f_list == 0.2)]]
summary(tar)
