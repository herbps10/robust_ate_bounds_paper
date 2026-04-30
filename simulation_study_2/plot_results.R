library(tidyverse)
results <- read_rds("simulation_study_2/results/simulation_study_2_summary.rds")

make_plot <- function(outcome = "width") {
  if (outcome == "width") {
    onestep_outcome <- "width_onestep_mean"
    bounds_outcome <- "width_uniform_mean"
    ylab <- "Mean CI Width"
  } else if (outcome == "coverage") {
    onestep_outcome <- "onestep_coverage"
    bounds_outcome <- "uniform_coverage"
    ylab <- "Empirical 95% Coverage"
  } else if (outcome == "power") {
    onestep_outcome <- "onestep_power"
    bounds_outcome <- "uniform_power"
    ylab <- "Power"
  }

  results |>
    mutate(gamma0 = paste0("gamma[0]==", gamma0)) |>
    ggplot(aes(x = N)) +
    geom_point(aes(y = .data[[onestep_outcome]], color = "DR Wald-type CI")) +
    geom_point(aes(
      y = .data[[bounds_outcome]],
      color = "Non-overlap Bounds CI"
    )) +
    geom_line(aes(y = .data[[onestep_outcome]], color = "DR Wald-type CI")) +
    geom_line(aes(
      y = .data[[bounds_outcome]],
      color = "Non-overlap Bounds CI"
    )) +
    facet_wrap(~gamma0, labeller = label_parsed) +
    scale_x_log10() +
    labs(
      x = expression(N),
      y = ylab,
      color = "",
      title = paste0("Simulation Study 2: ", ylab)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

make_plot("width")
ggsave("simulation_study_2/gamma0_ci_width.pdf", width = 8, height = 4)

make_plot("power")
ggsave("simulation_study_2/gamma0_ci_power.pdf", width = 8, height = 4)

make_plot("coverage")
ggsave("simulation_study_2/gamma0_ci_coverage.pdf", width = 8, height = 4)
