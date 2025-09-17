library(tidyverse)
library(kableExtra)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_1")
source(glue::glue("{basepath}/env.R"))
source(glue::glue("{basepath}/../R/manuscript.R"))

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds")) |>
  filter(index > 1)

q90 <- \(x) quantile(x, 0.90)

results_summarized <- simulations |>
  filter(smoothness > 1e-4) |>
  group_by(N, alpha, effect_size, smoothness, pscore_threshold) |>
  mutate(
    width_onestep = pmin(width_onestep, 2),
    width_uniform = pmin(width_uniform, 2)
  ) |>
  summarize_at(vars(width_onestep, width_uniform), list(mean = mean, median = median, max = max, sd = sd, q90 = q90))

results_summarized <- results_summarized |> left_join(
  simulations |> 
    group_by(N, alpha, effect_size, smoothness, pscore_threshold) |>
    summarize(
      mean_critical_value = mean(critical_value),
      onestep_coverage   = mean(lower_onestep <= true_ate & upper_onestep >= true_ate),
      uniform_coverage   = mean(lower_uniform <= true_ate & upper_uniform >= true_ate),
      uniform_power = mean(uniform_test == TRUE),
      onestep_power = mean(test_onestep == TRUE),
      power         = mean(test == TRUE),
      n             = n()
    ), 
  by = c("N", "alpha", "effect_size", "smoothness", "pscore_threshold")) |>
  mutate(
    width_diff = width_onestep_mean - width_uniform_mean,
    width_q90_diff = width_onestep_q90 - width_uniform_q90,
    width_sd_diff = width_onestep_sd - width_uniform_sd,
    power_diff = onestep_power - uniform_power,
    miscoverage_diff = abs(0.95 - onestep_coverage) - abs(0.95 - uniform_coverage)
  ) |>
  arrange(N)

width_comparison <- results_summarized |> 
  group_by(N, alpha, effect_size, smoothness, pscore_threshold) |>
  select(width_onestep_mean, width_uniform_mean, width_diff, width_onestep_q90, width_uniform_q90, width_q90_diff, width_onestep_sd, width_uniform_sd, width_sd_diff, n) |> 
  arrange(width_diff) |>
  ungroup()

power_comparison <- results_summarized |>
  group_by(N, alpha, effect_size, smoothness, pscore_threshold) |>
  select(onestep_power, uniform_power, power_diff) |>
  arrange(power_diff) |>
  ungroup()

coverage_comparison <- results_summarized |>
  group_by(N, alpha, effect_size, smoothness, pscore_threshold) |>
  select(onestep_coverage, uniform_coverage, miscoverage_diff) |>
  arrange(miscoverage_diff) |>
  ungroup()

table <- results_summarized |>
  ungroup() |>
  select(smoothness, pscore_threshold, effect_size, alpha, N, ends_with("mean"), ends_with("median"), ends_with("q90"), ends_with("sd"), ends_with("coverage"), ends_with("power")) |> 
  #filter(effect_size == 1, alpha %in% c(0, 3)) |>
  arrange(pscore_threshold, effect_size, alpha, N) |>
  mutate_at(vars(ends_with("mean")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(ends_with("median")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(ends_with("q90")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(ends_with("sd")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(ends_with("power")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(ends_with("coverage")), scales::percent_format(accuracy = 0.1)) 

make_width_table <- function(target_alpha = 5, target_pscore_threshold = Inf) {
  table |>
    filter(alpha == target_alpha, pscore_threshold == target_pscore_threshold) |>
    arrange(N, smoothness) |>
    select(N, smoothness, width_uniform_mean, width_onestep_mean, width_uniform_median, width_onestep_median, width_uniform_q90, width_onestep_q90, width_uniform_sd, width_onestep_sd ) |>
    mutate_at(vars(N), remove_dups) |>
    rowwise() |>
    bold_best(width_uniform_mean, width_onestep_mean, 0) |>
    bold_best(width_uniform_median, width_onestep_median, 0) |>
    bold_best(width_uniform_q90, width_onestep_q90, 0) |>
    bold_best(width_uniform_sd, width_onestep_sd, 0) |>
    mutate(smoothness = case_when(
      smoothness == 0.1 ~ "$10^{-1}$",
      smoothness == 0.01 ~ "$10^{-2}$",
      smoothness == 0.001 ~ "$10^{-3}$"
    )) |>
    knitr::kable(format = "latex", escape = FALSE) |> 
    kable_styling()
}

make_coverage_table <- function(target_alpha = 5, target_pscore_threshold = Inf) {
  table |>
    filter(alpha == target_alpha, pscore_threshold == target_pscore_threshold) |>
    arrange(N, smoothness) |>
    select(N, smoothness, uniform_coverage, onestep_coverage, uniform_power, onestep_power) |>
    mutate_at(vars(N), remove_dups) |>
    rowwise() |>
    bold_best(uniform_coverage, onestep_coverage, 95) |>
    bold_best(uniform_power, onestep_power, 1) |>
    mutate(smoothness = case_when(
      smoothness == 0.1 ~ "$10^{-1}$",
      smoothness == 0.01 ~ "$10^{-2}$",
      smoothness == 0.001 ~ "$10^{-3}$"
    )) |>

    knitr::kable(format = "latex", escape = FALSE) |> 
    kable_styling()
}
