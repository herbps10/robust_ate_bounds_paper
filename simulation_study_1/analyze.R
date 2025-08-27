library(tidyverse)
library(kableExtra)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_1")
source(glue::glue("{basepath}/env.R"))
source(glue::glue("{basepath}/../R/manuscript.R"))

results_path <- Sys.getenv("SIMULATION_RESULTS_PATH")
if(results_path == "") stop("Please set SIMULATION_RESULTS_PATH environment variable.")

simulations <- read_rds(glue::glue("{results_path}/simulation_results.rds"))

results_summarized <- simulations |>
  group_by(N, alpha, effect_size, smoothness) |>
  mutate(
    width    = pmin(width, 2),
    width_dr = pmin(width_dr, 2),
  ) |>
  summarize(
    mean_width_dr      = mean(width_dr), median_width_dr = median(width_dr),
    mean_width         = mean(width), median_width = median(width),
    mean_width_uniform = mean(width_uniform), median_width_uniform = median(width_uniform),

    shortest_threshold         = mean(shortest_threshold),
    shortest_uniform_threshold = mean(shortest_uniform_threshold),
    
    coverage_dr   = mean(lower_dr <= true_ate & upper_dr >= true_ate),
    coverage_uniform = mean(lower_uniform <= true_ate & upper_uniform >= true_ate),
    coverage      = mean(lower <= true_ate & upper >= true_ate),
    power_uniform = mean(uniform_test == TRUE),
    power_dr      = mean(dr_test == TRUE),
    power         = mean(test == TRUE),
    n             = n()
  ) |>
  mutate(
    width_diff = mean_width_dr - mean_width_uniform,
    power_diff = power_dr - power_uniform
  ) |>
  arrange(N)

width_comparison <- results_summarized |> 
  select(mean_width_dr, mean_width_uniform, width_diff) |> 
  arrange(width_diff)

power_comparison <- results_summarized |>
  select(power_dr, power_uniform, power_diff) |>
  arrange(power_diff)

table <- results_summarized |>
  ungroup() |>
  select(effect_size, alpha, N, starts_with("mean_width"), starts_with("coverage"), starts_with("power")) |> 
  select(-mean_width, -power, -coverage, -power_diff) |>
  arrange(effect_size, alpha, N) |>
  mutate_at(vars(starts_with("mean")), scales::number_format(accuracy = 0.001)) |>
  mutate_at(vars(starts_with("power")), scales::number_format(accuracy = 0.01)) |>
  mutate_at(vars(starts_with("coverage")), scales::percent_format(accuracy = 0.1)) |>
  #mutate(power = ifelse(effect_size == 0, "", power)) |>
  #mutate(power_dr = ifelse(effect_size == 0, "", power_dr)) |>
  #mutate(power_uniform = ifelse(effect_size == 0, "", power_uniform)) |>
  mutate_at(vars(effect_size, alpha), remove_dups)

table_latex <- table |>
  select(effect_size, alpha, N, mean_width_dr, mean_width_uniform, coverage_dr, coverage_uniform, power_dr, power_uniform) |>
  rowwise() |>
  bold_best(mean_width_uniform, mean_width_dr, 0) |>
  bold_best(coverage_uniform, coverage_dr, 95) |>
  bold_best(power_uniform, power_dr, 1) |>
  knitr::kable(format = "latex", escape = FALSE) |> 
  kable_styling()
