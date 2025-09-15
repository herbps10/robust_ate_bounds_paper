#
# Simulation Study
#

library(tidyverse)

root <- rprojroot::is_git_root
basepath <- root$find_file("simulation_study_1")

# Load simulation files
source(glue::glue("{basepath}/simulate.R"))
source(glue::glue("{basepath}/wrapper.R"))
source(glue::glue("{basepath}/env.R"))

cache_path <- Sys.getenv("SIMULATION_CACHE_PATH")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(cache_path == "") stop("Please set SIMULATION_CACHE_PATH environment variable.")
if(task_id == "") stop("Task id not set. Please set SLURM_ARRAY_TASK_ID, or run simulations through a Slurm job array, which will set this environment variable for you.")

task_id <- as.numeric(task_id)

N_simulations <- 50
simulations <- expand_grid(
  index = (task_id * N_simulations):((task_id + 1) * N_simulations - 1),
  N = c(100, 250, 500, 1e3),
  alpha = c(1, 5),
  beta = c(0.1),
  pscore_threshold = c(Inf, 0.01),
  effect_size = 1,
) |>
  filter(!(alpha == 1 & pscore_threshold == 0.01))

simulations <- simulations |>
  mutate(
    seed = index
  ) |>
  mutate(
    data = pmap(list(seed, N, alpha, beta, effect_size, pscore_threshold), simulate_data),
    path = glue::glue("{cache_path}/{task_id}.rds"),
    fits = pmap(list(index, N, alpha, beta, effect_size, pscore_threshold, data, path), wrapper)
  )
