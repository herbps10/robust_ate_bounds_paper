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

# Load R files
source(root$find_file("R/eifs.R"))
source(root$find_file("R/helpers.R"))
source(root$find_file("R/smooth_bounds.R"))

cache_path <- Sys.getenv("SIMULATION_CACHE_PATH")
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(cache_path == "") stop("Please set SIMULATION_CACHE_PATH environment variable.")
if(task_id == "") stop("Task id not set. Please set SLURM_ARRAY_TASK_ID, or run simulations through a Slurm job array, which will set this environment variable for you.")

N_simulations <- 1
simulations <- expand_grid(
  index = as.numeric(task_id),
  N = c(100, 250, 500, 1e3),
  alpha = 0:3,
  beta = c(0.1),
  effect_size = 0:3,
  smoothness = 0
) 

simulations <- simulations %>%
  mutate(
    seed = index
  ) %>%
  mutate(
    data = pmap(list(seed, N, alpha, beta, effect_size), simulate_data),
    path = glue::glue("{cache_path}/{N}/{alpha}_{beta}_{effect_size}_{smoothness}_{index}.rds"),
    fits = pmap(list(index, N, alpha, beta, effect_size, smoothness, data, path), wrapper)
  )
