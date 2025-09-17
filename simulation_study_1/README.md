# Simulation Study
This folder contains the simulation study replication materials. The simulation code is designed to be run on a SLURM parallel computing cluster.  The structure of the simulation study code is based on [this template](https://github.com/herbps10/simulation_study_template).

The results reported in the paper are saved as `results/simulation_results.rds`.

The main simulation files are:
- `simstudy.R`: generates the simulation datasets and runs the statistical analysis on each one.
- `wrapper.R`: implements the statistical analysis for each dataset.
- `collect.R`: gathers the cached results from all the workers and combines them into one main results file.
- `analyze.R`: loads the main results file, summarizes the results, and generates any tables and figures.

Configuration is done by setting environment variables in the `env.sh` file.

The `sim` shell script gathers several common tasks into one place. 
```
Usage: ./sim <command>

Available commands:
  start     Start batch job
  cancel    Cancel running batch job
  reset     Clear logs
  clean     Clear cache
  status    Print current job status
  watch     Watch current job status
  logs      Print all logs
  cache     List existing cache files
  setup     Create cache and log directories
  help      Show this help message
```

To rerun the simulation study, start the batch job by running:
```
./sim start
```

To generate the tables from the manuscript, without rerunning the full analysis and using the saved simulation results included in the repository, simply run the script `analyze.R` and use the functions `make_coverage_table` and `make_width_table`:

```
# Table 1
make_width_table(target_alpha = 5, target_pscore_threshold = 0.01)
make_width_table(target_alpha = 5, target_pscore_threshold = Inf)

# Table 2
make_coverage_table(target_alpha = 5, target_pscore_threshold = 0.01)
make_coverage_table(target_alpha = 5, target_pscore_threshold = Inf)

# Appendix Table 3
make_width_table(target_alpha = 1, target_pscore_threshold = Inf)

# Appendix Table 4
make_coverage_table(target_alpha = 1, target_pscore_threshold = Inf)
```
