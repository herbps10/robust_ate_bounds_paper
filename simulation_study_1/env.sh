#/bin/bash

# Path to main simulation study R script
export SIMULATION_STUDY_SCRIPT_PATH="/gpfs/home/susmah01/robust_ate_bounds_paper/simulation_study_1/simstudy.R"

# Path to cache folder for temporary results
export SIMULATION_CACHE_PATH="/gpfs/scratch/susmah01/robust_ate_bounds_paper/cache"

# Path to log folder
export SIMULATION_LOG_PATH="/gpfs/scratch/susmah01/robust_ate_bounds_paper/logs"

# Path to results folder for final results
export SIMULATION_RESULTS_PATH="/gpfs/home/susmah01/robust_ate_bounds_paper/simulation_study_1/results"

# Slurm job configuration

export JOB_NAME="atebounds" 

export PARTITION=cpu_short # SLURM partition

export MAIL_USER=herbert.susmann@nyulangone.org

export CPUS_PER_TASK=1 # CPU to assign to each task in the array

export MEM_PER_CPU=2GB # Memory assigned per CPU

export TIME=12:00:00 # Job time limit

export ARRAY=1-100 # Array indexes to run
