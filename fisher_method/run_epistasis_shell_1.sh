#!/bin/bash
#SBATCH -J epi_mem_test       ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics-himem                ## Queue
#SBATCH -t 48:00:00             ## Walltime/duration of the job
#SBATCH --mem=0           ## Total memory in GB needed for a job.
#SBATCH --cpus-per-task=1           ## Number of processors for your task

module load R/4.0.3


Rscript --vanilla Parallel_GI_Itot.R epistasis_flat_file.tsv 10
