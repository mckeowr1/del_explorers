#!/bin/bash
#SBATCH -J subset      ## Name of job
#SBATCH -A b1042               ## Allocation
#SBATCH -p genomics                ## Queue
#SBATCH -t 01:00:00             ## Walltime/duration of the job
#SBATCH --mem-per-cpu=16G               ## Memory per node in GB needed for a job.
#SBATCH --cpus-per-task=1            ## Number of processors for your task

module purge

conda activate bcftools 1_12

for region in III:2994514-3040846 X:14147700-14152070;
do
bcftools view --regions $region -Ov /projects/b1059/analysis/WI-20210121/isotype_only/WI.20210121.hard-filter.isotype.bcsq.20210401.vcf.gz> ${region}.vcf ;
done
